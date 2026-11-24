"""입력 준비의 읽기 QC와 대형 BAM 쓰기 예약을 분리한다."""
from pathlib import Path
import fcntl
import json
import os
import tempfile
import sys
from .stage_storage import wait_for_storage


def context():
    value=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if not value:return None
    path=Path(value);root=path.parent.parent;policy_path=root/'policy.json'
    policy=json.loads(policy_path.read_text()) if policy_path.exists() else {}
    if policy.get('prepare_phased') is not True:return None
    row=json.loads(path.read_text())
    if row.get('branch')!='WGS' or row.get('stage')!='01':return None
    return path,root


def mark_read_only_prepare(*, index_threads=None):
    item=context()
    if item is None:return
    path,root=item
    with (root/'pool.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);row=json.loads(path.read_text())
        row['prepare_phase']='read_only_qc'
        if index_threads is not None:
            # 마지막 쓰기가 끝난 색인/QC에서만 줄이며, 이후 큰 RAM 재예약은 없다.
            before=int(row['memory_bytes'])
            peak=4*1024**3+(max(0,int(index_threads))+1)*64*1024**2
            row['memory_bytes']=min(before,peak)
            row['prepare_finish_plan']={'phase':'final_index_qc','estimated_peak_bytes':peak,
                                        'prior_memory_bytes':before,'index_additional_threads':index_threads}
            print('[PREPARE-FINISH-PLAN] '+json.dumps(row['prepare_finish_plan']),file=sys.stderr,flush=True)
        fd,temp=tempfile.mkstemp(prefix='.'+path.name,dir=path.parent)
        try:
            with os.fdopen(fd,'w') as handle:json.dump(row,handle,ensure_ascii=False);handle.write('\n')
            os.replace(temp,path)
        finally:Path(temp).unlink(missing_ok=True)


def reserve_prepare_write():
    if context() is None:return None
    granted=wait_for_storage(os.environ.get('MFM_PREPARE_INPUT_CONTRACT'))
    if granted is not None:
        os.environ.update(MFM_SAMPLE_CPU_BUDGET=str(granted),MFM_CPU_POOL_GRANTED=str(granted))
    return granted


def prepare_entry(contract_path):
    os.environ['MFM_PREPARE_INPUT_CONTRACT']=str(contract_path)
    item=context()
    if item is not None:
        contract=json.loads(Path(contract_path).read_text());source=contract.get('source',{})
        files=source.get('files',[])
        if source.get('alignment_states')==['aligned'] and files and all(r.get('format','').lower()=='bam' for r in files):
            mark_read_only_prepare()
            return int(json.loads(item[0].read_text())['cpus'])
    return wait_for_storage(contract_path)


if __name__=='__main__':
    granted=prepare_entry(sys.argv[1])
    if granted is not None:print(granted)
