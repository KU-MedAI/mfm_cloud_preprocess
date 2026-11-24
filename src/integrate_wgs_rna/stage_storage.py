"""대형 BAM 준비 직전 디스크를 확인하고 대기 중 CPU·RAM 예약을 반환한다."""
from pathlib import Path
import fcntl
import json
import os
import tempfile
import time
import sys
from .memory_policy import GIB,usable_memory
from .pool_leases import regular_cpu_claim


def wait_for_storage(contract_path=None):
    value=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if not value:return None
    path=Path(value);root=path.parent.parent
    policy=json.loads((root/'policy.json').read_text())
    floor=int(policy.get('prepare_min_free_bytes',0));folder=policy.get('storage_path')
    if not floor or not folder:return int(json.loads(path.read_text())['cpus'])
    if contract_path:
        contract=json.loads(Path(contract_path).read_text())
        inputs=contract.get('source',{}).get('files',[])
        size=sum(Path(item['path']).stat().st_size for item in inputs)
        floor=max(floor,2*size+50_000_000_000)
    prior=json.loads(path.read_text());label=str(prior.get('branch',''))+':'+str(prior.get('stage',''))
    if label!='WGS:01':return int(prior['cpus'])
    waiting=False
    while True:
        with (root/'pool.lock').open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            policy=json.loads((root/'policy.json').read_text())
            row=json.loads(path.read_text());others=[]
            for p in (root/'leases').glob('*.json'):
                if p==path:continue
                try:others.append(json.loads(p.read_text()))
                except FileNotFoundError:pass
            stat=os.statvfs(folder);free=stat.f_bavail*stat.f_frsize
            active=sum(1 for r in others if r.get('branch')=='WGS' and r.get('stage')=='01'
                       and not r.get('stage_storage_wait') and not r.get('storage_paused')
                       and r.get('prepare_phase')!='read_only_qc')
            capacity=int(policy.get('stage_concurrency',{}).get('WGS:01',1))
            cpu=(int(row['total_cpus'])-sum(int(r.get('cpus',0)) for r in others if r.get('cpu_backfill'))
                 if prior.get('cpu_backfill') else int(row['total_cpus'])-sum(regular_cpu_claim(r) for r in others))
            ram=usable_memory(root)-sum(int(r.get('memory_bytes',0)) for r in others)
            ready=free>=floor and active<capacity
            grant=min(int(prior['cpus']),max(0,cpu))
            minimum=min(6 if os.environ.get('MFM_IO_BACKEND')=='nvcomp' else 4,int(prior['cpus']))
            if ready and (not waiting or grant>=minimum and ram>=int(prior['memory_bytes'])):
                row.update(cpus=int(row['cpus']) if not waiting else grant,memory_bytes=prior['memory_bytes'],
                           stage_storage_wait=False,prepare_phase='writing_bam')
            else:
                row.update(cpus=0,memory_bytes=4*GIB,stage_storage_wait=True,
                           storage_wait_reason='BAM 준비 전 디스크 여유·동시 작업·CPU/RAM 재확보 대기')
                if not waiting:print('[STORAGE-WAIT] 디스크 여유 확보까지 CPU/RAM을 반환하고 대기함',file=sys.stderr,flush=True)
                waiting=True
            fd,temp=tempfile.mkstemp(prefix='.'+path.name,dir=path.parent)
            try:
                with os.fdopen(fd,'w') as handle:json.dump(row,handle,ensure_ascii=False);handle.write('\n')
                os.replace(temp,path)
            finally:Path(temp).unlink(missing_ok=True)
            if ready and not row['stage_storage_wait']:return int(row['cpus'])
        time.sleep(2)


if __name__=='__main__':
    granted=wait_for_storage(sys.argv[1] if len(sys.argv)>1 else None)
    if granted is not None:print(granted)
