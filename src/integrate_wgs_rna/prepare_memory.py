"""aligned BAM 변환의 실제 스레드·sort 버퍼로 RAM 피크를 예약한다."""
from pathlib import Path
import fcntl
import json
import os
import time
from .memory_policy import GIB,usable_memory,limit_worker_address_space
from .pool_leases import read_pool_records
from .realign.align_resources import parse_memory


def transcode_plan(cpus,sort_threads,sort_memory):
    sort=(max(0,int(sort_threads))+1)*parse_memory(sort_memory)
    buffers=8*GIB+max(1,int(cpus))*64*1024**2
    required=max(32*GIB,(sort+buffers)*11//10)
    return {'schema':'mfm_aligned_bam_memory_v1','cpus':int(cpus),'sort_threads':int(sort_threads),
            'sort_buffer_bytes':sort,'pipeline_parent_thread_bytes':buffers,
            'margin_percent':10,'reservation_bytes':required,
            'floor_reason':'단일 long aligned BAM: legacy집계28GiB×1.1·index/QC 포함'}


def reserve_transcode(cpus,sort_threads,sort_memory):
    lease=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if not lease:return
    path=Path(lease);root=path.parent.parent;plan=transcode_plan(cpus,sort_threads,sort_memory)
    contract_path=os.environ.get('MFM_PREPARE_INPUT_CONTRACT')
    contract=json.loads(Path(contract_path).read_text()) if contract_path else {}
    source=contract.get('source',{});files=source.get('files',[])
    lower=(contract.get('platform') in {'hifi','ont','pacbio_hifi','nanopore'}
           and os.environ.get('MFM_IO_BACKEND','cpu')=='cpu'
           and source.get('alignment_states')==['aligned'] and len(files)==1
           and files[0].get('format','').lower()=='bam')
    while True:
        with (root/'pool.lock').open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX);row=json.loads(path.read_text())
            required=plan['reservation_bytes'] if lower else max(int(row['memory_bytes']),plan['reservation_bytes'])
            if required>usable_memory(root):raise ValueError('BAM 변환 예상 RAM이 운영 예산을 초과함')
            other=sum(int(r.get('memory_bytes',0)) for p,r in read_pool_records(root/'leases') if p!=path)
            if other+required<=usable_memory(root):
                row.update(memory_bytes=required,prepare_ram_plan=plan)
                temporary=path.with_name('.'+path.name+f'.{os.getpid()}.part')
                temporary.write_text(json.dumps(row,ensure_ascii=False)+'\n');os.replace(temporary,path)
                limit_worker_address_space(required)
                print('[ALIGNED-RAM-PLAN] '+json.dumps(plan,ensure_ascii=False),flush=True)
                return
        time.sleep(.2)
