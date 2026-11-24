"""STAR 입력 크기·정렬 버퍼·스레드별 RAM을 합산하여 예약을 정정한다."""
from pathlib import Path
import fcntl
import json
import os
import tempfile
import time

from .memory_policy import GIB, usable_memory


def estimate_star_memory(mode,threads,*,reference_bases=4_000_000_000,index_bytes=None,
                         genome_generate_ram=31_000_000_000,bam_sort_ram=0,io_buffer_bytes=80_000_000):
    threads=max(1,int(threads))
    buffers=8*GIB
    # STAR 기본 I/O 30+50MB/스레드에 junction·정렬 제어 버퍼까지 여유를 둔다.
    thread_buffers=threads*max(GIB//4,int(io_buffer_bytes))
    if mode=='genomeGenerate':
        native=max(2*genome_generate_ram,reference_bases*24)
        parts={'reference_and_suffix_buffers':native,'thread_buffers':thread_buffers,
               'other_buffers':buffers}
    elif mode=='alignReads':
        index_bytes=int(index_bytes or reference_bases*10)
        parts={'index_and_two_pass_buffers':2*index_bytes,
               'bam_sort_buffer':int(bam_sort_ram or index_bytes),
               'thread_buffers':thread_buffers,'other_buffers':buffers}
    else:
        raise ValueError('STAR 자원 계획 mode가 올바르지 않음')
    predicted=(sum(parts.values())*115+99)//100
    reserved=max(128*GIB,((predicted+GIB-1)//GIB)*GIB)
    return {'schema':'mfm_star_memory_v1','mode':mode,'threads':threads,'components':parts,
            'margin_percent':15,'memory_bytes':reserved,'reference_bases':reference_bases,
            'index_bytes':index_bytes}


def command_memory_plan(command):
    if not command or not Path(command[0]).name.startswith('STAR') or '--runMode' not in command:
        return None
    def value(option,default=None):
        return command[command.index(option)+1] if option in command else default
    mode=value('--runMode');threads=int(value('--runThreadN','1'))
    if mode=='alignReads':
        directory=Path(value('--genomeDir'))
        size=sum((directory/name).stat().st_size for name in ['Genome','SA','SAindex'])
        io_buffer_bytes = 80_000_000
        if '--limitIObufferSize' in command:
            begin = command.index('--limitIObufferSize') + 1
            io_buffer_bytes = sum(int(x) for x in command[begin:begin+2])
        return estimate_star_memory(mode,threads,index_bytes=size,
                                    bam_sort_ram=int(value('--limitBAMsortRAM','0')), io_buffer_bytes=io_buffer_bytes)
    if mode=='genomeGenerate':
        begin=command.index('--genomeFastaFiles')+1
        paths=[]
        for item in command[begin:]:
            if item.startswith('--'):break
            paths.append(Path(item))
        if not paths:raise ValueError('STAR index FASTA 입력이 없음')
        bases=0
        for path in paths:
            fai=Path(str(path)+'.fai')
            bases+=sum(int(line.split('\t')[1]) for line in fai.read_text().splitlines()) if fai.exists() else path.stat().st_size
        return estimate_star_memory(mode,threads,reference_bases=bases,
                                    genome_generate_ram=int(value('--limitGenomeGenerateRAM','31000000000')))
    return None


def apply_star_memory(lease_path,plan,*,running=False):
    """계산 결과는 바꾸지 않고 공유 RAM 예약만 갱신한다. 실행 중 증가는 거부한다."""
    path=Path(lease_path);root=path.parent.parent
    with (root/'pool.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        row=json.loads(path.read_text())
        if row.get('branch')!='RNA':return None
        required=int(plan['memory_bytes']);prior=int(row['memory_bytes'])
        if running and required>prior:
            raise ValueError('실행 중 STAR의 RAM 예약 증가는 허용하지 않음')
        others=sum(int(json.loads(p.read_text()).get('memory_bytes',0)) for p in (root/'leases').glob('*.json') if p!=path)
        if others+required>usable_memory(root):
            raise ValueError('STAR 입력 기반 RAM 예측이 공유 예산을 초과함. 계산 시작 전 중단')
        row.update(memory_bytes=required,star_memory_plan=plan,
                   star_previous_memory_bytes=prior,star_memory_updated_at=time.time())
        fd,name=tempfile.mkstemp(prefix='.'+path.name,dir=path.parent)
        try:
            with os.fdopen(fd,'w') as handle:json.dump(row,handle,ensure_ascii=False);handle.write('\n')
            os.replace(name,path)
        finally:Path(name).unlink(missing_ok=True)
        return row


def prepare_star_resources(command):
    lease=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if not lease or os.environ.get('MFM_RNA_MEMORY_PROFILE','input_based')!='input_based':
        return
    plan=command_memory_plan(command)
    if plan is not None:
        row=apply_star_memory(lease,plan)
        if row is not None:print('[STAR-RAM-PLAN] '+json.dumps(plan,ensure_ascii=False),flush=True)
