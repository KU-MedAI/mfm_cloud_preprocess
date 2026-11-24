"""FASTQ 복사/검사를 bounded native loop로 처리하고 통계를 재사용함."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import gzip,json,os,subprocess
from .record_stats_parallel import native_tool,NativeUnavailable


def process(source,destination=None,mate=None):
    if not Path(source).exists():return {'records':0,'max_read_length':0,'backend':'missing_optional'}
    try:tool=native_tool('rna_fastq')
    except NativeUnavailable:tool=None
    if tool:
        result=subprocess.run([str(tool),str(source),str(destination) if destination is not None else '-',str(mate or 0)],capture_output=True,text=True)
        if result.returncode:raise ValueError('FASTQ 손상 또는 복사/검사 실패: '+result.stderr[-2000:])
        return json.loads(result.stdout)
    count=maximum=0
    from contextlib import nullcontext
    with gzip.open(source,'rb') if str(source).endswith('.gz') else open(source,'rb') as stream:
        with open(destination,'wb') if destination is not None else nullcontext() as out:
            while first:=stream.readline():
                seq=stream.readline();plus=stream.readline();qual=stream.readline()
                if not first.startswith(b'@') or not plus.startswith(b'+') or not qual or len(seq.rstrip(b'\r\n'))!=len(qual.rstrip(b'\r\n')):raise ValueError('FASTQ record 손상')
                name=first.rstrip(b'\r\n')
                if mate is not None:
                    suffix=b'/'+str(mate).encode()
                    if not name.endswith(suffix):raise ValueError('FASTQ mate suffix 누락')
                    name=name[:-2]
                if out:out.write(name+b'\n'+seq+plus+qual)
                count+=1;maximum=max(maximum,len(seq.rstrip(b'\n')))
    return {'records':count,'max_read_length':maximum,'backend':'python_fallback'}


def reserve_conversion(plan):
    lease=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if not lease:return None
    row=json.loads(Path(lease).read_text())
    expected='11' if os.environ.get('MFM_CPU_POOL_LAYOUT')=='legacy' else '13'
    if row.get('branch')!='RNA' or str(row.get('stage'))!=expected or row.get('cpu_backfill'):return None
    from .cpu_executor import WorkerBudget
    budget=WorkerBudget(1,1,True,cpus_per_worker=plan['cpu_upper'])
    budget.base_memory=plan['predicted_peak_ram_bytes']
    if budget.capacity(1,grow=False)!=1:raise ValueError('FASTQ 변환 CPU/RAM 예약 부족')
    return {'reserved_cpus':plan['cpu_upper'],'reserved_memory_bytes':budget.base_memory}


def batch(tasks,cpus):
    if not tasks:return []
    workers=max(1,min(len(tasks),int(cpus)-1));budget=None
    lease=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if lease:
        row=json.loads(Path(lease).read_text())
        expected='11' if os.environ.get('MFM_CPU_POOL_LAYOUT')=='legacy' else '13'
        if row.get('branch')=='RNA' and str(row.get('stage'))==expected:
            from .cpu_executor import WorkerBudget
            budget=WorkerBudget(workers,workers,True,memory_per_worker=64<<20,parent_cpus=1)
            workers=budget.capacity(len(tasks),grow=False)
            if workers<1:raise ValueError('FASTQ 검사 worker CPU/RAM 예약 부족')
            print('[FASTQ-STATS-PLAN] '+json.dumps({'workers':workers,'parent_cpu':1,'cpu_upper':workers+1,'phase':'conversion_complete'}),flush=True)
    try:
        if int(cpus)<=1:return [process(*task) for task in tasks]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            return list(pool.map(lambda args:process(*args),tasks))
    finally:
        if budget is not None:budget.capacity(0,grow=False)
