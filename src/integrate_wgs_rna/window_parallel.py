"""염색체×HAP별 window 주석을 독립 계산하고 원래 순서로 반환함."""
import math
import os
import time
import json
from pathlib import Path
from multiprocessing import get_context
from .cpu_executor import cpu_results
from .input_parallel import cpu_budget
from .memory_policy import GIB,memory_limit,process_memory_snapshot


def _job(task):
    from .window_index import windows_for_fasta
    rows,blacklist,stats=windows_for_fasta(**task)
    return task['hap'],rows,blacklist,stats


def _sequence_bases(task):
    fasta=task['fasta'];index=Path(str(fasta)+'.fai')
    if index.exists():
        lengths={parts[0]:int(parts[1]) for line in index.read_text().splitlines()
                 if (parts:=line.split('\t')) and len(parts)>1}
        for name in [task['contig'],task['contig'].removeprefix('chr')]:
            if name in lengths:return lengths[name]
        if len(lengths)==1:return next(iter(lengths.values()))
        raise KeyError('window RAM예측 대상 contig가FAI에없음: '+task['contig'])
    if fasta.suffix.lower() in {'.gz','.bgz'}:
        raise ValueError('압축FASTA RAM예측에는FAI가필요함: '+str(fasta))
    return fasta.stat().st_size


def run(tasks, *, parallel):
    if not tasks:return
    if not parallel:
        for task in tasks:yield _job(task)
        return
    requested=os.environ.get('MFM_WINDOW_CONTIG_JOBS','auto')
    if requested!='auto' and (not requested.isdigit() or int(requested)<1):
        raise ValueError('MFM_WINDOW_CONTIG_JOBS는 auto 또는 양수여야 함')
    cap=len(tasks) if requested=='auto' else min(len(tasks),int(requested))
    if os.environ.get('MFM_CPU_POOL_STAGE_MAX'):
        cap=min(cap,max(1,int(os.environ['MFM_CPU_POOL_STAGE_MAX'])))
    estimates=[];output_bytes=0
    for task in tasks:
        size=_sequence_bases(task)
        chain=task['chain'].stat().st_size if task['chain'] is not None else 0
        rows=math.ceil(size/task['config'].window_size)
        output_bytes+=rows*4096
        estimates.append(max(4*GIB,GIB+size*8+chain*32+rows*4096))
    worker=max(estimates,default=4*GIB)
    parent=4*GIB+output_bytes
    usable=memory_limit()*9//10
    memory_cap=max(0,(usable-parent)//worker)
    cap=min(cap,memory_cap)
    if cap<1:raise ValueError('window 작업의 부모/worker RAM 예상량이 전체 예산을 초과함')
    address=max(8*GIB,int(process_memory_snapshot().get('vms_bytes') or 0)+worker+GIB)
    workers=max(1,min(cap,cpu_budget()))
    elastic=(requested!='1' and os.environ.get('MFM_CPU_UTILIZATION','maximum')=='maximum'
             and bool(os.environ.get('MFM_CPU_POOL_LEASE_FILE')))
    print('[WINDOW-PARALLEL] '+json.dumps({'phase':'plan','tasks':len(tasks),
        'initial_workers':workers,'worker_cap':cap,'worker_memory_bytes':worker,
        'parent_memory_bytes':parent,'worker_address_space_bytes':address,
        'elastic':elastic}),flush=True)
    started=time.monotonic();results={}
    names=['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']
    previous={name:os.environ.get(name) for name in names}
    try:
        # 각 Python worker는 한 CPU 계산이다. spawn import의 암묵적 BLAS 증식을 막음.
        os.environ.update({name:'1' for name in names})
        for index,result,error in cpu_results(_job,[(task,) for task in tasks],workers,
                elastic=elastic,max_workers=cap,cpus_per_worker=1,
                memory_per_worker_bytes=worker,address_space_per_worker_bytes=address,
                base_memory_bytes=parent,mp_context=get_context('spawn')):
            if error is not None:raise error
            results[index]=result
    finally:
        for name,value in previous.items():
            if value is None:os.environ.pop(name,None)
            else:os.environ[name]=value
    print('[WINDOW-PARALLEL] '+json.dumps({'phase':'complete','tasks':len(results),
                                         'wall_seconds':time.monotonic()-started}),flush=True)
    for index in range(len(tasks)):yield results.pop(index)
