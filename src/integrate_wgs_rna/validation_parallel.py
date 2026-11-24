"""최종 검증의 독립 작업을 현재 CPU/RAM 예약 안에서 실행함."""
from __future__ import annotations
import os
from contextlib import contextmanager
import time
from .cpu_executor import cpu_results, WorkerBudget

GIB=1<<30


def cpu_budget():
    caps=[len(os.sched_getaffinity(0))]
    for name in ['MFM_CPU_POOL_GRANTED','MFM_SAMPLE_CPU_BUDGET','MFM_FIXED_SAMPLE_CPUS','MFM_CPU_POOL_STAGE_MAX','SLURM_CPUS_PER_TASK']:
        if os.environ.get(name):caps.append(max(1,int(os.environ[name])))
    # 예약 없는 독립 검증은 기존 직렬 동작을 유지함.
    return min(caps) if len(caps)>1 else 1



@contextmanager
def hash_workers():
    # 앞선 병렬 검증에서 반환한 CPU를 실제 해시 읽기 전에 재확보함.
    limit=max(1,cpu_budget()-1)
    budget=WorkerBudget(limit,limit,False,memory_per_worker=16<<20,parent_cpus=int(cpu_budget()>1))
    workers=0
    try:
        while workers<1:
            workers=budget.capacity(limit)
            if not workers:time.sleep(.1)
        yield workers
    finally:
        budget.capacity(0,grow=False)


def _sequence_task(directory,source):
    from numcodecs import blosc
    from .sequence_only import validate_contig
    import json
    from pathlib import Path
    if json.loads((Path(directory)/'contig_metadata.json').read_text()).get('wgs_payload')=='pileup_feature':
        from .pileup_features.export import validate_contig
    blosc.set_nthreads(1)
    return validate_contig(directory,source=source)


def sequence_checks(directories, *, source=False, max_workers=None):
    directories = list(directories)
    cpus = cpu_budget()
    workers = min(len(directories), max(1, cpus - 1))
    if max_workers is not None:
        workers = min(workers, max(1, max_workers))
    if not directories:return []
    print(f'[VALIDATION-PLAN] tasks={len(directories)} workers={workers} parent_cpus=1 '
          f'worker_ram_bytes={2*GIB} peak_ram_bytes={(4+2*workers)*GIB}',flush=True)
    results=[None]*len(directories)
    # 부모 1CPU + worker별 1CPU/2GiB. 계산 전에 기존 실행기가 RAM을 재예약함.
    for index,value,error in cpu_results(_sequence_task,[(str(p),source) for p in directories],workers,
            max_workers=workers,memory_per_worker_bytes=2*GIB,parent_cpus=int(cpus>1),
            address_space_per_worker_bytes=4*GIB):
        results[index]=(value,str(error) if error is not None else None)
    return results


def uint32_stats(array):
    import hashlib
    import numpy as np
    if array.ndim!=1 or array.dtype!=np.dtype(np.uint32):
        raise ValueError(f'coverage dtype/shape 불일치: shape={array.shape} dtype={array.dtype}')
    chunk=int(array.chunks[0]) if array.chunks else 1<<20
    total=nonzero=0;digest=hashlib.sha256()
    for start in range(0,int(array.shape[0]),chunk):
        values=np.ascontiguousarray(array[start:start+chunk],dtype='<u4')
        total+=int(values.sum(dtype=np.uint64));nonzero+=int(np.count_nonzero(values))
        digest.update(memoryview(values).cast('B'))
    return total,nonzero,digest.hexdigest()


def _depth_task(store,contig):
    import zarr
    from numcodecs import blosc
    blosc.set_nthreads(1)
    return uint32_stats(zarr.open_group(store,mode='r')[contig])


def depth_checks(store,contigs):
    contigs=list(contigs);workers=min(len(contigs),max(1,cpu_budget()-1))
    results=[None]*len(contigs)
    if not contigs:return results
    print(f'[RNA-VALIDATION-PLAN] tasks={len(contigs)} workers={workers} '
          f'peak_ram_bytes={(4+2*workers)*GIB}',flush=True)
    for index,value,error in cpu_results(_depth_task,[(str(store),c) for c in contigs],workers,
            max_workers=workers,memory_per_worker_bytes=2*GIB,parent_cpus=int(cpu_budget()>1),
            address_space_per_worker_bytes=4*GIB):
        if error is not None:raise error
        results[index]=value
    return results
