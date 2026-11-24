"""읽기전용 CIGAR QC를 bounded NumPy chunk와 여유CPU thread로 수행함."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from pathlib import Path
import fcntl
import json
import os
import resource
import uuid

import numpy as np
from ..memory_policy import GIB,usable_memory
from ..pool_leases import prune_leases

MAX_OPS=8_388_608
RECORD_CHUNK=262_144


def record_ranges(arrays,lo,hi):
    offsets=np.asarray(arrays['cigar_offset'][lo:hi],dtype=np.int64)
    sizes=np.asarray(arrays['cigar_len'][lo:hi],dtype=np.int64)
    if hi>lo and int(np.max(offsets+sizes)-np.min(offsets))>MAX_OPS and hi-lo>1:
        middle=(lo+hi)//2
        yield from record_ranges(arrays,lo,middle)
        yield from record_ranges(arrays,middle,hi)
    else:
        yield lo,hi


def validate_chunk(arrays,cigar_ops,lo,hi):
    offsets=np.asarray(arrays['cigar_offset'][lo:hi],dtype=np.int64)
    sizes=np.asarray(arrays['cigar_len'][lo:hi],dtype=np.int64)
    if not sizes.size:return 0
    if np.any(sizes<=0):raise ValueError(f'full validation empty CIGAR for read_id={lo+int(np.flatnonzero(sizes<=0)[0])}')
    op_lo=int(offsets.min());op_hi=int(np.max(offsets+sizes))
    if op_hi-op_lo>MAX_OPS:
        from ..direct_compact_writer import _ReadValidationChunkCache as ReadValidationChunkCache
        cache=ReadValidationChunkCache(arrays,cigar_ops,int(arrays['lengths'].shape[0]),chunk_size=1,max_chunks=1)
        for read_id in range(lo,hi):cache.decode(read_id)
        return hi-lo
    ops=np.asarray(cigar_ops[op_lo:op_hi],dtype=np.uint32)
    if len(ops)!=op_hi-op_lo:raise ValueError('full validation truncated CIGAR')
    starts=offsets-op_lo;ends=starts+sizes
    if os.environ.get('MFM_QC_BACKEND','cpu')=='cupy':
        from ..gpu_arrays import cigar_summary
        return cigar_summary(starts,sizes,ops,np.asarray(arrays['lengths'][lo:hi],dtype=np.int64),
            np.asarray(arrays['ref_end'][lo:hi],dtype=np.int64)-np.asarray(arrays['ref_start'][lo:hi],dtype=np.int64),lo)
    codes=ops&15;counts=(ops>>4).astype(np.int64)
    valid=np.isin(codes,np.asarray([0,1,2,3,6,7,8],dtype=np.uint32))&(counts>0)
    prefix=np.concatenate([np.asarray([0],dtype=np.int64),np.cumsum(~valid,dtype=np.int64)])
    bad=prefix[ends]-prefix[starts]
    if np.any(bad):raise ValueError(f'full validation invalid stored CIGAR op read_id={lo+int(np.flatnonzero(bad)[0])}')
    del prefix,bad,valid
    prefix=np.concatenate([np.asarray([0],dtype=np.int64),np.cumsum(counts*np.isin(codes,[0,1,7,8]),dtype=np.int64)])
    query=prefix[ends]-prefix[starts];del prefix
    expected=np.asarray(arrays['lengths'][lo:hi],dtype=np.int64)
    if np.any(query!=expected):raise ValueError(f'full validation CIGAR/query mismatch read_id={lo+int(np.flatnonzero(query!=expected)[0])}')
    prefix=np.concatenate([np.asarray([0],dtype=np.int64),np.cumsum(counts*np.isin(codes,[0,2,3,7,8]),dtype=np.int64)])
    reference=prefix[ends]-prefix[starts]
    expected=np.asarray(arrays['ref_end'][lo:hi],dtype=np.int64)-np.asarray(arrays['ref_start'][lo:hi],dtype=np.int64)
    if np.any(reference!=expected):raise ValueError(f'full validation CIGAR/reference mismatch read_id={lo+int(np.flatnonzero(reference!=expected)[0])}')
    return hi-lo


@contextmanager
def borrowed_threads(task_count):
    requested=os.environ.get('MFM_QC_CHUNK_THREADS','auto')
    if requested!='auto' and (not requested.isdigit() or int(requested)<1):raise ValueError('MFM_QC_CHUNK_THREADS는 auto 또는 양수여야 함')
    limit=task_count if requested=='auto' else min(task_count,int(requested))
    from numcodecs import blosc
    if blosc.get_nthreads()!=1:limit=1

    soft,_=resource.getrlimit(resource.RLIMIT_AS)
    if soft!=resource.RLIM_INFINITY:
        status=Path('/proc/self/status').read_text();vm=int(next(line.split()[1] for line in status.splitlines() if line.startswith('VmSize:')))*1024
        limit=min(limit,1+max(0,(soft-vm-GIB)//(2*GIB)))
    parent_value=os.environ.get('MFM_CPU_POOL_LEASE_FILE');path=None;workers=1
    if parent_value and limit>1:
        parent=Path(parent_value);root=parent.parent.parent
        with (root/'pool.lock').open('a+') as handle:
            fcntl.flock(handle,fcntl.LOCK_EX);prune_leases(root);owner=json.loads(parent.read_text())
            rows=[json.loads(p.read_text()) for p in (root/'leases').glob('*.json')]
            cpu=max(0,int(owner['total_cpus'])-sum(int(r['cpus']) for r in rows))
            ram=max(0,usable_memory(root)-sum(int(r.get('memory_bytes',0)) for r in rows))
            extra=min(limit-1,cpu,ram//(2*GIB))
            if extra:
                path=root/'leases'/('qc_extra_'+uuid.uuid4().hex+'.json')
                # thread의RSS는parent stage guard가이미포함함. child_pid를넣지않음.
                path.write_text(json.dumps({'sample':owner['sample'],'branch':'QC','stage':'cigar_chunks',
                    'manager_pid':os.getpid(),'child_pid':None,'cpus':extra,'memory_bytes':extra*2*GIB,
                    'total_cpus':owner['total_cpus'],'accounting_parent_lease':str(parent)}))
                workers=1+extra
    elif not parent_value:
        # pool없는호출은알려진CPU예산에서만실행함. 다른sample예산을추측해빌리지않음.
        caps=[max(1,int(os.environ.get('MFM_SAMPLE_CPU_BUDGET','1')))];ram=usable_memory(None)
        for key in ['MFM_CPU_BUDGET','SLURM_CPUS_PER_TASK']:
            if os.environ.get(key):caps.append(max(1,int(os.environ[key])))
        for key,factor in [('SLURM_MEM_PER_NODE',1),('SLURM_MEM_PER_CPU',min(caps))]:
            if int(os.environ.get(key,'0'))>0:ram=min(ram,int(os.environ[key])*1024**2*factor*9//10)
        workers=min(max(1,limit),min(caps),max(1,ram//(2*GIB)))
    try:yield max(1,workers)
    finally:
        if path is not None:
            with (path.parent.parent/'pool.lock').open('a+') as handle:
                fcntl.flock(handle,fcntl.LOCK_EX);path.unlink(missing_ok=True)


def validate_cigars(arrays,cigar_ops,record_count):
    ranges=[]
    for lo in range(0,record_count,RECORD_CHUNK):ranges.extend(record_ranges(arrays,lo,min(record_count,lo+RECORD_CHUNK)))
    with borrowed_threads(len(ranges)) as workers:
        if workers==1:count=sum(validate_chunk(arrays,cigar_ops,lo,hi) for lo,hi in ranges)
        else:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                count=sum(pool.map(lambda bounds:validate_chunk(arrays,cigar_ops,*bounds),ranges))
    return {'mode':('cupy_bounded_chunks' if os.environ.get('MFM_QC_BACKEND','cpu')=='cupy' else 'numpy_bounded_chunks'),'records':count,'chunks':len(ranges),'threads':workers,'max_ops_per_chunk':MAX_OPS}
