"""CPU에서 BAM/CIGAR를 읽고 GPU에서 1bp 정수 RNA depth를 누적함."""
import hashlib,time
import numpy as np
import pysam,zarr
from . import gpu_arrays

def intervals(read,start,end,include_deletions=True):
    position=int(read.reference_start)
    for op,length in read.cigartuples or []:
        if op in {0,7,8} or (op==2 and include_deletions):
            lo=max(start,position);hi=min(end,position+length)
            if hi>lo:yield lo-start,hi-start
        if op in {0,2,3,7,8}:position+=length

def contig_depth(task):
    contig,length,bam,samtools,store=task
    from numcodecs import blosc
    blosc.set_nthreads(1)  # contig worker 수와 CPU 압축 thread가 곱으로 늘지 않게 함.
    root=zarr.open_group(store,mode='r+');array=root[contig];chunk=int(array.chunks[0])
    digest=hashlib.sha256();total=nonzero=0;started=time.monotonic()
    with pysam.AlignmentFile(bam,'rb') as source:
        for start in range(0,length,chunk):
            end=min(length,start+chunk);left=[];right=[]
            for read in source.fetch(contig,start,end):
                if read.flag&2820 or read.mapping_quality<255:continue
                # -J: deletion 포함; CIGAR N 제외; duplicate 포함; mate 겹침 두 번 계수함.
                for lo,hi in intervals(read,start,end,True):left.append(lo);right.append(hi)
            if gpu_arrays.enabled('qc') and getattr(gpu_arrays,'FUSED_DEPTH_QC',False):
                values,metrics=gpu_arrays.depth_from_intervals(left,right,end-start,return_metrics=True)
                nz=metrics['depth_nonzero_count']
            else:
                values=gpu_arrays.depth_from_intervals(left,right,end-start)
                from .window_quality import depth_metrics
                metrics=depth_metrics(values,len(values))
                nz=int(np.count_nonzero(values))
            array[start:end]=values
            total+=metrics['depth_sum'];nonzero+=nz
            digest.update(memoryview(values.astype('<u4',copy=False)).cast('B'))
    print(f'[cov:gpu] {contig} depth={total:,} >0={nonzero:,} elapsed={time.monotonic()-started:.1f}s',flush=True)
    return total,nonzero,total,digest.hexdigest()
