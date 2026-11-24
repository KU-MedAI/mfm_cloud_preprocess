"""Contig별 독립 RNA depth worker. 공유 Zarr chunk에 동시에 쓰지 않는다."""
from __future__ import annotations

import hashlib
import subprocess
import tempfile
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import zarr
from numcodecs import blosc

import rna_profile
from rna_resources import worker_count


def _contig_depth(task):
    if len(task)>6 and task[6]:
        from integrate_wgs_rna.rna_depth_native import contig_depth
        return contig_depth(task)
    return _legacy_contig_depth(task[:5])


def _legacy_contig_depth(task: tuple) -> tuple[int, int, int, str]:
    import os
    if os.environ.get('MFM_COVERAGE_BACKEND','cpu')=='cupy':
        from integrate_wgs_rna.gpu_rna_coverage import contig_depth
        return contig_depth(task)
    contig, length, bam, samtools, store = task
    started = time.monotonic()
    # Zarr 압축의 숨은 중첩 worker를 막는다. 코덱/압축률은 그대로다.
    blosc.set_nthreads(1)
    print(f"[cov:start] {contig} len={length:,}", flush=True)
    cov = np.zeros(length, dtype=np.uint32)
    command = [samtools, "depth", *rna_profile.COVERAGE_DEPTH_ARGS, "-r", contig, bam]
    # stderr PIPE를 먼저 못 읽어 외부 도구가 멈추는 위험을 피한다.
    with tempfile.TemporaryFile(mode="w+t") as stderr:
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=stderr,
                                text=True, bufsize=1 << 20)
        try:
            assert proc.stdout is not None
            for line in proc.stdout:
                fields = line.split("\t")
                pos0 = int(fields[1]) - 1
                if 0 <= pos0 < length:
                    cov[pos0] = int(fields[2])
            rc = proc.wait()
            if rc:
                stderr.seek(0)
                raise RuntimeError(f"samtools depth {contig} exit {rc}: {stderr.read()[-2000:]}")
        finally:
            if proc.poll() is None:
                proc.terminate()
                proc.wait()
            if proc.stdout is not None:
                proc.stdout.close()
    # Parent가 만들어 둔 서로 다른 contig 배열만 갱신한다. root attrs는 parent 전용이다.
    zarr.open_group(store, mode="r+")[contig][:] = cov
    mask = cov > 0
    digest = hashlib.sha256(memoryview(cov).cast("B")).hexdigest()
    stats = int(cov[mask].sum()), int(mask.sum()), int(cov.sum()), digest
    print(f"[cov:done] {contig} depth={stats[2]:,} >0={stats[1]:,} "
          f"elapsed={time.monotonic() - started:.1f}s", flush=True)
    return stats


def collect_coverage(contigs: list[str], lengths: dict[str, int], bam: str,
                     samtools: str, store: str, requested: str, snapshot_dir=None) -> tuple[int, int, int, dict[str, str]]:
    if not contigs:
        raise ValueError("coverage 대상 contig가 없습니다")
    import os
    from integrate_wgs_rna.record_stats_parallel import native_tool,NativeUnavailable
    mode=os.environ.get('RNA_DEPTH_BACKEND','auto')
    if mode not in {'auto','native','legacy'}:raise ValueError('RNA_DEPTH_BACKEND는 auto/native/legacy임')
    native=False
    if mode!='legacy':
        try:native_tool('rna_depth');native=True
        except NativeUnavailable:
            if mode=='native':raise
    cpu_cost=1 if native else 2
    jobs=worker_count(requested,len(contigs),cpu_cost=cpu_cost)
    print(f'[cov] contigs={len(contigs)} workers={jobs} cpu_cost={cpu_cost}',flush=True)
    subprocess.run([samtools,'quickcheck',bam],check=True)
    print(f'[cov:engine] native={native} CPU_cost={cpu_cost} peak_diff_bytes={8*(max(lengths[c] for c in contigs)+1)}',flush=True)
    tasks = [(c, lengths[c], bam, samtools, store,str(snapshot_dir) if snapshot_dir and native else None,native) for c in contigs]
    from integrate_wgs_rna.cpu_executor import cpu_results
    from integrate_wgs_rna.memory_policy import GIB,memory_limit
    maximum=max(lengths[c] for c in contigs)
    per_worker=16*maximum+GIB if native else 12*maximum+GIB
    base=4*GIB
    memory_upper=min(1_650_000_000_000,memory_limit())*9//10
    capacity=(memory_upper-base)//per_worker
    if capacity<1:raise ValueError('RNA depth worker의 예상 RAM이 예산을 초과함')
    from integrate_wgs_rna.gpu_compare_resources import cpu_limit
    upper=cpu_limit(int(os.environ.get('MFM_CPU_BUDGET',os.environ.get('MFM_SAMPLE_CPU_BUDGET','68'))))
    if os.environ.get('MFM_CPU_POOL_STAGE_MAX'):upper=min(upper,max(1,int(os.environ['MFM_CPU_POOL_STAGE_MAX'])))
    worker_upper=min(len(tasks),max(1,(upper-1)//cpu_cost),capacity)
    if requested!='auto':worker_upper=min(worker_upper,int(requested))
    jobs=min(jobs,worker_upper)
    scratch=4*sum(lengths[c] for c in contigs) if snapshot_dir and native else 0
    if native and os.environ.get('MFM_COVERAGE_BACKEND','cpu')=='cupy':scratch+=12*sum(sorted((lengths[c] for c in contigs),reverse=True)[:worker_upper])
    if native and shutil.disk_usage(snapshot_dir or tempfile.gettempdir()).free<scratch+GIB:raise ValueError('RNA depth 임시 buffer 공간 부족')
    # tmpfs도 cgroup RAM을 사용하므로 임시파일 상한을 부모 예약에 포함함.
    base+=scratch
    worker_upper=min(worker_upper,(memory_upper-base)//per_worker)
    if worker_upper<1:raise ValueError('RNA depth+임시 buffer 예상 RAM이 예산을 초과함')
    jobs=min(jobs,worker_upper)
    print(f'[cov:preflight] workers_upper={worker_upper} memory_worker={per_worker} memory_peak_upper={base+worker_upper*per_worker} scratch_upper={scratch}',flush=True)
    previous_base=os.environ.get('MFM_MEMORY_STAGE_BASE_BYTES')
    os.environ['MFM_MEMORY_STAGE_BASE_BYTES']=str(base)
    try:
        results={}
        for index,value,error in cpu_results(_contig_depth,[(task,) for task in tasks],jobs,elastic=True,
                max_workers=worker_upper,memory_per_worker_bytes=per_worker,cpus_per_worker=cpu_cost):
            if error:raise error
            results[contigs[index]]=value
            print(f'[cov:progress] {len(results)}/{len(tasks)} completed',flush=True)
    finally:
        if previous_base is None:os.environ.pop('MFM_MEMORY_STAGE_BASE_BYTES',None)
        else:os.environ['MFM_MEMORY_STAGE_BASE_BYTES']=previous_base
    totals = tuple(sum(results[c][i] for c in contigs) for i in range(3))
    return *totals, {c: results[c][3] for c in contigs}
