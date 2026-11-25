"""Picard heap·native buffer·인덱스/검사 스레드를 실제 lease 안에서 계획함."""
import json,os
from pathlib import Path
from .memory_policy import memory_limit
from .pool_leases import cgroup_cpu_limit
GIB=1024**3


def markdup_memory_bytes(heap_gb):
    if int(heap_gb)<1:raise ValueError('Picard heap은 양수여야 함')
    # JVM/native와 후속 HTSlib buffer는 2GiB를 예약하고10%여유를 더함.
    return ((int(heap_gb)+2)*GIB*11+9)//10


def markdup_plan(heap_gb,input_bam):
    lease={}
    if os.environ.get('MFM_CPU_POOL_LEASE_FILE'):
        lease=json.loads(Path(os.environ['MFM_CPU_POOL_LEASE_FILE']).read_text())
    caps=[len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else os.cpu_count() or 1]
    if cgroup_cpu_limit():caps.append(cgroup_cpu_limit())
    for name in ['MFM_CPU_BUDGET','SLURM_CPUS_PER_TASK']:
        if os.environ.get(name):caps.append(int(os.environ[name]))
    if lease.get('cpus'):caps.append(int(lease['cpus']))
    # HAP worker의 몫은 부모 lease보다 작을 수 있으므로 항상 함께 적용함.
    for name in ['MFM_CPU_POOL_GRANTED','MFM_SAMPLE_CPU_BUDGET']:
        if os.environ.get(name):caps.append(int(os.environ[name]))
    cpu=max(1,min(caps));ram=memory_limit()
    for name,multiplier in [('SLURM_MEM_PER_NODE',1),('SLURM_MEM_PER_CPU',min(caps))]:
        if int(os.environ.get(name,'0'))>0:ram=min(ram,int(os.environ[name])*1024**2*multiplier)
    ram=ram*9//10
    if lease.get('memory_bytes'):ram=min(ram,int(lease['memory_bytes']))
    estimate=markdup_memory_bytes(heap_gb)
    if estimate>ram:raise ValueError('Picard heap/native 예상 RAM이 예약 예산보다 큼')
    return {'cpu_grant':cpu,'java_active_processors':cpu,'htslib_aux_threads':max(0,cpu-1),
            'heap_bytes':int(heap_gb)*GIB,'native_buffer_bytes':2*GIB,'input_bytes':Path(input_bam).stat().st_size,
            'estimated_peak_rss_bytes':estimate,'ram_budget_bytes':ram}
