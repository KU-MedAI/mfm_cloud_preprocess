"""승인된 GPU 비교에서 예약 대기를 생략하고 실제 사용량만 기록함."""
from contextlib import contextmanager
from pathlib import Path
import fcntl,json,os,resource,signal,subprocess,time
from .gpu_pool import GIB
from .realign.phases import _atomic as _save
from .memory_policy import memory_limit,process_rss,protect_pool_memory,register_memory_observer
from .pool_leases import cgroup_cpu_limit

FLAG='MFM_GPU_COMPARE_UNRESERVED'

def enabled(environment=None):
    env=os.environ if environment is None else environment
    if env.get(FLAG)!='1':return False
    if not str(env.get('RUN_ID','')).startswith('gpu_compare_'):
        raise ValueError('예약 예외는 gpu_compare 실행에서만 허용함')
    cpu_control=(env.get('MFM_CPU_COMPARISON_CONTROL')=='1' and all(env.get(k,'cpu')=='cpu' for k in ['MFM_WGS_BACKEND','MFM_RNA_BACKEND','MFM_BAM_BACKEND']))
    if not cpu_control and not any(env.get(k)=='parabricks' for k in ['MFM_WGS_BACKEND','MFM_RNA_BACKEND','MFM_BAM_BACKEND']):
        raise ValueError('GPU 비교 backend 없이 예약 예외를 사용할 수 없음')
    return True

def config_environment(config):
    env=dict(config.get('pipeline_env',{}),RUN_ID=config.get('run_id',''))
    for key in ['wgs_backend','rna_backend','bam_backend']:
        if key in config.get('acceleration',{}):env['MFM_'+key.upper()]=config['acceleration'][key]
    return env

def cpu_limit(requested=68):
    caps=[requested,len(os.sched_getaffinity(0)),cgroup_cpu_limit() or requested]
    if os.environ.get('SLURM_CPUS_PER_TASK'):caps.append(int(os.environ['SLURM_CPUS_PER_TASK']))
    return max(1,min(caps))

@contextmanager
def sample_slot(config,sample):
    env=config_environment(config)
    if not enabled(env):yield;return
    # GPU 단계의 기존 독점 잠금만 사용하여 CPU 출력/QC 중 다음 샘플을 막지 않음.
    record=Path(config['state_root'])/'queue/gpu_compare_slot.json';record.parent.mkdir(parents=True,exist_ok=True)
    value={'sample':sample,'pid':os.getpid(),'status':'active','started_at':time.time(),
           'cpu_reservation':0,'ram_reservation':0,'serialization':'GPU 도구 단계에서만 적용함'}
    _save(record,value)
    try:yield
    finally:value.update(status='released',finished_at=time.time());_save(record,value)

def run_unreserved(args,command,adjust,audit):
    cpus=cpu_limit(args.total_cpus)
    if os.environ.get('MFM_CPU_POOL_STAGE_MAX'):cpus=min(cpus,max(1,int(os.environ['MFM_CPU_POOL_STAGE_MAX'])))
    if args.branch.upper()=='IO':cpus=min(cpus,6 if args.stage=='bam_roundtrip' else 4)
    effective=adjust(command,args.branch,args.stage,cpus)
    env=os.environ.copy();env.pop('MFM_CPU_POOL_LEASE_FILE',None)
    env.update(MFM_CPU_BUDGET=str(cpu_limit(args.total_cpus)),MFM_SAMPLE_CPU_BUDGET=str(cpus),
               MFM_CPU_POOL_GRANTED=str(cpus),MFM_CPU_POOL_STAGE=args.stage,MFM_DYNAMIC_CPU_POOL_ACTIVE='1',
               OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
    from .stage_cpu_config import apply as apply_stage_cpus
    apply_stage_cpus(env,args.branch,args.stage,cpus)
    root=Path(args.pool_dir);limit=memory_limit(root)
    start=time.time();before=resource.getrusage(resource.RUSAGE_CHILDREN);peak=0
    print(f'[GPU-COMPARE-UNRESERVED] sample={args.sample} stage={args.stage} thread_limit={cpus} CPU예약=0 RAM예약=0',flush=True)
    try:process=subprocess.Popen(effective,env=env,start_new_session=True)
    except PermissionError:
        if not effective[0].endswith('.sh'):raise
        process=subprocess.Popen(['bash',*effective],env=env,start_new_session=True)
    def forward(signum,frame):
        if process.poll() is None:os.killpg(process.pid,signum)
    handlers={s:signal.signal(s,forward) for s in [signal.SIGTERM,signal.SIGINT]}
    observation=None
    try:
        observation=register_memory_observer(root,{'manager_pid':os.getpid(),'child_pid':process.pid,
            'sample':args.sample,'branch':args.branch,'stage':args.stage,'started_at':start})
        while process.poll() is None:
            peak=max(peak,process_rss({'stage':process.pid})['stage'])
            with (root/'pool.lock').open('a+') as lock:
                fcntl.flock(lock,fcntl.LOCK_EX);guard=protect_pool_memory(root)
            limit=guard['limit_bytes']
            if guard.get('action')=='stop_largest_stage' and guard.get('stopped_pid')==process.pid:
                process.wait();raise MemoryError(f"전처리 전체 RSS 보호 종료: {guard['rss_bytes']}/{guard['trigger_bytes']}")
            time.sleep(.5)
        return process.returncode
    finally:
        if process.poll() is None:os.killpg(process.pid,signal.SIGTERM);process.wait()
        if observation is not None:
            with (root/'pool.lock').open('a+') as lock:
                fcntl.flock(lock,fcntl.LOCK_EX);observation.unlink(missing_ok=True)
        for s,h in handlers.items():signal.signal(s,h)
        after=resource.getrusage(resource.RUSAGE_CHILDREN);finished=time.time()
        audit({'schema':'mfm_cpu_lease_v1','sample':args.sample,'branch':args.branch,'stage':args.stage,
               'resource_mode':'unreserved_cpu_comparison_control' if env.get('MFM_CPU_COMPARISON_CONTROL')=='1' else 'unreserved_gpu_comparison','granted_cpus':cpus,'reserved_cpus':0,
               'reserved_memory_bytes':0,'cpu_seconds':max(0.,after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime),
               'started_at':start,'finished_at':finished,'wall_seconds':finished-start,'queue_wait_seconds':0,
               'peak_tree_rss_bytes':peak,'rss_limit_bytes':limit,'rss_limit_scope':'preprocessing_pool_total',
               'rss_trigger_bytes':limit*9//10,'tool_host_memory_bytes':int(env.get('MFM_GPU_HOST_MEMORY_BYTES',str(384*GIB))),
               'exit_code':process.returncode,'command':effective})
