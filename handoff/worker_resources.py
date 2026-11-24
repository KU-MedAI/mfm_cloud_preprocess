"""살아 있는 큐의 다음 worker만 명시된 풀 CPU 정책으로 전환함."""
from pathlib import Path
import copy,json,os,time
from . import resources


def apply(config, executor='local'):
    if executor != 'local' or os.environ.get('SLURM_JOB_ID'):
        return config
    if config.get('pipeline_env',{}).get('MFM_GPU_COMPARE_UNRESERVED')=='1':
        from integrate_wgs_rna.gpu_compare_resources import enabled,config_environment
        if enabled(config_environment(config)):return config
    pool=config.get('pipeline_env',{}).get('MFM_CPU_POOL_DIR')
    path=Path(pool)/'policy.json' if pool else None
    policy=json.loads(path.read_text()) if path and path.is_file() else {}
    requested=policy.get('future_sample_resources')
    if requested is None:return config
    if not isinstance(requested,dict) or set(requested)!={'cpu_allocation_mode','sample_cpus','local_jobs'}:
        raise ValueError('후속 worker CPU 정책의 필드가 잘못됨')
    if requested['cpu_allocation_mode']!='fixed_sample':
        raise ValueError('후속 worker 정책은 fixed_sample만 허용함')
    for key in ['sample_cpus','local_jobs']:
        if isinstance(requested[key],bool) or not isinstance(requested[key],int) or requested[key]<1:
            raise ValueError('후속 worker CPU/동시수는 양의 정수여야 함')
    updated=copy.deepcopy(config)
    if 'MFM_CPU_ALLOCATION_MODE' in updated.get('pipeline_env',{}):
        updated['pipeline_env']['MFM_CPU_ALLOCATION_MODE']='fixed_sample'
    updated=resources.resolve(updated,dict(requested,
        local_jobs=min(config['local_jobs'],requested['local_jobs'])))
    updated['_worker_resource_policy']={'source':str(path),'applied_at':time.time(),
        'cpu_allocation_mode':updated['cpu_allocation_mode'],'sample_cpus':updated['sample_cpus'],
        'local_jobs':updated['local_jobs'],'cpu_pool_limit':updated['cpu_pool_limit']}
    return updated
