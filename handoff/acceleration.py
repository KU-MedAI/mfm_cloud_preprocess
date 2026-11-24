"""단일 실행 진입점의 CPU/GPU 선택을 명시적인 환경 계약으로 전달함."""
from pathlib import Path
import copy,os

FIELDS={'io_backend':('MFM_IO_BACKEND',{'cpu','nvcomp'}),
        'wgs_backend':('MFM_WGS_BACKEND',{'cpu','parabricks'}),
        'rna_backend':('MFM_RNA_BACKEND',{'cpu','parabricks'}),
        'bam_backend':('MFM_BAM_BACKEND',{'cpu','parabricks'}),
        'coverage_backend':('MFM_COVERAGE_BACKEND',{'cpu','cupy'}),
        'qc_backend':('MFM_QC_BACKEND',{'cpu','cupy'}),
        'zarr_backend':('MFM_ZARR_BACKEND',{'cpu','cupy'})}
PARAMETERS={'gpu_devices':'MFM_GPU_DEVICES','gpu_runtime':'MFM_GPU_RUNTIME',
            'parabricks_image':'MFM_PARABRICKS_IMAGE','gpu_pool_dir':'MFM_GPU_POOL_DIR',
            'gpu_rna_star':'MFM_GPU_RNA_STAR','rna_gpu_policy':'MFM_GPU_RNA_POLICY','gpu_wait_seconds':'MFM_GPU_WAIT_SECONDS'}

def apply(config,overrides=None):
    c=copy.deepcopy(config);choices=dict(c.get('acceleration',{}));env=c.setdefault('pipeline_env',{})
    for k,v in (overrides or {}).items():
        if v is not None:choices[k]=v
    mode=choices.get('mode')
    if mode is not None:
        if mode not in {'cpu','gpu'}:raise ValueError('acceleration.mode는 cpu/gpu여야 함')
        preset={key:'cpu' for key in FIELDS}
        if mode=='gpu':preset.update(wgs_backend='parabricks',rna_backend='parabricks',bam_backend='parabricks',coverage_backend='cupy',qc_backend='cupy',zarr_backend='cupy')
        choices={**preset,**choices}
        if mode=='gpu' and 'MFM_GPU_WAIT_SECONDS' not in env:choices.setdefault('gpu_wait_seconds',86400)
    for key,(variable,allowed) in FIELDS.items():
        # 일반 CPU 실행에서는 기존 config에 새 키를 주입하지 않음. 운영 재시작 호환 보존함.
        default='parabricks' if key=='bam_backend' and any(env.get(k,os.environ.get(k))=='parabricks' for k in ['MFM_WGS_BACKEND','MFM_RNA_BACKEND']) else 'cpu'
        value=choices.get(key,env.get(variable,os.environ.get(variable,default)))
        if value not in allowed:raise ValueError(key+' 허용값: '+','.join(sorted(allowed)))
        if key in choices or variable in env or variable in os.environ or value!='cpu':env[variable]=value
    for key,variable in PARAMETERS.items():
        if key in choices:env[variable]=str(choices[key])
    if any(env.get(variable,'cpu')!='cpu' for variable,_ in FIELDS.values()):
        c.setdefault('slurm',{}).setdefault('gpus_per_task',1)
    if choices:c['acceleration']=choices
    return c

def cli_arguments(parser):
    parser.add_argument('--acceleration-mode',choices=['cpu','gpu'])
    for name,(_,allowed) in FIELDS.items():parser.add_argument('--'+name.replace('_','-'),choices=sorted(allowed))
    parser.add_argument('--gpu-devices',help='auto 또는 GPU 번호/UUID를 쉼표로 구분함')
    parser.add_argument('--gpu-runtime',choices=['auto','native','docker','apptainer'])
    parser.add_argument('--parabricks-image')
    parser.add_argument('--gpu-pool-dir')
    parser.add_argument('--rna-gpu-policy',choices=['strict','experimental'])
    parser.add_argument('--gpu-rna-star',help='GPU STAR와 호환되는 CPU STAR2.7.2a 실행 파일')
    parser.add_argument('--gpu-wait-seconds',type=int)

def overrides(args):return {**{key:getattr(args,key,None) for key in [*FIELDS,*PARAMETERS]},'mode':getattr(args,'acceleration_mode',None)}

def plan(config):
    c=apply(config);env=c.get('pipeline_env',{})
    return {'backends':{key:env.get(var,'cpu') for key,(var,_) in FIELDS.items()},
            'gpu_devices':env.get('MFM_GPU_DEVICES','auto'),
            'gpu_runtime':env.get('MFM_GPU_RUNTIME','auto'),
            'rna_gpu_policy':env.get('MFM_GPU_RNA_POLICY','strict'),
            'rna_gpu_compatibility':'GPU STAR2.7.2a; 적용 인자/현재 CPU profile 차이 기록함',
            'zarr_storage':'4bit ACGTNP/depth uint32/CPU Blosc-Zstd5 유지; GPU 패킹/QC 선택 가능함'}

def preflight(config,rows=None):
    c=apply(config);env=c.get('pipeline_env',{})
    if rows is not None and all(row.get('data_group')=='wgs' for row in rows):env['MFM_RNA_BACKEND']='cpu'
    if not any(env.get(variable,'cpu')!='cpu' for variable,_ in FIELDS.values()):return
    # 기존 환경을 잠시 바꾸는 작업은 단일 제어 큐의 시작 전 검증에서만 수행함.
    previous={key:os.environ.get(key) for key in env}
    try:
        os.environ.update({key:str(value) for key,value in env.items()})
        from integrate_wgs_rna.gpu_pool import inventory,allowed_devices
        if not allowed_devices(inventory()):raise RuntimeError('선택/할당 범위에 GPU가 없음')
        if env.get('MFM_IO_BACKEND')=='nvcomp':
            from integrate_wgs_rna.gpu_nvcomp import library_path
            library_path()
        if any(env.get(k)=='cupy' for k in ['MFM_COVERAGE_BACKEND','MFM_QC_BACKEND','MFM_ZARR_BACKEND']):
            import importlib.util
            if importlib.util.find_spec('cupy') is None:raise RuntimeError('00_environment.sh install-gpu로 CuPy 설치 필요함')
        if env.get('MFM_RNA_BACKEND')=='parabricks':
            from integrate_wgs_rna.gpu_rna import validate_policy
            validate_policy()
        if any(env.get(k)=='parabricks' for k in ['MFM_WGS_BACKEND','MFM_RNA_BACKEND','MFM_BAM_BACKEND']):
            from integrate_wgs_rna.parabricks_runtime import Runtime
            Runtime().validate()
        if env.get('MFM_RNA_BACKEND')=='parabricks':
            from integrate_wgs_rna.gpu_rna import index_star
            index_star()
    finally:
        for key,value in previous.items():
            if value is None:os.environ.pop(key,None)
            else:os.environ[key]=value
