"""샘플 시작/다운로드 자식에 GPU IO의 명시 설정만 전달함."""
import os


def apply(config):
    policy=config.get('pipeline_env',{})
    value=policy.get('MFM_IO_BACKEND',config.get('acceleration',{}).get('io_backend','cpu'))
    if value not in {'cpu','nvcomp'}:raise ValueError('io_backend는cpu/nvcomp여야함')
    if value=='cpu':
        if 'MFM_IO_BACKEND' in policy or 'io_backend' in config.get('acceleration',{}):os.environ['MFM_IO_BACKEND']='cpu'
        return
    os.environ['MFM_IO_BACKEND']='nvcomp'
    for field in ['CUDA_VISIBLE_DEVICES','MFM_GPU_DEVICES','MFM_GPU_POOL_DIR','MFM_GPU_WAIT_SECONDS','MFM_GPU_CRAM_MIN_BYTES','MFM_CPU_POOL_DIR','MFM_CPU_BUDGET','MFM_RAM_LIMIT_BYTES']:
        if field in policy:os.environ[field]=str(policy[field])
