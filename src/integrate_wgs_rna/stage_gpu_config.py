"""같은 실행의 effective config에서 다음 단계의 GPU 풀 설정만 갱신함."""
from pathlib import Path
import json

BACKENDS=('MFM_IO_BACKEND','MFM_WGS_BACKEND','MFM_RNA_BACKEND','MFM_COVERAGE_BACKEND','MFM_QC_BACKEND','MFM_ZARR_BACKEND','MFM_BAM_BACKEND')
SETTINGS=('MFM_GPU_STREAM_CPUS','MFM_GPU_FASTQ_BACKEND','MFM_GPU_PREPARE_OVERLAP','RNA_FASTQ_BACKEND','RNA_GPU_MARKDUP_MODE','MFM_GPU_DEVICES','MFM_GPU_POOL_DIR','MFM_GPU_HAP_JOBS','WGS_HAP_JOBS','MFM_PARABRICKS_MARKDUP_VALIDATION','MFM_GPU_ALIGNMENT_GPUS','MFM_GPU_SCRATCH_ROOT','MFM_GPU_COMPARE_UNRESERVED','MFM_GPU_HOST_MEMORY_BYTES','MFM_GPU_HIFI_AUTO_PROFILE','MFM_GPU_MINIMAP_ESTIMATED_READ_LENGTH','MFM_GPU_MINIMAP_ESTIMATED_ANCHOR_PER_READ','MFM_GPU_MINIMAP_LARGE_CHAINING_SIZE','MFM_GPU_MINIMAP_VERBOSE','MFM_GPU_INPUT_MODE','MFM_GPU_TAG_RESTORE','MFM_CPU_COMPARISON_CONTROL','MFM_COMPACT_HAP_JOBS','MFM_GPU_CIGAR_ON_GPU','MFM_GPU_BWA_PARTITION','MFM_TRACE_HASH_BACKEND','MFM_GPU_INPUT_CACHE_ROOT','MFM_GPU_BWA_INDEX_ROOT','MFM_GPU_BAM_INDEX_OVERLAP')


def apply_stage_environment(environment):
    """실행 중 계산을 건드리지 않고 새 자식 단계의 장치 목록/동시 상한만 적용함."""
    if not any(environment.get(key,'cpu')!='cpu' for key in BACKENDS):
        return None
    manifest=environment.get('PIPELINE_MANIFEST')
    if not manifest:
        return None
    parents=[p for p in Path(manifest).parents if p.name=='sample_manifests']
    if len(parents)!=1:
        return None
    path=parents[0].parent/'effective_config.json'
    if not path.is_file():
        return None
    config=json.loads(path.read_text())
    if config.get('run_id')!=environment.get('RUN_ID'):
        return None
    for field,variable in [('data_root','INPUT_ROOT'),('output_root','OUTPUT_ROOT')]:
        if not config.get(field) or not environment.get(variable) or Path(config[field]).resolve()!=Path(environment[variable]).resolve():
            return None
    choices=config.get('acceleration',{})
    configured=config.get('pipeline_env',{})
    settings={key:str(configured[key]) for key in SETTINGS if key in configured}
    for field,variable in [('gpu_devices','MFM_GPU_DEVICES'),('gpu_pool_dir','MFM_GPU_POOL_DIR')]:
        if field in choices:
            settings[variable]=str(choices[field])
    io_backend=choices.get('io_backend',configured.get('MFM_IO_BACKEND'))
    if io_backend is not None:
        if io_backend not in {'cpu','nvcomp'}:raise ValueError('io_backend는cpu/nvcomp여야함')
        settings['MFM_IO_BACKEND']=io_backend
    bam_backend=choices.get('bam_backend',configured.get('MFM_BAM_BACKEND'))
    if bam_backend is None and choices.get('mode')=='gpu':bam_backend='parabricks'
    if bam_backend is not None:
        if bam_backend not in {'cpu','parabricks'}:raise ValueError('bam_backend는 cpu/parabricks여야 함')
        settings['MFM_BAM_BACKEND']=bam_backend
    if 'RNA_FASTQ_BACKEND' in settings and settings['RNA_FASTQ_BACKEND'] not in {'auto','samtools','parabricks'}:
        raise ValueError('RNA_FASTQ_BACKEND는auto/samtools/parabricks이어야함')
    if 'RNA_GPU_MARKDUP_MODE' in settings and settings['RNA_GPU_MARKDUP_MODE'] not in {'integrated','separate'}:
        raise ValueError('RNA_GPU_MARKDUP_MODE는integrated/separate이어야함')
    if 'MFM_GPU_HAP_JOBS' in settings and settings['MFM_GPU_HAP_JOBS'] not in {'1','2'}:
        raise ValueError('MFM_GPU_HAP_JOBS는1또는2여야함')
    if 'WGS_HAP_JOBS' in settings:
        if settings['WGS_HAP_JOBS'] not in {'auto','1','2'}:
            raise ValueError('WGS_HAP_JOBS는auto/1/2여야함')
        settings['MFM_WGS_HAP_REQUESTED']=settings['WGS_HAP_JOBS']
    # 로컬 실행의 명시적인 CUDA 설정 갱신만 허용함. Slurm/컨테이너 범위는 넓히지 않음.
    if 'CUDA_VISIBLE_DEVICES' in configured and not environment.get('SLURM_JOB_ID'):
        settings['CUDA_VISIBLE_DEVICES']=str(configured['CUDA_VISIBLE_DEVICES'])
    environment.update(settings)
    return {'source':str(path),'settings':settings,
            'slurm_cuda_visibility_preserved':bool(environment.get('SLURM_JOB_ID')),
            'nvidia_visible_devices':environment.get('NVIDIA_VISIBLE_DEVICES')}
