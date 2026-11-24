"""처리 배치·버퍼 설정을 검증하여 기존 도구의 환경변수에 연결함."""
import re

FIELDS = {
    'minimap_batch_bases': 'MFM_MINIMAP_BATCH_BASES',
    'windows_per_chunk': 'WINDOWS_PER_CHUNK',
    'trace_prefetch_records': 'MFM_TRACE_PREFETCH_RECORDS',
    'trace_prefetch_packed_mib': 'MFM_TRACE_PREFETCH_PACKED_MIB',
    'ase_batch_fragments': 'RNA_ALLELIC_BATCH_SIZE',
    'sort_memory_per_thread': 'WGS_SORT_MEM',
    'star_read_buffer_bytes': 'MFM_STAR_READ_BUFFER_BYTES',
    'star_bam_sort_memory_bytes': 'MFM_STAR_BAM_SORT_MEMORY_BYTES',
}


def positive(value, name):
    if isinstance(value, bool) or not str(value).isdigit() or int(value) < 1:
        raise ValueError(f'{name}은 양의 정수여야 함')
    return int(value)


def apply(config):
    result = dict(config)
    values = config.get('processing_buffers', {})
    if not isinstance(values, dict) or set(values) - set(FIELDS):
        raise ValueError('processing_buffers의 지원하지 않는 항목 또는 형식')
    environment = dict(config.get('pipeline_env', {}))
    for name, raw in values.items():
        if name == 'sort_memory_per_thread':
            if not isinstance(raw, str) or not re.fullmatch(r'[1-9][0-9]*[KMG]?', raw.upper()):
                raise ValueError('sort_memory_per_thread는 512M, 1G 등 양수여야 함')
            value = raw.upper()
            from .realign.align_resources import parse_memory
            if parse_memory(value) < 1024**2:
                raise ValueError('samtools sort 버퍼는 최소1MiB여야 함')
        else:
            value = str(positive(raw, name))
            if name == 'minimap_batch_bases' and int(value) > 500_000_000:
                raise ValueError('메모리 감소용 minimap_batch_bases는 기본500M 이하여야 함')
        key = FIELDS[name]
        if key in environment and str(environment[key]) != value:
            raise ValueError(f'processing_buffers와 pipeline_env의 {key} 설정이 충돌함')
        environment[key] = value
        if name == 'sort_memory_per_thread':
            if 'STAGE04_SORT_MEM' in environment and str(environment['STAGE04_SORT_MEM']) != value:
                raise ValueError('processing_buffers와 STAGE04_SORT_MEM 설정이 충돌함')
            environment['STAGE04_SORT_MEM'] = value
    result['pipeline_env'] = environment
    return result


def minimap_arguments(environment):
    if not environment.get('MFM_MINIMAP_BATCH_BASES'):
        return ''
    value = positive(environment['MFM_MINIMAP_BATCH_BASES'], 'MFM_MINIMAP_BATCH_BASES')
    if value > 500_000_000:
        raise ValueError('메모리 감소용 MFM_MINIMAP_BATCH_BASES는 기본500M 이하여야 함')
    return ' -K ' + str(value)


def star_arguments(environment):
    args = []
    if environment.get('MFM_STAR_READ_BUFFER_BYTES'):
        value = str(positive(environment['MFM_STAR_READ_BUFFER_BYTES'], 'MFM_STAR_READ_BUFFER_BYTES'))
        args.extend(['--limitIObufferSize', value, value])
    if environment.get('MFM_STAR_BAM_SORT_MEMORY_BYTES'):
        value = str(positive(environment['MFM_STAR_BAM_SORT_MEMORY_BYTES'], 'MFM_STAR_BAM_SORT_MEMORY_BYTES'))
        args.extend(['--limitBAMsortRAM', value])
    return args
