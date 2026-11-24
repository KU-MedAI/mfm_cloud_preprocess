"""샘플별 재개 단계와 RNA 대체 경로를 적용하고 필수 blacklist를 먼저 검사한다."""
from pathlib import Path
from .downloads import safe_path

def sample_environment(config, sample):
    environment = dict(config.get('pipeline_env', {}))
    overrides = config.get('sample_pipeline_env', {}).get(sample, {})
    if not isinstance(overrides, dict):
        raise ValueError('sample_pipeline_env의 샘플별 값은 객체여야 합니다.')
    environment.update(overrides)
    return environment

def rna_relative(config, row):
    relative = config.get('sample_rna_paths', {}).get(row['sample_id'], row['rna_path'])
    if relative:
        safe_path(config['data_root'], relative)
    if not row['rna_path'] and relative:
        raise ValueError('WGS-only 샘플에 RNA 경로를 추가할 수 없습니다.')
    if row['rna_path'] and not relative:
        raise ValueError('pair 샘플의 RNA 입력을 제거할 수 없습니다.')
    return relative

def alignment_relative(config, row):
    relative = config.get('sample_wgs_paths', {}).get(row['sample_id'], row.get('alignment_path',''))
    if relative:
        safe_path(config['data_root'], relative)
    if not row.get('alignment_path') and relative:
        raise ValueError('read_files 입력 샘플에 alignment 대체 경로를 추가할 수 없습니다.')
    if row.get('alignment_path') and not relative:
        raise ValueError('WGS alignment 입력을 제거할 수 없습니다.')
    return relative

def sample_downloads(config, row, resources):
    from .sequence_mode import needs_wgs,defer_rna
    if not needs_wgs(config,row):
        resources=[r for r in resources if r['kind'] not in {'wgs_alignment','wgs_index','wgs_read'}]
    if defer_rna(config,row):
        resources=[r for r in resources if r['kind'] not in {'rna_alignment','rna_index'}]
    relative = rna_relative(config, row)
    wgs = alignment_relative(config,row)
    if relative == row['rna_path'] and wgs == row.get('alignment_path',''):
        return resources
    result = [dict(resource) for resource in resources]
    for resource in result:
        if resource['kind'] == 'wgs_alignment' and wgs:
            resource['path'] = wgs
        elif resource['kind'] == 'wgs_index' and wgs:
            resource['path'] = wgs + ('.crai' if wgs.endswith('.cram') else '.bai')
        elif resource['kind'] == 'rna_alignment':
            resource['path'] = relative
        elif resource['kind'] == 'rna_index':
            suffix = '.crai' if relative.endswith('.cram') else '.bai'
            resource['path'] = relative + suffix
    # URL, sample identity 및 크기/checksum은 동결 catalog를 그대로 사용한다.
    return result

def require_blacklist(config, root, rows):
    for row in rows:
        env = sample_environment(config, row['sample_id'])
        # shell config와 같은 기본 경로를 검사하여 대용량 작업 전에 누락을 검출한다.
        value = env.get('BLACKLIST_BED') or str(Path(root) / 'references/hg38_autosomes_blacklist_union.bed')
        path = Path(value)
        if not path.is_absolute():
            path = Path(root) / path
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError('필수 blacklist BED 없음: ' + str(path))
