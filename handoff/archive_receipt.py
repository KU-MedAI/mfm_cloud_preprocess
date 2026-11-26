"""NAS 검증 영수증으로 로컬 자료 회수 후의 완료 상태와 실제 저장 위치를 연결함."""
from pathlib import Path
import hashlib
import json

SCHEMA = 'mfm_nas_archive_v1'
RECEIPT = 'metadata/nas_archive.json'
RETAINED = {'qc/final_package_validation.json', 'metadata/cleanup_completed.json',
            'metadata/window_policy.json'}


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            value.update(block)
    return value.hexdigest()


def read(root):
    root = Path(root).resolve()
    path = root / RECEIPT
    if path.is_symlink():raise ValueError('NAS 영수증 symlink 거부')
    if not path.exists():
        return None
    value = json.loads(path.read_text())
    if (value.get('schema') != SCHEMA or value.get('status') not in {'nas_verified', 'archived'}
            or value.get('sample') != root.name or value.get('source_root') != str(root)
            or not Path(value.get('destination', '')).is_absolute()
            or value.get('verification') != 'sha256_all_files'
            or len(value.get('inventory_sha256', '')) != 64
            or set(value.get('retained_hashes', {})) != RETAINED):
        raise ValueError('NAS 이관 영수증 계약 불일치: ' + str(path))
    batch_path = Path(value['batch_file'])
    if batch_path.name != 'batch.json' or batch_path.parent.name != value['batch_id']:
        raise ValueError('NAS 이관 배치 경로 불일치')
    batch = json.loads(batch_path.read_text())
    item = next((row for row in batch['samples'] if row['sample'] == root.name), None)
    if (item is None or batch.get('schema') != 'mfm_nas_batch_v1'
            or any(value.get(key) != item.get(key) for key in
                   ['source_root', 'destination', 'inventory_sha256', 'retained_hashes'])
            or any(value.get(key) != batch.get(key) for key in ['run_id', 'selection_fingerprint'])):
        raise ValueError('NAS 이관 영수증과 원본 배치가 다름')
    for name, expected in value['retained_hashes'].items():
        if (root/name).is_symlink() or not (root/name).resolve().is_relative_to(root) or digest(root / name) != expected:
            raise ValueError('NAS 이관 뒤 보존 QC/정책 변경: ' + name)
    return value


def matches(config, row, root):
    value = read(root)
    if value is None:
        return None
    from .sequence_mode import completion_mode, payload_matches
    from .resume import sample_environment
    if (value['run_id'] != config['run_id']
            or value['selection_fingerprint'] != config['selection_fingerprint']
            or value['mode'] != completion_mode(config, row)
            or not payload_matches(config, row, root)):
        raise ValueError('다른 실행/출력 정책의 NAS 완료 기록임: ' + row['sample_id'])
    expected = sample_environment(config, row['sample_id']).get('WINDOW_STORAGE_POLICY')
    if expected is not None:
        policy = json.loads((Path(root) / 'metadata/window_policy.json').read_text())
        if policy.get('storage_policy') != expected:
            raise ValueError('NAS 이관 window 정책 불일치')
    return value


def catalog_row(root):
    value = read(root)
    if value is None:
        return None
    row = dict(value['catalog_snapshot'])
    prefix = value['source_root'] + '/'
    for name, item in row.items():
        if isinstance(item, str) and item.startswith(prefix):
            row[name] = value['destination'] + '/' + item[len(prefix):]
    row.update(data_location='nas', archived_root=value['destination'], processing_status='archived',
               archive_training_eligible=row.get('training_eligible', 0),
               qc_all_paths_exist='False', training_eligible='False', missing='archived_to_nas')
    return row
