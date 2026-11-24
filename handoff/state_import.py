"""과거 wave별 장부를 내장 큐의 데이터 연결 정보로 이관함."""
from pathlib import Path
import fcntl
import hashlib
import json
import time

from . import download_only, planner


def import_states(config, rows, source, api):
    """검증/중복 제어를 마친 뒤 데이터 장부만 기록함."""
    state = Path(config['state_root'])
    queue = state / 'queue'
    queue.mkdir(parents=True, exist_ok=True)
    with (queue / 'control.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return _import(config, rows, source, api)


def _import(config, rows, source, api):
    source = Path(source).resolve()
    controller = source / '연속제어상태.json'
    if controller.is_file():
        old = json.loads(controller.read_text())
        pid = int(old.get('pid', 0))
        if pid and (Path('/proc') / str(pid)).exists() and old.get('status') != 'draining':
            raise ValueError('이전 제어 큐를 신규 투입 중지(draining)한 뒤 이관해야 합니다.')
    known = {r['sample_id']: r for r in rows}
    items = {}
    selections = list(source.glob('묶음*.tsv')) + list(source.glob('long*_*.tsv'))
    for path in selections:
        cfg_path = path.with_suffix('.json')
        if not cfg_path.is_file():
            continue
        c = json.loads(cfg_path.read_text())
        selected = planner.read_tsv(path)
        for index, row in enumerate(selected):
            sample = row['sample_id']
            if sample not in known:
                continue
            if sample in items:
                raise ValueError('이전 선택 장부에 샘플 중복: ' + sample)
            if any(row.get(k) != known[sample].get(k) for k in ['project_id', 'platform', 'data_group']):
                raise ValueError('이전 샘플의 데이터 계약이 다름: ' + sample)
            for key in ['data_root', 'output_root']:
                if Path(c[key]).resolve() != Path(config[key]).resolve():
                    raise ValueError('입력/출력 경로가 다른 실행을 이관할 수 없음: ' + sample)
            items[sample] = {'config': str(cfg_path), 'selection': str(path),
                'selection_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                'task_index': index, 'state_root': c['state_root'],
                'download_state_root': download_only.configure(c)['state_root']}
    if set(items) != set(known):
        raise ValueError('이전 장부에서 찾지 못한 샘플: ' + ','.join(sorted(set(known) - set(items))))
    state = Path(config['state_root'])
    target = state / 'imported_samples.json'
    value = {'schema': 'mfm_imported_sample_states_v1', 'imported_at': time.time(),
             'source': str(source), 'samples': items}
    if target.exists() and json.loads(target.read_text()).get('samples') != items:
        raise ValueError('기존 이관 장부가 다름. 덮어쓰지 않습니다.')
    queue = state / 'queue'
    pending = {}
    for name in ['재시도상태.json', '작업재개요청.json']:
        path = source / name
        if not path.is_file():
            continue
        old = json.loads(path.read_text())
        if name == '작업재개요청.json':
            old['samples'] = {s:r for s,r in old.get('samples', {}).items() if s in known}
        else:
            for key in ['downloads', 'workers', 'applied_recoveries', 'ready_since', 'purged_inputs']:
                if isinstance(old.get(key), dict):
                    old[key] = {s:r for s,r in old[key].items() if s in known}
        if (queue / name).exists():
            if json.loads((queue / name).read_text()) != old:
                raise ValueError('내장 큐 장부가 이미 존재하며 내용이 다름: ' + str(queue / name))
        pending[name] = old
    # 모든 대상과 기존 장부를 확인한 다음 기록함. 모델 파일은 이동하지 않음.
    for name, old in pending.items():
        api.atomic_json(queue / name, old)
    api.atomic_json(target, value)
    from .state_bindings import StateBindings
    bindings = StateBindings(config, rows)
    bindings.mirror(api)
    result = {'status': 'imported', 'samples': len(items), 'state_root': str(state),
              'model_data_moved': False, 'outside_python_imported': False}
    api.atomic_json(state / 'import_report.json', result)
    return result
