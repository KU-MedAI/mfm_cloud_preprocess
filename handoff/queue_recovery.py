"""내장 CLI의 명시적 단계 재개 요청과 분리 실행 중복 방지를 관리함."""
from pathlib import Path
import fcntl
import json
import time


def stages(original=True):
    return {'01','02','06','10','12','23'} if original else {'01','02','06','08','10','23'}


def request_recovery(config, sample, request):
    env = config.setdefault('sample_pipeline_env', {}).setdefault(sample, {})
    env['WGS_START_AT'] = request['stage']
    if request.get('overwrite_prepared'):
        env['STAGE04_OVERWRITE'] = '1'
    if request.get('rebuild_haplotag_partial'):
        env['HAPLOTAG_REBUILD_PARTIAL'] = '1'


def recovery_stage_running(request):
    stage = request.get('wait_for_stage')
    if not stage:
        return False
    pid = int(stage['pid'])
    proc = Path('/proc') / str(pid)
    try:
        stat = (proc/'stat').read_text()
        fields = stat[stat.rfind(')')+2:].split()
        if int(fields[19]) != int(stage['started']):
            raise RuntimeError('분리 단계 PID가 재사용됨. 복구 기록 확인 필요함.')
        args = [x.decode() for x in (proc/'cmdline').read_bytes().split(b'\0') if x]
        if args != stage['argv'] and fields[0] not in {'Z','X'}:
            raise RuntimeError('분리 단계 명령이 복구 기록과 다름')
        if fields[0] not in {'Z','X'}:
            return True
    except FileNotFoundError:
        pass
    from integrate_wgs_rna.pool_leases import group_alive
    return group_alive(pid)


def enqueue(config, rows, samples, stage, overwrite, rebuild, save, original=True):
    if not samples or len(set(samples)) != len(samples):
        raise ValueError('--sample ID를 중복 없이 지정해야 함')
    stage = str(stage).zfill(2)
    if stage not in stages(original):
        raise ValueError('검토된 재개 단계만 허용함: '+','.join(sorted(stages(original))))
    known = {r['sample_id'] for r in rows}
    if set(samples) - known:
        raise ValueError('선택 범위 밖의 재개 요청임')
    state = Path(config['state_root'])
    effective = state/'effective_config.json'
    if not effective.is_file():
        raise ValueError('기존 실행 장부가 없음. 먼저 run으로 실행해야 함')
    from .local_queue import locked
    from .state_bindings import StateBindings
    from .planner import read_tsv
    running_config = json.loads(effective.read_text())
    selected = read_tsv(state/'selected_samples.tsv')
    if set(samples) - {r['sample_id'] for r in selected}:
        raise ValueError('기존 큐의 선택 범위 밖임')
    from .archive_receipt import matches as archived_match
    for row in selected:
        if row['sample_id'] in samples:
            root=Path(running_config['output_root'])/row['data_group']/row['project_id']/row['sample_id']
            if archived_match(running_config,row,root) is not None:
                raise ValueError('NAS 이관 완료 샘플은 재개하지 않음. 재처리는 새 output/state를 사용하세요.')
    bindings = StateBindings(running_config, selected)
    if any(locked(bindings.state(s), s) for s in samples):
        raise ValueError('실행 중인 샘플은 재개 요청할 수 없음. 완료/실패 후 요청해야 함')
    queue = state/'queue'
    path = queue/'작업재개요청.json'
    with (queue/'recovery.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        value = json.loads(path.read_text()) if path.exists() else {'samples': {}}
        stamp = time.time_ns()/1e9
        for sample in samples:
            value.setdefault('samples', {})[sample] = dict(
                stage=stage, requested_at=stamp, overwrite_prepared=bool(overwrite),
                rebuild_haplotag_partial=bool(rebuild), retry_download=stage == '01',
                reason='run_pipeline.sh resume 명령으로 요청함')
        save(path, value)
    return {'status':'recovery_requested','samples':samples,'from_stage':stage,
            'queue_running_action':'실행 중인 내장 큐가 다음 확인 주기에 적용함. 중지 상태이면 run으로 재개함.'}
