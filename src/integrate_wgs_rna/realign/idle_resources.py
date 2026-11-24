"""재정렬의 계산 경계에서만 실측 유휴 CPU와 RAM을 함께 대여한다."""
import json
import math
import os
import time

from ..pool_backfill import cgroup_sample, in_cpu_group
from ..cpu_loan_policy import available_loan_cpus
from ..memory_policy import available_memory


def loan_capacity(root, total, row, others, memory_bytes):
    """pool.lock을 가진 호출자만 사용하며 반환한 한도 안에서 계획을 확정한다."""
    if os.environ.get('MFM_FIXED_SAMPLE_TOKEN'):
        return None
    policy = json.loads((root / 'policy.json').read_text()) if (root / 'policy.json').exists() else {}
    setting = policy.get('cpu_backfill', {})
    if (setting.get('enabled') is not True or setting.get('realign_phases') is not True
            or row.get('branch') != 'WGS' or row.get('stage') not in {'08', '10'}
            or os.environ.get('MFM_CPU_UTILIZATION', 'maximum') != 'maximum'
            or (os.environ.get('MFM_WGS_BACKEND','cpu')!='cpu' and not (os.environ.get('MFM_WGS_BACKEND')=='parabricks' and setting.get('gpu_realign_phases') is True))):
        return None
    try:
        state = json.loads((root / 'watchdog.json').read_text())
        observed = state['cpu_backfill']
        if (time.time() - float(state['updated']) > 3 or time.time() < float(state['updated'])
                or state.get('guard_action') != 'continue' or not observed.get('ready')
                or abs(float(observed['quota']) - total) > .001
                or float(observed['observation_seconds']) < 10):
            return None
        events = root / 'cpu_backfill_events.jsonl'
        if events.exists() and time.time() - events.stat().st_mtime < 10:
            return None
        current = cgroup_sample()
        if (current is None or abs(current['quota'] - total) > .001
                or not in_cpu_group(os.getpid(), current['cpu_cgroup'])):
            return None
        nice = os.getpriority(os.PRIO_PROCESS, 0)
        if nice < max(0, min(19, int(setting.get('nice', 19)))):
            return None
        idle = int(observed['idle_cpus'])
        if setting.get('use_average_idle') is True:
            idle = max(0, math.floor(total - float(observed['used_cpus']) - 1))
        physical = int(idle * max(1, min(2, float(setting.get('grant_multiplier', 1)))))
        # 다른 프로젝트의 CPU 사용량으로 전처리 미사용 예산을 차감하지 않음.
        # quota/동일 cgroup/RAM 보호와 전체 논리 대여 상한은 별도로 유지함.
        from ..pool_cpu_observer import reuse_capacity
        reuse = reuse_capacity(total, setting, observed.get('reservation_reuse'))
        basis = 'reservation_reuse' if reuse > physical else 'cgroup_idle'
        cpus = min(total, max(physical, reuse), available_loan_cpus(total, setting, others))
        ram = min(max(0, int(memory_bytes)),available_memory(root,others,observed=current))
        return cpus, ram, {'schema': 'mfm_idle_cpu_loan_v1', 'cpus': cpus, 'nice': nice,
                          'allocation_basis': basis, 'reservation_reuse_cpus': reuse,
                          'kind': 'gpu_realign_phase' if os.environ.get('MFM_WGS_BACKEND')=='parabricks' else 'realign_phase', 'observed_used_cpus': observed['used_cpus'],
                          'observed_peak_cpus': observed['peak_used_cpus'],
                          'cpu_cgroup': current['cpu_cgroup']}
    except (OSError, ValueError, KeyError):
        return None


def commit_loan(root, row, loan):
    """대기/완료 경계의 최종 예약만 장부에 남긴다. 실행 중 thread는 바꾸지 않는다."""
    if loan is None:
        row.pop('cpu_backfill', None)
        return
    row['cpu_backfill'] = dict(loan, cpus=int(row['cpus']))
    with (root / 'cpu_backfill_events.jsonl').open('a') as handle:
        handle.write(json.dumps(dict(row, phase_loan_at=time.time()), ensure_ascii=False) + '\n')
