"""단일 샘플 haplotag의 모든 염색체 종료 후 CPU/RAM을 단계별 재예약함."""
from pathlib import Path
import fcntl
import json
import os
import time

from .cpu_executor import WorkerBudget
from .memory_policy import GIB
from .pipeline.checkpoints import atomic_json


def reserve_tail(config, sample, phase, maximum=1, memory_per_cpu=0):
    value = os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if not value or config.jobs != 1:
        return None
    path = Path(value)
    root = path.parent.parent
    with (root/'pool.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        row = json.loads(path.read_text())
        if (row.get('branch'), row.get('stage'), row.get('sample')) != ('WGS', '06', sample):
            return None
        row.setdefault('haplotag_initial_memory_bytes', row['memory_bytes'])
        policy_path = root/'policy.json'
        policy = json.loads(policy_path.read_text()) if policy_path.exists() else {}
        loan = row.get('cpu_backfill')
        if (loan and policy.get('cpu_backfill', {}).get('enabled') is True
                and policy.get('cpu_backfill', {}).get('haplotag_phases') is True):
            loan['phase_reacquire_cpus'] = int(loan['cpus'])
        row['haplotag_phase'] = phase + ':waiting'
        atomic_json(path, row)
    started = time.monotonic()
    budget = WorkerBudget(maximum, maximum, True, memory_per_cpu)
    budget.base_memory = 4*GIB
    count = budget.capacity(maximum)
    while count < 1:
        time.sleep(.1)
        count = budget.capacity(maximum)
    with (root/'pool.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        row = json.loads(path.read_text())
        event = {'phase': phase, 'cpus': count, 'memory_bytes': row['memory_bytes'],
                 'memory_per_cpu': memory_per_cpu, 'queue_wait_seconds': time.monotonic()-started}
        row['haplotag_phase'] = phase
        row.setdefault('haplotag_resource_events', []).append(event)
        atomic_json(path, row)
    print('[HAPLOTAG-PHASE] '+json.dumps(event, ensure_ascii=False), flush=True)
    return count
