"""완료 경계에서 오래 기다린 최소 CPU를 모음. 실행 중 예약은 회수하지 않음."""
import os
import time
from . import memory_policy
from .pool_leases import regular_cpu_claim, live_waiters

WAIT_SECONDS = 30


def request_key(row):
    return (row.get('manager_pid', row.get('pid')), row.get('ticket'),
            row.get('sample'), row.get('fixed_sample_token', ''))


def priority_request(root, total, leases, *, waiters=None, now=None):
    """CPU/RAM 최소량이 전체 예산에 맞는 오래된 요청 하나만 보호함."""
    now = time.time() if now is None else now
    token = os.environ.get('MFM_FIXED_SAMPLE_TOKEN', '')
    memory = memory_policy.usable_memory(root) - sum(max(0, int(r.get('memory_bytes', 0))) for r in leases)
    candidates = []
    for row in leases:
        if row.get('sample_cpu_owner') or str(row.get('fixed_sample_token', '')) != token:
            continue
        pending = row.get('realign_wait') or row.get('resource_wait')
        if not pending:
            worker = row.get('worker_allocation', {})
            if worker.get('requested_workers', 0) and not worker.get('granted_workers', 0) and not worker.get('running_workers', 0):
                pending = dict(minimum_cpus=worker.get('minimum_cpus', worker.get('cpus_per_worker', 1)),
                               memory_bytes=worker.get('minimum_memory_bytes', row.get('memory_bytes', 0)),
                               queued_at=worker.get('queued_at', row.get('acquired_at')))
        if not pending:
            continue
        candidates.append(('internal', row, pending, int(row.get('memory_bytes', 0))))
    waiters = live_waiters(root/'waiters') if waiters is None else waiters
    for _, row in waiters:
        if str(row.get('fixed_sample_token', '')) == token:
            candidates.append(('stage', row, dict(minimum_cpus=row.get('min_cpus', 1),
                              memory_bytes=row.get('memory_bytes', 0),queued_at=row.get('queued_at')), 0))
    eligible = []
    for kind, row, pending, held_memory in candidates:
        queued = pending.get('queued_at', row.get('acquired_at'))
        if queued is None or not 0 <= float(queued) <= now-WAIT_SECONDS:
            continue
        cpus = int(pending['minimum_cpus'])
        needed_memory = int(pending.get('memory_bytes', 0))
        if not 1 <= cpus <= total or needed_memory > memory+held_memory:
            continue
        extra = max(0, cpus-(regular_cpu_claim(row) if kind=='internal' else 0))
        if extra:
            eligible.append(dict(kind=kind,key=request_key(row),minimum_cpus=cpus,
                additional_cpus=extra,memory_bytes=needed_memory,queued_at=float(queued),
                sample=row.get('sample'),stage=row.get('stage')))
    return min(eligible,key=lambda r:(r['queued_at'],str(r['key']))) if eligible else None


def fair_cpu_capacity(root, total, leases, available, *, request=None, waiters=None):
    """다른 작업의 신규 확대만 제한함. 현재 running worker 하한은 호출자가 보존함."""
    priority = priority_request(root,total,leases,waiters=waiters)
    if priority is None or request is not None and priority['key']==request_key(request):
        return available
    return max(0, available-priority['additional_cpus'])
