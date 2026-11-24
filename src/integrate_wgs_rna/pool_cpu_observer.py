"""전처리 제어기 트리의 실제 CPU를 관측함. 다른 프로젝트는 포함하지 않음."""
from pathlib import Path
import math
import os
import time

from .pool_backfill import CpuObserver
from .pool_leases import descendants


def tree_cpu_sample(roots, *, proc=Path('/proc')):
    """회수된 자식 CPU까지 합산하고 PID 시작시각으로 관측 범위를 고정함."""
    rows = {}
    for name in os.listdir(proc):
        if not name.isdigit():
            continue
        try:
            raw = (proc / name / 'stat').read_text()
            f = raw[raw.rfind(')') + 2:].split()
            rows[int(name)] = (int(f[1]), int(f[19]), sum(int(f[i]) for i in (11, 12, 13, 14)))
        except (OSError, ValueError, IndexError):
            continue
    selected = {(int(r['pid']), int(r['started'])) for r in roots
                if int(r['pid']) in rows and rows[int(r['pid'])][1] == int(r['started'])}
    if not selected:
        return None
    children = {}
    for pid, row in rows.items():
        children.setdefault(row[0], []).append(pid)
    members = descendants(children, [pid for pid, _ in selected])
    ticks = sum(rows[pid][2] for pid in members if pid in rows)
    return {'scope': sorted(selected), 'usage_usec': ticks * 1_000_000 // os.sysconf('SC_CLK_TCK'),
            'process_count': len(members)}


class PreprocessingCpuObserver:
    def __init__(self):
        self.observer = CpuObserver()
        self.scope = None
        self.last_time = None
        self.state = {'ready': False}

    def observe(self, total, setting, *, now=None, reader=tree_cpu_sample):
        """/proc 순회는 풀 잠금 밖에서 최대 2초에 한 번 수행함."""
        explicit_time = now is not None
        now = time.monotonic() if now is None else now
        started = now
        reuse = setting.get('reservation_reuse', {})
        if setting.get('enabled') is not True or reuse.get('enabled') is not True:
            self.observer.history.clear()
            self.last_time = None
            return {'ready': False, 'enabled': False}
        if self.last_time is not None and now - self.last_time < 2:
            return self.state
        try:
            sample = reader(reuse.get('roots', []))
        except (OSError, ValueError, KeyError, TypeError):
            sample = None
        if not explicit_time:
            now = time.monotonic()
        self.last_time = now
        if sample is None:
            self.observer.history.clear()
            self.state = {'ready': False, 'enabled': True, 'reason': '전처리 관측 PID 확인 필요'}
            return self.state
        if sample['scope'] != self.scope:
            self.observer.history.clear()
            self.scope = sample['scope']
        state = self.observer.observe(total, now=now, reading=dict(sample, quota=total))
        # 비원자적 /proc 순회 중 자식 회수가 겹친 이상 표본은 대여 근거로 사용하지 않음.
        if state.get('ready') and state['peak_used_cpus'] > total + 1:
            self.observer.history.clear()
            state.update(ready=False, reason='전처리 CPU 관측 재안정화 중')
        state.update(enabled=True, cpu_scope='preprocessing_trees', sampled_at=now,
                     collection_seconds=max(0, now-started))
        self.state = state
        return state


def reuse_capacity(total, setting, observed, *, now=None):
    """커널 유휴와 구분한 전처리 예산 미사용분임. RAM 예산은 늘리지 않음."""
    now = time.monotonic() if now is None else now
    if (setting.get('reservation_reuse', {}).get('enabled') is not True
            or not observed or not observed.get('ready')
            or not 0 <= now - observed.get('sampled_at', -math.inf) <= 3):
        return 0
    return max(0, math.floor(total - float(observed['peak_used_cpus']) - 2))
