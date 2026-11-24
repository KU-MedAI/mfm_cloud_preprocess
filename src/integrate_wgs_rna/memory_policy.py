"""전처리 전체 RAM 예약, 프로세스 주소 공간 제한과 실사용 보호를 관리한다."""
from __future__ import annotations

import json
import fcntl
import os
import resource
import signal
import time
from pathlib import Path
from .pool_leases import lease_alive, descendants

GIB = 1024 ** 3
DEFAULT_LIMIT = 1_800_000_000_000


def allocation_memory_limit():
    """물리 메모리와 현재 Slurm/cgroup 할당 중 작은 한도를 반환함."""
    limits = []
    try:
        limits.append(os.sysconf('SC_PHYS_PAGES') * os.sysconf('SC_PAGE_SIZE'))
    except (OSError, ValueError):
        pass
    if os.environ.get('SLURM_MEM_PER_NODE'):
        limits.append(int(os.environ['SLURM_MEM_PER_NODE']) * 1024**2)
    elif os.environ.get('SLURM_MEM_PER_CPU'):
        cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', os.environ.get('MFM_CPU_BUDGET', '1')))
        limits.append(int(os.environ['SLURM_MEM_PER_CPU']) * 1024**2 * cpus)
    try:
        relative = next(line.split('::', 1)[1].strip().lstrip('/')
                        for line in Path('/proc/self/cgroup').read_text().splitlines() if '::' in line)
        mount = Path('/sys/fs/cgroup')
        leaf = mount / relative
        for parent in [leaf, *leaf.parents]:
            if parent != mount and mount not in parent.parents:
                break
            path = parent / 'memory.max'
            raw = path.read_text().strip() if path.exists() else 'max'
            if raw != 'max':
                limits.append(int(raw))
    except (OSError, ValueError, StopIteration):
        pass
    return min(limits) if limits else None


def memory_limit(root: Path | None = None) -> int:
    value = min(DEFAULT_LIMIT, int(os.environ.get('MFM_RAM_LIMIT_BYTES', DEFAULT_LIMIT)))
    if root is not None:
        try:
            policy = json.loads((root / 'policy.json').read_text())
            configured = int(policy.get('memory_limit_bytes', value))
            value = min(DEFAULT_LIMIT, configured) if policy.get('live_limits') is True else min(value, configured)
        except FileNotFoundError:
            pass
    allocated = allocation_memory_limit()
    if allocated is not None:
        value = min(value, allocated)
    if value < 1:
        raise ValueError('MFM_RAM_LIMIT_BYTES는 양의 정수여야 합니다.')
    return value


def usable_memory(root: Path) -> int:
    # Python 부모·비계산 프로세스와 예측 오차를 위해 전체의 10%를 남긴다.
    return memory_limit(root) * 9 // 10


def available_memory(root, leases, *, ceiling_bytes=None, observed=None, rss_bytes=None):
    """기존 상한과 확보한 관측으로 RAM 잔여를 계산함. 새 프로세스 실측은 하지 않음."""
    budget = usable_memory(root)
    if ceiling_bytes is not None:
        budget = min(budget, max(0, int(ceiling_bytes)))
    capacity = max(0, budget - sum(max(0, int(row.get('memory_bytes', 0))) for row in leases))
    if observed and observed.get('memory_max'):
        pressure = observed.get('memory_working_bytes')
        if pressure is None:
            pressure = observed.get('memory_current')
        if pressure is None:
            return 0
        capacity = min(capacity, max(0, int(observed['memory_max']) * 9 // 10 - int(pressure)))
    if rss_bytes is not None:
        capacity = min(capacity, max(0, budget - max(0, int(rss_bytes))))
    return capacity


def memory_active(root: Path | None) -> bool:
    if 'MFM_RAM_LIMIT_BYTES' in os.environ:
        return True
    if root is None:
        return False
    try:
        return 'memory_limit_bytes' in json.loads((root / 'policy.json').read_text())
    except FileNotFoundError:
        return False


def process_memory_snapshot():
    """현재 프로세스의 RSS/PSS와 주소공간을 구분한다. 계측 실패는 계산을 막지 않음."""
    result = {'pid': os.getpid(), 'time': time.time(), 'scope': '현재 worker만; 자식 제외',
              'rss_bytes': None, 'pss_bytes': None, 'vms_bytes': None,
              'process_lifetime_peak_rss_bytes': None, 'process_lifetime_peak_vms_bytes': None}
    names = {'VmRSS': 'rss_bytes', 'VmSize': 'vms_bytes',
             'VmHWM': 'process_lifetime_peak_rss_bytes', 'VmPeak': 'process_lifetime_peak_vms_bytes'}
    try:
        for line in Path('/proc/self/status').read_text().splitlines():
            key, _, value = line.partition(':')
            if key in names:
                result[names[key]] = int(value.split()[0]) * 1024
    except (OSError, ValueError, IndexError):
        pass
    try:
        for line in Path('/proc/self/smaps_rollup').read_text().splitlines():
            if line.startswith('Pss:'):
                result['pss_bytes'] = int(line.split()[1]) * 1024
                break
    except (OSError, ValueError, IndexError):
        pass
    return result


def limit_worker_address_space(limit: int) -> None:
    if limit:
        _, previous = resource.getrlimit(resource.RLIMIT_AS)
        target = limit if previous == resource.RLIM_INFINITY else min(previous, limit)
        resource.setrlimit(resource.RLIMIT_AS, (target, previous))


def process_rss(roots: dict[str, int], *, unique_roots=False) -> dict[str, int]:
    processes = {}
    # PID 목록을 먼저 고정한다. Path.glob의 지연 순회 중 PID가 종료되면
    # try 블록에 들어가기 전 ProcessLookupError가 발생할 수 있다.
    for name in os.listdir('/proc'):
        if not name.isdigit():
            continue
        try:
            text = (Path('/proc') / name / 'stat').read_text()
            fields = text[text.rfind(')') + 2:].split()
            processes[int(name)] = (int(fields[1]), int(fields[21]) * os.sysconf('SC_PAGE_SIZE'))
        except (OSError, ValueError, IndexError):
            pass
    children = {}
    for pid, (parent, _) in processes.items():
        children.setdefault(parent, []).append(pid)
    trees = {key: descendants(children, [root]) for key, root in roots.items()}
    usage = {}
    seen = set()
    for key, members in sorted(trees.items(), key=lambda item: -len(item[1])):
        counted = members - seen if unique_roots else members
        usage[key] = sum(processes[pid][1] for pid in counted if pid in processes)
        seen.update(members)
    return usage


def register_memory_observer(root: Path, row: dict) -> Path:
    """예약·입장을 변경하지 않고 승인된 GPU 비교를 전체 RSS 감시에 포함함."""
    path = root / 'memory_observers' / f"{row['manager_pid']}_{row['child_pid']}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    with (root / 'pool.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        candidate = path.with_suffix(f'.{os.getpid()}.tmp')
        candidate.write_text(json.dumps(dict(row, observation_only=True)) + '\n')
        os.replace(candidate, path)
    return path


def protect_pool_memory(root: Path, *, force=False) -> dict:
    """pool.lock 아래 호출한다. 다른 전처리·사용자의 프로세스는 건드리지 않는다."""
    guard = root / 'memory_guard.json'
    try:
        previous = json.loads(guard.read_text())
    except FileNotFoundError:
        previous = {}
    if not force and time.time() - previous.get('time', 0) < 2:
        return previous
    leases = {}
    for folder in ['leases', 'memory_observers']:
        for path in (root / folder).glob('*.json'):
            try:
                row = json.loads(path.read_text())
                if folder == 'memory_observers' and not lease_alive(row):
                    path.unlink(missing_ok=True)
                    continue
                if row.get('child_pid'):
                    key = path.name if folder == 'leases' else folder + '/' + path.name
                    leases[key] = row
            except FileNotFoundError:
                pass
    usage = process_rss({key: row['child_pid'] for key, row in leases.items()}, unique_roots=True)
    total = sum(usage.values())
    limit = memory_limit(root)
    snapshot = {'time': time.time(), 'rss_bytes': total, 'limit_bytes': limit,
                'trigger_bytes': limit * 9 // 10, 'rss_scope': 'preprocessing_including_gpu_comparisons',
                'usage_by_lease': usage, 'action': 'continue'}
    if total >= limit * 9 // 10 and usage and time.time() - previous.get('last_stop', 0) >= 5:
        key = max(usage, key=usage.get)
        row = leases[key]
        snapshot.update(action='stop_largest_stage', sample=row['sample'], stage=row['stage'],
                        last_stop=time.time(), stopped_pid=row['child_pid'])
        try:
            os.killpg(row['child_pid'], signal.SIGTERM)
        except ProcessLookupError:
            pass
    elif previous.get('last_stop'):
        snapshot['last_stop'] = previous['last_stop']
    candidate = guard.with_name(guard.name + f'.{os.getpid()}.tmp')
    candidate.write_text(json.dumps(snapshot) + '\n')
    os.replace(candidate, guard)
    return snapshot
