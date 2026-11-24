"""종료된 CPU 예약을 정리하되 분리 세션의 생존 계산은 보존함."""
from pathlib import Path
import json
import os
import socket
import time


def cgroup_cpu_limit():
    values = []
    try:
        relative = next(line.split('::', 1)[1].strip().lstrip('/')
                        for line in Path('/proc/self/cgroup').read_text().splitlines() if '::' in line)
        mount = Path('/sys/fs/cgroup')
        leaf = mount / relative
        for parent in [leaf, *leaf.parents]:
            if parent != mount and mount not in parent.parents:
                break
            path = parent / 'cpu.max'
            raw = path.read_text().split() if path.exists() else ['max']
            if raw[0] != 'max':
                values.append(max(1, int(raw[0]) // int(raw[1])))
    except (OSError, ValueError, StopIteration):
        pass
    return min(values) if values else None


def pool_cpu_limit(root, requested):
    """승인된 공유 CPU 정책을 실제 affinity·Slurm 범위 안에서 반영함."""
    value = int(requested)
    try:
        policy = json.loads((Path(root) / 'policy.json').read_text())
        if policy.get('live_limits') is True:
            value = int(policy.get('cpu_pool_limit', value))
    except FileNotFoundError:
        pass
    affinity = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else os.cpu_count() or 1
    value = min(value, affinity)
    if os.environ.get('SLURM_CPUS_PER_TASK'):
        value = min(value, int(os.environ['SLURM_CPUS_PER_TASK']))
    quota = cgroup_cpu_limit()
    if quota is not None:
        value = min(value, quota)
    if value < 1:
        raise ValueError('공유 CPU 상한은 양수여야 함')
    return value


def fitting_waiters(rows, cpus, memory_bytes, *, stage_limits=None, active=(),prepare_read_only=False):
    """최소 CPU와 RAM이 동시에 들어가는 요청만 FIFO 순서로 선택함."""
    selected = []
    token = os.environ.get('MFM_FIXED_SAMPLE_TOKEN', '')
    counts={}
    for row in active:
        if row.get('stage_storage_wait') or row.get('storage_paused') or row.get('prepare_phase')=='read_only_qc':
            continue
        label=str(row.get('branch','')).upper()+':'+str(row.get('stage',''))
        counts[label]=counts.get(label,0)+1
    for path, row in rows:
        if str(row.get('fixed_sample_token', '')) != token:
            continue
        label=str(row.get('branch','')).upper()+':'+str(row.get('stage',''))
        read_entry=prepare_read_only and label=='WGS:01'
        if not read_entry and stage_limits and label in stage_limits and counts.get(label,0)>=int(stage_limits[label]):
            continue
        need_cpu = int(row.get('min_cpus', 1))
        need_memory = int(row.get('memory_bytes', 0))
        if need_cpu <= cpus and need_memory <= memory_bytes:
            selected.append((path, row))
            cpus -= need_cpu
            memory_bytes -= need_memory
            if not read_entry:counts[label]=counts.get(label,0)+1
    return selected


def logical_cpu_claim(row, *, token=None):
    """전역에서는 샘플 소유권만, 고정 샘플 내부에서는 같은 샘플의 단계만 셈."""
    token = os.environ.get('MFM_FIXED_SAMPLE_TOKEN', '') if token is None else token
    sample = str(row.get('fixed_sample_token', ''))
    if token:
        if sample != token or row.get('sample_cpu_owner'):
            return 0
    elif sample and not row.get('sample_cpu_owner'):
        return 0
    return max(0, int(row.get('cpus', 0)))


def native_allocator_ready(root):
    """동일 호스트의 살아 있는 내장 감시기에 정상 입장을 맡김. 고정 하위 풀은 제외함."""
    if os.environ.get('MFM_FIXED_SAMPLE_TOKEN'):
        return False
    try:
        state=json.loads((Path(root)/'watchdog.json').read_text())
        owner=state['watchdog_owner']
        if (state.get('schema')!='mfm_resource_watchdog_v1' or owner['host']!=socket.gethostname()
                or not 0 <= time.time()-float(state['updated']) <= 10):
            return False
        raw=Path(f"/proc/{int(owner['pid'])}/stat").read_text()
        fields=raw[raw.rfind(')')+2:].split()
        return fields[0] not in {'Z','X'} and int(fields[19])==int(owner['started'])
    except (OSError,ValueError,KeyError,TypeError,IndexError):
        return False


def regular_cpu_claim(row):
    """대여·재정렬 대기는 정상 예약에서 제외하고 샘플 소유권은 한 번만 셈."""
    if row.get('cpu_backfill') or str(row.get('realign_phase','')).endswith(':waiting'):
        return 0
    return logical_cpu_claim(row)


def descendants(children, roots):
    """부모별 자식 색인으로 각 PID를 한 번만 방문함. 중복/순환도 허용함."""
    members = set(roots)
    pending = list(members)
    while pending:
        for child in children.get(pending.pop(), ()):
            if child not in members:
                members.add(child)
                pending.append(child)
    return members


def process_alive(pid):
    if not pid:
        return False
    try:
        value = (Path('/proc') / str(int(pid)) / 'stat').read_text()
        return value[value.rfind(')') + 2:].split()[0] not in {'Z', 'X'}
    except FileNotFoundError:
        return False
    except PermissionError:
        return True


def group_alive(group):
    if not group:
        return False
    for path in Path('/proc').iterdir():
        if not path.name.isdigit():
            continue
        try:
            value = (path / 'stat').read_text()
            fields = value[value.rfind(')') + 2:].split()
            if int(fields[2]) == int(group) and fields[0] not in {'Z', 'X'}:
                return True
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    return False


def lease_alive(row):
    if row.get('sample_cpu_owner'):
        from .sample_cpu import owner_alive
        return owner_alive(row)
    # 외부 수동 예약처럼 소유 PID 정보가 없으면 임의 삭제하지 않음.
    if 'manager_pid' not in row and 'child_pid' not in row:
        return True
    return (process_alive(row.get('manager_pid')) or process_alive(row.get('child_pid'))
            or group_alive(row.get('child_pid')))


def prune_leases(root):
    """호출자가 pool.lock을 보유한 경우에만 사용함."""
    removed = []
    for path in (Path(root) / 'leases').glob('*.json'):
        try:
            row = json.loads(path.read_text())
            if not lease_alive(row):
                path.unlink(missing_ok=True)
                removed.append(path.name)
        except FileNotFoundError:
            continue
    return removed


def read_pool_records(directory):
    """이전 manager의 잠금 없는 반환과 겹쳐 사라진 다른 장부만 건너뜀."""
    for path in Path(directory).glob('*.json'):
        try:
            yield path, json.loads(path.read_text())
        except FileNotFoundError:
            continue


def live_waiters(directory):
    """종료된 요청은 우선순위에서 제외하고 PID 없는 수동 요청은 보존함."""
    for path,row in read_pool_records(directory):
        try:
            if 'pid' in row and not process_alive(row['pid']):
                continue
        except (TypeError,ValueError):
            continue
        yield path,row
