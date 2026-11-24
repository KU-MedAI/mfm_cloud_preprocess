"""샘플 전체 CPU 소유권을 유지하고 단계들은 같은 공유 RAM 장부를 사용함."""
from contextlib import contextmanager
from pathlib import Path
import fcntl
import json
import os
import signal
import subprocess
import time
import uuid

TOKEN = 'MFM_FIXED_SAMPLE_TOKEN'


def read(path):
    try:
        return json.loads(Path(path).read_text())
    except FileNotFoundError:
        return {}


def save(path, value):
    temporary = path.with_name(path.name + f'.{os.getpid()}.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False) + '\n')
    temporary.replace(path)


def process_start(pid):
    try:
        raw = Path(f'/proc/{int(pid)}/stat').read_text()
        fields = raw[raw.rfind(')') + 2:].split()
        return int(fields[19]) if fields[0] not in {'Z', 'X'} else None
    except (FileNotFoundError, ProcessLookupError):
        return None


def token():
    return os.environ.get(TOKEN, '')


def owner(root, value=None):
    value = value or token()
    if not value or not value.isalnum():
        raise ValueError('고정 샘플 CPU 소유권 식별자가 없거나 잘못됨')
    path = Path(root) / 'leases' / f'sample_{value}.json'
    row = read(path)
    if not row.get('sample_cpu_owner') or row.get('fixed_sample_token') != value:
        raise ValueError('살아 있는 샘플 CPU 소유권 기록이 없음')
    return path, row


def capacity(root):
    return int(owner(root)[1]['cpus']) if token() else None


def global_claim(row):
    if row.get('fixed_sample_token') and not row.get('sample_cpu_owner'):
        return 0  # 단계 CPU는 부모 샘플 소유권에 이미 포함됨.
    if row.get('cpu_backfill') or str(row.get('realign_phase', '')).endswith(':waiting'):
        return 0
    return max(0, int(row.get('cpus', 0)))


def owner_alive(row):
    started=process_start(row.get('manager_pid', 0))
    if started is not None and started == row.get('manager_start'):
        return True
    from .pool_leases import lease_alive
    for path in (Path(row['pool_dir']) / 'leases').glob('*.json'):
        child = read(path)
        if (not child.get('sample_cpu_owner') and child.get('fixed_sample_token') == row['fixed_sample_token']
                and lease_alive(child)):
            return True
    return False


def available(root, total):
    from .pool_leases import prune_leases, read_pool_records
    prune_leases(Path(root))
    from .cpu_loan_policy import available_regular_cpus
    rows = [row for _, row in read_pool_records(Path(root) / 'leases')]
    return available_regular_cpus(total, read(Path(root) / 'policy.json').get('cpu_backfill', {}), rows)


@contextmanager
def sample_slot(config, sample):
    if config.get('cpu_allocation_mode', 'dynamic') != 'fixed_sample':
        yield {}
        return
    if (os.environ.get('MFM_GPU_COMPARE_UNRESERVED') == '1'
            or str(config.get('pipeline_env', {}).get('MFM_GPU_COMPARE_UNRESERVED', '0')) == '1'):
        raise ValueError('고정 CPU 모드와 CPU/RAM 예약 우회 비교 모드는 함께 사용할 수 없음')
    root = Path(config.get('pipeline_env', {}).get('MFM_CPU_POOL_DIR',
                str(Path(config['state_root']) / 'cpu_pool' / config['run_id']))).resolve()
    cpus = int(config['sample_cpus'])
    if token():
        _, row = owner(root)
        if row['sample'] != sample or int(row['cpus']) != cpus:
            raise ValueError('상속한 샘플 CPU 소유권과 요청 샘플/CPU가 다름')
        yield {TOKEN: token(), 'MFM_CPU_BUDGET': str(cpus), 'MFM_SAMPLE_CPU_BUDGET': str(cpus)}
        return
    from .pool_leases import pool_cpu_limit, prune_leases, read_pool_records
    from .memory_policy import GIB, available_memory
    from .cpu_loan_policy import available_regular_cpus
    root.mkdir(parents=True, exist_ok=True)
    (root / 'leases').mkdir(exist_ok=True)
    (root / 'waiters').mkdir(exist_ok=True)
    requested = (config.get('effective_cpu_budget') or
                 config.get('resource_resolution', {}).get('effective_total_cpu_budget') or
                 config.get('cpu_pool_limit', 'auto'))
    ceiling = len(os.sched_getaffinity(0)) if requested == 'auto' else int(requested)
    total = min(ceiling, pool_cpu_limit(root, ceiling))
    allowed = sorted(os.sched_getaffinity(0))[:total]
    if cpus < 1 or cpus > total:
        raise ValueError('고정 sample_cpus가 실제 CPU 풀 한도를 초과함')
    value = uuid.uuid4().hex
    path = root / 'leases' / f'sample_{value}.json'
    queued = time.monotonic()
    base_memory = 4 * GIB
    waiting_path=root/'leases'/f'samplewait_{value}.json'
    try:
        while True:
            with (root / 'pool.lock').open('a+') as lock:
                fcntl.flock(lock, fcntl.LOCK_EX)
                policy = read(root / 'policy.json')
                if not policy:
                    policy = dict(live_limits=True, cpu_pool_limit=total,
                                  memory_limit_bytes=config.get('memory_limit_bytes', 1_650_000_000_000),
                                  profile='throughput', cpu_backfill={'enabled': False})
                    save(root / 'policy.json', policy)
                elif (policy.get('memory_limit_bytes', config.get('memory_limit_bytes', 1_650_000_000_000))
                      > config.get('memory_limit_bytes', 1_650_000_000_000)):
                    raise ValueError('기존 pool RAM 정책이 요청보다 큼. 새 run_id/state 또는 일치하는 정책이 필요함')
                elif 'memory_limit_bytes' not in policy:
                    policy['memory_limit_bytes'] = config.get('memory_limit_bytes', 1_650_000_000_000)
                    save(root / 'policy.json', policy)
                prune_leases(root)
                total=min(ceiling,pool_cpu_limit(root,ceiling))
                allowed=sorted(os.sched_getaffinity(0))[:total]
                if cpus>total:raise ValueError('대기 중 CPU 풀이 고정 샘플 예산보다 작아짐')
                rows = [row for record_path,row in read_pool_records(root / 'leases') if record_path!=waiting_path]
                pending=read(waiting_path)
                if not pending and available_memory(root,rows)>=512*(1<<20):
                    pending=dict(sample=sample,branch='SAMPLE',stage='sample_wait',cpus=0,
                        memory_bytes=512*(1<<20),manager_pid=os.getpid(),manager_start=process_start(os.getpid()),
                        resource_wait=dict(minimum_cpus=cpus,memory_bytes=base_memory,queued_at=time.time()))
                    save(waiting_path,pending)
                used = {cpu for row in rows if row.get('sample_cpu_owner') for cpu in row.get('cpu_ids', [])}
                ids = [cpu for cpu in allowed if cpu not in used]
                from .resource_fairness import fair_cpu_capacity
                claims=[*rows,pending] if pending else rows
                free=fair_cpu_capacity(root,total,claims,available_regular_cpus(total,policy.get('cpu_backfill',{}),rows),request=pending or None)
                if (free >= cpus and len(ids) >= cpus
                        and available_memory(root, rows) >= base_memory):
                    ids = ids[:cpus]
                    row = dict(sample=sample, branch='SAMPLE', stage='sample', sample_cpu_owner=True,
                               fixed_sample_token=value, cpus=cpus, cpu_ids=ids, total_cpus=total,
                               memory_bytes=base_memory, manager_pid=os.getpid(), manager_start=process_start(os.getpid()),
                               pool_dir=str(root), acquired_at=time.time())
                    waiting_path.unlink(missing_ok=True)
                    save(path, row)
                    break
            time.sleep(.1)
    finally:
        with (root/'pool.lock').open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            waiting_path.unlink(missing_ok=True)
    before = os.sched_getaffinity(0)
    env = {'MFM_CPU_ALLOCATION_MODE': 'fixed_sample', TOKEN: value,
           'MFM_FIXED_SAMPLE_CPUS': str(cpus), 'MFM_FIXED_SAMPLE_ID': sample,
           'MFM_CPU_POOL_DIR': str(root), 'MFM_CPU_BUDGET': str(cpus),
           'MFM_SAMPLE_CPU_BUDGET': str(cpus)}
    print('[SAMPLE-CPU] ' + json.dumps(dict(sample=sample, cpus=cpus, cpu_ids=ids,
          queue_wait_seconds=time.monotonic()-queued, token=value), ensure_ascii=False), flush=True)
    try:
        os.sched_setaffinity(0, ids)
        yield env
    finally:
        os.sched_setaffinity(0, before)
        with (root / 'pool.lock').open('a+') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            from .pool_leases import lease_alive
            live = any(not child.get('sample_cpu_owner') and child.get('fixed_sample_token') == value and lease_alive(child)
                       for _, child in read_pool_records(root / 'leases'))
            if not live:
                path.unlink(missing_ok=True)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--sample-cpus', type=int, required=True)
    parser.add_argument('--total-cpus', type=int, required=True)
    parser.add_argument('--pool-dir', type=Path, required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('실행 명령이 필요함')
    config = dict(cpu_allocation_mode='fixed_sample', sample_cpus=args.sample_cpus,
                  effective_cpu_budget=args.total_cpus, state_root=str(args.pool_dir.parent), run_id='direct',
                  memory_limit_bytes=int(os.environ.get('MFM_RAM_LIMIT_BYTES', '1650000000000')),
                  pipeline_env={'MFM_CPU_POOL_DIR': str(args.pool_dir)})
    with sample_slot(config, args.sample) as values:
        proc = subprocess.Popen(command, env=dict(os.environ, **values), start_new_session=True)
        def forward(sig, _frame):
            try:
                os.killpg(proc.pid, sig)
            except ProcessLookupError:
                pass
        for sig in [signal.SIGTERM, signal.SIGINT]:
            signal.signal(sig, forward)
        return proc.wait()


if __name__ == '__main__':
    raise SystemExit(main())
