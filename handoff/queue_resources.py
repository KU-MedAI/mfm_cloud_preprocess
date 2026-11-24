"""전체 CPU 공유 풀과 RAM 예약을 계산하며 디스크 투입 차단은 하지 않음."""
from pathlib import Path
import json
import os
import time

from .resources import cpu_budget, memory_budget
from .stage_resources import longread_memory_plan, starting_memory
from .lease_identity import lease_matches_config
from .download_budget import memory_plan
import sys
try:
    from integrate_wgs_rna.pool_leases import lease_alive
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
    from integrate_wgs_rna.pool_leases import lease_alive

GIB = 1024 ** 3


from .host_limits import cgroup_limits
from integrate_wgs_rna.cgroup_memory import pressure_usage


def effective_memory(config):
    limit = memory_budget(config)
    _, cgroup, _ = cgroup_limits()
    return min(limit, cgroup) if cgroup else limit


def cpu_limit(config):
    limit = cpu_budget(config, 'local')
    _, _, quota = cgroup_limits()
    return min(limit, quota) if quota else limit


def full_ram(rows, config=None):
    if not rows:
        return 0
    total = 16 * GIB
    for row in rows:
        size = int(row.get('known_input_total_bytes') or row.get('wgs_input_bytes') or 0)
        if row['platform'] == 'illumina':
            total += 192 * GIB + max(8 * GIB, size * 3 // 5)
        else:
            plan = longread_memory_plan(row, config)
            total += plan['admission_bytes'] if plan else max(4 * size + 160_000_000_000, 192 * GIB)
    return total


def live_leases(pool):
    root=Path(pool)
    if not (root/'leases').is_dir():return {}
    result = {}
    # 큐의 읽기 판정은 원자 JSON 스냅샷만 사용한다. 실제 예약은 할당기가 잠금 안에서 재검사한다.
    for path in (root/'leases').glob('*.json'):
        try:
            row = json.loads(path.read_text())
            if lease_alive(row):
                result.setdefault(row['sample'], []).append(row)
        except (OSError, ValueError, KeyError):
            continue
    return result


def active_ram(pool, rows, candidate=None, final_stage='12', config=None):
    all_rows = [*rows, *([candidate] if candidate else [])]
    if not all_rows:
        return 0
    current = live_leases(pool)
    total = 16 * GIB
    for row in all_rows:
        full = full_ram([row], config) - 16 * GIB
        entry=starting_memory(row,config)
        fallback=entry if entry is not None else full
        leases = current.get(row['sample_id'], []) if row is not candidate else []
        # 0/미검증 예약은 전체 모델로 되돌리며, 살아있는 단계의 실제 예약은 중복하지 않음.
        owned = [lease for lease in leases if lease.get('branch') in {'WGS','RNA','QC'}
                 and lease_matches_config(lease, config)]
        validated = owned and all(int(lease.get('memory_bytes',0)) > 0 for lease in owned)
        total += sum(int(lease['memory_bytes']) for lease in owned) if validated else fallback
    return total


def output_space(config):
    """입력 tmpfs와 독립적으로 실제 출력 저장소 여유를 확인함."""
    minimum = int(config.get('min_output_free_bytes', 0))
    if minimum < 0: raise ValueError('min_output_free_bytes는 0 이상이어야 함')
    if not minimum: return True, ''
    try:
        path = Path(config['output_root'])
        while not path.exists(): path = path.parent
        usage = os.statvfs(path)
        free = usage.f_bavail * usage.f_frsize
    except OSError:
        return False, '최종 출력 디스크 여유 확인 실패'
    return (True, '') if free >= minimum else (False, f'최종 출력 디스크 공간 대기: {free}/{minimum}B')


def admit(config, pool, rows, candidate=None, prefetch=(), download_target=3):
    allowed, reason = output_space(config)
    if not allowed: return False, reason
    from integrate_wgs_rna.gpu_compare_resources import enabled,config_environment
    if enabled(config_environment(config)):return True,'GPU 비교 CPU/RAM 예약 예외'
    if candidate and config.get('cpu_allocation_mode') == 'fixed_sample':
        from integrate_wgs_rna.sample_cpu import available
        Path(pool).mkdir(parents=True, exist_ok=True)
        with (Path(pool)/'pool.lock').open('a+') as handle:
            import fcntl
            fcntl.flock(handle, fcntl.LOCK_EX)
            if available(pool, cpu_limit(config)) < config['sample_cpus']:
                return False, '고정 샘플 CPU 확보 대기'
    current, cgroup, _ = cgroup_limits()
    current = pressure_usage(current, cgroup)
    limit = effective_memory(config) * 9 // 10
    download_ram = sum(int(lease.get('memory_bytes', 0))
                       for group in live_leases(pool).values() for lease in group
                       if lease.get('branch') == 'DOWNLOAD')
    peak = active_ram(pool, rows, candidate, final_stage=config.get('_wgs_compact_stage', '12'), config=config) + download_ram
    if peak > limit:
        return False, f'RAM 예약 초과: {peak}/{limit}'
    increment = active_ram(pool,[],candidate,config=config)-16*GIB if candidate else 0
    if current is not None and cgroup and current + increment + download_ram > cgroup * 9 // 10:
        return False, 'cgroup 실측+신규 예약 여유 부족'
    try:
        guard = json.loads((Path(pool) / 'memory_guard.json').read_text())
        if time.time() - float(guard['time']) <= 30 and int(guard['rss_bytes']) + increment + download_ram > limit:
            return False, '전처리 RSS+신규 예약 여유 부족'
    except (OSError, ValueError, KeyError):
        pass
    return True, ''


def download_admit(config, pool, *, local_only=False):
    """다음 샘플의 전장 RAM 대신 실제 전송 버퍼 피크만 추가 예약함."""
    from integrate_wgs_rna.gpu_compare_resources import enabled,config_environment
    if enabled(config_environment(config)):return True,'명시 비교 입력 검증·전송 예약 예외'
    current, cgroup, _ = cgroup_limits()
    current = pressure_usage(current, cgroup)
    limit = effective_memory(config) * 9 // 10
    increment = memory_plan(config)['reservation_bytes']
    leases = live_leases(pool)
    try:
        policy=json.loads((Path(pool)/'policy.json').read_text())
    except FileNotFoundError:
        policy={}
    limit_jobs=policy.get('max_download_jobs')
    active_downloads=sum(1 for group in leases.values() for row in group if row.get('branch')=='DOWNLOAD')
    if not local_only and limit_jobs is not None and active_downloads>=int(limit_jobs):
        return False,'전처리 계산 우선: 선행 다운로드 동시성 제한'
    floor=int(policy.get('download_min_free_bytes',0))
    if floor and not local_only:
        try:
            stat=os.statvfs(config['data_root'])
            if stat.f_bavail*stat.f_frsize<floor:
                return False,'디스크 여유 확보까지 선행 다운로드 대기'
        except OSError:
            return False,'다운로드 전 디스크 여유 확인 실패'
    reserved = sum(int(lease.get('memory_bytes', 0)) for group in leases.values() for lease in group)
    if reserved + increment > limit:
        return False, f'공유 RAM 예약+전송 {increment}B 상한 초과'
    if current is not None and cgroup and current + increment > cgroup * 9 // 10:
        return False, f'cgroup 실측+전송 {increment}B 여유 부족'
    try:
        guard = json.loads((Path(pool) / 'memory_guard.json').read_text())
        if time.time() - float(guard['time']) <= 30 and int(guard['rss_bytes']) + increment > limit:
            return False, f'전처리 RSS+전송 {increment}B 여유 부족'
    except (OSError, ValueError, KeyError):
        pass
    return True, ''
