"""전송 중 입력 크기 합계와 디스크 안전 여유를 다운로드 전에 확인한다."""
from pathlib import Path
import json
import os


def memory_plan(config):
    workers = config.get('http_range_workers')
    if workers is None:
        workers = config.get('pipeline_env', {}).get(
            'MFM_HTTP_RANGE_WORKERS', os.environ.get('MFM_HTTP_RANGE_WORKERS'))
    if workers is None:
        pool = config.get('pipeline_env', {}).get('MFM_CPU_POOL_DIR')
        policy = Path(pool)/'policy.json' if pool else None
        settings = json.loads(policy.read_text()) if policy and policy.is_file() else {}
        workers = settings.get('download_range_workers', 1)
    if isinstance(workers, bool) or not str(workers).isdigit() or not 1 <= int(workers) <= 16:
        raise ValueError('HTTP 범위 worker 수는 1~16이어야 함')
    workers = int(workers)
    # Python/스레드512MiB, 범위별 응답·검증·쓰기32MiB 세 벌, 20%여유.
    peak = ((512 << 20) + workers * 3 * (32 << 20)) * 12 // 10
    return {'range_workers':workers, 'estimated_peak_bytes':peak,
            'reservation_bytes':max(1_000_000_000, peak)}


def storage_admit(config, pool, pending, candidate):
    maximum = config.get('prefetch_max_bytes')
    if maximum is None:
        return True, ''
    if isinstance(maximum, bool) or not isinstance(maximum, int) or maximum < 1:
        raise ValueError('prefetch_max_bytes는 양의 바이트 정수여야 함')
    rows = {row['sample_id']: row for row in [*pending, candidate]}
    sizes = [int(row.get('known_input_total_bytes') or 0) for row in rows.values()]
    if any(size <= 0 for size in sizes):
        return False, '선행 다운로드 입력 크기 사전 산정 필요'
    committed = sum(sizes)
    if committed > maximum:
        return False, f'선행 전송 입력 합계 상한: {committed}/{maximum}'
    path = Path(pool) / 'policy.json'
    policy = json.loads(path.read_text()) if path.exists() else {}
    floor = int(policy.get('download_min_free_bytes', 0))
    stat = os.statvfs(config['data_root'])
    free = stat.f_bavail * stat.f_frsize
    if committed + floor > free:
        return False, f'남은 전송 원본과 디스크 안전 여유 대기: {committed}+{floor}/{free}'
    return True, ''
