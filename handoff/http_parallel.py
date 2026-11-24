"""검증된 HTTP 구간을 병렬 수신하고 연속된 prefix만 순서대로 저장함."""
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from contextvars import copy_context
from .download_telemetry import progress


def transfer_ranges(remote, partial, offset, chunk_bytes, attempts, workers):
    from .http_download import read_range, error_payload, InvalidResponse, CHUNK_BYTES
    if not 1 <= workers <= 16 or not 1 <= chunk_bytes <= CHUNK_BYTES:
        raise ValueError('병렬 전송의 작업 수/버퍼 범위 오류')
    # 최대16×32MiB 결과와 수신 중 임시 복사본을 제한하여 기존4GB 예약 안에서 수행함.
    print(f'[download-parallel] workers={workers} max_pending_bytes={workers*chunk_bytes}', flush=True)
    pending = deque()
    next_start = offset
    def submit(pool):
        nonlocal next_start
        if next_start >= remote.size:
            return
        start = next_start
        end = min(remote.size-1, start+chunk_bytes-1)
        future = pool.submit(copy_context().run, read_range, remote, start, end, attempts)
        pending.append((start, end, future))
        next_start = end+1
    with ThreadPoolExecutor(max_workers=workers) as pool:
        try:
            for _ in range(workers):
                submit(pool)
            while pending:
                start, end, future = pending.popleft()
                data = future.result()
                if start != offset:
                    raise ValueError('검증 prefix의 구간 순서 불일치')
                if partial.exists():
                    with partial.open('rb') as previous:
                        previous.seek(max(0, offset-16384))
                        tail = previous.read()
                    if error_payload(tail+data[:16384]):
                        raise InvalidResponse('구간 경계에서 XML 오류 응답 혼입')
                with partial.open('ab') as stream:
                    stream.write(data)
                offset = end+1
                progress(offset)
                print(f'[download-range] verified_bytes={offset}/{remote.size} path={partial}', flush=True)
                submit(pool)
        finally:
            for _, _, future in pending:
                future.cancel()
