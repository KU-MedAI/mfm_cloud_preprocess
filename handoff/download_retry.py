"""끊긴 전송마다 현재 부분 파일 크기를 다시 읽어 Range 다운로드를 재개한다."""
from pathlib import Path
import subprocess
import time
from .download_telemetry import NetworkRequest, retry_wait, progress

RETRY_CODES = {5, 6, 7, 18, 28, 35, 52, 55, 56, 92}
RETRY_HTTP = {408, 429, 500, 502, 503, 504}

def transfer(url, partial, complete, attempts=12):
    partial = Path(partial)
    from urllib.error import HTTPError
    from .http_download import RangeUnsupported, transfer as ranged_transfer, error_offsets, InvalidResponse
    if url.startswith(('https://','http://')):
        try:
            return ranged_transfer(url,partial,complete,attempts)
        except RangeUnsupported:
            pass
        except HTTPError as error:
            raise subprocess.CalledProcessError(22,['HTTP',url],output=str(error.code)) from error
    if partial.exists() and error_offsets(partial):
        raise InvalidResponse('기존 부분 파일에 XML 오류 응답이 있습니다. 별도 복구 경로가 필요합니다.')
    for attempt in range(1, attempts + 1):
        if complete(partial):
            return
        offset = partial.stat().st_size if partial.exists() else 0
        print(f'[download] attempt={attempt}/{attempts} resume_bytes={offset} path={partial}', flush=True)
        command = ['curl', '--fail', '--location', '--silent', '--show-error',
                   '--connect-timeout', '30', '--speed-limit', '1024', '--speed-time', '120',
                   '--continue-at', '-', '--output', str(partial), '--write-out',
                   '%{size_download}\t%{time_total}\t%{http_code}', url]
        with NetworkRequest() as observed:
            result = subprocess.run(command, stdout=subprocess.PIPE, text=True)
            fields = result.stdout.strip().split('\t')
            try:
                observed.body_bytes = int(float(fields[-3]))
                observed.seconds = float(fields[-2])
            except (ValueError, IndexError):
                observed.body_bytes = max(0, partial.stat().st_size-offset) if partial.exists() else 0
        if partial.exists():
            progress(partial.stat().st_size)
        if partial.exists() and error_offsets(partial):
            raise InvalidResponse('전송 파일에 XML 오류 응답 혼입. 손상 부분 파일 보존.')
        if result.returncode == 0:
            return
        # EOF 직후 연결이 끊겨도 전체 크기/checksum이 맞으면 완료된 전송이다.
        if complete(partial):
            return
        try:
            status = int(result.stdout.strip()[-3:])
        except (ValueError, TypeError):
            status = 0
        retryable = result.returncode in RETRY_CODES or (result.returncode == 22 and status in RETRY_HTTP)
        if not retryable or attempt == attempts:
            raise subprocess.CalledProcessError(result.returncode, command, output=result.stdout)
        delay = min(30, 3 * 2 ** (attempt - 1))
        print(f'[download] retry curl={result.returncode} http={status} delay_seconds={delay}', flush=True)
        with retry_wait():
            time.sleep(delay)
