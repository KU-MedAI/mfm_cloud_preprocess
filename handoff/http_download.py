"""원격 버전을 고정하고 검증된 작은 HTTP 구간만 다운로드 파일에 저장한다."""
from __future__ import annotations
import http.client
import os
import json
import re
import socket
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
from urllib.parse import urlsplit, urlunsplit
from .download_telemetry import NetworkRequest, retry_wait, progress

CHUNK_BYTES = 32 << 20
ERROR_XML = re.compile(rb'<Error>\s*<Code>[^<]{1,256}</Code>\s*<Message>.{0,8192}?</Message>.{0,4096}?</Error>', re.S)
RETRY_HTTP = {408, 429, 500, 502, 503, 504}
_EBI_HTTPS_MIRRORS = {}


def _open(request, timeout):
    """ENA HTTPS403만 공식 EBI 동일 경로로 재시도하며 원본 식별자는 유지함."""
    parts = urlsplit(request.full_url)
    eligible = (parts.scheme == 'https' and parts.netloc == 'ftp.sra.ebi.ac.uk'
                and parts.path.startswith('/vol1/'))
    def mirror_request():
        url = urlunsplit(parts._replace(netloc='ftp.ebi.ac.uk'))
        return Request(url, headers=dict(request.header_items()), method=request.get_method())
    if eligible and parts.netloc in _EBI_HTTPS_MIRRORS:
        try:
            return urlopen(mirror_request(), timeout=timeout)
        except HTTPError as error:
            if error.code != 403:
                raise
            error.close()
            _EBI_HTTPS_MIRRORS.pop(parts.netloc, None)
    try:
        return urlopen(request, timeout=timeout)
    except HTTPError as error:
        if not eligible or error.code != 403:
            raise
        error.close()
        alternative = mirror_request()
        response = urlopen(alternative, timeout=timeout)
        _EBI_HTTPS_MIRRORS[parts.netloc] = alternative.host
        print(f'[download-range] HTTPS403 official_mirror={alternative.full_url}', flush=True)
        return response


def _retry_http(error, url):
    parts = urlsplit(url)
    return error.code in RETRY_HTTP or (error.code == 403 and parts.scheme == 'https'
            and parts.netloc == 'ftp.sra.ebi.ac.uk' and parts.path.startswith('/vol1/'))

class RangeUnsupported(ValueError):
    """Range 미지원 서버는 기존 전송 방식으로 처리한다."""

class InvalidResponse(ValueError):
    """길이/구간/오류 응답 검증 실패는 아직 파일에 저장하지 않는다."""

class RemoteChanged(ValueError):
    """원격 버전 변경 시 서로 다른 버전의 데이터를 혼합하지 않는다."""

@dataclass(frozen=True)
class RemoteObject:
    url: str
    size: int
    etag: str | None
    modified: str | None

    def headers(self):
        values = {'Accept-Encoding': 'identity'}
        if self.etag and not self.etag.startswith('W/'):
            values['If-Match'] = self.etag
        elif self.modified:
            values['If-Unmodified-Since'] = self.modified
        return values

def error_payload(data):
    return ERROR_XML.search(data)

def error_offsets(path, chunk_bytes=8 << 20):
    """구간 경계를 걸친 XML도 찾고 전체 파일을 메모리에 올리지 않는다."""
    result = []
    tail = b''
    consumed = 0
    with Path(path).open('rb') as stream:
        while chunk := stream.read(chunk_bytes):
            data = tail + chunk
            base = consumed - len(tail)
            for found in ERROR_XML.finditer(data):
                offset = base + found.start()
                if not result or offset > result[-1]:
                    result.append(offset)
            consumed += len(chunk)
            tail = data[-16384:]
    return result

def _identity(response):
    return response.headers.get('ETag'), response.headers.get('Last-Modified')

def _identify_once(url, timeout=45):
    request = Request(url, headers={'Range':'bytes=0-0', 'Accept-Encoding':'identity'})
    with NetworkRequest() as observed:
        with _open(request, timeout=timeout) as response:
            if response.status != 206:
                raise RangeUnsupported('서버가 HTTP Range를 지원하지 않습니다.')
            match = re.fullmatch(r'bytes 0-0/(\d+)', response.headers.get('Content-Range',''))
            if not match or len(observed.read(response, 2)) != 1:
                raise InvalidResponse('원격 파일의 크기/Range 식별 실패')
            etag, modified = _identity(response)
            if not etag and not modified:
                raise RangeUnsupported('원격 버전을 고정할 식별 정보가 없습니다.')
            return RemoteObject(url, int(match.group(1)), etag, modified)

def identify(url, timeout=45, attempts=12):
    for attempt in range(1,attempts+1):
        try:
            return _identify_once(url,timeout)
        except HTTPError as error:
            error.close()
            if not _retry_http(error, url) or attempt==attempts:raise
        except (InvalidResponse,URLError,http.client.HTTPException,socket.timeout,ConnectionError):
            if attempt==attempts:raise
        delay=min(30,3*2**(attempt-1))
        print(f'[download-range] identify_retry={attempt}/{attempts} delay={delay}',flush=True)
        with retry_wait():
            time.sleep(delay)

def read_range(remote, start, end, attempts=12, timeout=120):
    """응답 전체를 임시 메모리에서 검사하며 실패 응답은 호출자에게 전달하지 않는다."""
    if not 0 <= start <= end < remote.size or end-start+1 > CHUNK_BYTES:
        raise ValueError('원격 Range 범위가 올바르지 않습니다.')
    wanted = end-start+1
    for attempt in range(1, attempts+1):
        try:
            headers = dict(remote.headers(), Range=f'bytes={start}-{end}')
            with NetworkRequest() as observed:
                with _open(Request(remote.url,headers=headers),timeout=timeout) as response:
                    if response.status != 206:
                        raise InvalidResponse('부분 요청에 대한206 응답이 아닙니다.')
                    if response.headers.get('Content-Range') != f'bytes {start}-{end}/{remote.size}':
                        raise InvalidResponse('요청 범위와 Content-Range가 다릅니다.')
                    etag, modified = _identity(response)
                    if (remote.etag and etag != remote.etag) or (remote.modified and modified != remote.modified):
                        raise RemoteChanged('원격 파일 버전이 바뀌었습니다.')
                    encoding = response.headers.get('Content-Encoding','identity')
                    if encoding not in {'identity',''}:
                        raise InvalidResponse('HTTP 전송 중 데이터 인코딩이 변경됐습니다.')
                    data = observed.read(response, wanted+1)
            if len(data) != wanted:
                raise InvalidResponse(f'부분 응답 크기 불일치: {len(data)}/{wanted}')
            if error_payload(data):
                raise InvalidResponse('바이너리 내부 XML 오류 응답 혼입')
            return data
        except HTTPError as error:
            error.close()
            if error.code == 412:
                raise RemoteChanged('원격 파일 조건부 요청 실패:412') from error
            if not _retry_http(error, remote.url) or attempt == attempts:
                raise
            detail = f'HTTP{error.code}'
        except (InvalidResponse, URLError, http.client.HTTPException, socket.timeout, ConnectionError) as error:
            if attempt == attempts:
                raise
            detail = type(error).__name__ + ': ' + str(error)
        delay = min(30, 3*2**(attempt-1))
        print(f'[download-range] retry={attempt}/{attempts} range={start}-{end} reason={detail} delay={delay}',flush=True)
        with retry_wait():
            time.sleep(delay)

def transfer(url, partial, complete, attempts=12, chunk_bytes=CHUNK_BYTES):
    partial = Path(partial)
    remote = identify(url)
    metadata = Path(str(partial)+'.http.json')
    identity = asdict(remote)
    # 본문이 없으면 혼합할 데이터도 없으므로 오래된 표지만으로 새 전송을 막지 않는다.
    has_partial_data = partial.exists() and partial.stat().st_size > 0
    if has_partial_data and metadata.exists() and json.loads(metadata.read_text()) != identity:
        raise RemoteChanged('기존 부분 다운로드와 원격 버전이 다릅니다. 기존 파일을 보존합니다.')
    # 이전 버전 downloader가 남긴 XML을 새 전송으로 이어붙이지 않는다.
    if partial.exists():
        offsets = error_offsets(partial)
        if offsets:
            raise InvalidResponse('기존 부분 파일에 오류 응답이 있습니다. 별도 복구 경로가 필요합니다.')
        if not metadata.exists() and partial.stat().st_size:
            # 구버전 부분 파일의 시작/중간/마지막 구간을 원격과 비교한다.
            size=partial.stat().st_size
            if size>remote.size:raise InvalidResponse('부분 파일이 원격보다 큽니다.')
            starts=sorted({0,max(0,size//2-32768),max(0,size-65536)})
            with partial.open('rb') as source:
                for start in starts:
                    end=min(size-1,start+65535)
                    source.seek(start)
                    if source.read(end-start+1)!=read_range(remote,start,end,attempts):
                        raise InvalidResponse('기존 부분 파일과 원격 구간이 다릅니다. 기존 파일을 보존합니다.')
    metadata.write_text(json.dumps(identity,ensure_ascii=False,indent=2)+'\n')
    offset = partial.stat().st_size if partial.exists() else 0
    if offset > remote.size:
        raise InvalidResponse('부분 파일이 원격 파일보다 큽니다.')
    workers = int(os.environ.get('MFM_HTTP_RANGE_WORKERS', '1'))
    if not 1 <= workers <= 16:
        raise ValueError('MFM_HTTP_RANGE_WORKERS는 1~16이어야 함')
    if workers > 1 and offset < remote.size:
        from .http_parallel import transfer_ranges
        return transfer_ranges(remote, partial, offset, chunk_bytes, attempts, workers)
    while offset < remote.size:
        end = min(remote.size-1,offset+chunk_bytes-1)
        data = read_range(remote,offset,end,attempts)
        # 경계 양쪽에 나뉜 XML 응답을 확인한 뒤 검증된 구간만 추가한다.
        if partial.exists():
            with partial.open('rb') as previous:
                previous.seek(max(0,offset-16384));tail=previous.read()
            if error_payload(tail+data[:16384]):
                raise InvalidResponse('구간 경계에서 XML 오류 응답 혼입')
        with partial.open('ab') as stream:
            stream.write(data)
        offset=end+1
        progress(offset)
        print(f'[download-range] verified_bytes={offset}/{remote.size} path={partial}',flush=True)
    # 전체 크기와 checksum의 최종 판정은 fetch에서 유지한다.
