"""다운로드 본문 수신률·재개 바이트·검증 시간을 서로 구분해 기록한다."""
from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from datetime import datetime, timezone
import fcntl
import http.client
import json
from pathlib import Path
import time
import threading
from urllib.parse import urlsplit, urlunsplit
import uuid

_active = ContextVar('mfm_download_transfer', default=None)


def public_url(url):
    """로그에는 사용자 정보·서명 query·fragment를 기록하지 않는다."""
    parts = urlsplit(url)
    host = parts.hostname or ''
    if ':' in host:
        host = '[' + host + ']'
    if parts.port:
        host += ':' + str(parts.port)
    return urlunsplit((parts.scheme, host, parts.path, '', ''))


class NetworkRequest:
    def __enter__(self):
        self.owner = _active.get()
        self.started = time.perf_counter()
        self.body_bytes = 0
        self.seconds = None
        return self

    def read(self, response, size):
        try:
            data = response.read(size)
        except http.client.IncompleteRead as error:
            self.body_bytes += len(error.partial)
            raise
        self.body_bytes += len(data)
        return data

    def __exit__(self, kind, error, traceback):
        if self.owner is not None:
            with self.owner._lock:
                self.owner.network_seconds += self.seconds if self.seconds is not None else time.perf_counter() - self.started
                self.owner.received_body_bytes += self.body_bytes
                self.owner.network_requests += 1
                self.owner.failed_requests += kind is not None


@contextmanager
def retry_wait():
    owner = _active.get()
    if owner is not None:
        with owner._lock:
            owner.retries += 1
    started = time.perf_counter()
    try:
        yield
    finally:
        if owner is not None:
            with owner._lock:
                owner.retry_wait_seconds += time.perf_counter() - started


def progress(file_bytes):
    owner = _active.get()
    if owner is not None:
        owner.progress(file_bytes)


class TransferLog:
    def __init__(self, url, target, partial, state_root, data_root=None, *, kind='', purpose='download'):
        self._lock = threading.RLock()
        self.url = public_url(url)
        self.source_host = urlsplit(url).hostname or ''
        self.target = Path(target)
        self.partial = Path(partial)
        self.initial_bytes = self.partial.stat().st_size if self.partial.exists() else 0
        self.new_file_bytes = 0
        self.network_seconds = 0.0
        self.received_body_bytes = 0
        self.network_requests = 0
        self.failed_requests = 0
        self.retry_wait_seconds = 0.0
        self.retries = 0
        self.phases = {'transfer': 0.0, 'verification': 0.0}
        self.kind = kind
        self.purpose = purpose
        self.transfer_id = uuid.uuid4().hex
        self.status = 'running'
        self.failure_type = None
        self.last_progress = None
        self.finished = None
        self.finished_clock = None
        self.log_errors = []
        paths = [Path(state_root) / 'stats/download_transfers.jsonl']
        if data_root is not None:
            paths.append(Path(data_root) / '.download_records/transfer_rates.jsonl')
        self.log_paths = list(dict.fromkeys(p.resolve() for p in paths))

    def __enter__(self):
        self.started = time.time()
        self.started_clock = time.perf_counter()
        self.token = _active.set(self)
        self.emit('start')
        return self

    @contextmanager
    def phase(self, name):
        started = time.perf_counter()
        try:
            yield
        finally:
            self.phases[name] += time.perf_counter() - started

    def progress(self, file_bytes):
        self.new_file_bytes = max(0, file_bytes - self.initial_bytes)
        now = time.perf_counter()
        if self.last_progress is None or now - self.last_progress >= 30:
            self.last_progress = now
            self.emit('progress')

    def summary(self):
        now = self.finished_clock if self.finished_clock is not None else time.perf_counter()
        network_rate = self.received_body_bytes / self.network_seconds if self.network_seconds else None
        transfer_seconds = self.phases['transfer']
        useful_rate = self.new_file_bytes / transfer_seconds if transfer_seconds else None
        return {
            'schema': 'mfm_download_transfer_v1', 'transfer_id': self.transfer_id,
            'url': self.url, 'source_host': self.source_host,
            'target': str(self.target), 'kind': self.kind, 'purpose': self.purpose,
            'status': self.status, 'started_at': datetime.fromtimestamp(self.started, timezone.utc).isoformat(),
            'finished_at': datetime.fromtimestamp(self.finished, timezone.utc).isoformat() if self.finished is not None else None,
            'elapsed_seconds': now - self.started_clock,
            'resume_offset_bytes': self.initial_bytes, 'new_file_bytes': self.new_file_bytes,
            'received_body_bytes': self.received_body_bytes,
            'network_seconds': self.network_seconds, 'network_requests': self.network_requests,
            'failed_requests': self.failed_requests, 'retry_count': self.retries,
            'retry_wait_seconds': self.retry_wait_seconds,
            'transfer_seconds': transfer_seconds, 'verification_seconds': self.phases['verification'],
            'network_MBps': network_rate / 1_000_000 if network_rate is not None else None,
            'network_Gbps': network_rate * 8 / 1_000_000_000 if network_rate is not None else None,
            'useful_transfer_MBps': useful_rate / 1_000_000 if useful_rate is not None else None,
            'failure_type': self.failure_type, 'log_errors': list(self.log_errors),
            'measurement_scope': 'HTTP/curl 수신 본문 기준(헤더/TLS 제외). network는 연결·요청·본문수신 시간이며 retry 대기·디스크쓰기·최종검증 제외. useful_transfer는 새 파일바이트/전체전송단계시간.',
        }

    def emit(self, event):
        payload = self.summary()
        payload.update(event=event, recorded_at=datetime.now(timezone.utc).isoformat())
        line = json.dumps(payload, ensure_ascii=False) + '\n'
        for path in self.log_paths:
            try:
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open('a') as handle:
                    fcntl.flock(handle, fcntl.LOCK_EX)
                    handle.write(line)
            except OSError as error:
                self.log_errors.append({'path': str(path), 'error_type': type(error).__name__})
                print(f'[download-rate] log_error={type(error).__name__} path={path}', flush=True)
        rate = payload['network_MBps']
        shown = 'NA' if rate is None else f'{rate:.3f}'
        gbps = payload['network_Gbps']
        shown_gbps = 'NA' if gbps is None else f'{gbps:.6f}'
        print(f'[download-rate] event={event} source={self.source_host} new_bytes={self.new_file_bytes} '
              f'received_bytes={self.received_body_bytes} network_seconds={self.network_seconds:.3f} '
              f'MBps={shown} Gbps={shown_gbps} target={self.target}', flush=True)

    def __exit__(self, kind, error, traceback):
        try:
            if self.partial.exists():
                self.new_file_bytes = max(0, self.partial.stat().st_size - self.initial_bytes)
            self.finished = time.time()
            self.finished_clock = time.perf_counter()
            self.status = 'failed' if kind is not None else 'complete'
            self.failure_type = kind.__name__ if kind is not None else None
            self.emit(self.status)
        finally:
            _active.reset(self.token)
