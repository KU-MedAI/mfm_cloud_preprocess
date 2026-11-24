"""실제 HTTP 병렬 수신의 순서·재개·오류 차단·수신량을 검증함."""
from pathlib import Path
from unittest.mock import patch
import os
import tempfile
import threading
import unittest
from handoff import http_download as http
from handoff.download_telemetry import TransferLog
from test_http_integrity import LocalRemote


class ParallelRangeTests(unittest.TestCase):
    def test_실제HTTP가동시에수신되어도파일과계수가동일함(self):
        data = bytes(range(256))*64
        original = http.read_range
        lock = threading.Lock()
        gate = threading.Event()
        active = peak = 0
        def read(*args):
            nonlocal active, peak
            with lock:
                active += 1
                peak = max(peak, active)
                if active >= 4:
                    gate.set()
            gate.wait(2)
            try:
                return original(*args)
            finally:
                with lock:
                    active -= 1
        with tempfile.TemporaryDirectory() as raw, LocalRemote(data) as remote:
            root = Path(raw)
            part = root/'input.part'
            with patch.dict(os.environ, {'MFM_HTTP_RANGE_WORKERS': '4'}), patch.object(http, 'read_range', side_effect=read):
                with TransferLog(remote.url, root/'input', part, root/'state') as metrics:
                    with metrics.phase('transfer'):
                        http.transfer(remote.url, part, lambda _: False, chunk_bytes=1024)
            self.assertEqual(part.read_bytes(), data)
            self.assertEqual(metrics.received_body_bytes, len(data)+1)
            self.assertEqual(metrics.network_requests, 17)
            self.assertEqual(peak, 4)

    def test_검증된연속prefix에서재개함(self):
        data = bytes(range(256))*32
        with tempfile.TemporaryDirectory() as raw, LocalRemote(data) as remote, patch.dict(os.environ, {'MFM_HTTP_RANGE_WORKERS': '4'}):
            part = Path(raw)/'input.part'
            part.write_bytes(data[:1536])
            http.transfer(remote.url, part, lambda _: False, chunk_bytes=1024)
            self.assertEqual(part.read_bytes(), data)
            self.assertIn((1536, 2559), remote.requests)

    def test_XML재시도도수신량에반영하고정상구간만저장함(self):
        data = bytes(range(256))*32
        with tempfile.TemporaryDirectory() as raw, LocalRemote(data, 'xml') as remote:
            root = Path(raw)
            part = root/'input.part'
            with patch.dict(os.environ, {'MFM_HTTP_RANGE_WORKERS': '4'}), patch.object(http.time, 'sleep'):
                with TransferLog(remote.url, root/'input', part, root/'state') as metrics:
                    http.transfer(remote.url, part, lambda _: False, chunk_bytes=1024)
            self.assertEqual(part.read_bytes(), data)
            self.assertEqual(metrics.received_body_bytes, len(data)+1024+1)
            self.assertEqual(metrics.retries, 1)

    def test_원격버전변경시저장된prefix를섞지않음(self):
        data = bytes(range(256))*32
        with tempfile.TemporaryDirectory() as raw, LocalRemote(data, 'etag') as remote, patch.dict(os.environ, {'MFM_HTTP_RANGE_WORKERS': '4'}):
            part = Path(raw)/'input.part'
            with self.assertRaises(http.RemoteChanged):
                http.transfer(remote.url, part, lambda _: False, chunk_bytes=1024)
            if part.exists():
                self.assertEqual(part.read_bytes(), data[:part.stat().st_size])


if __name__ == '__main__':
    unittest.main()
