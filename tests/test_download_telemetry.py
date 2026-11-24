"""실제 HTTP의 재개·재시도·검증 실패·병렬 전송률 기록을 검사한다."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import time
import unittest
from unittest.mock import patch

from handoff import downloads, download_retry, http_download
from test_http_integrity import LocalRemote


class DownloadTelemetryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.data = self.root / 'data'
        self.state = self.root / 'state'

    def tearDown(self):
        self.tmp.cleanup()

    def fetch(self, url, name='input.bin', **kwargs):
        return downloads.fetch({'path': name, 'url': url, 'kind': 'test_binary', **kwargs}, self.data, self.state)

    def events(self):
        return [json.loads(line) for line in (self.state / 'stats/download_transfers.jsonl').read_text().splitlines()]

    def test_http_complete_records_body_and_skips_cached_files(self):
        data = b'ACGT' * 4096
        with LocalRemote(data) as remote:
            url = remote.url + '?X-Amz-Signature=SECRET_TOKEN#private'
            target = self.fetch(url, size_bytes=len(data))
            self.assertEqual(target.read_bytes(), data)
            events = self.events()
            self.assertEqual([r['event'] for r in events], ['start', 'progress', 'complete'])
            last = events[-1]
            self.assertEqual(last['new_file_bytes'], len(data))
            self.assertEqual(last['received_body_bytes'], len(data) + 1)
            self.assertEqual(last['network_requests'], 2)
            self.assertGreater(last['network_MBps'], 0)
            self.assertAlmostEqual(last['network_Gbps'], last['network_MBps'] * .008)
            self.assertNotIn('SECRET_TOKEN', json.dumps(events))
            self.assertNotIn('private', last['url'])
            self.fetch(url, size_bytes=len(data))
            self.assertEqual(self.events(), events)
        receipts = list((self.data / '.download_records').glob('*.json'))
        receipt = next(json.loads(p.read_text()) for p in receipts if not p.name.endswith('.http.json'))
        self.assertEqual(receipt['transfer_metrics']['new_file_bytes'], len(data))
        shared = (self.data / '.download_records/transfer_rates.jsonl').read_text().splitlines()
        self.assertEqual(len(shared), len(events))

    def test_resume_counts_only_new_file_bytes(self):
        data = bytes(range(256)) * 8
        self.data.mkdir()
        (self.data / 'input.bin.part').write_bytes(data[:512])
        with LocalRemote(data) as remote:
            self.fetch(remote.url, size_bytes=len(data))
        last = self.events()[-1]
        self.assertEqual(last['resume_offset_bytes'], 512)
        self.assertEqual(last['new_file_bytes'], len(data) - 512)
        self.assertGreater(last['received_body_bytes'], last['new_file_bytes'])
        self.assertEqual((self.data / 'input.bin').read_bytes(), data)

    def test_retried_payload_is_received_twice_but_saved_once(self):
        data = bytes(range(256)) * 4
        with LocalRemote(data, 'xml') as remote, patch.object(http_download.time, 'sleep'):
            self.fetch(remote.url, size_bytes=len(data))
        last = self.events()[-1]
        self.assertEqual(last['new_file_bytes'], len(data))
        self.assertEqual(last['received_body_bytes'], 2 * len(data) + 1)
        self.assertEqual(last['retry_count'], 1)

    def test_remote_change_records_failure_without_publishing(self):
        with LocalRemote(b'x' * 1024, 'etag') as remote:
            with self.assertRaises(http_download.RemoteChanged):
                self.fetch(remote.url)
        last = self.events()[-1]
        self.assertEqual(last['event'], 'failed')
        self.assertEqual(last['failure_type'], 'RemoteChanged')
        self.assertEqual(last['new_file_bytes'], 0)
        self.assertFalse((self.data / 'input.bin').exists())

    def test_checksum_failure_preserves_partial_and_records_failure(self):
        data = b'x' * 1024
        with LocalRemote(data) as remote:
            with self.assertRaises(ValueError):
                self.fetch(remote.url, size_bytes=len(data), checksum='0' * 64, checksum_type='sha256')
        last = self.events()[-1]
        self.assertEqual(last['event'], 'failed')
        self.assertEqual(last['new_file_bytes'], len(data))
        self.assertEqual((self.data / 'input.bin.part').read_bytes(), data)
        self.assertFalse((self.data / 'input.bin').exists())

    def test_verification_delay_is_separate_from_network(self):
        data = b'x' * 1024
        original = downloads.verify
        def delayed(path, *args):
            if str(path).endswith('.part'):
                time.sleep(.05)
            return original(path, *args)
        with LocalRemote(data) as remote, patch.object(downloads, 'verify', side_effect=delayed):
            self.fetch(remote.url, size_bytes=len(data))
        last = self.events()[-1]
        self.assertGreaterEqual(last['verification_seconds'], .05)
        self.assertGreaterEqual(last['elapsed_seconds'] - last['network_seconds'], .05)

    def test_threads_keep_independent_counters_and_valid_json_lines(self):
        data = b'x' * 4096
        with LocalRemote(data) as remote, ThreadPoolExecutor(2) as pool:
            paths = list(pool.map(lambda name: self.fetch(remote.url, name, size_bytes=len(data)), ['a.bin', 'b.bin']))
        self.assertEqual(len(paths), 2)
        events = [r for r in self.events() if r['event'] == 'complete']
        self.assertEqual(len(events), 2)
        self.assertEqual(len({r['transfer_id'] for r in events}), 2)
        self.assertTrue(all(r['new_file_bytes'] == len(data) and r['received_body_bytes'] == len(data)+1 for r in events))

    def test_curl_retry_uses_reported_bytes_time_and_current_offset(self):
        partial = self.data / 'input.bin.part'
        seen = []
        def run(command, **kwargs):
            seen.append(partial.stat().st_size if partial.exists() else 0)
            with partial.open('ab') as handle:
                handle.write(b'abc' if len(seen) == 1 else b'def')
            return subprocess.CompletedProcess(command, 56 if len(seen) == 1 else 0, '3\t1.5\t206')
        with patch.object(download_retry.subprocess, 'run', side_effect=run), patch.object(download_retry.time, 'sleep'):
            self.fetch('ftp://example.test/input', size_bytes=6)
        last = self.events()[-1]
        self.assertEqual(seen, [0, 3])
        self.assertEqual(last['new_file_bytes'], 6)
        self.assertEqual(last['received_body_bytes'], 6)
        self.assertEqual(last['network_seconds'], 3)
        self.assertEqual(last['retry_count'], 1)


if __name__ == '__main__':
    unittest.main()
