"""이미 있는 개인 입력의 검증·재사용과 검증 조건 변경을 검사한다."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from handoff import downloads


class InputReuseTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / 'remote.bin'
        self.source.write_bytes(b'ACGT')
        self.data = self.root / 'data'
        self.state = self.root / 'state'
        self.resource = {'path': 'wgs/input.bin', 'url': self.source.as_uri(),
                         'size_bytes': 4, 'checksum_type': 'sha256',
                         'checksum': hashlib.sha256(b'ACGT').hexdigest()}

    def fetch(self):
        return downloads.fetch(self.resource, self.data, self.state)

    def receipt(self, path):
        key = hashlib.sha256(str(path).encode()).hexdigest()
        return self.data / '.download_records' / (key + '.json')

    @staticmethod
    def signature(path):
        stat = path.stat()
        return {'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns,
                'ctime_ns': stat.st_ctime_ns}

    def test_missing_downloads_once_and_cached_input_skips_transfer(self):
        target = self.fetch()
        self.assertEqual(target.read_bytes(), b'ACGT')
        with patch.object(downloads.shutil, 'copyfile', side_effect=AssertionError('중복 전송')):
            self.assertEqual(self.fetch(), target)

    def test_preexisting_valid_input_is_cached_and_stays_unowned(self):
        target = self.data / self.resource['path']
        target.parent.mkdir(parents=True)
        target.write_bytes(b'ACGT')
        with patch.object(downloads.shutil, 'copyfile', side_effect=AssertionError('불필요한 전송')):
            self.assertEqual(self.fetch(), target)
        report = json.loads(self.receipt(target).read_text())
        self.assertFalse(report['owned'])
        self.assertEqual(report['status'], 'verified')
        with patch.object(downloads, 'verify', side_effect=AssertionError('변경 없는 원본 재해시')):
            self.assertEqual(self.fetch(), target)

    def test_invalid_preexisting_input_is_preserved(self):
        target = self.data / self.resource['path']
        target.parent.mkdir(parents=True)
        target.write_bytes(b'BAD!')
        with patch.object(downloads.shutil, 'copyfile') as transfer:
            with self.assertRaisesRegex(ValueError, '기존 파일 검증 실패'):
                self.fetch()
        transfer.assert_not_called()
        self.assertEqual(target.read_bytes(), b'BAD!')

    def test_changed_checksum_invalidates_old_receipt(self):
        target = self.fetch()
        self.source.write_bytes(b'TGCA')
        self.resource['checksum'] = hashlib.sha256(b'TGCA').hexdigest()
        self.fetch()
        self.assertEqual(target.read_bytes(), b'TGCA')

    def test_changed_size_invalidates_old_receipt_without_checksum(self):
        self.resource['checksum'] = None
        target = self.fetch()
        self.source.write_bytes(b'ACGTX')
        self.resource['size_bytes'] = 5
        self.fetch()
        self.assertEqual(target.read_bytes(), b'ACGTX')

    def test_generated_index_is_refreshed_when_parent_changes(self):
        target = self.fetch()
        parent = self.root / 'parent.bam'
        parent.write_bytes(b'parent')
        receipt = self.receipt(target)
        report = json.loads(receipt.read_text())
        report.update(status='generated_index', generated_from=str(parent),
                      generated_from_signature=self.signature(parent))
        receipt.write_text(json.dumps(report))
        self.resource.update(size_bytes=None, checksum=None)
        self.source.write_bytes(b'new-index')
        parent.write_bytes(b'changed-parent')
        self.fetch()
        self.assertEqual(target.read_bytes(), b'new-index')

    def test_unchanged_generated_index_keeps_local_result(self):
        target = self.fetch()
        parent = self.root / 'parent.bam'
        parent.write_bytes(b'parent')
        receipt = self.receipt(target)
        report = json.loads(receipt.read_text())
        report.update(status='generated_index', generated_from=str(parent),
                      generated_from_signature=self.signature(parent))
        receipt.write_text(json.dumps(report))
        self.resource.update(size_bytes=100, checksum='0' * 64)
        with patch.object(downloads.shutil, 'copyfile', side_effect=AssertionError('정상 로컬 index 재수신')):
            self.assertEqual(self.fetch(), target)


if __name__ == '__main__':
    unittest.main()
