"""실제 BAM의 직렬/병렬·기존 samtools 일치와 손상·자원 보호를 검사한다."""
from array import array
from pathlib import Path
import json
import os
import tempfile
import unittest
from unittest.mock import patch
import pysam

from integrate_wgs_rna import alignment_record_stats as stats
from integrate_wgs_rna.record_stats_parallel import native_tool, run_scan, scan_plan, NativeUnavailable, RecordSizeLimit


class RecordStatsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tool = native_tool()

    def fixture(self, root):
        path = root / 'unsorted.bam'
        header = {'HD': {'VN': '1.6', 'SO': 'unsorted'}, 'SQ': [{'SN': 'chr1', 'LN': 1000000}]}
        flags = [0, 4, 256, 2048, 2304, 1024, 1536, 512, 3584, 260, 2052, 99, 147]
        expected = dict.fromkeys(['records', 'mapped_records', 'primary_records', 'secondary_records',
                                 'supplementary_records', 'duplicate_records', 'primary_duplicate_records',
                                 'hp_tagged_records', 'ps_tagged_records'], 0)
        with pysam.AlignmentFile(path, 'wb', header=header) as out:
            for i, flag in enumerate(flags * 9):
                read = pysam.AlignedSegment()
                read.query_name = f'read_{i}'
                read.query_sequence = 'ACGT' * (50 + i)
                read.query_qualities = array('B', [30] * len(read.query_sequence))
                read.flag = flag
                read.reference_id = -1 if i % 5 == 0 and flag & 4 else 0
                read.reference_start = -1 if read.reference_id == -1 else 10000 - i
                if not flag & 4: read.cigarstring = str(len(read.query_sequence)) + 'M'
                # HP/PS 뒤의 aux도 읽어 배열·실수·문자열을 모두 보존하는지 검사한다.
                if i % 2 == 0: read.set_tag('HP', i % 3)
                if i % 3 == 0: read.set_tag('PS', 'value_' + str(i), value_type='Z')
                read.set_tag('XA', array('i', [1, -2, 3]))
                read.set_tag('XB', array('f', [1.25, -2.5]))
                read.set_tag('XC', 'a', value_type='A')
                read.set_tag('XD', 1.25, value_type='d')
                read.set_tag('XE', '0123ABCD', value_type='H')
                out.write(read)
                expected['records'] += 1
                expected['mapped_records'] += not flag & 4
                expected['secondary_records'] += bool(flag & 256)
                expected['supplementary_records'] += bool(flag & 2048) and not flag & 256
                expected['primary_records'] += not flag & (256 | 2048)
                expected['duplicate_records'] += bool(flag & 1024)
                expected['primary_duplicate_records'] += bool(flag & 1024) and not flag & (256 | 2048)
                expected['hp_tagged_records'] += i % 2 == 0
                expected['ps_tagged_records'] += i % 3 == 0
        return path, expected

    def test_인덱스없이미배치_QC실패_복합flag_모든태그집계일치(self):
        with tempfile.TemporaryDirectory() as raw:
            path, expected = self.fixture(Path(raw))
            for cpus in [1, 2, 3, 4, 8, 16]:
                result = run_scan(path, cpus, tool=self.tool)
                self.assertEqual(result['counts'], expected)
                self.assertLessEqual(result['execution']['native_threads'], cpus)
                self.assertLessEqual(result['execution']['peak_rss_bytes'], scan_plan(cpus)['estimated_peak_rss_bytes'])
            with patch.dict(os.environ, {'MFM_CPU_BUDGET': '8', 'MFM_CPU_POOL_GRANTED': '8'}):
                self.assertEqual(stats._legacy_record_stats(path), expected)

    def test_여러BGZF블록과배치경계에걸친레코드도정확히한번(self):
        with tempfile.TemporaryDirectory() as raw:
            path, expected = self.fixture(Path(raw))
            plan = scan_plan(4); plan['batch_bytes'] = 97
            self.assertEqual(run_scan(path, 4, tool=self.tool, plan=plan)['counts'], expected)

    def test_잘린BGZF본문과정상압축된잘린레코드거부(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw); path, _ = self.fixture(root)
            broken = root / 'broken.bam'; broken.write_bytes(path.read_bytes()[:-100])
            with self.assertRaises(RuntimeError): run_scan(broken, 4, tool=self.tool)
            with pysam.BGZFile(str(path), 'rb') as handle: decoded = handle.read()
            with pysam.BGZFile(str(broken), 'wb') as handle: handle.write(decoded[:-7])
            with self.assertRaisesRegex(RuntimeError, '레코드'): run_scan(broken, 4, tool=self.tool)

    def test_손상된aux를정상카운트로확정하지않음(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw); path, _ = self.fixture(root)
            with pysam.BGZFile(str(path), 'rb') as handle: decoded = bytearray(handle.read())
            index = decoded.index(b'XA' + b'B' + b'i')
            decoded[index + 2] = ord('?')
            with pysam.BGZFile(str(path), 'wb') as handle: handle.write(bytes(decoded))
            with self.assertRaisesRegex(RuntimeError, 'aux'): run_scan(path, 4, tool=self.tool)

    def test_거대레코드는사전RAM상한오류로구분(self):
        with tempfile.TemporaryDirectory() as raw:
            path, _ = self.fixture(Path(raw)); plan = scan_plan(4); plan['max_record_bytes'] = 64
            with self.assertRaises(RecordSizeLimit): run_scan(path, 4, tool=self.tool, plan=plan)

    def test_전수집계중입력변경거부(self):
        with tempfile.TemporaryDirectory() as raw:
            path, expected = self.fixture(Path(raw))
            def changed(*args, **kwargs):
                value = path.stat(); os.utime(path, ns=(value.st_atime_ns, value.st_mtime_ns + 1))
                return {'counts': expected, 'execution': {}}
            with patch.dict(os.environ, {'MFM_RECORD_STATS_SCAN_MODE': 'single_pass'}), \
                 patch('integrate_wgs_rna.record_stats_parallel.run_scan', side_effect=changed):
                with self.assertRaisesRegex(RuntimeError, '변경'): stats.native_record_stats(path)

    def test_컴파일도구없는환경은auto에서만기존검증경로(self):
        with tempfile.TemporaryDirectory() as raw:
            path, expected = self.fixture(Path(raw))
            with patch('integrate_wgs_rna.record_stats_parallel.native_tool', side_effect=NativeUnavailable('검사')), \
                 patch.dict(os.environ, {'MFM_RECORD_STATS_SCAN_MODE': 'auto'}):
                self.assertEqual(stats.native_record_stats(path), expected)
            with patch('integrate_wgs_rna.record_stats_parallel.native_tool', side_effect=NativeUnavailable('검사')), \
                 patch.dict(os.environ, {'MFM_RECORD_STATS_SCAN_MODE': 'single_pass'}):
                with self.assertRaises(NativeUnavailable): stats.native_record_stats(path)

    def test_현재leaseCPU보다많은스레드로시작하지않음(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw); path, expected = self.fixture(root)
            (root / 'leases').mkdir(); (root / 'waiters').mkdir()
            lease = root / 'leases/self.json'
            lease.write_text(json.dumps(dict(manager_pid=os.getpid(), child_pid=None, cpus=2,
                total_cpus=8, memory_bytes=4 * 1024**3)))
            with patch.dict(os.environ, {'MFM_CPU_POOL_LEASE_FILE': str(lease), 'MFM_CPU_BUDGET': '8',
                 'MFM_CPU_POOL_GRANTED': '8', 'MFM_RECORD_STATS_SCAN_MODE': 'single_pass',
                 'MFM_RAM_LIMIT_BYTES': '1650000000000'}):
                self.assertEqual(stats.native_record_stats(path), expected)
                self.assertLessEqual(json.loads(lease.read_text())['cpus'], 2)

    def test_모든CPU예산에서스레드와RAM계획상한(self):
        for cpus in range(1, 69):
            plan = scan_plan(cpus)
            self.assertLessEqual(plan['native_thread_limit'], cpus)
            self.assertLess(plan['estimated_peak_rss_bytes'], 8 * 1024**3)


if __name__ == '__main__': unittest.main()
