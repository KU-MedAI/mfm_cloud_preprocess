"""pair 계산 우선과 다중 선행 전송·저장공간 보호를 함께 검사한다."""
from pathlib import Path
from types import SimpleNamespace
import copy, json, tempfile, unittest
from unittest.mock import patch
from test_ready_queue import Harness
from handoff.download_budget import storage_admit, memory_plan
from handoff.local_queue import execution_contract


class PrefetchTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)

    def harness(self, names):
        h = Harness(self.temp.name, names)
        for row in h.rows:
            row.update(read_type='long_read', data_group='wgs', known_input_total_bytes='100')
        h.queue.PRIORITY_POLICY = dict(mode='long_read_only', pair_first=True,
            blocked_samples={n: 'pair 계산 대기' for n in names if n != 'pair'},
            download_allowed_block_reasons=['pair 계산 대기'])
        h.queue.DOWNLOAD_TARGET = 4
        h.queue.PREFETCH_LIMIT = 6
        return h

    def test_계산보류샘플도4개까지다운로드하고pair계산만진입(self):
        h = self.harness(['wgs1', 'pair', 'wgs2', 'wgs3', 'wgs4'])
        h.rows[1]['data_group'] = 'wgs_rna_pair'
        s = h.scheduler(); s.step()
        self.assertEqual(len(s.downloads), 4)
        self.assertIn('pair', s.downloads)
        for name, (_, proc) in s.downloads.items():
            h.success(name); proc.result = 0
        s.step()
        self.assertEqual(set(s.active), {'pair'})
        self.assertFalse(s.deferred(s.entries[0], worker=False))
        self.assertTrue(s.deferred(s.entries[0]))

    def test_다른보류사유는다운로드예외로우회하지않음(self):
        h = self.harness(['wgs1'])
        h.queue.PRIORITY_POLICY['blocked_samples']['wgs1'] = '참조 불일치'
        s = h.scheduler(); s.step()
        self.assertFalse(h.started)

    def test_이미확보한WGS버퍼가미준비pair를막지않음(self):
        h = self.harness(['wgs1', 'wgs2', 'pair'])
        h.rows[2]['data_group'] = 'wgs_rna_pair'
        h.queue.PREFETCH_LIMIT = 2
        h.success('wgs1'); h.success('wgs2')
        s = h.scheduler(); s.step()
        self.assertIn('pair', s.downloads)
        self.assertFalse(s.active)

    def test_전송제어변경은계산계약과입력변경을구분(self):
        before = dict(sample_cpus=68, local_jobs=6, download_jobs=1, pipeline_env={'MFM_WGS_BACKEND':'cpu'})
        after = copy.deepcopy(before)
        after.update(download_jobs=4, prefetch_samples=6, prefetch_max_bytes=1000)
        self.assertEqual(execution_contract(before), execution_contract(after))
        after['pipeline_env']['MFM_WGS_BACKEND'] = 'parabricks'
        self.assertNotEqual(execution_contract(before), execution_contract(after))


class StorageTests(unittest.TestCase):
    def test_기본전송1GB와다중범위버퍼증가를동일계산(self):
        self.assertEqual(memory_plan({})['reservation_bytes'], 1_000_000_000)
        larger = memory_plan({'http_range_workers':16})
        self.assertGreater(larger['reservation_bytes'], 2_000_000_000)
        self.assertGreaterEqual(larger['reservation_bytes'], larger['estimated_peak_bytes'])
        with self.assertRaises(ValueError):memory_plan({'http_range_workers':17})

    def test_전송총량과디스크여유를각각보호하고중복합산방지(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            (root/'policy.json').write_text(json.dumps({'download_min_free_bytes':100}))
            config = dict(data_root=raw, prefetch_max_bytes=1000)
            a = dict(sample_id='a', known_input_total_bytes='400')
            b = dict(sample_id='b', known_input_total_bytes='400')
            with patch('handoff.download_budget.os.statvfs', return_value=SimpleNamespace(f_bavail=900, f_frsize=1)):
                self.assertTrue(storage_admit(config, root, [a], b)[0])
                self.assertTrue(storage_admit(config, root, [a], a)[0])
                self.assertFalse(storage_admit(config, root, [a], dict(sample_id='big', known_input_total_bytes='700'))[0])
            with patch('handoff.download_budget.os.statvfs', return_value=SimpleNamespace(f_bavail=899, f_frsize=1)):
                self.assertFalse(storage_admit(config, root, [a], b)[0])
            self.assertFalse(storage_admit(config, root, [], dict(sample_id='unknown'))[0])
            self.assertTrue(storage_admit({}, root, [], dict(sample_id='legacy'))[0])


if __name__ == '__main__':
    unittest.main()
