"""작은 실제 BAM으로 준비 캐시 재사용·무효화·동시 요청·실패를 검증함."""
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path
import json
import os
import random
import sys
import tempfile
import unittest
from unittest.mock import patch
import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'vendor/MFM-DNARNAFM')]
from integrate_wgs_rna import personal_realign as module
from integrate_wgs_rna.realign import prepared_cache as cache


class PreparedCacheTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'source.bam'
        self.work = self.root / 'work'
        self.work.mkdir()
        sequence = ''.join(random.Random(17).choices('ACGT', k=2000))
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'chr1', 'LN': 2000}]}
        with pysam.AlignmentFile(self.source, 'wb', header=header) as handle:
            for i, (name, flag, pos, mate) in enumerate([
                ('dup', 99 | 1024, 100, 300), ('dup', 147 | 1024, 300, 100),
                ('normal', 99, 500, 700), ('normal', 147, 700, 500)]):
                r = pysam.AlignedSegment()
                r.query_name, r.flag, r.reference_id, r.reference_start = name, flag, 0, pos
                r.query_sequence = sequence[pos:pos + 100]
                r.query_qualities = pysam.qualitystring_to_array('I' * 100)
                r.mapping_quality = 60
                r.cigar = [(0, 100)]
                r.next_reference_id, r.next_reference_start = 0, mate
                r.set_tag('HP', 1); r.set_tag('PS', 100)
                handle.write(r)
        pysam.index(str(self.source))
        self.config = module.PersonalRealignConfig(
            cohort='test', output_root=self.root, sample_list=None, samples=['S'],
            chroms=['chr1'], sequence_variant_policy='snv_indel', dataset_role='pre_training',
            platform='illumina', hap_read_mode='all_to_both', read_routing_policy='all_to_both',
            threads=1, sort_threads=1, route_jobs=1, sort_mem='32M')
        self.real = module._prepare_illumina_paired_input_uncached

    def prepare(self, config=None):
        return module.prepare_illumina_paired_input(
            config or self.config, 'S', 1, self.source, ['chr1'], self.work, None, 'all_to_both')

    def test_실제BAM은한번만생성하고재개때재사용(self):
        with patch.object(module, '_prepare_illumina_paired_input_uncached', wraps=self.real) as build:
            first = self.prepare()
            before = cache.fingerprint(first.path)
            second = self.prepare(replace(self.config, threads=3, route_jobs=2))
            self.assertEqual(build.call_count, 1)
        self.assertEqual(cache.fingerprint(second.path), before)
        self.assertEqual(first.input_counts, second.input_counts)
        self.assertEqual(second.stats['prepared_cache']['status'], 'reused')
        self.assertEqual(second.stats['input_metrics_execution']['mode'], 'cached')
        self.assertEqual(second.stats['alignment_input_primary'], 4)

    def test_입력지문이바뀌면다시생성(self):
        with patch.object(module, '_prepare_illumina_paired_input_uncached', wraps=self.real) as build:
            self.prepare(); s = self.source.stat()
            os.utime(self.source, ns=(s.st_atime_ns, s.st_mtime_ns + 1000000))
            result = self.prepare()
        self.assertEqual(build.call_count, 2)
        self.assertEqual(result.stats['prepared_cache']['status'], 'created')

    def test_정책이바뀌면다시생성(self):
        with patch.object(module, '_prepare_illumina_paired_input_uncached', wraps=self.real) as build:
            self.prepare(); self.prepare(replace(self.config, include_unmapped=False))
        self.assertEqual(build.call_count, 2)

    def test_출력지문이나증명손상은재사용하지않음(self):
        with patch.object(module, '_prepare_illumina_paired_input_uncached', wraps=self.real) as build:
            first = self.prepare(); first.path.touch()
            self.prepare(); cache.cache_path(first.path).write_text('{broken')
            self.prepare()
        self.assertEqual(build.call_count, 3)

    def test_기존무증명BAM은다시검증생성(self):
        first = self.prepare(); cache.discard(first.path)
        with patch.object(module, '_prepare_illumina_paired_input_uncached', wraps=self.real) as build:
            self.prepare()
        build.assert_called_once()

    def test_동시요청도한번만생성(self):
        with patch.object(module, '_prepare_illumina_paired_input_uncached', wraps=self.real) as build:
            with ThreadPoolExecutor(max_workers=2) as pool:
                results = list(pool.map(lambda _: self.prepare(), range(2)))
        self.assertEqual(build.call_count, 1)
        self.assertEqual({r.stats['prepared_cache']['status'] for r in results}, {'created', 'reused'})

    def test_실패한부분출력은성공증명이없고재개때다시생성(self):
        output = module.illumina_prepared_path(self.work, 'S', 'all_to_both')
        def broken(*args):
            output.write_bytes(b'partial')
            raise RuntimeError('의도한 실패')
        with patch.object(module, '_prepare_illumina_paired_input_uncached', side_effect=broken):
            with self.assertRaisesRegex(RuntimeError, '의도한'):
                self.prepare()
        self.assertFalse(cache.cache_path(output).exists())
        self.assertEqual(self.prepare().stats['prepared_cache']['status'], 'created')

    def test_생성중원본변경은성공증명을남기지않음(self):
        def changed(*args):
            result = self.real(*args)
            s = self.source.stat(); os.utime(self.source, ns=(s.st_atime_ns, s.st_mtime_ns + 1000000))
            return result
        with patch.object(module, '_prepare_illumina_paired_input_uncached', side_effect=changed):
            with self.assertRaisesRegex(RuntimeError, '생성 도중'):
                self.prepare()
        output = module.illumina_prepared_path(self.work, 'S', 'all_to_both')
        self.assertFalse(cache.cache_path(output).exists())

    def test_불완전metrics증명은다시생성(self):
        first=self.prepare();p=cache.cache_path(first.path)
        d=json.loads(p.read_text());d['stats']={};p.write_text(json.dumps(d))
        with patch.object(module, '_prepare_illumina_paired_input_uncached', wraps=self.real) as build:
            self.prepare()
        build.assert_called_once()

    def test_캐시옵션오타는명시오류(self):
        with patch.dict(os.environ,MFM_PREPARED_BAM_CACHE='typo'):
            with self.assertRaisesRegex(ValueError,'0 또는1'):
                self.prepare()

    def test_명시적재생성과캐시비활성화존중(self):
        self.prepare()
        with patch.object(module, '_prepare_illumina_paired_input_uncached', wraps=self.real) as build:
            self.prepare(replace(self.config, overwrite=True))
            with patch.dict(os.environ, MFM_PREPARED_BAM_CACHE='0'):
                result = self.prepare()
        self.assertEqual(build.call_count, 2)
        self.assertNotIn('prepared_cache', result.stats)
        self.assertFalse(cache.cache_path(result.path).exists())

if __name__ == '__main__':
    unittest.main()
