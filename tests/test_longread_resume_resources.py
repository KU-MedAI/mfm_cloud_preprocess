"""구실행 재개 auto 보존 및 긴 read 메모리 보호 회귀 검사함."""
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pysam
from test_hap_boundary import BoundaryTests
from integrate_wgs_rna.realign.boundary_upgrade import prepare_boundary
from integrate_wgs_rna.compact.memory_preflight import estimate_peak
from integrate_wgs_rna.memory_policy import GIB
from pre_training.preprocess_raw_readset.filters import rejection_reason
from pre_training.preprocess_raw_readset.config import ExportConfig

class LongreadResumeTests(unittest.TestCase):
    def test_legacy_missing_auto_expands_and_restores_environment(self):
        with tempfile.TemporaryDirectory() as raw:
            _, _, script, report, _, _, _, env = BoundaryTests().fixture(raw)
            env.pop('MFM_WGS_ALIGN_REQUESTED')
            env.pop('MFM_WGS_SORT_REQUESTED')
            with patch.dict(os.environ, env, clear=True):
                _, info = prepare_boundary(['bash', str(script)], report)
                self.assertGreater(info['threads'], 2)
                self.assertLessEqual(info['phase_plan']['reserved_cpus'] + 4, 16)
                self.assertNotIn('MFM_WGS_ALIGN_REQUESTED', os.environ)
                self.assertNotIn('MFM_WGS_SORT_REQUESTED', os.environ)
                self.assertNotIn('MFM_WGS_HAP_REQUESTED', os.environ)

    def test_legacy_explicit_thread_request_is_not_overridden(self):
        with tempfile.TemporaryDirectory() as raw:
            _, _, script, report, _, _, _, env = BoundaryTests().fixture(raw)
            before = script.read_bytes()
            env['MFM_WGS_ALIGN_REQUESTED'] = '2'
            with patch.dict(os.environ, env, clear=True):
                _, info = prepare_boundary(['bash', str(script)], report)
            self.assertIsNone(info)
            self.assertEqual(before, script.read_bytes())

    def test_long_reads_reserve_the_same_larger_worker_ceiling(self):
        short = dict(contig_lengths={'chr4': 190_000_000}, total_reads=9_000_000,
                     total_reference_bases=2_900_000_000, read_length=150)
        long = dict(short, read_length=36_000)
        a = estimate_peak(short, 'chr4', 2, 1_048_576)
        b = estimate_peak(long, 'chr4', 2, 1_048_576)
        self.assertEqual(a['worker_ceiling_bytes'], 8 * GIB)
        self.assertGreaterEqual(b['worker_ceiling_bytes'], 32 * GIB)
        self.assertGreaterEqual(b['worker_ceiling_bytes'], b['estimated_peak_bytes'])

    def test_sequence_filter_retains_policy_without_allocating_sequence_cache(self):
        class LengthOnlyRead:
            is_unmapped = is_secondary = is_supplementary = is_qcfail = is_duplicate = is_paired = False
            mapping_quality = 60
            query_length = 10
            @property
            def query_sequence(self):
                raise AssertionError('SEQ 문자열 캐시를 만들면 안 됨')
        config = ExportConfig()
        self.assertIsNone(rejection_reason(LengthOnlyRead(), config))
        read = LengthOnlyRead(); read.query_length = 0
        self.assertEqual(rejection_reason(read, config), 'missing_sequence')
        for sequence in ['ACGTN', None]:
            real = pysam.AlignedSegment(); real.query_name = 'r'
            real.query_sequence = sequence; real.flag = 0; real.mapping_quality = 60
            expected = None if sequence else 'missing_sequence'
            self.assertEqual(rejection_reason(real, config), expected)
