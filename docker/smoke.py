#!/usr/bin/env python3
"""이미지의 native 빌드와 소형 WGS/RNA 계산을 검증함. donor는 받지 않음."""
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
import io
import json
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'src'), str(ROOT / 'vendor/MFM-DNARNAFM'), str(ROOT / 'tests')]


def main():
    from integrate_wgs_rna.record_stats_parallel import native_tool
    names = ['record_stats', 'gpu_fastq', 'restore_tags', 'trace_hash', 'rna_depth', 'rna_fastq',
             'rna_ase_partition', 'rna_primary', 'rna_variant_support', 'rna_cnv_sum',
             'rna_ase_read', 'pileup_oracle', 'pileup_counts']
    built = {name: str(native_tool(name)) for name in names}
    suite = unittest.TestLoader().loadTestsFromNames([
        'test_pileup_native',
        'test_rna_native.NativeTests.test_fastq_copy_counts_gzip_and_corruption',
        'test_rna_window_package.RnaWindowPackageTest.test_legacy_coordinates_produce_identical_depth',
    ])
    stream = io.StringIO()
    with redirect_stdout(stream), redirect_stderr(stream):
        result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    passed = result.wasSuccessful() and not result.skipped
    print(json.dumps({'status': 'pass' if passed else 'fail', 'native_tools': built,
                      'tests': result.testsRun, 'failures': len(result.failures),
                      'errors': len(result.errors), 'skipped': len(result.skipped),
                      'actual_donor_processing': False}, ensure_ascii=False, indent=2))
    if not passed:
        print(stream.getvalue(), file=sys.stderr)
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
