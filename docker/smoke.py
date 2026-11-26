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
    names = ['record_stats', 'alignment_scan', 'gpu_fastq', 'restore_tags', 'trace_hash', 'rna_depth', 'rna_fastq',
             'rna_ase_partition', 'rna_primary', 'rna_variant_support', 'rna_cnv_sum',
             'rna_ase_read', 'pileup_oracle', 'pileup_counts', 'pileup_counts_complete', 'pileup_oracle_complete']
    built = {name: str(native_tool(name)) for name in names}
    suite = unittest.TestLoader().loadTestsFromNames([
        'test_pileup_native',
        'test_pileup_complete',
        'test_package_inventory',
        'test_file_inventory',
        'test_prepare_fast_input',
        'test_record_stats_parallel.RecordStatsTests.test_htslib_single_pass_cram_counts_and_real_coordinate_order',
        'test_real_markdup.MarkdupResourceTest.test_child_hap_budget_caps_parent_lease',
        'test_real_markdup.MarkdupResourceTest.test_index_and_independent_counts_parallel_match_serial',
        'test_cram_output',
        # 전체 단계 실행은 실제 작업 RAM 예약이 필요하므로 작은 계산/패키지만 검사함.
        'test_pileup_feature.CountsTest',
        'test_pileup_feature.PileupPackageTest.test_export_preserves_windows_insertion_sequence_trace_reader_resume_cleanup',
        'test_pileup_feature.PileupPackageTest.test_pileup_reader_keeps_rna_window_signal',
        'test_pileup_feature.PileupPackageTest.test_interrupted_publication_is_resumable',
        'test_rna_native.NativeTests.test_fastq_copy_counts_gzip_and_corruption',
        'test_rna_window_package.RnaWindowPackageTest.test_legacy_coordinates_produce_identical_depth',
        'test_align_resources.SortBufferPlans',
        'test_fixed_realign_cpu',
        'test_stage01_rna_defer',
        'test_nas_transfer',
        'test_fixed_sample_cpu.FixedSampleCpuTests.test_sort메모리_자동요청과명시값을_자식단계에도보존',
        'test_docker_package.DockerPackageTests.test_docker_config_preserves_biology_and_limits_fixed_cpu',
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
