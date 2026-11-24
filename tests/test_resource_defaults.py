#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from unittest.mock import patch
from contextlib import ExitStack
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from handoff import resources  # noqa: E402
from handoff import cli  # noqa: E402


def base() -> dict:
    return {
        "resource_profile": "auto",
        "sample_cpus": "auto",
        "local_jobs": "auto",
        "download_jobs": "auto",
        "slurm": {"cpus": 18, "concurrency": 4},
    }


class ResourceDefaultsTest(unittest.TestCase):
    def setUp(self):
        probe = patch.object(resources, 'cgroup_limits', return_value=(None, None, None))
        probe.start(); self.addCleanup(probe.stop)

    def test_worker_launches_with_full_cpu_budget_even_when_input_resources_is_a_list(self):
        # 실제 worker의 입력 목록 변수 resources가 자원 정책 모듈을 가리지 않는지 검사한다.
        for executor, mode, expected in [('local', 'maximum', 72), ('local', 'bounded', 16),
                                         ('slurm', 'maximum', 12)]:
            with self.subTest(executor=executor, mode=mode), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                state = root / 'state'; state.mkdir()
                reference = root / 'reference.fa'; reference.write_text('>chr1\nA\n')
                refs = {key: str(reference) for key in ['ncbi', 'autosomes', 'broad', 'gtf', 'rna_extra']}
                (state / 'references.json').write_text(json.dumps(refs))
                config = {**base(), 'state_root': str(state), 'data_root': str(root / 'data'),
                          'output_root': str(root / 'output'), 'run_id': 'cpu-test',
                          'resource_profile': 'balanced', 'local_jobs': 2, 'sample_cpus': 8,
                          'cpu_utilization': mode, 'resource_resolution': {},
                          'remove_owned_raw_inputs': False}
                row = {'sample_id': 'S', 'data_group': 'wgs_only', 'platform': 'hifi', 'rna_path': ''}
                def launch(*args, **kwargs):
                    env = kwargs['env']
                    self.assertEqual(str(expected), env['MFM_CPU_BUDGET'])
                    self.assertEqual(mode, env['MFM_CPU_UTILIZATION'])
                with ExitStack() as stack:
                    environment = {'SLURM_JOB_ID': 'test', 'SLURM_CPUS_PER_TASK': '12'} if executor == 'slurm' else {}
                    stack.enter_context(patch.dict('os.environ', environment, clear=True))
                    stack.enter_context(patch.object(cli.os, 'sched_getaffinity', return_value=set(range(72))))
                    for name, value in [('claim', True), ('sample_resources', []), ('sample_downloads', []),
                                        ('sample_root', root / 'output/S'), ('prepare_rna', (config, [])),
                                        ('sample_environment', {}), ('processing_complete', True), ('model_bytes', 0),
                                        ('bind_manifest', (root / 'manifest.tsv',
                                         {'common_vcf_dir': str(root), 'source_decode_reference_fasta': str(reference)}))]:
                        stack.enter_context(patch.object(cli, name, return_value=value))
                    runner = stack.enter_context(patch.object(cli.subprocess, 'run', side_effect=launch))
                    self.assertEqual(0, cli.worker(config, row, executor))
                    runner.assert_called_once()

    def test_maximum_is_default_and_does_not_truncate_odd_cpu_count(self):
        result = resources.resolve(base(), cpu_count=73,
                                   available_memory_bytes=1024 * resources.GIB,
                                   free_storage_bytes=2 * resources.TIB)
        self.assertEqual('maximum', result['cpu_utilization'])
        self.assertEqual(73, result['resource_resolution']['effective_total_cpu_budget'])

    def test_fewer_workers_share_entire_cpu_and_explicit_bounded_is_preserved(self):
        config = {**base(), 'local_jobs': 2, 'sample_cpus': 8}
        result = resources.resolve(config, cpu_count=72,
                                   available_memory_bytes=1024 * resources.GIB,
                                   free_storage_bytes=2 * resources.TIB)
        self.assertEqual(72, resources.cpu_budget(result, affinity=72))
        limited = resources.resolve({**config, 'cpu_utilization': 'bounded'}, cpu_count=72,
                                    available_memory_bytes=1024 * resources.GIB,
                                    free_storage_bytes=2 * resources.TIB)
        self.assertEqual(16, resources.cpu_budget(limited, affinity=72))

    def test_slurm_maximum_cannot_exceed_actual_allocation(self):
        with patch.dict('os.environ', {'SLURM_JOB_ID': '123', 'SLURM_CPUS_PER_TASK': '12'}):
            self.assertEqual(12, resources.cpu_budget({**base(), 'sample_cpus': 18}, 'slurm', affinity=72))
            self.assertEqual(8, resources.cpu_budget({**base(), 'sample_cpus': 18}, 'slurm', affinity=8))

    def test_invalid_cpu_utilization_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'cpu_utilization'):
            resources.resolve({**base(), 'cpu_utilization': 'unknown'}, cpu_count=72,
                              available_memory_bytes=1024 * resources.GIB,
                              free_storage_bytes=2 * resources.TIB)

    def test_72_cores_resolve_to_four_by_eighteen(self) -> None:
        result = resources.resolve(base(), cpu_count=72,
                                   available_memory_bytes=1024 * resources.GIB,
                                   free_storage_bytes=2 * resources.TIB)
        self.assertEqual((result["local_jobs"], result["sample_cpus"]), (4, 18))
        self.assertEqual(result["download_jobs"], 4)
        self.assertEqual(result["resource_resolution"]["effective_total_cpu_budget"], 72)
        self.assertFalse(result["resource_resolution"]["elastic_enabled"])

    def test_ample_storage_enables_six_elastic_workers(self) -> None:
        result = resources.resolve(base(), cpu_count=72,
                                   available_memory_bytes=1024 * resources.GIB,
                                   free_storage_bytes=20 * resources.TIB)
        self.assertEqual((result["local_jobs"], result["sample_cpus"]), (6, 18))
        self.assertEqual(result["resource_resolution"]["effective_total_cpu_budget"], 72)
        self.assertEqual(result["resource_resolution"]["worker_cpu_overcommit_ratio"], 1.5)
        self.assertTrue(result["resource_resolution"]["elastic_enabled"])

    def test_32_cores_resolve_to_two_by_sixteen(self) -> None:
        result = resources.resolve(base(), cpu_count=32,
                                   available_memory_bytes=512 * resources.GIB,
                                   free_storage_bytes=20 * resources.TIB)
        self.assertEqual((result["local_jobs"], result["sample_cpus"]), (2, 16))

    def test_memory_and_profile_limit_parallel_samples(self) -> None:
        result = resources.resolve(base(), {"resource_profile": "conservative"},
                                   cpu_count=72, available_memory_bytes=64 * resources.GIB,
                                   free_storage_bytes=20 * resources.TIB)
        self.assertEqual(result["local_jobs"], 1)
        self.assertEqual(result["sample_cpus"], 72)

    def test_explicit_worker_slots_may_exceed_cpu_share_but_pool_does_not(self) -> None:
        result = resources.resolve(base(), {"local_jobs": 5, "sample_cpus": 18,
                                            "download_jobs": 3},
                                   cpu_count=72, available_memory_bytes=1024 * resources.GIB,
                                   free_storage_bytes=2 * resources.TIB)
        self.assertEqual((result["local_jobs"], result["sample_cpus"]), (5, 18))
        self.assertEqual(result["download_jobs"], 3)
        self.assertEqual(result["resource_resolution"]["effective_total_cpu_budget"], 72)
        self.assertTrue(result["resource_resolution"]["adjustments"])

    def test_throughput_profile_forces_six_with_storage_warning(self) -> None:
        result = resources.resolve(base(), {"resource_profile": "throughput"},
                                   cpu_count=72, available_memory_bytes=1024 * resources.GIB,
                                   free_storage_bytes=2 * resources.TIB)
        self.assertEqual((result["local_jobs"], result["sample_cpus"]), (6, 18))
        self.assertTrue(result["resource_resolution"]["elastic_enabled"])
        self.assertTrue(result["resource_resolution"]["adjustments"])

    def test_slurm_cli_overrides(self) -> None:
        result = resources.resolve(base(), {"slurm_cpus": 24, "slurm_concurrency": 3},
                                   cpu_count=72, available_memory_bytes=1024 * resources.GIB,
                                   free_storage_bytes=2 * resources.TIB)
        self.assertEqual(result["slurm"]["cpus"], 24)
        self.assertEqual(result["slurm"]["concurrency"], 3)

    def test_invalid_values_are_rejected(self) -> None:
        for override in ({"local_jobs": 0}, {"sample_cpus": -1}, {"download_jobs": 0},
                         {"elastic_max_jobs": 0}):
            with self.subTest(override=override), self.assertRaises(ValueError):
                resources.resolve(base(), override, cpu_count=72,
                                  available_memory_bytes=1024 * resources.GIB,
                                  free_storage_bytes=2 * resources.TIB)

    def test_effective_config_is_not_recomputed(self) -> None:
        first = resources.resolve(base(), cpu_count=72,
                                  available_memory_bytes=1024 * resources.GIB,
                                  free_storage_bytes=2 * resources.TIB)
        second = resources.resolve(first, cpu_count=8,
                                   available_memory_bytes=64 * resources.GIB,
                                   free_storage_bytes=20 * resources.TIB)
        self.assertEqual((second["local_jobs"], second["sample_cpus"]), (4, 18))

    def test_cpu_pool_summary_aggregates_stage_audit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            state = Path(tmp)
            audit = state / "stats" / "cpu_pool_usage.jsonl"
            audit.parent.mkdir(parents=True)
            rows = [
                {"branch": "WGS", "stage": "10", "granted_cpus": 18,
                 "maximum_cpus": 18, "queue_wait_sec": 0.5, "elapsed_sec": 10},
                {"branch": "RNA", "stage": "14", "granted_cpus": 12,
                 "maximum_cpus": 18, "queue_wait_sec": 2.0, "elapsed_sec": 20},
            ]
            audit.write_text("\n".join(json.dumps(row) for row in rows) + "\n")
            summary = cli.cpu_pool_summary(
                state, {"resource_resolution": {"elastic_enabled": True}})
            self.assertEqual(summary["event_count"], 2)
            self.assertEqual(summary["average_granted_cpus"], 15)
            self.assertEqual(summary["max_granted_cpus"], 18)
            self.assertEqual(summary["queue_wait_over_1s_fraction"], 0.5)
            self.assertEqual(summary["below_stage_cap_fraction"], 0.5)
            self.assertEqual(summary["leased_cpu_seconds"], 420)
            self.assertEqual(summary["stage_event_counts"],
                             {"RNA:14": 1, "WGS:10": 1})


if __name__ == "__main__":
    unittest.main()
