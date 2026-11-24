"""실제 lease 파일의 단계 전환·대기·RAM 제한과 반환을 검증함."""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from integrate_wgs_rna.memory_policy import GIB
from integrate_wgs_rna.realign.config import PersonalRealignConfig
from integrate_wgs_rna.realign.phases import RealignPhases, _atomic, reference_bases


class PhaseTests(unittest.TestCase):
    def setUp(self):
        quota=patch('integrate_wgs_rna.pool_leases.cgroup_cpu_limit',return_value=None)
        quota.start();self.addCleanup(quota.stop)

    def fixture(self, directory, cpus=72):
        root = Path(directory)
        (root / 'leases').mkdir()
        (root / 'waiters').mkdir()
        path = root / 'leases/owner.json'
        _atomic(path, {'cpus': cpus, 'total_cpus': 72, 'memory_bytes': 192 * GIB,
                      'manager_pid': os.getpid(), 'child_pid': None})
        return root, path

    def test_index_releases_70_cpus_for_real_pool_waiter(self):
        with tempfile.TemporaryDirectory() as directory:
            root, path = self.fixture(directory)
            with patch.dict(os.environ, {'MFM_CPU_POOL_LEASE_FILE': str(path)}, clear=True):
                phases = RealignPhases()
                self.assertEqual(phases.index_jobs({1: 3_000_000_000, 2: 3_000_000_000}, 'auto'), 2)
                _atomic(root / 'waiters/next.json', {'pid': os.getpid(), 'queued_at': time.time() - 1,
                    'min_cpus': 8, 'max_cpus': 72, 'memory_bytes': 192 * GIB,
                    'branch': 'WGS', 'stage': '10', 'sample': 'NEXT'})
                spec = importlib.util.spec_from_file_location('phase_pool', ROOT / 'scripts/cpu_pool.py')
                pool = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(pool)
                pool.allocate(root, 72)
                self.assertEqual(json.loads((root / 'leases/next.json').read_text())['cpus'], 70)
                self.assertEqual(json.loads(path.read_text())['cpus'], 2)
                phases.release()
                self.assertEqual(json.loads(path.read_text())['cpus'], 1)

    def test_ram_allows_one_index_and_never_reserves_two(self):
        with tempfile.TemporaryDirectory() as directory:
            root, path = self.fixture(directory)
            with patch.dict(os.environ, {'MFM_CPU_POOL_LEASE_FILE': str(path),
                                        'MFM_RAM_LIMIT_BYTES': str(200 * GIB)}, clear=True):
                phases = RealignPhases()
                self.assertEqual(phases.index_jobs({1: 3_000_000_000, 2: 3_000_000_000}, 'auto'), 1)
                row = json.loads(path.read_text())
                self.assertLessEqual(row['memory_bytes'], 180 * GIB)
                self.assertEqual(row['cpus'], 1)

    def test_memory_wait_returns_old_ram_and_cpu(self):
        with tempfile.TemporaryDirectory() as directory:
            root, path = self.fixture(directory)
            other = root / 'leases/other.json'
            _atomic(other, {'cpus': 8, 'memory_bytes': 240 * GIB})
            errors = []
            with patch.dict(os.environ, {'MFM_CPU_POOL_LEASE_FILE': str(path),
                                        'MFM_RAM_LIMIT_BYTES': str(300 * GIB)}, clear=True):
                phases = RealignPhases()
                def run():
                    try:
                        phases.reserve('large', max_units=16, min_units=8, base_memory=100 * GIB)
                    except BaseException as error:
                        errors.append(error)
                thread = threading.Thread(target=run)
                thread.start()
                try:
                    deadline = time.monotonic() + 3
                    while time.monotonic() < deadline:
                        row = json.loads(path.read_text())
                        if row.get('realign_phase') == 'large:waiting':
                            break
                        time.sleep(.01)
                    self.assertEqual(row['cpus'], 1)
                    self.assertEqual(row['memory_bytes'], 4 * GIB)
                    self.assertTrue(thread.is_alive())
                finally:
                    other.unlink()
                    thread.join(3)
                self.assertFalse(thread.is_alive())
                self.assertFalse(errors)
                self.assertEqual(json.loads(path.read_text())['memory_bytes'], 100 * GIB)

    def test_alignment_preserves_bwa_threads_and_counts_auxiliary_processes(self):
        with tempfile.TemporaryDirectory() as directory:
            root, path = self.fixture(directory)
            _atomic(root / 'leases/other.json', {'cpus': 2, 'memory_bytes': 4 * GIB})
            config = PersonalRealignConfig(cohort='test', output_root=root, sample_list=None,
                samples=['S'], chroms=['chr1'], sequence_variant_policy='snv_indel',
                dataset_role='pre_training', threads=32, sort_threads=4, hap_jobs=2)
            with patch.dict(os.environ, {'MFM_CPU_POOL_LEASE_FILE': str(path)}, clear=True):
                phases = RealignPhases()
                actual = phases.align_config(config, {1: 3_000_000_000, 2: 3_000_000_000})
                self.assertEqual(actual.threads, 32)
                self.assertEqual(actual.sort_threads, 0)
                self.assertEqual(actual.fastq_threads, 0)
                row = json.loads(path.read_text())
                self.assertEqual(row['cpus'], 70)
                self.assertEqual(2 * (actual.threads + actual.sort_threads + 3), row['cpus'])
                phases.release()

    def test_impossible_minimum_fails_before_commands(self):
        with patch.dict(os.environ, {'MFM_CPU_BUDGET': '4'}, clear=True):
            with self.assertRaisesRegex(ValueError, '최소'):
                RealignPhases().reserve('oversized', max_units=8, min_units=8)

    def test_slurm_smaller_cpu_and_ram_take_precedence(self):
        with patch.dict(os.environ, {'MFM_CPU_BUDGET': '72', 'SLURM_CPUS_PER_TASK': '8',
                                    'SLURM_MEM_PER_NODE': str(32 * 1024)}, clear=True):
            phases = RealignPhases()
            self.assertEqual(phases.total, 8)
            self.assertEqual(phases.ram_limit, 32 * GIB)
            with self.assertRaisesRegex(ValueError, 'RAM'):
                phases.reserve('too_large', max_units=1, memory_per_unit=32 * GIB)

    def test_index_address_space_limit_is_inherited_by_tool(self):
        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / 'resource.json'
            limit = 128 * 1024 ** 2
            code = 'import resource; assert resource.getrlimit(resource.RLIMIT_AS)[0]==' + str(limit)
            completed = subprocess.run([sys.executable, str(ROOT / 'scripts/resource_runner.py'),
                '--output', str(report), '--address-space-bytes', str(limit), '--',
                sys.executable, '-c', code], check=False)
            self.assertEqual(completed.returncode, 0)
            self.assertEqual(json.loads(report.read_text())['address_space_limit_bytes'], limit)

    def test_fai_uses_sequence_length_and_missing_fai_is_conservative(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / 'chr1.fa'
            fasta.write_text('>chr1\nAAAA\n')
            self.assertEqual(reference_bases({1: root}, ['chr1']), {1: fasta.stat().st_size})
            Path(str(fasta) + '.fai').write_text('chr1\t4\t6\t4\t5\n')
            self.assertEqual(reference_bases({1: root}, ['chr1']), {1: 4})


if __name__ == '__main__':
    unittest.main()
