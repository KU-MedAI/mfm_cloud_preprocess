"""입력 RAM 예측, 전역 입장 제어와 제한 동작을 안전한 작은 입력으로 검증한다."""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from handoff.resources import cpu_budget, memory_budget
from integrate_wgs_rna import memory_policy as policy
from integrate_wgs_rna.cpu_executor import WorkerBudget
from integrate_wgs_rna.compact.memory_preflight import estimate_peak


class MemoryPolicyTest(unittest.TestCase):
    def test_global_limit_and_slurm_allocation_preserve_cpu_pool(self):
        with patch.dict(os.environ, {}, clear=True), patch('handoff.resources.cgroup_limits',return_value=(None,None,None)):
            self.assertEqual(1_800_000_000_000, memory_budget({}))
            self.assertEqual(1_800_000_000_000, memory_budget({'memory_limit_bytes': 2_000_000_000_000}))
            self.assertEqual(72, cpu_budget({'local_jobs': 3, 'sample_cpus': 18}, affinity=72))
            with patch.dict(os.environ, SLURM_MEM_PER_NODE='262144'):
                self.assertEqual(256 * policy.GIB, memory_budget({}, 'slurm'))
            with self.assertRaises(ValueError):
                memory_budget({'memory_limit_bytes': True})

    def test_memory_admission_and_release_respect_other_stage(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            (root / 'leases').mkdir(); (root / 'waiters').mkdir()
            own = root / 'leases/own.json'
            other = root / 'leases/other.json'
            own.write_text(json.dumps({'cpus': 72, 'total_cpus': 72}))
            other.write_text(json.dumps({'cpus': 4, 'memory_bytes': 400 * policy.GIB}))
            (root / 'policy.json').write_text(json.dumps({'memory_limit_bytes': 1_000_000_000_000}))
            with patch.dict(os.environ, MFM_CPU_POOL_LEASE_FILE=str(own)):
                budget = WorkerBudget(72, 72, True, memory_per_worker=20 * policy.GIB, cpus_per_worker=2)
                cap = budget.capacity(72)
                self.assertLess(cap, 36)
                self.assertEqual(cap * 2, json.loads(own.read_text())['cpus'])
                self.assertLessEqual(json.loads(own.read_text())['memory_bytes'] + 400 * policy.GIB, policy.usable_memory(root))
                self.assertEqual(0, budget.capacity(0, grow=False))
                self.assertEqual(4 * policy.GIB, json.loads(own.read_text())['memory_bytes'])

    def test_guard_only_signals_largest_tracked_stage(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw); (root / 'leases').mkdir()
            for sample, pid in [('A', 90001), ('B', 90002)]:
                (root / f'leases/{sample}.json').write_text(json.dumps({'sample': sample, 'stage': '12', 'child_pid': pid}))
            with patch.dict(os.environ, MFM_RAM_LIMIT_BYTES=str(100 * policy.GIB)), patch.object(policy, 'process_rss', return_value={'A.json': 30 * policy.GIB, 'B.json': 65 * policy.GIB}), patch.object(policy.os, 'killpg') as kill:
                result = policy.protect_pool_memory(root)
                self.assertEqual('B', result['sample'])
                self.assertEqual('stop_largest_stage', result['action'])
                kill.assert_called_once_with(90002, policy.signal.SIGTERM)
                policy.protect_pool_memory(root)
                self.assertEqual(1, kill.call_count)

    def test_rss_scan_survives_process_exit_during_proc_iteration(self):
        fields = ['S', '1'] + ['0'] * 19 + ['4']
        live_stat = '100 (worker) ' + ' '.join(fields)
        def read_stat(path, *args, **kwargs):
            if path.parent.name == '101':
                raise ProcessLookupError('process exited')
            return live_stat
        with patch.object(policy.os, 'listdir', return_value=['100', '101', 'self']):
            with patch.object(Path, 'read_text', autospec=True, side_effect=read_stat):
                with patch.object(Path, 'glob', side_effect=ProcessLookupError('lazy iteration race')):
                    self.assertEqual({'stage': 4 * os.sysconf('SC_PAGE_SIZE')},
                                     policy.process_rss({'stage': 100}))

    def test_worker_address_space_fails_allocation_before_host_oom(self):
        command = 'from integrate_wgs_rna.memory_policy import limit_worker_address_space\nlimit_worker_address_space(256 * 1024**2)\ntry:\n x=bytearray(384 * 1024**2)\nexcept MemoryError:\n print("bounded")\nelse:\n raise SystemExit(2)'
        result = subprocess.run([sys.executable, '-c', command], capture_output=True, text=True, timeout=15)
        self.assertEqual(0, result.returncode, result.stderr)
        self.assertEqual('bounded', result.stdout.strip())

    def test_nested_observations_count_each_pid_once(self):
        def stat(pid, parent, pages):
            return f'{pid} (worker) ' + ' '.join(['S',str(parent)]+['0']*19+[str(pages)])
        rows={'100':stat(100,1,2),'101':stat(101,100,3),'102':stat(102,101,5)}
        with patch.object(policy.os,'listdir',return_value=list(rows)), patch.object(Path,'read_text',autospec=True,side_effect=lambda p:rows[p.parent.name]):
            usage=policy.process_rss({'inner':101,'outer':100,'duplicate':100},unique_roots=True)
        self.assertEqual(sum(usage.values()),10*os.sysconf('SC_PAGE_SIZE'))
        self.assertEqual(max(usage.values()),10*os.sysconf('SC_PAGE_SIZE'))

    def test_unreserved_gpu_is_included_and_dead_observers_are_pruned(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir()
            (root/'leases/cpu.json').write_text(json.dumps({'sample':'CPU','stage':'10','child_pid':90001}))
            live=policy.register_memory_observer(root,{'sample':'GPU','stage':'10','manager_pid':90002,'child_pid':90003})
            dead=policy.register_memory_observer(root,{'sample':'DEAD','stage':'10','manager_pid':90004,'child_pid':90005})
            usage={'cpu.json':40*policy.GIB,'memory_observers/'+live.name:55*policy.GIB}
            with patch.dict(os.environ,MFM_RAM_LIMIT_BYTES=str(100*policy.GIB)), patch.object(policy,'lease_alive',side_effect=lambda row:row['sample']=='GPU'), patch.object(policy,'process_rss',return_value=usage) as rss, patch.object(policy.os,'killpg') as kill:
                result=policy.protect_pool_memory(root,force=True)
            self.assertEqual(result['sample'],'GPU');self.assertEqual(result['rss_bytes'],95*policy.GIB)
            self.assertEqual(result['trigger_bytes'],90*policy.GIB)
            kill.assert_called_once_with(90003,policy.signal.SIGTERM)
            self.assertEqual(set(rss.call_args.args[0]),set(usage))
            self.assertFalse(dead.exists());self.assertEqual(len(list((root/'leases').glob('*.json'))),1)

    def test_input_peak_reduced_by_batch_and_bounded_index(self):
        summary = {'contig_lengths': {'chr1': 249_000_000}, 'total_reads': 740_000_000,
                   'total_reference_bases': 2_900_000_000, 'read_length': 151}
        old = estimate_peak(summary, 'chr1', 32, 1_000_000, False)
        new = estimate_peak(summary, 'chr1', 2, 1_000_000, True)
        self.assertLess(new['estimated_peak_bytes'], old['estimated_peak_bytes'] // 5)
        self.assertGreaterEqual(new['worker_ceiling_bytes'], new['estimated_peak_bytes'])


if __name__ == '__main__':
    unittest.main()
