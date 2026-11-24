"""내부 유휴 대여·RAM 상한·정렬 버퍼와 실제 sort 결과를 검증한다."""
import json
import os
import random
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from integrate_wgs_rna.memory_policy import GIB
from integrate_wgs_rna.realign.config import PersonalRealignConfig
from integrate_wgs_rna.realign.phases import RealignPhases, _atomic
from integrate_wgs_rna.realign.idle_resources import loan_capacity
from integrate_wgs_rna.realign.hap_resources import HaplotypeResources
from integrate_wgs_rna import pool_backfill


class IdleTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / 'leases').mkdir()
        (self.root / 'waiters').mkdir()
        self.path = self.root / 'leases/owner.json'
        self.row = dict(cpus=1, total_cpus=68, memory_bytes=4*GIB,
                        manager_pid=os.getpid(), branch='WGS', stage='10', sample='S')
        _atomic(self.path, self.row)
        self.policy = {'live_limits': True, 'cpu_pool_limit': 68, 'memory_limit_bytes': 1650000000000,
                       'realign_sort_memory_bytes': 16*GIB,
                       'cpu_backfill': {'enabled': True, 'realign_phases': True, 'nice': 0,
                                        'grant_multiplier': 2, 'max_claim_pools': 3,
                                        'max_loan_pools': 2, 'use_average_idle': True}}
        _atomic(self.root / 'policy.json', self.policy)
        self.watch = {'updated': time.time(), 'guard_action': 'continue',
                      'cpu_backfill': {'ready': True, 'quota': 68, 'observation_seconds': 10,
                                       'used_cpus': 22, 'peak_used_cpus': 25, 'idle_cpus': 42}}
        _atomic(self.root / 'watchdog.json', self.watch)
        self.current = dict(quota=68, cpu_cgroup='/test', memory_max=1717986918400,
                            memory_current=500_000_000_000)
        for target, value in [('integrate_wgs_rna.realign.idle_resources.cgroup_sample', self.current),
                              ('integrate_wgs_rna.realign.idle_resources.in_cpu_group', True),
                              ('integrate_wgs_rna.realign.idle_resources.os.getpriority', 0),
                              ('integrate_wgs_rna.pool_leases.cgroup_cpu_limit', 68)]:
            mock = patch(target, return_value=value)
            mock.start()
            self.addCleanup(mock.stop)
        self.env = patch.dict(os.environ, {'MFM_CPU_POOL_LEASE_FILE': str(self.path),
            'MFM_WGS_ALIGN_REQUESTED': 'auto', 'MFM_WGS_SORT_REQUESTED': 'auto',
            'MFM_WGS_HAP_REQUESTED': 'auto', 'MFM_CPU_UTILIZATION': 'maximum',
            'MFM_WGS_CPU_BALANCE':'legacy'}, clear=True)
        self.env.start()
        self.addCleanup(self.env.stop)
        self.config = PersonalRealignConfig(cohort='test', output_root=self.root, sample_list=None,
            samples=['S'], chroms=['chr1'], sequence_variant_policy='snv_indel',
            dataset_role='pre_training', platform='hifi', threads=4, sort_threads=1,
            hap_jobs=2, sort_mem='2G')

    def capacity(self, others=(), memory=500*GIB):
        return loan_capacity(self.root, 68, self.row, others, memory)

    def test_explicit_flag_stale_observation_guard_and_priority_block(self):
        self.policy['cpu_backfill']['realign_phases'] = False
        _atomic(self.root / 'policy.json', self.policy)
        self.assertIsNone(self.capacity())
        self.policy['cpu_backfill']['realign_phases'] = True
        _atomic(self.root / 'policy.json', self.policy)
        self.watch['updated'] = time.time()-4
        _atomic(self.root / 'watchdog.json', self.watch)
        self.assertIsNone(self.capacity())
        self.watch.update(updated=time.time(), guard_action='pause')
        _atomic(self.root / 'watchdog.json', self.watch)
        self.assertIsNone(self.capacity())
        self.watch['guard_action'] = 'continue'
        _atomic(self.root / 'watchdog.json', self.watch)
        self.policy['cpu_backfill']['nice'] = 19
        _atomic(self.root / 'policy.json', self.policy)
        self.assertIsNone(self.capacity())

    def test_explicit_four_claim_three_loan_uses_same_limits(self):
        others=[dict(cpus=68),dict(cpus=136,cpu_backfill={'cpus':136})]
        self.assertEqual(self.capacity(others)[0],0)
        self.policy['cpu_backfill'].update(max_claim_pools=4,max_loan_pools=3)
        _atomic(self.root/'policy.json',self.policy)
        self.assertEqual(self.capacity(others)[0],68)
        self.current.update(quota=64)
        self.assertIsNone(self.capacity(others))

    def test_quota_memory_and_both_logical_cpu_limits(self):
        self.assertEqual(self.capacity()[0], 68)
        self.assertEqual(self.capacity([dict(cpus=128, cpu_backfill={'cpus':128})])[0], 8)
        self.assertEqual(self.capacity([dict(cpus=200)])[0], 4)
        self.current.update(memory_max=20*GIB, memory_current=19*GIB)
        self.assertEqual(self.capacity()[1], 0)
        self.current['quota'] = 4
        self.assertIsNone(self.capacity())

    def test_concurrent_gpu_policy_requires_separate_explicit_opt_in(self):
        with patch.dict(os.environ,MFM_WGS_BACKEND='parabricks'):
            self.assertIsNone(self.capacity())
            self.policy['cpu_backfill']['gpu_realign_phases']=True
            _atomic(self.root/'policy.json',self.policy)
            capacity=self.capacity()
            self.assertEqual(capacity[0],68)
            self.assertEqual(capacity[2]['kind'],'gpu_realign_phase')

    def test_preprocessing_budget_reuse_ignores_other_project_load(self):
        self.policy['cpu_backfill'].update(
            max_extra_cpus=16, reservation_reuse={'enabled': True})
        _atomic(self.root/'policy.json', self.policy)
        self.watch['cpu_backfill'].update(used_cpus=68, peak_used_cpus=68, idle_cpus=0,
            reservation_reuse={'ready': True, 'sampled_at': time.monotonic(),
                               'peak_used_cpus': 20})
        _atomic(self.root/'watchdog.json', self.watch)
        result = self.capacity([dict(cpus=68)])
        self.assertEqual(result[0], 16)
        self.assertEqual(result[2]['allocation_basis'], 'reservation_reuse')
        self.assertEqual(self.capacity([dict(cpus=68),
            dict(cpus=16, cpu_backfill={'cpus':16})])[0], 0)

    def test_stale_preprocessing_observation_cannot_reuse_budget(self):
        self.policy['cpu_backfill']['reservation_reuse'] = {'enabled': True}
        _atomic(self.root/'policy.json', self.policy)
        self.watch['cpu_backfill'].update(used_cpus=68, peak_used_cpus=68, idle_cpus=0,
            reservation_reuse={'ready': True, 'sampled_at': time.monotonic()-4,
                               'peak_used_cpus': 0})
        _atomic(self.root/'watchdog.json', self.watch)
        self.assertEqual(self.capacity()[0], 0)
        with patch.dict(os.environ, MFM_FIXED_SAMPLE_TOKEN='fixed-owner'):
            self.assertIsNone(self.capacity())

    def test_shared_cooldown_blocks_another_phase(self):
        (self.root / 'cpu_backfill_events.jsonl').write_text('{}\n')
        self.assertIsNone(self.capacity())
        stamp = time.time()-11
        os.utime(self.root / 'cpu_backfill_events.jsonl', (stamp, stamp))
        self.assertIsNotNone(self.capacity())

    def test_only_declared_dynamic_four_gib_entry_is_eligible(self):
        plan={'adaptive_realign':True,'stage_estimates_bytes':{'personal_realign_entry':4*GIB}}
        row={'pid':os.getpid(),'branch':'WGS','stage':'10','layout':'original','memory_bytes':4*GIB}
        def check(value,extra=''):
            data=('MFM_LONGREAD_MEMORY_PLAN='+json.dumps(value)+'\0'+extra).encode()
            with patch.object(Path,'read_bytes',return_value=data):
                return pool_backfill.eligible_waiter(row)
        self.assertTrue(check(plan))
        self.assertFalse(check({}))
        self.assertFalse(check([]))
        self.assertFalse(check(plan,'MFM_REALIGN_PHASE_RESOURCES=0\0'))
        self.assertFalse(check(plan,'MFM_WGS_BACKEND=parabricks\0'))
        row['memory_bytes']=128*GIB
        self.assertFalse(check(plan))

    def test_dynamic_entry_borrows_only_two_cpus_and_reacquires_peak_later(self):
        row={'pid':os.getpid(),'branch':'WGS','stage':'10','layout':'original','sample':'ENTRY',
             'queued_at':time.time()-1,'min_cpus':2,'max_cpus':68,'memory_bytes':4*GIB}
        _atomic(self.root/'waiters/entry.json',row)
        observed=dict(self.watch['cpu_backfill'],**self.current)
        plan={'adaptive_realign':True,'stage_estimates_bytes':{'personal_realign_entry':4*GIB}}
        raw=('MFM_LONGREAD_MEMORY_PLAN='+json.dumps(plan)+'\0').encode()
        with patch.object(Path,'read_bytes',return_value=raw),patch.object(pool_backfill,'in_cpu_group',return_value=True),patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=500*GIB),patch.object(os,'setpriority'),patch.object(os,'getpriority',return_value=0):
            result=pool_backfill.lend_idle_cpus(self.root,68,observed,pool_backfill.CpuObserver(),{'action':'continue','rss_bytes':0})
        self.assertEqual(result['launched'][0]['cpus'],2)
        entry=json.loads((self.root/'leases/entry.json').read_text())
        self.assertEqual(entry['memory_bytes'],4*GIB)
        self.assertEqual(entry['cpu_backfill']['kind'],'realign_entry')

    def test_waiting_alignment_borrows_and_reserves_larger_sort_before_start(self):
        _atomic(self.root / 'leases/other.json', dict(cpus=64, memory_bytes=4*GIB))
        phases = RealignPhases()
        config = phases.align_config(self.config, {1:3_000_000_000, 2:3_000_000_000})
        row = json.loads(self.path.read_text())
        self.assertEqual(row['cpu_backfill']['kind'], 'realign_phase')
        self.assertLessEqual(row['cpus']+64, 204)
        self.assertLessEqual(row['cpus'], 68)
        self.assertEqual(int(config.sort_mem), 16*GIB)
        expected = 16*GIB*max(1,config.sort_threads)*config.hap_jobs
        self.assertEqual(row['effective_alignment']['memory_model']['sort_buffers_bytes'], expected)
        self.assertLessEqual(row['memory_bytes']+4*GIB, 1485000000000)
        self.assertIsNone(self.capacity())
        phases.release()
        self.assertNotIn('cpu_backfill', json.loads(self.path.read_text()))

    def test_generic_boundary_loan_and_release(self):
        _atomic(self.root / 'leases/other.json', dict(cpus=68, memory_bytes=4*GIB))
        phases = RealignPhases()
        self.assertEqual(phases.reserve('output_qc', max_units=23, min_units=2,
                                      memory_per_unit=2*GIB), 23)
        self.assertEqual(json.loads(self.path.read_text())['memory_bytes'], 50*GIB)
        (self.root / 'leases/other.json').unlink()
        phases.release()
        self.assertNotIn('cpu_backfill', json.loads(self.path.read_text()))

    def test_partial_normal_grant_grows_with_measured_idle(self):
        _atomic(self.root/'leases/other.json',dict(cpus=64,memory_bytes=4*GIB))
        phases=RealignPhases()
        self.assertEqual(phases.reserve('input_counts',max_units=68,memory_per_unit=2*GIB),68)
        row=json.loads(self.path.read_text())
        self.assertEqual(row['cpu_backfill']['cpus'],68)
        self.assertEqual(row['memory_bytes'],140*GIB)

    def test_sort_falls_back_to_small_buffer_and_keeps_explicit_larger_value(self):
        phases = RealignPhases()
        lengths = {1:3_000_000_000}
        small = phases.alignment_plan(self.config, lengths, 9, 88*GIB)
        self.assertIsNotNone(small)
        self.assertEqual(small['sort_mem'], '2G')
        self.config.sort_mem = '32G'
        explicit = phases.alignment_plan(self.config, lengths, 9, 200*GIB)
        self.assertEqual(explicit['sort_mem'], '32G')

    def test_larger_buffer_never_serializes_haps_to_fit_memory(self):
        phases = RealignPhases()
        plan = phases.alignment_plan(self.config, {1:3_000_000_000,2:3_000_000_000},68,220*GIB)
        self.assertEqual(plan['hap_jobs'],2)
        from integrate_wgs_rna.realign.align_resources import parse_memory
        # 현재68CPU의 보조8개/HAP는16개 버퍼를 보유해4GiB가 RAM220GiB를 초과함.
        self.assertEqual(parse_memory(plan['sort_mem']),2*GIB)
        self.assertEqual(plan['reserved_cpus'],68)
        self.assertLessEqual(plan['reserved_memory_bytes'],220*GIB)

    def test_large_operating_buffer_steps_down_with_ram_and_preserves_two_haps(self):
        self.policy['realign_sort_memory_bytes']=64*GIB
        _atomic(self.root/'policy.json',self.policy)
        phases=RealignPhases()
        lengths={1:3_000_000_000,2:3_000_000_000}
        large=phases.alignment_plan(self.config,lengths,68,800*GIB)
        self.assertEqual(large['hap_jobs'],2)
        self.assertEqual(int(large['sort_mem']),32*GIB)
        self.assertEqual(large['reserved_cpus'],68)
        self.assertEqual(large['sort_threads'],8)
        self.assertLessEqual(large['reserved_memory_bytes'],800*GIB)
        smaller=phases.alignment_plan(self.config,lengths,68,320*GIB)
        self.assertEqual(smaller['hap_jobs'],2)
        self.assertEqual(int(smaller['sort_mem']),4*GIB)
        self.assertEqual(smaller['reserved_cpus'],68)
        # 한 단계 큰 버퍼는 같은68CPU/2HAP에서 RAM 예산을 넘는지 독립 계산함.
        from dataclasses import replace
        from integrate_wgs_rna.realign.align_memory import alignment_memory
        for current,ram in [(large,800*GIB),(smaller,320*GIB)]:
            bigger=replace(self.config,sort_mem=str(int(current['sort_mem'])*2))
            self.assertGreater(alignment_memory(bigger,lengths,2,current['threads'],current['sort_threads'],
                dict(os.environ,MFM_MINIMAP_DECODE_THREADS_EFFECTIVE=str(current['input_decode_threads'])))['reserved_memory_bytes'],ram)
        self.assertLessEqual(smaller['reserved_memory_bytes'],320*GIB)

    def test_phase_boundary_refreshes_live_memory_policy(self):
        self.policy['memory_limit_bytes']=500*GIB
        _atomic(self.root/'policy.json',self.policy)
        phases=RealignPhases()
        self.assertEqual(phases.ceiling(),450*GIB)
        self.policy['memory_limit_bytes']=800*GIB
        _atomic(self.root/'policy.json',self.policy)
        self.assertEqual(phases.ceiling(),720*GIB)

    def test_hap_wait_keeps_other_running_hap_loan_until_normal_capacity_returns(self):
        self.row.update(cpus=20, memory_bytes=40*GIB, cpu_backfill={'cpus':20})
        _atomic(self.path, self.row)
        other = self.root / 'leases/other.json'
        _atomic(other, dict(cpus=68, memory_bytes=4*GIB))
        phases = RealignPhases()
        haps = HaplotypeResources(phases, self.config, [1,2])
        errors = []
        def reserve():
            try:
                haps.reserve(1, 'output_index', 4, 1, GIB, 4*GIB)
            except BaseException as exc:
                errors.append(exc)
        with patch('integrate_wgs_rna.realign.hap_resources.loan_capacity', return_value=None):
            thread = threading.Thread(target=reserve)
            thread.start()
            try:
                deadline = time.monotonic()+3
                while time.monotonic()<deadline:
                    row = json.loads(self.path.read_text())
                    if row['hap_resources']['1']['phase'].endswith(':waiting'):
                        break
                    time.sleep(.01)
                self.assertTrue(thread.is_alive())
                self.assertIn('cpu_backfill', row)
                self.assertEqual(row['hap_resources']['2']['cpus'], 10)
                self.assertEqual(row['cpus'], 10)
                self.assertEqual(row['hap_resources']['1']['cpus'], 0)
            finally:
                other.unlink()
                thread.join(3)
        self.assertFalse(thread.is_alive())
        self.assertFalse(errors)
        self.assertNotIn('cpu_backfill', json.loads(self.path.read_text()))
        self.assertEqual(json.loads(self.path.read_text())['hap_resources']['2']['cpus'], 10)

    def test_actual_samtools_sort_spill_and_larger_buffer_preserve_records_tags_order(self):
        samtools = Path(sys.executable).parent / 'samtools'
        source = self.root / 'input.sam'
        randomizer = random.Random(17)
        with source.open('w') as handle:
            handle.write('@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n')
            for index in range(6000):
                position = randomizer.randrange(1,90000)
                handle.write(f'r{index}\t0\tchr1\t{position}\t60\t1000M\t*\t0\t0\t'+
                             'A'*1000+'\t'+'I'*1000+f'\tHP:i:{index%2+1}\tPS:i:12\n')
        outputs = []
        for memory in ['1M','16G']:
            output = self.root / (memory+'.bam')
            subprocess.run([str(samtools),'sort','-@','1','-m',memory,'-T',str(self.root/memory),
                            '-o',str(output),str(source)],check=True,capture_output=True)
            outputs.append(subprocess.run([str(samtools),'view',str(output)],check=True,
                                          capture_output=True).stdout)
        self.assertEqual(outputs[0],outputs[1])


if __name__ == '__main__':
    unittest.main()
