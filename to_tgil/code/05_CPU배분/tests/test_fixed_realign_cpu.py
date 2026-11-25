"""고정 소유 CPU를 검증한 HAP 겹침과 RAM/명시 설정의 보존을 검사함."""
from pathlib import Path
from dataclasses import replace
import json
import os
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'src')]
from integrate_wgs_rna.memory_policy import GIB
from integrate_wgs_rna.sample_cpu import process_start
from integrate_wgs_rna.realign.align_resources import alignment_plan
from integrate_wgs_rna.realign.config import PersonalRealignConfig
from integrate_wgs_rna.realign.hap_resources import HaplotypeResources
from integrate_wgs_rna.realign.phases import RealignPhases, _atomic
from integrate_wgs_rna.realign.stream_resources import mapper_finished


class FixedRealignCpuTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / 'leases').mkdir()
        (self.root / 'waiters').mkdir()
        self.ids = list(range(16))
        self.owner = self.root / 'leases/sample_A.json'
        _atomic(self.owner, dict(sample='S', sample_cpu_owner=True, fixed_sample_token='A',
            cpus=16, cpu_ids=self.ids, memory_bytes=4*GIB, pool_dir=str(self.root),
            manager_pid=os.getpid(), manager_start=process_start(os.getpid())))
        self.env = dict(MFM_CPU_ALLOCATION_MODE='fixed_sample', MFM_FIXED_SAMPLE_TOKEN='A',
            MFM_CPU_POOL_DIR=str(self.root), MFM_CPU_BUDGET='16', MFM_RAM_LIMIT_BYTES=str(512*GIB),
            MFM_WGS_ALIGN_REQUESTED='auto', MFM_WGS_SORT_REQUESTED='auto', MFM_WGS_HAP_REQUESTED='auto')
        self.config = PersonalRealignConfig(cohort='T', output_root=self.root, sample_list=None,
            samples=['S'], chroms=['chr1'], sequence_variant_policy='snv_indel',
            dataset_role='pre_training', platform='hifi', threads=2, sort_threads=1,
            hap_jobs=2, sort_mem='2G')
        self.lengths = {1: 3_000_000_000, 2: 3_000_000_000}
        mask = patch('os.sched_getaffinity', return_value=set(self.ids))
        mask.start(); self.addCleanup(mask.stop)

    def plan(self, env=None, ram=900*GIB, cpus=16, config=None, lengths=None):
        return alignment_plan(config or self.config, lengths or self.lengths,
                              cpus, ram, self.env if env is None else env)

    def test_owned_affinity_overlaps_two_haps_with_one_physical_budget(self):
        plan = self.plan()
        self.assertEqual((plan['hap_jobs'], plan['threads'], plan['sort_threads'],
                          plan['input_decode_threads']), (2, 8, 1, 1))
        self.assertEqual(plan['reserved_cpus'], 16)
        self.assertEqual(plan['logical_pipeline_cpus'], 30)
        self.assertEqual(plan['cpu_affinity_ids'], self.ids)
        self.assertEqual(plan['memory_model']['active_haps'], 2)
        self.assertEqual(plan['memory_model']['reference_index_bytes'], 96_000_000_000)

    def test_missing_dead_or_unpinned_owner_uses_strict_thread_budget(self):
        for overrides in [dict(MFM_CPU_ALLOCATION_MODE='dynamic'), dict(MFM_FIXED_SAMPLE_TOKEN=''),
                          dict(MFM_FIXED_SAMPLE_TOKEN='missing')]:
            plan = self.plan(dict(self.env, **overrides))
            self.assertEqual(plan['hap_jobs'], 1)
            self.assertNotIn('cpu_execution_policy', plan)
        with patch('os.sched_getaffinity', return_value=set(range(32))):
            self.assertEqual(self.plan()['hap_jobs'], 1)
        with patch('os.sched_getaffinity', return_value=set(range(8))):
            self.assertEqual(self.plan()['hap_jobs'], 1)
        row=json.loads(self.owner.read_text()); row['manager_start'] = -1; _atomic(self.owner,row)
        self.assertEqual(self.plan()['hap_jobs'], 1)

    def test_explicit_settings_other_platforms_and_budgets_keep_existing_policy(self):
        cases = [dict(WGS_AUTO_THREADS='0'), dict(MFM_WGS_ALIGN_REQUESTED='6'),
                 dict(MFM_WGS_SORT_REQUESTED='2'), dict(MFM_WGS_HAP_REQUESTED='1'),
                 dict(MFM_HAP_PARALLEL_ALLOWED='0'), dict(MFM_WGS_CPU_BALANCE='mapping'),
                 dict(MFM_MINIMAP_DECODE_THREADS='0'), dict(MFM_WGS_BACKEND='parabricks')]
        for overrides in cases:
            plan = self.plan(dict(self.env, **overrides))
            self.assertTrue(plan is None or 'cpu_execution_policy' not in plan, overrides)
        for platform in ['illumina', 'ont']:
            plan = self.plan(config=replace(self.config, platform=platform))
            self.assertTrue(plan is None or 'cpu_execution_policy' not in plan)
        for cpus in [8, 12, 32, 68]:
            plan = self.plan(cpus=cpus)
            self.assertTrue(plan is None or 'cpu_execution_policy' not in plan)
        self.assertEqual(self.plan(lengths={1:3_000_000_000})['hap_jobs'],1)

    def test_ram_shortage_falls_back_to_single_hap_before_launch(self):
        old = self.plan(dict(self.env, MFM_CPU_ALLOCATION_MODE='dynamic'))
        plan = self.plan(ram=old['reserved_memory_bytes'])
        self.assertEqual(plan['hap_jobs'], 1)
        self.assertLessEqual(plan['reserved_memory_bytes'],old['reserved_memory_bytes'])
        self.assertIsNone(self.plan(ram=8*GIB))

    def test_mapper_tail_and_hap_completion_preserve_sibling_memory_and_parent(self):
        lease=self.root/'leases/stage.json';plan=self.plan()
        _atomic(lease,dict(plan,cpus=16,total_cpus=16,memory_bytes=plan['reserved_memory_bytes'],
            fixed_sample_token='A',manager_pid=os.getpid(),realign_phase='align_restore_sort_qc'))
        env=dict(self.env,MFM_CPU_POOL_LEASE_FILE=str(lease))
        # CPU 장부의 가상 RAM을 작은 CI 호스트의 실제 RAM과 분리함. 도구는 실행하지 않음.
        with (patch.dict(os.environ,env,clear=True),
              patch('integrate_wgs_rna.pool_leases.cgroup_cpu_limit',return_value=None),
              patch('integrate_wgs_rna.realign.phases.RealignPhases.ceiling',return_value=450*GIB),
              patch('integrate_wgs_rna.realign.hap_resources.available_memory',return_value=450*GIB)):
            phase=RealignPhases();group=HaplotypeResources(phase,self.config,[1,2])
            before=json.loads(lease.read_text());sibling=before['hap_resources']['2'].copy()
            mapper_finished(1,1,lease)
            after=json.loads(lease.read_text())
            self.assertEqual(after['hap_resources']['2'],sibling)
            self.assertEqual(after['memory_bytes'],before['memory_bytes'])
            self.assertLessEqual(after['cpus'],16)
            group.controller(1).reserve('output_index',max_units=16,base_memory=4*GIB)
            self.assertLessEqual(json.loads(lease.read_text())['cpus'],16)
            group.release(1);group.release(2)
        self.assertEqual(json.loads(self.owner.read_text())['cpus'],16)
        self.assertEqual(json.loads(lease.read_text())['cpus'],1)


if __name__=='__main__':unittest.main()
