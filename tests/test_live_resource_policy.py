"""승인된1.5TB·실제할당 우선·들어갈 수 있는 대기만 CPU분배하는지 검증함."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch

from integrate_wgs_rna import memory_policy
from integrate_wgs_rna.pool_leases import fitting_waiters, pool_cpu_limit
from integrate_wgs_rna.realign.align_resources import alignment_plan
from integrate_wgs_rna.realign.config import PersonalRealignConfig
from handoff.local_queue import capacity_increase_only
from test_hap_boundary import BoundaryTests
from integrate_wgs_rna.realign.boundary_upgrade import prepare_boundary

ROOT = Path(__file__).resolve().parents[1]
AUTO = {'WGS_AUTO_THREADS':'1', 'MFM_WGS_ALIGN_REQUESTED':'auto',
        'MFM_WGS_SORT_REQUESTED':'auto', 'MFM_WGS_HAP_REQUESTED':'auto'}


class LiveResourceTests(unittest.TestCase):
    def setUp(self):
        quota=patch('integrate_wgs_rna.pool_leases.cgroup_cpu_limit',return_value=None)
        quota.start();self.addCleanup(quota.stop)

    def test_승인공유정책은구환경1TB를대체하고실제할당은우선(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'policy.json').write_text(json.dumps({'live_limits':True,'memory_limit_bytes':1_500_000_000_000,'cpu_pool_limit':72}))
            with patch.dict(os.environ,{'MFM_RAM_LIMIT_BYTES':'1000000000000'},clear=True), patch.object(memory_policy,'allocation_memory_limit',return_value=2_000_000_000_000):
                self.assertEqual(memory_policy.memory_limit(root),1_500_000_000_000)
                self.assertEqual(memory_policy.usable_memory(root),1_350_000_000_000)
                self.assertEqual(memory_policy.memory_limit(),1_000_000_000_000)
            with patch.object(memory_policy,'allocation_memory_limit',return_value=256*1024**3):
                self.assertEqual(memory_policy.memory_limit(root),256*1024**3)

    def test_명시승인없는기존상한은유지(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'policy.json').write_text(json.dumps({'memory_limit_bytes':1_500_000_000_000}))
            with patch.dict(os.environ,{'MFM_RAM_LIMIT_BYTES':'1000000000000'},clear=True):
                self.assertEqual(memory_policy.memory_limit(root),1_000_000_000_000)

    def test_CPU공유72와작은Slurm할당(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'policy.json').write_text(json.dumps({'live_limits':True,'cpu_pool_limit':72}))
            with patch.dict(os.environ,{},clear=True),patch('os.sched_getaffinity',return_value=set(range(72))):
                self.assertEqual(pool_cpu_limit(root,68),72)
                with patch.dict(os.environ,{'SLURM_CPUS_PER_TASK':'12'}):
                    self.assertEqual(pool_cpu_limit(root,68),12)
                with patch('integrate_wgs_rna.pool_leases.cgroup_cpu_limit',return_value=8):
                    self.assertEqual(pool_cpu_limit(root,68),8)

    def test_RAM막힌대기와최소CPU부족대기는분배에서제외(self):
        rows=[('large',{'min_cpus':8,'memory_bytes':400}),('a',{'min_cpus':4,'memory_bytes':60}),('b',{'min_cpus':4,'memory_bytes':60})]
        self.assertEqual([p for p,_ in fitting_waiters(rows,24,150)],['a','b'])
        self.assertEqual([p for p,_ in fitting_waiters(rows,7,500)],['a'])

    def test_동시에RAM이맞는한작업만모든가용CPU배정(self):
        spec=importlib.util.spec_from_file_location('live_pool',ROOT/'scripts/cpu_pool.py');pool=importlib.util.module_from_spec(spec);spec.loader.exec_module(pool)
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir();(root/'waiters').mkdir()
            for i in range(3):
                (root/f'waiters/{i}.json').write_text(json.dumps({'pid':os.getpid(),'queued_at':time.time()-5,'min_cpus':8,'max_cpus':72,'memory_bytes':400,'branch':'RNA','stage':'16','sample':str(i)}))
            with patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=600):pool.allocate(root,24)
            leases=list((root/'leases').glob('*.json'));self.assertEqual(len(leases),1)
            self.assertEqual(json.loads(leases[0].read_text())['cpus'],24)

    def test_구관리자도매할당때최신CPU상한을읽고기존계산은보존(self):
        spec=importlib.util.spec_from_file_location('resize_pool',ROOT/'scripts/cpu_pool.py')
        pool=importlib.util.module_from_spec(spec);spec.loader.exec_module(pool)
        with tempfile.TemporaryDirectory() as raw, patch.dict(os.environ,{},clear=True), \
                patch('os.sched_getaffinity',return_value=set(range(72))), \
                patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=1000):
            root=Path(raw);(root/'leases').mkdir();(root/'waiters').mkdir()
            active=root/'leases/active.json'
            pool.atomic_json(active,{'manager_pid':os.getpid(),'cpus':44,'memory_bytes':10})
            original=active.read_bytes()
            waiter={'pid':os.getpid(),'queued_at':time.time()-5,'min_cpus':2,
                    'max_cpus':68,'memory_bytes':10,'branch':'WGS','stage':'03','sample':'S'}
            pool.atomic_json(root/'waiters/first.json',waiter)
            pool.atomic_json(root/'policy.json',{'live_limits':True,'cpu_pool_limit':40})
            pool.allocate(root,68)
            self.assertEqual(active.read_bytes(),original)
            self.assertTrue((root/'waiters/first.json').exists())
            self.assertFalse((root/'leases/first.json').exists())
            pool.atomic_json(active,{'manager_pid':os.getpid(),'cpus':36,'memory_bytes':10})
            pool.allocate(root,68)
            first=pool.read_json(root/'leases/first.json')
            self.assertEqual((first['cpus'],first['total_cpus']),(4,40))
            pool.atomic_json(root/'policy.json',{'live_limits':True,'cpu_pool_limit':60})
            pool.atomic_json(root/'waiters/second.json',waiter)
            pool.allocate(root,68)
            second=pool.read_json(root/'leases/second.json')
            self.assertEqual((second['cpus'],second['total_cpus']),(20,60))

    def test_auto정렬은최소4계산스레드를확보(self):
        cfg=PersonalRealignConfig(cohort='t',output_root=Path('/tmp'),sample_list=None,samples=['S'],chroms=['chr1'],sequence_variant_policy='snv_indel',dataset_role='pre_training',platform='hifi')
        self.assertIsNone(alignment_plan(cfg,{1:3_000_000_000},6,1_350_000_000_000,AUTO))
        plan=alignment_plan(cfg,{1:3_000_000_000},16,1_350_000_000_000,AUTO)
        self.assertGreaterEqual(plan['threads'],4)
        fixed=alignment_plan(cfg,{1:3_000_000_000},6,1_350_000_000_000,dict(AUTO,MFM_WGS_ALIGN_REQUESTED='1'))
        self.assertEqual(fixed['threads'],1)

    def test_공유감시기가새상한과RSS기록을갱신(self):
        spec=importlib.util.spec_from_file_location('test_watchdog',ROOT/'scripts/resource_watchdog.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir();(root/'waiters').mkdir()
            (root/'policy.json').write_text(json.dumps({'live_limits':True,'memory_limit_bytes':1_500_000_000_000,'cpu_pool_limit':68}))
            with patch.dict(os.environ,{'MFM_RAM_LIMIT_BYTES':'1000000000000'},clear=True):
                state=module.tick(root,68,os.getpid())
            self.assertEqual(state['memory_limit_bytes'],1_500_000_000_000)
            self.assertEqual(state['reservation_bytes'],1_350_000_000_000)
            self.assertEqual(json.loads((root/'memory_guard.json').read_text())['limit_bytes'],1_500_000_000_000)

    def test_활성데이터계약은상한증가외변경을거부(self):
        before={'memory_limit_bytes':1_000_000_000_000,'cpu_pool_limit':68,'sample_cpus':17,'local_jobs':4,'pipeline_env':{'MFM_RAM_LIMIT_BYTES':'1000000000000','WINDOW_STORAGE_POLICY':'keep'}}
        after=copy.deepcopy(before);after.update(memory_limit_bytes=1_500_000_000_000,cpu_pool_limit=72);after['pipeline_env']['MFM_RAM_LIMIT_BYTES']='1500000000000'
        with patch('handoff.local_queue.effective_memory',side_effect=lambda c:c['memory_limit_bytes']),patch('handoff.local_queue.cpu_limit',side_effect=lambda c:c['cpu_pool_limit']):
            self.assertTrue(capacity_increase_only(before,after))
            bad=copy.deepcopy(after);bad['pipeline_env']['WINDOW_STORAGE_POLICY']='changed';self.assertFalse(capacity_increase_only(before,bad))
            bad=copy.deepcopy(after);bad['memory_limit_bytes']=500_000_000_000;self.assertFalse(capacity_increase_only(before,bad))

    def test_minimap2_HAP1보존후HAP2경계만증설(self):
        with tempfile.TemporaryDirectory() as raw:
            pool,root,script,report,lease,requests,metadata,env=BoundaryTests().fixture(raw)
            script.write_text('#!/bin/bash\ngrp() { "$SAMTOOLS" view -u "$INPUT" | "$SAMTOOLS" fastq -; }\ngrp | "$MINIMAP2" -ax map-hifi -y -Y --MD -t 1 -R RG target - | "$SAMTOOLS" sort -@ 0 -m 2G -o out -\n')
            before=metadata.read_bytes()
            with patch.dict(os.environ,env,clear=True):_,info=prepare_boundary(['bash',str(script)],report)
            self.assertEqual(info['platform'],'hifi');self.assertGreaterEqual(info['threads'],4)
            self.assertEqual(metadata.read_bytes(),before);self.assertIn('-t '+str(info['threads']),script.read_text())


if __name__=='__main__':unittest.main()
