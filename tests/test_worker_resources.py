"""신규 worker 경계 전환·Slurm 제외·입력 계약 보존을 확인함."""
import copy,json,os,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from handoff.worker_resources import apply
from integrate_wgs_rna.vcf_execution import compression_threads
from types import SimpleNamespace

class WorkerPolicyTests(unittest.TestCase):
    def test_local_policy_changes_only_future_resource_config(self):
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{},clear=True):
            root=Path(raw);(root/'policy.json').write_text(json.dumps({'future_sample_resources':{
                'cpu_allocation_mode':'fixed_sample','sample_cpus':16,'local_jobs':4}}))
            c={'pipeline_env':{'MFM_CPU_POOL_DIR':raw,'WINDOW_DROP_N':'0'},'local_jobs':6,
               'sample_cpus':68,'cpu_allocation_mode':'dynamic','run_id':'same','output_root':'same-output'}
            before=copy.deepcopy(c)
            def resolve(value,overrides):return {**value,**overrides,'cpu_pool_limit':68}
            with patch('handoff.worker_resources.resources.resolve',side_effect=resolve):
                updated=apply(c)
                self.assertEqual((updated['cpu_allocation_mode'],updated['sample_cpus'],updated['local_jobs']),('fixed_sample',16,4))
                self.assertEqual((updated['run_id'],updated['output_root']),('same','same-output'))
                self.assertEqual(c,before)
                c['local_jobs']=2
                self.assertEqual(apply(c)['local_jobs'],2)
                self.assertIs(apply(c,'slurm'),c)
                gpu=copy.deepcopy(c);gpu['run_id']='gpu_compare_test';gpu['pipeline_env'].update(MFM_GPU_COMPARE_UNRESERVED='1',MFM_WGS_BACKEND='parabricks')
                self.assertIs(apply(gpu),gpu)
                with patch.dict(os.environ,SLURM_JOB_ID='job'):self.assertIs(apply(c),c)
    def test_missing_and_invalid_policy(self):
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{},clear=True):
            c={'pipeline_env':{'MFM_CPU_POOL_DIR':raw}}
            self.assertIs(apply(c),c)
            p=Path(raw)/'policy.json'
            for value in [False,{},dict(cpu_allocation_mode='fixed_sample',sample_cpus=True,local_jobs=4)]:
                p.write_text(json.dumps({'future_sample_resources':value}))
                with self.assertRaises(ValueError):apply(c)
    def test_vcf_selection_and_conversion_fit_same_cpu_and_memory(self):
        from integrate_wgs_rna.vcf_execution import resource_plan
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{},clear=True):
            lease=Path(raw)/'lease.json';lease.write_text(json.dumps({'cpus':16,'memory_bytes':20<<30}))
            with patch.dict(os.environ,MFM_CPU_BUDGET='16',MFM_RAM_LIMIT_BYTES=str(128<<30),MFM_CPU_POOL_LEASE_FILE=str(lease)):
                p=resource_plan();self.assertEqual((p['selection_workers'],p['chromosome_workers']),(15,5))
                self.assertLessEqual(p['estimated_peak_bytes'],p['ram_budget_bytes'])
                p=resource_plan(4);self.assertLessEqual(p['estimated_peak_bytes']*4,20<<30)
                lease.write_text(json.dumps({'cpus':16,'memory_bytes':1800<<20}))
                with self.assertRaisesRegex(ValueError,'RAM'):resource_plan()
    def test_compression_threads_respects_main_and_explicit_limits(self):
        self.assertEqual(compression_threads(SimpleNamespace(concat_threads=0),16),15)
        self.assertEqual(compression_threads(SimpleNamespace(concat_threads=4),16),4)
        self.assertEqual(compression_threads(SimpleNamespace(concat_threads=24),16),15)
        self.assertEqual(compression_threads(SimpleNamespace(concat_threads=0),1),0)
        with self.assertRaises(ValueError):compression_threads(SimpleNamespace(concat_threads=-1),16)

class FixedWaitTests(unittest.TestCase):
    def setup_pool(self,root):
        import os
        (root/'leases').mkdir();(root/'waiters').mkdir()
        (root/'policy.json').write_text(json.dumps({'live_limits':True,'cpu_pool_limit':4,'memory_limit_bytes':1650000000000,'cpu_backfill':{'enabled':False}}))
        (root/'leases/block.json').write_text(json.dumps({'cpus':4,'memory_bytes':0,'manager_pid':os.getpid(),'sample':'OTHER'}))
        return dict(cpu_allocation_mode='fixed_sample',sample_cpus=4,cpu_pool_limit=4,
                    state_root=str(root/'state'),run_id='test',memory_limit_bytes=1650000000000,
                    pipeline_env={'MFM_CPU_POOL_DIR':str(root)})
    def test_waiter_reserves_minimum_fairly_and_cleans_on_grant(self):
        import time
        from integrate_wgs_rna.sample_cpu import sample_slot
        from integrate_wgs_rna.resource_fairness import priority_request,fair_cpu_capacity
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{},clear=True):
            root=Path(raw);c=self.setup_pool(root);seen=[]
            def release(_):
                path=next((root/'leases').glob('samplewait_*.json'));row=json.loads(path.read_text())
                self.assertEqual(row['cpus'],0);self.assertEqual(row['memory_bytes'],512<<20)
                row['resource_wait']['queued_at']=time.time()-31;path.write_text(json.dumps(row))
                self.assertEqual(priority_request(root,4,[row])['sample'],'NEW')
                self.assertEqual(fair_cpu_capacity(root,4,[row],4,request={'sample':'OTHER'}),0)
                seen.append(True);(root/'leases/block.json').unlink()
            with patch('integrate_wgs_rna.sample_cpu.time.sleep',side_effect=release):
                with sample_slot(c,'NEW'):
                    self.assertFalse(list((root/'leases').glob('samplewait_*.json')))
                    self.assertEqual(len(os.sched_getaffinity(0)),4)
            self.assertTrue(seen);self.assertFalse(list((root/'leases').glob('sample*.json')))
    def test_waiter_cleanup_on_exception(self):
        from integrate_wgs_rna.sample_cpu import sample_slot
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{},clear=True):
            root=Path(raw);c=self.setup_pool(root)
            with patch('integrate_wgs_rna.sample_cpu.time.sleep',side_effect=RuntimeError('중단')):
                with self.assertRaisesRegex(RuntimeError,'중단'):
                    with sample_slot(c,'NEW'):self.fail('할당 불가')
            self.assertFalse(list((root/'leases').glob('samplewait_*.json')))
