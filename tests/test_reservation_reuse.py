"""전처리 예약 재사용의 CPU 관측·한도·RAM 보호·기존 정책 보존을 검증함."""
from pathlib import Path
import importlib.util
import json
import os
import tempfile
import time
import unittest
from unittest.mock import patch

from integrate_wgs_rna import pool_backfill as backfill
from integrate_wgs_rna.cpu_loan_policy import loan_limits
from integrate_wgs_rna.pool_cpu_observer import PreprocessingCpuObserver, reuse_capacity, tree_cpu_sample

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('reuse_test_pool', ROOT/'scripts/cpu_pool.py')
POOL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(POOL)


class ReservationReuseTests(unittest.TestCase):
    def test_기존한도보존과추가16CPU한도(self):
        self.assertEqual(loan_limits(68, {}), (136, 68))
        self.assertEqual(loan_limits(68, {'max_extra_cpus': 16}), (84, 16))
        self.assertEqual(loan_limits(68, {'max_extra_cpus': 0}), (68, 0))

    def test_PID시작시각과회수된자식CPU및중복트리(self):
        with tempfile.TemporaryDirectory() as raw:
            proc=Path(raw)
            def record(pid, parent, start, ticks, child_ticks=0):
                fields=['0']*22
                fields[0]='S';fields[1]=str(parent);fields[11]=str(ticks)
                fields[13]=str(child_ticks);fields[19]=str(start)
                p=proc/str(pid);p.mkdir(exist_ok=True)
                (p/'stat').write_text(str(pid)+' (test worker) '+' '.join(fields))
            record(10, 1, 100, 5);record(11, 10, 200, 7)
            roots=[{'pid':10,'started':100},{'pid':11,'started':200}]
            before=tree_cpu_sample(roots,proc=proc)
            self.assertEqual(before['usage_usec'],12*1_000_000//os.sysconf('SC_CLK_TCK'))
            (proc/'11/stat').unlink();(proc/'11').rmdir()
            record(10,1,100,5,7)
            self.assertEqual(tree_cpu_sample(roots,proc=proc)['usage_usec'],before['usage_usec'])
            self.assertIsNone(tree_cpu_sample([{'pid':10,'started':999}],proc=proc))

    def test_10초관측과범위변경및실종시재관측(self):
        observer=PreprocessingCpuObserver()
        setting={'enabled':True,'reservation_reuse':{'enabled':True,'roots':[]}}
        for moment in range(0,12,2):
            state=observer.observe(8,setting,now=moment,reader=lambda roots:dict(scope=[(1,100)],usage_usec=moment*1_000_000))
            self.assertEqual(state['ready'],moment>=10)
        self.assertEqual(reuse_capacity(8,setting,state,now=10),5)
        self.assertEqual(reuse_capacity(8,setting,state,now=14),0)
        state=observer.observe(8,setting,now=12,reader=lambda roots:dict(scope=[(2,100)],usage_usec=12_000_000))
        self.assertFalse(state['ready'])
        state=observer.observe(8,setting,now=14,reader=lambda roots:None)
        self.assertFalse(state['ready'])

    def test_비활성및누계감소시대여금지(self):
        observer=PreprocessingCpuObserver()
        self.assertFalse(observer.observe(8,{},reader=lambda roots:self.fail())['ready'])
        setting={'enabled':True,'reservation_reuse':{'enabled':True}}
        for moment in range(0,12,2):
            observer.observe(8,setting,now=moment,reader=lambda roots:dict(scope=[(1,1)],usage_usec=moment*1_000_000))
        state=observer.observe(8,setting,now=12,reader=lambda roots:dict(scope=[(1,1)],usage_usec=0))
        self.assertFalse(state['ready'])

    def fixture(self, root, *, memory=100):
        for directory in ['leases','waiters']:(root/directory).mkdir()
        setting={'enabled':True,'nice':0,'max_extra_cpus':4,'reservation_reuse':{'enabled':True}}
        (root/'policy.json').write_text(json.dumps({'cpu_backfill':setting}))
        (root/'leases/held.json').write_text(json.dumps({'cpus':8,'memory_bytes':memory}))
        for name in ['a','b']:
            (root/f'waiters/{name}.json').write_text(json.dumps(dict(pid=os.getpid(),sample=name,branch='RNA',stage='15',
                min_cpus=4,max_cpus=8,memory_bytes=200,queued_at=time.time()-1)))
        shared=dict(ready=True,idle_cpus=0,used_cpus=8,peak_used_cpus=8,
                    memory_current=100,memory_max=10000,cpu_cgroup='fake',quota=8)
        own=dict(ready=True,used_cpus=1,peak_used_cpus=1,sampled_at=time.monotonic())
        return setting,shared,own

    def lend(self, root, shared, own):
        with patch.object(backfill,'eligible_waiter',return_value=True),patch.object(backfill,'in_cpu_group',return_value=True),\
                patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=1000),patch.object(os,'setpriority') as priority,\
                patch.object(os,'getpriority',return_value=0):
            result=backfill.lend_idle_cpus(root,8,shared,backfill.CpuObserver(),{'action':'continue','rss_bytes':100},preprocessing=own)
            if result['launched']:
                priority.assert_called_once_with(os.PRIO_PROCESS,os.getpid(),0)
            return result

    def test_공유CPU포화여도자체미사용예약만4CPU재사용하고RAM전량유지(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);setting,shared,own=self.fixture(root)
            result=self.lend(root,shared,own)
            self.assertEqual(result['allocation_basis'],'reservation_reuse')
            self.assertEqual([r['cpus'] for r in result['launched']],[4])
            rows=[json.loads(p.read_text()) for p in (root/'leases').glob('*.json')]
            self.assertEqual(sum(r['cpus'] for r in rows),12)
            self.assertEqual(sum(r['memory_bytes'] for r in rows),300)
            self.assertTrue((root/'waiters/b.json').exists())

    def test_관측미완료및자체포화및RAM부족은차단(self):
        for case in ['unready','busy','ram','disabled']:
            with self.subTest(case=case),tempfile.TemporaryDirectory() as raw:
                root=Path(raw);setting,shared,own=self.fixture(root,memory=900 if case=='ram' else 100)
                if case=='unready':own['ready']=False
                if case=='busy':own['peak_used_cpus']=8
                if case=='disabled':
                    setting.pop('reservation_reuse')
                    (root/'policy.json').write_text(json.dumps({'cpu_backfill':setting}))
                self.assertEqual(self.lend(root,shared,own)['launched'],[])

    def test_서열전용검증정리만4GiB예약(self):
        with patch.dict(os.environ,{'MFM_WGS_PAYLOAD':'sequence_only'},clear=True):
            self.assertEqual(POOL.stage_memory('WGS','23'),4*1024**3)
            self.assertEqual(POOL.stage_memory('WGS','24'),4*1024**3)
            self.assertEqual(POOL.stage_memory('RNA','23'),64*1024**3)
        with patch.dict(os.environ,{},clear=True):
            self.assertEqual(POOL.stage_memory('WGS','23'),64*1024**3)


if __name__=='__main__':unittest.main()
