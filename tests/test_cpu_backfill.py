"""실제 유휴 관측·엄격 RAM·정상 작업 우선·자식 nice 상속을 검증한다."""
from pathlib import Path
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

from integrate_wgs_rna import pool_backfill as module
from integrate_wgs_rna.pool_leases import regular_cpu_claim
from integrate_wgs_rna.realign.phases import RealignPhases

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('test_backfill_pool',ROOT/'scripts/cpu_pool.py')
POOL=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(POOL)


class BackfillTests(unittest.TestCase):
    def fixture(self,raw,*,memory=400,request=200):
        root=Path(raw)
        for name in ['leases','waiters']:(root/name).mkdir()
        (root/'policy.json').write_text(json.dumps({'cpu_backfill':{'enabled':True}}))
        (root/'leases/held.json').write_text(json.dumps({'cpus':8,'memory_bytes':memory}))
        (root/'waiters/a.json').write_text(json.dumps({'pid':os.getpid(),'sample':'S','branch':'WGS',
            'stage':'01','min_cpus':4,'max_cpus':8,'memory_bytes':request,'queued_at':time.time()-1}))
        observed=dict(ready=True,idle_cpus=6,used_cpus=1,peak_used_cpus=1,
                      memory_current=100,memory_max=10000,cpu_cgroup='fake',quota=8)
        return root,observed,module.CpuObserver(),{'action':'continue','rss_bytes':100}

    def lend(self,root,observed,observer,guard):
        with patch.object(module,'eligible_waiter',return_value=True),patch.object(module,'in_cpu_group',return_value=True),patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=1000),patch.object(os,'setpriority'),patch.object(os,'getpriority',return_value=19):
            return module.lend_idle_cpus(root,8,observed,observer,guard)

    def test_읽기QC대여전환은CPU수와RAM을보존하고정상예약만반환(self):
        with tempfile.TemporaryDirectory() as raw:
            root,obs,observer,guard=self.fixture(raw)
            p=root/'leases/read.json'
            row={'ticket':'read','manager_pid':os.getpid(),'sample':'R','branch':'WGS','stage':'01',
                 'prepare_phase':'read_only_qc','cpus':6,'memory_bytes':400}
            p.write_text(json.dumps(row))
            with patch.object(module,'in_cpu_group',return_value=True),patch.object(module,'lower_process_tree_priority',return_value=True):
                result=module.relend_fixed_prepare(root,[row],6,6,obs)
            self.assertEqual(result[0]['cpus'],6)
            converted=json.loads(p.read_text())
            self.assertEqual((converted['cpus'],converted['memory_bytes']),(6,400))
            self.assertEqual(regular_cpu_claim(converted),0)
            self.assertTrue(converted['cpu_backfill']['converted_fixed_prepare'])

    def test_다른단계전환금지와유휴대여상한및우선순위확인(self):
        for field,value in [('stage','10'),('cpus',7),('manager_pid',0)]:
            with self.subTest(field=field),tempfile.TemporaryDirectory() as raw:
                root,obs,observer,guard=self.fixture(raw)
                row={'ticket':'read','manager_pid':os.getpid(),'sample':'R','branch':'WGS','stage':'01',
                     'prepare_phase':'read_only_qc','cpus':6,'memory_bytes':400};row[field]=value
                p=root/'leases/read.json';p.write_text(json.dumps(row))
                with patch.object(module,'in_cpu_group',return_value=True),patch.object(module,'lower_process_tree_priority',return_value=True):
                    self.assertEqual(module.relend_fixed_prepare(root,[row],6,6,obs),[])
                self.assertNotIn('cpu_backfill',json.loads(p.read_text()))
        with tempfile.TemporaryDirectory() as raw:
            root,obs,observer,guard=self.fixture(raw)
            row={'ticket':'read','manager_pid':os.getpid(),'sample':'R','branch':'WGS','stage':'01',
                 'prepare_phase':'read_only_qc','cpus':6,'memory_bytes':400}
            (root/'leases/read.json').write_text(json.dumps(row))
            with patch.object(module,'in_cpu_group',return_value=True),patch.object(module,'lower_process_tree_priority',return_value=False):
                self.assertEqual(module.relend_fixed_prepare(root,[row],6,6,obs),[])

    def test_기존자식스레드의nice19및호출부모보존(self):
        before=os.getpriority(os.PRIO_PROCESS,0)
        proc=subprocess.Popen([sys.executable,'-c','import threading,time;threading.Thread(target=time.sleep,args=(30,)).start();time.sleep(30)'])
        try:
            deadline=time.monotonic()+3
            while len(list((Path('/proc')/str(proc.pid)/'task').glob('*')))<2 and time.monotonic()<deadline:time.sleep(.02)
            self.assertTrue(module.lower_process_tree_priority(proc.pid))
            tasks=list((Path('/proc')/str(proc.pid)/'task').glob('*'))
            self.assertGreaterEqual(len(tasks),2)
            self.assertTrue(all(os.getpriority(os.PRIO_PROCESS,int(p.name))==19 for p in tasks))
            self.assertEqual(os.getpriority(os.PRIO_PROCESS,0),before)
        finally:proc.terminate();proc.wait()

    def test_전체관측기간과실제quota가필요(self):
        observer=module.CpuObserver()
        for moment in range(0,12,2):
            state=observer.observe(8,now=moment,reading={'quota':8,'usage_usec':moment*2_000_000})
            self.assertEqual(state['ready'],moment>=10)
        self.assertEqual(state['idle_cpus'],5)
        self.assertFalse(observer.observe(8,now=12,reading={'quota':16,'usage_usec':1})['ready'])

    def test_순간평균보다최근최대사용량을우선(self):
        observer=module.CpuObserver()
        for moment in range(0,10,2):
            observer.observe(8,now=moment,reading={'quota':8,'usage_usec':moment*1_000_000})
        state=observer.observe(8,now=10,reading={'quota':8,'usage_usec':22_000_000})
        self.assertEqual(state['peak_used_cpus'],7)
        self.assertEqual(state['idle_cpus'],0)

    def test_CPU여유와RAM예산을모두지키며대여(self):
        with tempfile.TemporaryDirectory() as raw:
            root,obs,observer,guard=self.fixture(raw)
            result=self.lend(root,obs,observer,guard)
            self.assertEqual(result['launched'][0]['cpus'],6)
            row=json.loads((root/'leases/a.json').read_text())
            self.assertEqual(row['cpu_backfill']['nice'],19)
            self.assertEqual(regular_cpu_claim(row),0)
            self.assertEqual(sum(json.loads(p.read_text())['memory_bytes'] for p in (root/'leases').glob('*.json')),600)
            self.assertGreater(observer.cooldown_until,time.monotonic())

    def test_RAM예약이부족하면CPU가비어도대여하지않음(self):
        with tempfile.TemporaryDirectory() as raw:
            root,obs,observer,guard=self.fixture(raw,memory=900)
            self.assertEqual(self.lend(root,obs,observer,guard)['launched'],[])

    def test_cgroup실측과RSS도신규RAM예측에포함(self):
        for field,value in [('memory_current',8900),('rss_bytes',900)]:
            with self.subTest(field=field),tempfile.TemporaryDirectory() as raw:
                root,obs,observer,guard=self.fixture(raw)
                (guard if field=='rss_bytes' else obs)[field]=value
                self.assertEqual(self.lend(root,obs,observer,guard)['launched'],[])

    def test_기본예약과대여합계는최대2배(self):
        with tempfile.TemporaryDirectory() as raw:
            root,obs,observer,guard=self.fixture(raw)
            (root/'leases/borrowed.json').write_text(json.dumps({'cpus':8,'memory_bytes':100,'cpu_backfill':{'nice':19}}))
            self.assertEqual(self.lend(root,obs,observer,guard)['launched'],[])

    def test_IO대기감안대여는늘려도전체논리2배와낮은우선순위유지(self):
        with tempfile.TemporaryDirectory() as raw:
            root,obs,observer,guard=self.fixture(raw)
            (root/'policy.json').write_text(json.dumps({'cpu_backfill':{'enabled':True,'grant_multiplier':2}}))
            result=self.lend(root,obs,observer,guard)
            self.assertEqual(result['launched'][0]['cpus'],8)
            self.assertEqual(result['launched'][0]['nice'],19)

    def test_일시peak만으로평균유휴CPU대여를막지않음(self):
        with tempfile.TemporaryDirectory() as raw:
            root,obs,observer,guard=self.fixture(raw);obs.update(idle_cpus=0,used_cpus=2,peak_used_cpus=8)
            (root/'policy.json').write_text(json.dumps({'cpu_backfill':{'enabled':True,'grant_multiplier':2,'use_average_idle':True}}))
            result=self.lend(root,obs,observer,guard)
            self.assertEqual(result['launched'][0]['cpus'],8)
            self.assertEqual(result['launched'][0]['nice'],19)

    def test_명시한IO처리량정책은추가읽기를허용하고상한유지(self):
        with tempfile.TemporaryDirectory() as raw:
            root,obs,observer,guard=self.fixture(raw)
            (root/'leases/old_loan.json').write_text(json.dumps({'cpus':8,'memory_bytes':100,'cpu_backfill':{'nice':19}}))
            (root/'policy.json').write_text(json.dumps({'cpu_backfill':{'enabled':True,'grant_multiplier':2,'max_claim_pools':3,'max_loan_pools':2}}))
            result=self.lend(root,obs,observer,guard)
            self.assertEqual(result['launched'][0]['cpus'],8)
            self.assertEqual(result['logical_cpu_ceiling'],24)
            self.assertEqual(result['loan_cpu_ceiling'],16)
            self.assertEqual(result['launched'][0]['nice'],19)

    def test_명시확장도실측유휴_RAM과논리상한을유지(self):
        self.assertEqual(module.loan_limits(68,{}),(136,68))
        self.assertEqual(module.loan_limits(68,{'max_claim_pools':99,'max_loan_pools':99}),(272,204))
        for used,memory,expected in [(1,400,True),(8,400,False),(1,900,False)]:
            with self.subTest(used=used,memory=memory),tempfile.TemporaryDirectory() as raw:
                root,obs,observer,guard=self.fixture(raw,memory=memory)
                obs.update(used_cpus=used)
                (root/'leases/loans.json').write_text(json.dumps({'cpus':16,'memory_bytes':100,'cpu_backfill':{'nice':19}}))
                (root/'policy.json').write_text(json.dumps({'cpu_backfill':{'enabled':True,'grant_multiplier':2,
                    'max_claim_pools':4,'max_loan_pools':3,'use_average_idle':True}}))
                result=self.lend(root,obs,observer,guard)
                self.assertEqual(bool(result['launched']),expected)
                self.assertEqual(result['logical_cpu_ceiling'],32)
                self.assertEqual(result['loan_cpu_ceiling'],24)
                if expected:
                    self.assertLessEqual(result['launched'][0]['cpus'],8)
                    self.assertLessEqual(result['launched'][0]['cpus']+24,32)

    def test_명시한pair대여우선과기본FIFO를보존(self):
        for preferred,expected in [([], 'S'),(['P'],'P')]:
            with self.subTest(preferred=preferred),tempfile.TemporaryDirectory() as raw:
                root,obs,observer,guard=self.fixture(raw)
                second={'pid':os.getpid(),'sample':'P','branch':'WGS','stage':'02',
                        'min_cpus':4,'max_cpus':8,'memory_bytes':200,'queued_at':time.time()-1}
                (root/'waiters/b.json').write_text(json.dumps(second))
                (root/'policy.json').write_text(json.dumps({'cpu_backfill':{'enabled':True,
                        'validation_phases':True,'preferred_samples':preferred}}))
                result=self.lend(root,obs,observer,guard)
                self.assertEqual(result['launched'][0]['sample'],expected)
                self.assertEqual(len(result['launched']),1)

    def test_비활성정책과우선순위하향실패는대여없음(self):
        with tempfile.TemporaryDirectory() as raw:
            root,obs,observer,guard=self.fixture(raw)
            (root/'policy.json').write_text('{}')
            self.assertEqual(self.lend(root,obs,observer,guard)['launched'],[])
            (root/'policy.json').write_text(json.dumps({'cpu_backfill':{'enabled':True}}))
            with patch.object(module,'eligible_waiter',return_value=True),patch.object(module,'in_cpu_group',return_value=True),patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=1000),patch.object(os,'setpriority',side_effect=PermissionError):
                self.assertEqual(module.lend_idle_cpus(root,8,obs,observer,guard)['launched'],[])

    def test_정상예약이적어도대여자체는한풀까지만(self):
        with tempfile.TemporaryDirectory() as raw:
            root,obs,observer,guard=self.fixture(raw)
            (root/'leases/held.json').write_text(json.dumps({'cpus':0,'memory_bytes':100}))
            (root/'leases/borrowed.json').write_text(json.dumps({'cpus':8,'memory_bytes':100,'cpu_backfill':{'nice':19}}))
            self.assertEqual(self.lend(root,obs,observer,guard)['launched'],[])

    def test_재정렬단계는고정자원대여에서제외(self):
        row={'pid':os.getpid(),'branch':'WGS','stage':'10','cpu_utilization':'maximum'}
        self.assertFalse(module.eligible_waiter(row))
        row.update(stage='01',cpu_utilization='bounded')
        self.assertFalse(module.eligible_waiter(row))
        row['cpu_utilization']='maximum'
        self.assertTrue(module.eligible_waiter(row))

    def test_고정heap마크듀프는대여하고가변coverage는제외(self):
        row={'pid':os.getpid(),'branch':'RNA','stage':'17','layout':'original','cpu_utilization':'maximum'}
        self.assertTrue(module.eligible_waiter(row))
        row['stage']='18';self.assertFalse(module.eligible_waiter(row))
        row.update(stage='15',layout='legacy');self.assertTrue(module.eligible_waiter(row))

    def test_정상할당은대여CPU에막히지않고RAM은합산(self):
        with tempfile.TemporaryDirectory() as raw:
            root,obs,observer,guard=self.fixture(raw,memory=200)
            (root/'leases/held.json').write_text(json.dumps({'cpus':2,'memory_bytes':200}))
            (root/'leases/loan.json').write_text(json.dumps({'cpus':6,'memory_bytes':200,'cpu_backfill':{'nice':19}}))
            with patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=1000):POOL.allocate(root,8)
            row=json.loads((root/'leases/a.json').read_text())
            self.assertEqual(row['cpus'],6)
            self.assertEqual(sum(regular_cpu_claim(json.loads(p.read_text())) for p in (root/'leases').glob('*.json')),8)

    def test_HAP경계도대여CPU때문에막히지않음(self):
        with tempfile.TemporaryDirectory() as raw:
            root,obs,observer,guard=self.fixture(raw,memory=10*1024**3)
            (root/'waiters/a.json').unlink()
            (root/'leases/held.json').write_text(json.dumps({'cpus':3,'memory_bytes':10*1024**3}))
            (root/'leases/loan.json').write_text(json.dumps({'cpus':6,'memory_bytes':10*1024**3,'cpu_backfill':{'nice':19}}))
            own=root/'leases/own.json';own.write_text(json.dumps({'cpus':1,'total_cpus':8,'memory_bytes':4*1024**3}))
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(own),'MFM_CPU_BUDGET':'8'},clear=True):
                self.assertEqual(RealignPhases().reserve('test',max_units=8),5)

    def test_실제대기manager에서자식으로nice19를상속(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir();(root/'waiters').mkdir()
            (root/'policy.json').write_text(json.dumps({'memory_limit_bytes':500_000_000_000,'cpu_backfill':{'enabled':True}}))
            (root/'leases/held.json').write_text(json.dumps({'cpus':8,'memory_bytes':64*1024**3}))
            result=root/'child.json';audit=root/'audit.jsonl'
            code="import os,json,pathlib;pathlib.Path("+repr(str(result))+").write_text(json.dumps({'nice':os.getpriority(os.PRIO_PROCESS,0),'cpus':os.environ['MFM_CPU_POOL_GRANTED']}))"
            env={**os.environ,'MFM_CPU_UTILIZATION':'maximum','MFM_RAM_LIMIT_BYTES':'500000000000'}
            with (root/'output.log').open('w') as out:
                process=subprocess.Popen([sys.executable,str(ROOT/'scripts/cpu_pool.py'),'--pool-dir',str(root),
                    '--total-cpus','8','--branch','WGS','--stage','01','--sample','B','--audit',str(audit),
                    '--',sys.executable,'-c',code],env=env,stdout=out,stderr=subprocess.STDOUT)
                try:
                    deadline=time.monotonic()+5
                    while not list((root/'waiters').glob('*.json')) and time.monotonic()<deadline:time.sleep(.05)
                    time.sleep(.4)
                    obs=dict(ready=True,idle_cpus=6,used_cpus=1,peak_used_cpus=1,memory_current=0,
                             memory_max=2_000_000_000_000,cpu_cgroup=module.cgroup_sample()['cpu_cgroup'],quota=8)
                    with (root/'pool.lock').open('a+') as lock:
                        import fcntl
                        fcntl.flock(lock,fcntl.LOCK_EX)
                        loan=module.lend_idle_cpus(root,8,obs,module.CpuObserver(),{'action':'continue','rss_bytes':0})
                    self.assertEqual(len(loan['launched']),1)
                    self.assertEqual(process.wait(timeout=6),0,(root/'output.log').read_text())
                    self.assertEqual(json.loads(result.read_text()),{'nice':19,'cpus':'6'})
                    self.assertEqual(json.loads(audit.read_text())['cpu_backfill']['nice'],19)
                finally:
                    if process.poll() is None:process.kill();process.wait()


if __name__=='__main__':unittest.main()
