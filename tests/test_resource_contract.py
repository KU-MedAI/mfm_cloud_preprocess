"""일반 입장·내부 재확보·대여·RAM 잔여량의 공통 계약을 검증함."""
from pathlib import Path
from unittest.mock import patch
import importlib.util
import json
import os
import tempfile
import unittest

from integrate_wgs_rna.cpu_executor import WorkerBudget

GIB=1024**3
ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('contract_pool',ROOT/'scripts/cpu_pool.py')
POOL=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(POOL)


class ResourceContractTests(unittest.TestCase):
    def fixture(self, raw, *, held=1):
        root=Path(raw)
        for name in ['leases','waiters']:(root/name).mkdir()
        setting=dict(enabled=True,validation_phases=True,haplotag_phases=True,
                     max_claim_pools=2,max_loan_pools=1,max_extra_cpus=16)
        policy=dict(memory_limit_bytes=200*GIB,cpu_backfill=setting)
        (root/'policy.json').write_text(json.dumps(policy))
        for name, row in [('regular',dict(cpus=68,memory_bytes=4*GIB)),
                          ('loan',dict(cpus=12,memory_bytes=4*GIB,cpu_backfill={'cpus':12}))]:
            row['manager_pid']=os.getpid()
            (root/'leases'/(name+'.json')).write_text(json.dumps(row))
        path=root/'leases/own.json'
        row=dict(manager_pid=os.getpid(),child_pid=None,cpus=held,total_cpus=68,
                 memory_bytes=16*GIB,branch='WGS',stage='02',
                 cpu_backfill=dict(cpus=8,phase_reacquire_cpus=8,nice=0))
        path.write_text(json.dumps(row))
        return root,path,policy

    def capacity(self,path,**kwargs):
        with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path),'MFM_RAM_LIMIT_BYTES':str(200*GIB)}),\
             patch('integrate_wgs_rna.cpu_executor.pool_cpu_limit',return_value=68):
            return WorkerBudget(8,8,True,GIB,**kwargs).capacity(8)

    def test_68정상12대여중최초8을요청해도추가16중남은4만재확보(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,policy=self.fixture(raw)
            self.assertEqual(self.capacity(path),4)
            rows=[json.loads(p.read_text()) for p in (root/'leases').glob('*.json')]
            self.assertEqual(sum(r['cpus'] for r in rows),84)

    def test_전역해제및단계해제시현재보유량보다확대하지않음(self):
        for field in ['enabled','validation_phases']:
            with self.subTest(field=field),tempfile.TemporaryDirectory() as raw:
                root,path,policy=self.fixture(raw)
                policy['cpu_backfill'][field]=False
                (root/'policy.json').write_text(json.dumps(policy))
                self.assertEqual(self.capacity(path),1)

    def test_부모와worker비용을합산하고최초grant도지킴(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,policy=self.fixture(raw,held=0)
            self.assertEqual(self.capacity(path,cpus_per_worker=2,parent_cpus=1),1)
            self.assertEqual(json.loads(path.read_text())['cpus'],3)

    def test_자식시작직전에도최신추가상한을재검사(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,policy=self.fixture(raw,held=8)
            with patch.object(POOL,'pool_cpu_limit',return_value=68):
                row=POOL.reconcile_launch_memory(root,path,16*GIB,2)
            self.assertEqual(row['cpus'],4)

    def test_추가상한이없으면기존배수계약과최초grant유지(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,policy=self.fixture(raw)
            policy['cpu_backfill'].pop('max_extra_cpus')
            (root/'policy.json').write_text(json.dumps(policy))
            self.assertEqual(self.capacity(path),8)

    def test_대여정책축소시실행중worker는보존함(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,policy=self.fixture(raw,held=4)
            policy['cpu_backfill']['max_extra_cpus']=0
            (root/'policy.json').write_text(json.dumps(policy))
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path),'MFM_RAM_LIMIT_BYTES':str(200*GIB)}),patch('integrate_wgs_rna.cpu_executor.pool_cpu_limit',return_value=68):
                budget=WorkerBudget(8,8,True,GIB,cpus_per_worker=2)
                self.assertEqual(budget.capacity(8,running=2),2)
                self.assertEqual(json.loads(path.read_text())['cpus'],4)

    def test_시작전최소CPU미달은반환뒤정책잔여로재시도함(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,policy=self.fixture(raw,held=8)
            def release(_):
                row=json.loads(path.read_text())
                self.assertEqual((row['cpus'],row['memory_bytes']),(0,4*GIB))
                other=root/'leases/loan.json';row=json.loads(other.read_text());row['cpus']=10
                POOL.atomic_json(other,row)
            with patch.object(POOL,'pool_cpu_limit',return_value=68),patch.object(POOL.time,'sleep',side_effect=release):
                row=POOL.reconcile_launch_memory(root,path,16*GIB,6)
            self.assertEqual(row['cpus'],6)

    def test_정상풀도대기부모를포함한전체논리상한을지킴(self):
        from integrate_wgs_rna.cpu_loan_policy import available_regular_cpus,available_loan_cpus
        setting={'enabled':True,'max_extra_cpus':16}
        rows=[dict(cpus=67),dict(cpus=16,cpu_backfill={'cpus':16}),dict(cpus=1,realign_phase='align:waiting')]
        self.assertEqual(available_regular_cpus(68,setting,rows),0)
        self.assertEqual(available_loan_cpus(68,setting,rows),0)
        rows=[dict(cpus=40),dict(cpus=12,cpu_backfill={'cpus':12})]
        self.assertEqual(available_loan_cpus(68,setting,iter(rows)),4)
        self.assertEqual(available_regular_cpus(68,setting,iter(rows)),28)

    def test_계산이없는순수대기는부모CPU도반환함(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,policy=self.fixture(raw,held=0)
            policy['cpu_backfill']['enabled']=False
            (root/'policy.json').write_text(json.dumps(policy))
            self.assertEqual(self.capacity(path,parent_cpus=1),0)
            self.assertEqual(json.loads(path.read_text())['cpus'],0)

    def test_고정샘플은CPU소유권을한번만세고RAM은전역합산함(self):
        from integrate_wgs_rna.cpu_loan_policy import available_regular_cpus,available_loan_cpus
        from integrate_wgs_rna.pool_leases import fitting_waiters
        from integrate_wgs_rna.memory_policy import available_memory
        rows=[dict(cpus=4,memory_bytes=4,fixed_sample_token='A',sample_cpu_owner=True),
              dict(cpus=2,memory_bytes=10,fixed_sample_token='A'),
              dict(cpus=4,memory_bytes=4,fixed_sample_token='B',sample_cpu_owner=True),
              dict(cpus=3,memory_bytes=20,fixed_sample_token='B')]
        with patch.dict(os.environ,{'MFM_FIXED_SAMPLE_TOKEN':''}):
            self.assertEqual(available_regular_cpus(8,{},rows),0)
        with patch.dict(os.environ,{'MFM_FIXED_SAMPLE_TOKEN':'A'}):
            self.assertEqual(available_regular_cpus(4,{},rows),2)
            self.assertEqual(available_loan_cpus(4,{'enabled':True},rows),0)
            self.assertEqual(fitting_waiters([(Path('B'),dict(fixed_sample_token='B',min_cpus=1))],2,100),[])
            with patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=100):
                self.assertEqual(available_memory(Path('.'),rows),62)

    def test_RAM잔여계산은예약_관측_더작은상한을함께반영(self):
        from integrate_wgs_rna.memory_policy import available_memory
        with patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=1000):
            others=[dict(memory_bytes=200)]
            self.assertEqual(available_memory(Path('.'),others),800)
            self.assertEqual(available_memory(Path('.'),others,ceiling_bytes=700),500)
            observed=dict(memory_max=1000,memory_current=850,memory_working_bytes=600)
            self.assertEqual(available_memory(Path('.'),others,observed=observed),300)
            observed.pop('memory_working_bytes')
            self.assertEqual(available_memory(Path('.'),others,observed=observed),50)
            self.assertEqual(available_memory(Path('.'),others,rss_bytes=900),100)

    def test_재정렬도공통캐시메모리판정을사용(self):
        from tests.test_realign_idle import IdleTests
        case=IdleTests()
        case.setUp()
        try:
            case.current.update(memory_current=1_500_000_000_000,memory_working_bytes=450_000_000_000)
            self.assertEqual(case.capacity(memory=400_000_000_000)[1],400_000_000_000)
        finally:
            case.doCleanups()


if __name__=='__main__':unittest.main()
