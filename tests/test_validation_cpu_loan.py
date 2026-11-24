"""검증 하위 묶음의 CPU 반환·재진입과 대여·RAM 상한을 검사한다."""
from pathlib import Path
import json
import os
import tempfile
import time
import unittest
from unittest.mock import patch
from integrate_wgs_rna.cpu_executor import WorkerBudget
from integrate_wgs_rna.input_parallel import parallel_calls
from integrate_wgs_rna import pool_backfill as module

GIB=1024**3


class ValidationLoanTests(unittest.TestCase):
    def fixture(self, folder):
        root=Path(folder)
        for name in ['leases','waiters']:(root/name).mkdir()
        (root/'policy.json').write_text(json.dumps({'memory_limit_bytes':200*GIB,
            'cpu_backfill':{'enabled':True,'validation_phases':True,'max_claim_pools':3,'max_loan_pools':2}}))
        (root/'leases/other.json').write_text(json.dumps({'manager_pid':os.getpid(),'cpus':8,'memory_bytes':4*GIB}))
        own=root/'leases/own.json'
        own.write_text(json.dumps({'manager_pid':os.getpid(),'cpus':8,'total_cpus':8,'memory_bytes':16*GIB,
            'branch':'WGS','stage':'02','cpu_backfill':{'cpus':8,'nice':19,'phase_reacquire_cpus':8}}))
        return root,own

    def test_정상풀이가득차도두묶음의병렬검증이반환후재확보(self):
        with tempfile.TemporaryDirectory() as raw:
            root,own=self.fixture(raw)
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(own),'MFM_RAM_LIMIT_BYTES':str(200*GIB)}), \
                 patch('integrate_wgs_rna.cpu_executor.pool_cpu_limit',return_value=8):
                def task(value):time.sleep(.04);return value*value
                for _ in range(2):
                    self.assertEqual(parallel_calls(task,[(i,) for i in range(8)],4,memory_per_task=GIB),[i*i for i in range(8)])
                    row=json.loads(own.read_text());self.assertEqual(row['cpus'],1)
                    self.assertLessEqual(row['memory_bytes'],16*GIB)

    def test_대여상한과각작업최초grant를넘지않고RAM도제한(self):
        with tempfile.TemporaryDirectory() as raw:
            root,own=self.fixture(raw)
            (root/'leases/loan.json').write_text(json.dumps({'manager_pid':os.getpid(),'cpus':12,'memory_bytes':4*GIB,'cpu_backfill':{'nice':19}}))
            row=json.loads(own.read_text());row['cpus']=1;own.write_text(json.dumps(row))
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(own),'MFM_RAM_LIMIT_BYTES':str(200*GIB)}), \
                 patch('integrate_wgs_rna.cpu_executor.pool_cpu_limit',return_value=8):
                budget=WorkerBudget(4,4,True,GIB,cpus_per_worker=2)
                self.assertEqual(budget.capacity(4),2)
                self.assertEqual(json.loads(own.read_text())['cpus'],4)
                with patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=8*GIB):
                    self.assertEqual(budget.capacity(4),0)

    def test_02진입은명시정책과nice19및관측을요구(self):
        with tempfile.TemporaryDirectory() as raw:
            root,own=self.fixture(raw);own.unlink()
            row={'pid':os.getpid(),'sample':'S','branch':'WGS','stage':'02','memory_bytes':4*GIB,
                 'queued_at':time.time()-1,'min_cpus':2,'max_cpus':4,'cpu_utilization':'maximum'}
            p=root/'waiters/a.json';p.write_text(json.dumps(row))
            obs={'ready':True,'idle_cpus':4,'used_cpus':3,'peak_used_cpus':3,'memory_current':0,
                 'memory_max':200*GIB,'cpu_cgroup':'test'}
            observer=module.CpuObserver();guard={'action':'continue','rss_bytes':0}
            with patch.object(module,'in_cpu_group',return_value=True),patch.object(os,'setpriority'),patch.object(os,'getpriority',return_value=19):
                policy=json.loads((root/'policy.json').read_text());policy['cpu_backfill']['validation_phases']=False
                (root/'policy.json').write_text(json.dumps(policy))
                self.assertEqual(module.lend_idle_cpus(root,8,obs,observer,guard)['launched'],[])
                policy['cpu_backfill']['validation_phases']=True;(root/'policy.json').write_text(json.dumps(policy))
                value=module.lend_idle_cpus(root,8,obs,observer,guard)
                self.assertEqual(value['launched'][0]['cpus'],4)
                loan=json.loads((root/'leases/a.json').read_text())['cpu_backfill']
                self.assertEqual(loan['phase_reacquire_cpus'],4)
                self.assertEqual(loan['nice'],19)


if __name__=='__main__':unittest.main()
