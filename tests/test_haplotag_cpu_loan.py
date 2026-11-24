"""haplotag의 최초 대여 grant·RAM 보존·실제 WhatsHap 동등성을 검증한다."""
from pathlib import Path
import json
import os
import tempfile
import time
import unittest
from unittest.mock import patch
from integrate_wgs_rna.haplotag_parallel import ChromosomeBudget
from integrate_wgs_rna import pool_backfill as module


class HaplotagLoanTests(unittest.TestCase):
    def pool(self,root,cpus=3):
        for name in ['leases','waiters']:(root/name).mkdir()
        (root/'policy.json').write_text(json.dumps({'memory_limit_bytes':500*1024**3,
            'cpu_backfill':{'enabled':True,'haplotag_phases':True}}))
        path=root/'leases/own.json'
        path.write_text(json.dumps({'manager_pid':os.getpid(),'branch':'WGS','stage':'06','cpus':cpus,
            'total_cpus':8,'memory_bytes':128*1024**3,'cpu_backfill':{'cpus':cpus,'nice':19}}))
        return path

    def test_대기없어도최초grant초과증설금지_RAM유지(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);path=self.pool(root)
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)}):
                budget=ChromosomeBudget(3,22)
                self.assertEqual(budget.capacity(True),3)
                budget.remaining=1;self.assertEqual(budget.capacity(False),1)
            self.assertEqual(json.loads(path.read_text())['memory_bytes'],128*1024**3)

    def test_실제WhatsHap가대여중에도직렬과모든record_tag순서동일(self):
        from test_haplotag_parallel import ParallelHaplotagTest
        helper=ParallelHaplotagTest();helper.setUp()
        try:
            path=self.pool(helper.root)
            (helper.root/'leases/normal.json').write_text(json.dumps({'manager_pid':os.getpid(),'cpus':8,'memory_bytes':4*1024**3}))
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)}):helper.compare(helper.bam)
            self.assertEqual(json.loads(path.read_text())['cpus'],1)
        finally:helper.tearDown()

    def test_06입장은명시정책에서만최초grant를대여(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);path=self.pool(root);path.unlink()
            row={'pid':os.getpid(),'sample':'S','branch':'WGS','stage':'06','memory_bytes':128*1024**3,
                'queued_at':time.time()-1,'min_cpus':2,'max_cpus':4,'cpu_utilization':'maximum'}
            (root/'waiters/a.json').write_text(json.dumps(row));obs={'ready':True,'idle_cpus':4,'used_cpus':3,
                'peak_used_cpus':3,'memory_current':0,'memory_max':500*1024**3,'cpu_cgroup':'test'}
            observer=module.CpuObserver();guard={'action':'continue','rss_bytes':0}
            with patch.object(module,'in_cpu_group',return_value=True),patch.object(os,'setpriority'),patch.object(os,'getpriority',return_value=19):
                policy=json.loads((root/'policy.json').read_text());policy['cpu_backfill']['haplotag_phases']=False;(root/'policy.json').write_text(json.dumps(policy))
                self.assertFalse(module.lend_idle_cpus(root,8,obs,observer,guard)['launched'])
                policy['cpu_backfill']['haplotag_phases']=True;(root/'policy.json').write_text(json.dumps(policy))
                self.assertEqual(module.lend_idle_cpus(root,8,obs,observer,guard)['launched'][0]['cpus'],4)


if __name__=='__main__':unittest.main()
