"""후반 RAM 반환·index 재확보·대여 상한·예산 대기를 검증함."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import json
import os
import tempfile
import threading
import time
import unittest

from integrate_wgs_rna.haplotag_resources import reserve_tail
from integrate_wgs_rna.pipeline.checkpoints import atomic_json

GIB=1024**3


class HaplotagResourcesTest(unittest.TestCase):
    def fixture(self, raw, *, loan=False):
        root=Path(raw);(root/'leases').mkdir();(root/'waiters').mkdir()
        policy={'live_limits':True,'cpu_pool_limit':8,'memory_limit_bytes':256*GIB,
                'cpu_backfill':{'enabled':True,'haplotag_phases':True,'max_loan_pools':1,'max_claim_pools':2}}
        atomic_json(root/'policy.json',policy)
        row={'branch':'WGS','stage':'06','sample':'S','cpus':1,'total_cpus':8,
             'memory_bytes':128*GIB,'manager_pid':os.getpid()}
        if loan:row['cpu_backfill']={'cpus':6,'nice':19}
        path=root/'leases/ours.json';atomic_json(path,row)
        return root,path,policy

    def test_결합후RAM반환_index재확보_QC반환(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,_=self.fixture(raw)
            atomic_json(root/'leases/other.json',{'cpus':6,'memory_bytes':64*GIB})
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)}):
                self.assertEqual(reserve_tail(SimpleNamespace(jobs=1),'S','combine_and_tags'),1)
                row=json.loads(path.read_text());self.assertEqual(row['memory_bytes'],4*GIB)
                self.assertEqual(row['haplotag_initial_memory_bytes'],128*GIB)
                self.assertEqual(reserve_tail(SimpleNamespace(jobs=1),'S','output_index',5,64<<20),2)
                row=json.loads(path.read_text());self.assertEqual(row['memory_bytes'],4*GIB+2*(64<<20))
                self.assertEqual(reserve_tail(SimpleNamespace(jobs=1),'S','tag_counts'),1)
                self.assertEqual(json.loads(path.read_text())['memory_bytes'],4*GIB)

    def test_최초대여와전체논리상한으로만재확보(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,_=self.fixture(raw,loan=True)
            atomic_json(root/'leases/regular.json',{'cpus':8,'memory_bytes':32*GIB})
            atomic_json(root/'leases/loan.json',{'cpus':5,'memory_bytes':4*GIB,'cpu_backfill':{'cpus':5,'nice':19}})
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)}):
                self.assertEqual(reserve_tail(SimpleNamespace(jobs=1),'S','output_index',7,64<<20),3)
            row=json.loads(path.read_text());self.assertLessEqual(row['cpus'],row['cpu_backfill']['cpus'])
            self.assertEqual(row['cpu_backfill']['nice'],19)
            self.assertLessEqual(row['cpus']+5,8)
            self.assertLessEqual(row['cpus']+5+8,16)

    def test_대여정책해제시기존한CPU이상늘리지않음(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,policy=self.fixture(raw,loan=True)
            policy['cpu_backfill']['haplotag_phases']=False;atomic_json(root/'policy.json',policy)
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)}):
                self.assertEqual(reserve_tail(SimpleNamespace(jobs=1),'S','output_index',5,64<<20),1)

    def test_RAM부족시계산없이기다렸다재확보(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,policy=self.fixture(raw);policy['memory_limit_bytes']=100*GIB;atomic_json(root/'policy.json',policy)
            row=json.loads(path.read_text());row['memory_bytes']=4*GIB;atomic_json(path,row)
            other=root/'leases/other.json';atomic_json(other,{'cpus':2,'memory_bytes':86*GIB})
            result=[]
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)}):
                thread=threading.Thread(target=lambda:result.append(reserve_tail(SimpleNamespace(jobs=1),'S','output_index',5,64<<20)))
                thread.start()
                try:
                    time.sleep(.15);self.assertTrue(thread.is_alive());self.assertEqual(json.loads(path.read_text())['memory_bytes'],4*GIB)
                finally:
                    atomic_json(other,{'cpus':2,'memory_bytes':85*GIB});thread.join(3)
                self.assertFalse(thread.is_alive());self.assertEqual(result,[5])
                self.assertLessEqual(json.loads(path.read_text())['memory_bytes']+85*GIB,90*GIB)

    def test_다중샘플공유예약과다른샘플은변경하지않음(self):
        with tempfile.TemporaryDirectory() as raw:
            _,path,_=self.fixture(raw);before=path.read_bytes()
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)}):
                self.assertIsNone(reserve_tail(SimpleNamespace(jobs=2),'S','combine_and_tags'))
                self.assertIsNone(reserve_tail(SimpleNamespace(jobs=1),'OTHER','combine_and_tags'))
            self.assertEqual(path.read_bytes(),before)

    def test_본인예약누락은오류임(self):
        with tempfile.TemporaryDirectory() as raw:
            _,path,_=self.fixture(raw);path.unlink()
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)}):
                with self.assertRaises(FileNotFoundError):reserve_tail(SimpleNamespace(jobs=1),'S','output_index',5)


if __name__=='__main__':unittest.main()
