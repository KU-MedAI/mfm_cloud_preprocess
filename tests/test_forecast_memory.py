"""예측만 하는 대기와 실제 계산·보존 중 메모리 계약을 구분해 검증함."""
import importlib.util,json,os,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from integrate_wgs_rna.cpu_executor import WorkerBudget
from integrate_wgs_rna.memory_policy import GIB

class ForecastMemoryTest(unittest.TestCase):
    def test_wait_releases_forecast_and_reacquires_without_releasing_running_workers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'leases').mkdir();(root/'waiters').mkdir()
            own=root/'leases/own.json';other=root/'leases/other.json'
            own.write_text(json.dumps(dict(cpus=1,total_cpus=8,memory_bytes=64*GIB)))
            other.write_text(json.dumps(dict(cpus=1,memory_bytes=48*GIB)))
            with patch.dict(os.environ,MFM_CPU_POOL_LEASE_FILE=str(own),MFM_RAM_LIMIT_BYTES=str(80*GIB)):
                legacy=WorkerBudget(1,8,True);legacy.base_memory=64*GIB
                self.assertEqual(0,legacy.capacity(8))
                self.assertEqual(64*GIB,json.loads(own.read_text())['memory_bytes'])
                budget=WorkerBudget(1,8,True,memory_per_worker=GIB,wait_memory_bytes=4*GIB);budget.base_memory=64*GIB
                self.assertEqual(0,budget.capacity(8))
                self.assertEqual(4*GIB,json.loads(own.read_text())['memory_bytes'])
                self.assertEqual(2,budget.capacity(8,running=2))
                self.assertGreaterEqual(json.loads(own.read_text())['memory_bytes'],66*GIB)
                other.unlink();self.assertEqual(8,budget.capacity(8))
                self.assertEqual(72*GIB,json.loads(own.read_text())['memory_bytes'])
                self.assertEqual(0,budget.capacity(0))
                self.assertEqual(4*GIB,json.loads(own.read_text())['memory_bytes'])

    def test_wait_preserves_actual_parent_rss_and_explicit_buffer(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'leases').mkdir();(root/'waiters').mkdir();own=root/'leases/own.json'
            own.write_text(json.dumps(dict(cpus=1,total_cpus=8,memory_bytes=64*GIB)))
            with patch.dict(os.environ,MFM_CPU_POOL_LEASE_FILE=str(own),MFM_RAM_LIMIT_BYTES=str(32*GIB)):
                budget=WorkerBudget(1,8,True,wait_memory_bytes=4*GIB);budget.base_memory=64*GIB
                read=Path.read_text
                def statm(path,*a,**kw):
                    return '0 '+str(8*GIB//os.sysconf('SC_PAGE_SIZE')) if str(path)=='/proc/self/statm' else read(path,*a,**kw)
                with patch.object(Path,'read_text',statm):self.assertEqual(0,budget.capacity(8))
                self.assertEqual((8*GIB*11+9)//10,json.loads(own.read_text())['memory_bytes'])
                budget.wait_memory_bytes=12*GIB;budget.capacity(8)
                self.assertEqual(12*GIB,json.loads(own.read_text())['memory_bytes'])
                with self.assertRaises(ValueError):WorkerBudget(1,1,True,wait_memory_bytes=-1)

    def test_bench_entry_override_never_lowers_real_pipeline_entry(self):
        import integrate_wgs_rna
        root=Path(integrate_wgs_rna.__file__).resolve().parents[2]
        spec=importlib.util.spec_from_file_location('forecast_pool',root/'scripts/cpu_pool.py')
        pool=importlib.util.module_from_spec(spec);spec.loader.exec_module(pool)
        with patch.dict(os.environ):
            os.environ.pop('MFM_LONGREAD_MEMORY_PLAN',None);os.environ.pop('MFM_BENCH_ENTRY_MEMORY_BYTES',None)
            expected=pool.stage_memory('RNA','22');self.assertEqual(64*GIB,pool.stage_memory('BENCH','rna_fused'))
            os.environ['MFM_BENCH_ENTRY_MEMORY_BYTES']=str(4*GIB)
            self.assertEqual(4*GIB,pool.stage_memory('BENCH','rna_fused'))
            self.assertEqual(expected,pool.stage_memory('RNA','22'))
            os.environ['MFM_BENCH_ENTRY_MEMORY_BYTES']='1'
            with self.assertRaises(ValueError):pool.stage_memory('BENCH','rna_fused')
            self.assertEqual(expected,pool.stage_memory('RNA','22'))

if __name__=='__main__':unittest.main()
