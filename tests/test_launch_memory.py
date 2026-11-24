"""공유풀의이전버전RAM예약을실행직전에확정하는보호를검증함."""
from pathlib import Path
from unittest.mock import patch
import importlib.util,json,os,tempfile,threading,time,unittest
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('launch_pool',ROOT/'scripts/cpu_pool.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
GIB=1024**3
class GuardTests(unittest.TestCase):
    def fixture(self,raw):
        root=Path(raw);(root/'leases').mkdir();(root/'waiters').mkdir();owner=root/'leases/owner.json';module.atomic_json(owner,{'cpus':8,'total_cpus':16,'memory_bytes':8*GIB,'manager_pid':os.getpid(),'child_pid':None});return root,owner
    def test_동일예약은그대로사용(self):
        with tempfile.TemporaryDirectory() as raw:
            root,owner=self.fixture(raw);before=owner.read_bytes();row=module.reconcile_launch_memory(root,owner,8*GIB,1);self.assertEqual(row['cpus'],8);self.assertEqual(owner.read_bytes(),before)
    def test_과소예약을계산시작전수정(self):
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{'MFM_RAM_LIMIT_BYTES':str(100*GIB)}):
            root,owner=self.fixture(raw);row=module.reconcile_launch_memory(root,owner,40*GIB,2);self.assertEqual(row['memory_bytes'],40*GIB);self.assertEqual(row['initial_allocator_memory_bytes'],8*GIB);self.assertTrue(row['memory_plan_reconciled'])
    def test_여유없으면예약반환후기다리고재확정(self):
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{'MFM_RAM_LIMIT_BYTES':str(100*GIB)}),patch.object(module,'POLL_SECONDS',.01):
            root,owner=self.fixture(raw);other=root/'leases/other.json';module.atomic_json(other,{'cpus':8,'memory_bytes':80*GIB,'manager_pid':os.getpid(),'child_pid':None});result=[]
            t=threading.Thread(target=lambda:result.append(module.reconcile_launch_memory(root,owner,40*GIB,2)));t.start()
            try:
                deadline=time.monotonic()+2
                while time.monotonic()<deadline:
                    row=module.read_json(owner)
                    if row.get('memory_plan_reconciled') is False:break
                    time.sleep(.01)
                self.assertEqual(row['cpus'],0);self.assertEqual(row['memory_bytes'],4*GIB);self.assertTrue(t.is_alive())
            finally:other.unlink();t.join(3)
            self.assertFalse(t.is_alive());self.assertEqual(result[0]['cpus'],8);self.assertEqual(result[0]['memory_bytes'],40*GIB)
    def test_이미시작한도구는보정거부(self):
        with tempfile.TemporaryDirectory() as raw:
            root,owner=self.fixture(raw);row=module.read_json(owner);row['child_pid']=os.getpid();module.atomic_json(owner,row)
            with self.assertRaisesRegex(ValueError,'실행중'):module.reconcile_launch_memory(root,owner,40*GIB,2)
    def test_RAM동일해도구할당CPU를시작전에재확정(self):
        with tempfile.TemporaryDirectory() as raw,patch('integrate_wgs_rna.pool_leases.cgroup_cpu_limit',return_value=None),patch('os.sched_getaffinity',return_value=set(range(72))):
            root,owner=self.fixture(raw);row=module.read_json(owner);row['total_cpus']=68;module.atomic_json(owner,row)
            module.atomic_json(root/'policy.json',{'live_limits':True,'cpu_pool_limit':40})
            other=root/'leases/other.json';module.atomic_json(other,{'cpus':36,'memory_bytes':8*GIB,'manager_pid':os.getpid()})
            before=other.read_bytes();row=module.reconcile_launch_memory(root,owner,8*GIB,2)
            self.assertEqual((row['cpus'],row['total_cpus']),(4,40));self.assertEqual(other.read_bytes(),before)
    def test_CPU대기는0코어를반환해재확보교착을막음(self):
        with tempfile.TemporaryDirectory() as raw,patch('integrate_wgs_rna.pool_leases.cgroup_cpu_limit',return_value=None),patch('os.sched_getaffinity',return_value=set(range(72))):
            root,owner=self.fixture(raw);row=module.read_json(owner);row['total_cpus']=68;module.atomic_json(owner,row)
            module.atomic_json(root/'policy.json',{'live_limits':True,'cpu_pool_limit':40})
            other=root/'leases/other.json';module.atomic_json(other,{'cpus':40,'memory_bytes':8*GIB,'manager_pid':os.getpid()})
            def release(_):
                waiting=module.read_json(owner);self.assertEqual(waiting['cpus'],0)
                self.assertEqual(waiting['total_cpus'],40);self.assertEqual(waiting['memory_bytes'],4*GIB)
                value=module.read_json(other);value['cpus']=36;module.atomic_json(other,value)
            with patch.object(module.time,'sleep',side_effect=release) as sleep:
                row=module.reconcile_launch_memory(root,owner,8*GIB,2)
            sleep.assert_called_once();self.assertEqual(row['cpus'],4)
    def test_시작전취소된유휴대여는정상40코어안에서재확정(self):
        with tempfile.TemporaryDirectory() as raw,patch('integrate_wgs_rna.pool_leases.cgroup_cpu_limit',return_value=None),patch('os.sched_getaffinity',return_value=set(range(72))):
            root,owner=self.fixture(raw);row=module.read_json(owner);row.update(total_cpus=68,cpu_backfill={'cpus':8});module.atomic_json(owner,row)
            module.atomic_json(root/'policy.json',{'live_limits':True,'cpu_pool_limit':40,'cpu_backfill':{'enabled':False}})
            module.atomic_json(root/'leases/other.json',{'cpus':38,'memory_bytes':8*GIB,'manager_pid':os.getpid()})
            row=module.reconcile_launch_memory(root,owner,8*GIB,2)
            self.assertEqual((row['cpus'],row['total_cpus']),(2,40));self.assertNotIn('cpu_backfill',row)
            self.assertEqual(row['revoked_cpu_backfill'],{'cpus':8})
if __name__=='__main__':unittest.main()
