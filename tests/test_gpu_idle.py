"""GPU 분산·외부 작업 보호·해독 예산·작업 종료시 context 회수를 검증함."""
from pathlib import Path
from unittest.mock import patch
import json,os,tempfile,unittest,multiprocessing
from integrate_wgs_rna import gpu_pool
from integrate_wgs_rna.compact.resources import plan
from integrate_wgs_rna.cpu_executor import cpu_results


def pid_task(value):return os.getpid(),value


class IdleTests(unittest.TestCase):
    def test_22개worker의순차HAP해독을66CPU로배분함(self):
        with patch.dict(os.environ,{'MFM_CPU_BUDGET':'68','MFM_ZARR_BACKEND':'cupy','MFM_BLOSC_THREADS':'1'},clear=True),patch('os.sched_getaffinity',return_value=set(range(72))),patch('integrate_wgs_rna.compact.resources.cgroup_cpu_limit',return_value=68):
            p=plan(22,22);self.assertEqual(p['workers'],22);self.assertEqual(p['alignment_threads'],3);self.assertEqual(p['cpu_per_worker']*p['workers'],66);self.assertEqual(p['extra_decode_buffers_bytes'],64*(1<<20));self.assertEqual(p['recycle_tasks'],1)
    def test_CPU도해독병렬화하고worker재사용은유지함(self):
        with patch.dict(os.environ,{'MFM_CPU_BUDGET':'68','MFM_BLOSC_THREADS':'1'},clear=True),patch('os.sched_getaffinity',return_value=set(range(72))),patch('integrate_wgs_rna.compact.resources.cgroup_cpu_limit',return_value=68):
            p=plan(22,22);self.assertEqual(p['alignment_threads'],3);self.assertEqual(p['workers']*p['cpu_per_worker'],66);self.assertIsNone(p['recycle_tasks'])
    def test_작은grant와명시thread도CPU한도를지킴(self):
        for total in [1,2,4,16,68]:
            for requested in [1,2,22,44]:
                with self.subTest(total=total,workers=requested),patch.dict(os.environ,{'MFM_CPU_POOL_GRANTED':str(total),'MFM_CPU_BUDGET':'68'},clear=True),patch('os.sched_getaffinity',return_value=set(range(72))),patch('integrate_wgs_rna.compact.resources.cgroup_cpu_limit',return_value=68):
                    p=plan(requested,44);self.assertLessEqual(p['workers']*p['cpu_per_worker'],total)
        with patch.dict(os.environ,{'MFM_CPU_BUDGET':'68','MFM_COMPACT_ALIGNMENT_THREADS':'4','MFM_BLOSC_THREADS':'2'},clear=True),patch('os.sched_getaffinity',return_value=set(range(72))),patch('integrate_wgs_rna.compact.resources.cgroup_cpu_limit',return_value=68):
            p=plan(22,22);self.assertEqual(p['alignment_threads'],4);self.assertEqual(p['cpu_per_worker'],5);self.assertEqual(p['workers'],13)
        with patch.dict(os.environ,{'MFM_CPU_POOL_GRANTED':'2','MFM_COMPACT_ALIGNMENT_THREADS':'4'},clear=True):
            with self.assertRaises(ValueError):plan(1,1)
    def test_작은공유풀은호스트보다우선하며worker를굶기지않음(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'leases/own.json';path.parent.mkdir();path.write_text(json.dumps({'cpus':1,'total_cpus':4}))
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path),'MFM_CPU_BUDGET':'68'},clear=True),patch('os.sched_getaffinity',return_value=set(range(72))),patch('integrate_wgs_rna.compact.resources.cgroup_cpu_limit',return_value=68):
                p=plan(2,2);self.assertEqual(p['cpu_limit'],4);self.assertEqual(p['workers'],2);self.assertEqual(p['cpu_per_worker'],2)
    def test_드라이버잔여5MiB는PID없는경우만사용함(self):
        with tempfile.TemporaryDirectory() as tmp,patch.dict(os.environ,{'MFM_GPU_POOL_DIR':tmp,'MFM_GPU_DEVICES':'0,1'},clear=True):
            devices=[dict(index=i,uuid=f'GPU-{i}',pci=f'0000:0{i}:00.0',total=32*gpu_pool.GIB,used=5*(1<<20),utilization=0,compute_pids=[999] if i==0 else []) for i in range(2)]
            with patch.object(gpu_pool,'inventory',return_value=devices):
                lease=gpu_pool.acquire(gpu_pool.GIB,purpose='test',timeout=0);self.assertEqual(lease.device['index'],1);lease.close()
                devices[1]['compute_pids']=[888]
                devices[1]['used']=0  # 메모리 표본이 0이어도 다른 PID가 있으면 보호함.
                with self.assertRaises(TimeoutError):gpu_pool.acquire(gpu_pool.GIB,purpose='test',timeout=0)
    def test_공유worker4개를GPU4개에분산함(self):
        with tempfile.TemporaryDirectory() as tmp,patch.dict(os.environ,{'MFM_GPU_POOL_DIR':tmp,'MFM_GPU_DEVICES':'0,1,2,3'},clear=True):
            devices=[dict(index=i,uuid=f'GPU-{i}',pci=f'0000:0{i}:00.0',total=32*gpu_pool.GIB,used=0,utilization=0) for i in range(4)]
            with patch.object(gpu_pool,'inventory',return_value=devices):
                leases=[gpu_pool.acquire(gpu_pool.GIB,purpose='test',timeout=0) for _ in range(4)]
                try:self.assertEqual([x.device['index'] for x in leases],[0,1,2,3])
                finally:
                    for x in leases:x.close()
    def test_spawn작업완료후같은PID를재사용하지않음(self):
        with patch.dict(os.environ,{},clear=True):
            rows=list(cpu_results(pid_task,[(i,) for i in range(3)],1,mp_context=multiprocessing.get_context('spawn'),max_tasks_per_child=1))
        self.assertTrue(all(e is None for _,_,e in rows));self.assertEqual(len({r[0] for _,r,_ in rows}),3)


if __name__=='__main__':unittest.main()
