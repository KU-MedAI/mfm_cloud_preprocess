"""실행 가능한 작업을 예약 때문에 놀리지 않는 공유 CPU/RAM 배분 검증."""
import json,os,tempfile,threading,time,unittest
from pathlib import Path
from unittest.mock import patch
from integrate_wgs_rna.memory_policy import GIB
from integrate_wgs_rna.compact.resources import reserve_shared,fit_memory
from integrate_wgs_rna.cpu_executor import WorkerBudget

class SharedCPUTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        for name in ['leases','waiters']:(self.root/name).mkdir()
        self.own=self.root/'leases/own.json'
        self.own.write_text(json.dumps(dict(cpus=1,total_cpus=68,memory_bytes=4*GIB)))
        (self.root/'policy.json').write_text(json.dumps(dict(live_limits=True,cpu_pool_limit=68,memory_limit_bytes=512*GIB)))
        self.env={'MFM_CPU_POOL_LEASE_FILE':str(self.own),'MFM_CPU_BUDGET':'68',
                  'MFM_CPU_POOL_GRANTED':'1','MFM_RAM_LIMIT_BYTES':str(512*GIB),
                  'MFM_CPU_UTILIZATION':'maximum','MFM_COMPACT_JOBS_REQUESTED':'auto',
                  'MFM_BLOSC_THREADS':'1'}
        for context in [patch.dict(os.environ,self.env,clear=True),
                        patch('os.sched_getaffinity',return_value=set(range(72))),
                        patch('integrate_wgs_rna.pool_leases.cgroup_cpu_limit',return_value=68),
                        patch('integrate_wgs_rna.compact.resources.cgroup_cpu_limit',return_value=68),
                        patch('integrate_wgs_rna.memory_policy.allocation_memory_limit',return_value=None)]:
            context.start();self.addCleanup(context.stop)

    def forecast(self):
        return dict(worker_ceiling_bytes=32*GIB,worker_reservation_bytes=8*GIB,
                    estimated_peak_per_worker_bytes=4*GIB,memory_limit_bytes=512*GIB)

    def test_hifi_mapping_balance_assigns_more_real_mapping_threads_within_68(self):
        from types import SimpleNamespace
        from integrate_wgs_rna.realign.align_resources import alignment_plan
        config=SimpleNamespace(platform='hifi',threads=30,sort_threads=4,hap_jobs=2,sort_mem='2G')
        env={'MFM_WGS_ALIGN_REQUESTED':'auto','MFM_WGS_SORT_REQUESTED':'auto','MFM_WGS_HAP_REQUESTED':'auto'}
        lengths={1:2_876_676_610,2:2_876_520_048}
        old=alignment_plan(config,lengths,68,1485000000000,dict(env,MFM_WGS_CPU_BALANCE='legacy'))
        new=alignment_plan(config,lengths,68,1485000000000,env)
        self.assertEqual((old['threads'],old['sort_threads'],old['input_decode_threads']),(17,8,4))
        self.assertEqual((new['threads'],new['sort_threads'],new['input_decode_threads']),(26,2,1))
        previous=alignment_plan(config,lengths,68,1485000000000,dict(env,MFM_WGS_CPU_BALANCE='mapping'))
        self.assertEqual((previous['threads'],previous['sort_threads'],previous['input_decode_threads']),(23,4,2))
        self.assertEqual(new['reserved_cpus'],68)
        self.assertGreater(new['threads']*new['hap_jobs'],old['threads']*old['hap_jobs'])
        config.platform='ont'
        self.assertEqual(alignment_plan(config,lengths,68,1485000000000,env)['cpu_balance'],'legacy')
        with self.assertRaises(ValueError):
            alignment_plan(config,lengths,68,1485000000000,dict(env,MFM_WGS_CPU_BALANCE='bad'))

    def test_small_initial_grant_does_not_cap_actual_reacquired_pool(self):
        forecast=self.forecast();execution=reserve_shared(1,22,forecast,hap_split=True)
        own=json.loads(self.own.read_text())
        self.assertEqual((execution['workers'],execution['cpu_per_worker']),(17,4))
        self.assertEqual(own['cpus'],68)
        self.assertEqual(own['memory_bytes'],4*GIB+17*forecast['worker_reservation_bytes'])

    def test_smaller_slurm_allocation_stays_above_stale_grant_but_below_six(self):
        with patch.dict(os.environ,SLURM_CPUS_PER_TASK='6'):
            execution=reserve_shared(1,22,self.forecast(),hap_split=True)
        self.assertLessEqual(json.loads(self.own.read_text())['cpus'],6)
        self.assertGreater(execution['workers']*execution['cpu_per_worker'],1)

    def test_same_stage_waiter_blocked_by_concurrency_does_not_deadlock_active_owner(self):
        policy=json.loads((self.root/'policy.json').read_text());policy['stage_concurrency']={'WGS:12':1}
        (self.root/'policy.json').write_text(json.dumps(policy))
        own=json.loads(self.own.read_text());own.update(branch='WGS',stage='12');self.own.write_text(json.dumps(own))
        (self.root/'waiters/next.json').write_text(json.dumps(dict(branch='WGS',stage='12',min_cpus=68,memory_bytes=100*GIB)))
        reserve_shared(1,22,self.forecast(),hap_split=True)
        self.assertEqual(json.loads(self.own.read_text())['cpus'],68)

    def test_dead_waiter_does_not_block_entry_or_executor_growth(self):
        waiter=self.root/'waiters/dead.json';waiter.write_text(json.dumps(dict(pid=999999999,min_cpus=68,memory_bytes=100*GIB)))
        reserve_shared(1,22,self.forecast(),hap_split=True)
        self.assertEqual(json.loads(self.own.read_text())['cpus'],68)
        own=json.loads(self.own.read_text());own['cpus']=1;self.own.write_text(json.dumps(own))
        self.assertEqual(WorkerBudget(1,17,True,8*GIB,cpus_per_worker=4).capacity(22),17)
        self.assertTrue(waiter.exists())

    def test_other_leases_limit_cpu_and_ram_in_the_same_atomic_plan(self):
        other=self.root/'leases/other.json'
        original=json.dumps(dict(cpus=40,memory_bytes=400*GIB));other.write_text(original)
        forecast=self.forecast();execution=reserve_shared(1,22,forecast,hap_split=True)
        own=json.loads(self.own.read_text())
        self.assertLessEqual(own['cpus']+40,68)
        self.assertLessEqual(own['memory_bytes']+400*GIB,512*GIB*9//10)
        self.assertEqual(execution['workers'],3)
        self.assertGreater(execution['alignment_threads'],2)
        self.assertEqual(other.read_text(),original)

    def test_fitting_waiter_keeps_its_cpu_and_ram_room(self):
        (self.root/'waiters/first.json').write_text(json.dumps(dict(min_cpus=60,memory_bytes=100*GIB)))
        reserve_shared(1,22,self.forecast(),hap_split=True)
        own=json.loads(self.own.read_text())
        self.assertLessEqual(own['cpus'],8)
        self.assertLessEqual(own['memory_bytes']+100*GIB,512*GIB*9//10)
        self.assertEqual(own['compact_entry_allocation']['priority_waiters'],1)

    def test_unfittable_waiter_does_not_block_ready_compact_work(self):
        (self.root/'waiters/large.json').write_text(json.dumps(dict(min_cpus=68,memory_bytes=1024*GIB)))
        reserve_shared(1,22,self.forecast(),hap_split=True)
        self.assertEqual(json.loads(self.own.read_text())['cpus'],68)

    def test_no_compute_wait_returns_cpu_and_replans_after_release(self):
        other=self.root/'leases/other.json';other.write_text(json.dumps(dict(cpus=68,memory_bytes=100*GIB)))
        results=[];errors=[]
        def run():
            try:results.append(reserve_shared(1,22,self.forecast(),hap_split=True))
            except BaseException as exc:errors.append(exc)
        worker=threading.Thread(target=run);worker.start()
        try:
            until=time.monotonic()+3
            while time.monotonic()<until:
                own=json.loads(self.own.read_text())
                if own.get('compact_entry_allocation',{}).get('waiting'):break
                time.sleep(.01)
            self.assertTrue(worker.is_alive())
            self.assertEqual(own['cpus'],0)
            self.assertGreaterEqual(own['memory_bytes'],4*GIB)
        finally:
            other.unlink(missing_ok=True);worker.join(5)
        self.assertFalse(worker.is_alive());self.assertFalse(errors)
        self.assertEqual(json.loads(self.own.read_text())['cpus'],68)

    def test_memory_pressure_uses_sequential_haps_only_for_auto(self):
        forecast=self.forecast();available=18*GIB
        self.assertEqual(fit_memory(22,22,forecast,hap_split=True,
                                   available_cpus=68,available_memory=available)['hap_jobs'],1)
        self.assertLessEqual(forecast['worker_reservation_bytes']+4*GIB,available)
        with patch.dict(os.environ,MFM_COMPACT_HAP_JOBS='2'):
            with self.assertRaises(MemoryError):
                fit_memory(22,22,self.forecast(),hap_split=True,
                           available_cpus=68,available_memory=available)

    def test_ram_waiting_worker_releases_cpu_for_another_ready_stage(self):
        (self.root/'policy.json').write_text(json.dumps(dict(live_limits=True,cpu_pool_limit=8,memory_limit_bytes=64*GIB)))
        self.own.write_text(json.dumps(dict(cpus=8,total_cpus=8,memory_bytes=4*GIB)))
        (self.root/'leases/holder.json').write_text(json.dumps(dict(cpus=0,memory_bytes=40*GIB)))
        ready=self.root/'leases/ready.json';ready.write_text(json.dumps(dict(cpus=0,total_cpus=8,memory_bytes=4*GIB)))
        blocked=WorkerBudget(8,8,True,16*GIB)
        self.assertEqual(blocked.capacity(8),0)
        self.assertEqual(json.loads(self.own.read_text())['cpus'],0)
        with patch.dict(os.environ,MFM_CPU_POOL_LEASE_FILE=str(ready)):
            compute=WorkerBudget(2,2,True,GIB,cpus_per_worker=2)
            self.assertEqual(compute.capacity(2),2)
        self.assertEqual(json.loads(ready.read_text())['cpus'],4)
        # 살아있는 계산과 명시적 부모 CPU는 반납하지 않는다.
        protected=WorkerBudget(2,2,True,16*GIB,cpus_per_worker=2,parent_cpus=1)
        self.assertEqual(protected.capacity(2,running=1),1)
        self.assertEqual(json.loads(self.own.read_text())['cpus'],3)

if __name__=='__main__':unittest.main()
