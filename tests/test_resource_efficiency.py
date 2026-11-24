"""큰 요청의 기아 방지·RAM 안전성·고정 CPU 범위와 관측 시각을 검증함."""
from pathlib import Path
from unittest.mock import patch
import importlib.util,json,os,tempfile,time,unittest
from integrate_wgs_rna.resource_fairness import priority_request,fair_cpu_capacity,request_key
from integrate_wgs_rna.pool_cpu_observer import PreprocessingCpuObserver
from integrate_wgs_rna.pool_leases import descendants

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('efficiency_pool',ROOT/'scripts/cpu_pool.py')
POOL=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(POOL)

class EfficiencyTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        for d in ['leases','waiters']:(self.root/d).mkdir()
        self.addCleanup(patch.stopall)
        patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=1000).start()
        patch.dict(os.environ,{'MFM_FIXED_SAMPLE_TOKEN':''}).start()
        self.now=time.time()
    def pending(self,sample='A',queued=None,cpus=9,ram=100):
        return dict(manager_pid=os.getpid(),sample=sample,stage='10',cpus=0,memory_bytes=4,
            acquired_at=self.now-100,realign_wait=dict(minimum_cpus=cpus,memory_bytes=ram,
                queued_at=self.now-100 if queued is None else queued))
    def test_old_internal_collects_cpus_and_own_request_can_start(self):
        own=self.pending();other=dict(cpus=60,memory_bytes=100)
        self.assertEqual(0,fair_cpu_capacity(self.root,68,[own,other],8))
        self.assertEqual(8,fair_cpu_capacity(self.root,68,[own,other],8,request=own))
        self.assertEqual(3,fair_cpu_capacity(self.root,68,[own,other],12))
    def test_insufficient_ram_does_not_block_small_work(self):
        self.assertIsNone(priority_request(self.root,68,[self.pending(ram=950),dict(cpus=60,memory_bytes=100)]))
        self.assertEqual(8,fair_cpu_capacity(self.root,68,[self.pending(ram=950),dict(cpus=60,memory_bytes=100)],8))
    def test_oldest_of_two_internal_requests_avoids_mutual_deadlock(self):
        a=self.pending('A',queued=self.now-100);b=self.pending('B',queued=self.now-90)
        self.assertEqual(10,fair_cpu_capacity(self.root,68,[a,b],10,request=a))
        self.assertEqual(1,fair_cpu_capacity(self.root,68,[a,b],10,request=b))
    def test_new_or_oversized_request_does_not_hold_cpus(self):
        for own in [self.pending(queued=self.now-1),self.pending(cpus=69)]:
            self.assertIsNone(priority_request(self.root,68,[own]))
    def test_fixed_sample_does_not_wait_for_other_sample_cpu(self):
        own=self.pending();own['fixed_sample_token']='B'
        with patch.dict(os.environ,{'MFM_FIXED_SAMPLE_TOKEN':'A'}):
            self.assertEqual(4,fair_cpu_capacity(self.root,4,[own],4))
    def test_old_regular_stage_also_collects_minimum_cpus(self):
        row=dict(pid=os.getpid(),sample='A',min_cpus=9,memory_bytes=100,queued_at=self.now-100)
        p=priority_request(self.root,68,[],waiters=[(Path('A'),row)])
        self.assertEqual('stage',p['kind']);self.assertEqual(0,fair_cpu_capacity(self.root,68,[],4,waiters=[(Path('A'),row)]))
    def test_allocator_stops_small_arrivals_until_internal_minimum_is_free(self):
        (self.root/'policy.json').write_text(json.dumps({'cpu_backfill':{'max_extra_cpus':16}}))
        own=self.pending();(self.root/'leases/long.json').write_text(json.dumps(own))
        (self.root/'leases/busy.json').write_text(json.dumps(dict(cpus=62,memory_bytes=100)))
        small=self.root/'waiters/small.json'
        small.write_text(json.dumps(dict(pid=os.getpid(),sample='S',branch='WGS',stage='03',min_cpus=2,max_cpus=2,memory_bytes=10,queued_at=self.now-1)))
        with patch.object(POOL,'pool_cpu_limit',return_value=68):
            POOL.allocate(self.root,68);self.assertTrue(small.exists())
            busy=self.root/'leases/busy.json';busy.write_text(json.dumps(dict(cpus=58,memory_bytes=100)))
            POOL.allocate(self.root,68);self.assertTrue(small.exists())
        self.assertEqual(10,fair_cpu_capacity(self.root,68,[own,dict(cpus=58,memory_bytes=100)],10,request=own))
    def test_mapper_exit_releases_only_ended_cpu_and_keeps_sort_and_other_hap_ram(self):
        from integrate_wgs_rna.realign.stream_resources import mapper_finished
        path=self.root/'leases/own.json'
        row=dict(cpus=44,memory_bytes=200,realign_phase='align_restore_sort_qc',
            hap_resources={'1':dict(cpus=22,memory_bytes=100,phase='align_restore_sort_qc'),
                           '2':dict(cpus=22,memory_bytes=100,phase='align_restore_sort_qc')})
        path.write_text(json.dumps(row));event=mapper_finished(1,2,path)
        after=json.loads(path.read_text())
        self.assertEqual((event['cpus_before'],event['cpus_after']),(44,26))
        self.assertEqual(after['memory_bytes'],200)
        self.assertEqual(after['hap_resources']['2'],row['hap_resources']['2'])
        self.assertEqual(after['hap_resources']['1']['memory_bytes'],100)
        self.assertEqual(mapper_finished(1,2,path)['cpus_after'],26)

    def test_idle_live_workers_keep_memory_until_join(self):
        from integrate_wgs_rna.cpu_executor import WorkerBudget
        gib=1024**3;path=self.root/'leases/own.json'
        path.write_text(json.dumps(dict(cpus=4,total_cpus=4,memory_bytes=8*gib)))
        with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path),'MFM_RAM_LIMIT_BYTES':str(200*gib)}),patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=180*gib):
            budget=WorkerBudget(4,4,True,memory_per_worker=gib)
            budget.resident_workers=3
            self.assertEqual(1,budget.capacity(1,grow=False))
            self.assertEqual(7*gib,json.loads(path.read_text())['memory_bytes'])
            budget.resident_workers=0
            budget.capacity(0,grow=False)
            self.assertEqual(4*gib,json.loads(path.read_text())['memory_bytes'])

    def test_only_live_same_host_native_allocator_suppresses_redundant_allocation(self):
        import socket
        from integrate_wgs_rna.pool_leases import native_allocator_ready
        raw=Path(f'/proc/{os.getpid()}/stat').read_text();start=int(raw[raw.rfind(')')+2:].split()[19])
        record=dict(schema='mfm_resource_watchdog_v1',updated=time.time(),
                    watchdog_owner=dict(host=socket.gethostname(),pid=os.getpid(),started=start))
        path=self.root/'watchdog.json';path.write_text(json.dumps(record))
        self.assertTrue(native_allocator_ready(self.root))
        with patch.dict(os.environ,{'MFM_FIXED_SAMPLE_TOKEN':'A'}):self.assertFalse(native_allocator_ready(self.root))
        record['watchdog_owner']['started']+=1;path.write_text(json.dumps(record));self.assertFalse(native_allocator_ready(self.root))
        record['watchdog_owner']['started']=start;record['updated']=time.time()-11
        path.write_text(json.dumps(record));self.assertFalse(native_allocator_ready(self.root))

    def test_identical_worker_claim_is_not_written_again_but_changes_are_immediate(self):
        from integrate_wgs_rna.cpu_executor import WorkerBudget
        gib=1024**3;path=self.root/'leases/own.json'
        path.write_text(json.dumps(dict(cpus=4,total_cpus=4,memory_bytes=4*gib)))
        with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path),'MFM_RAM_LIMIT_BYTES':str(200*gib)}),patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=180*gib):
            budget=WorkerBudget(4,4,True)
            with patch('integrate_wgs_rna.cpu_executor.os.replace',wraps=os.replace) as write:
                budget.capacity(4);budget.capacity(4)
                self.assertEqual(1,write.call_count)
                budget.capacity(2)
                self.assertEqual(2,write.call_count)
                self.assertEqual(2,json.loads(path.read_text())['cpus'])

    def test_unknown_running_count_keeps_held_cpus_but_confirmed_idle_can_yield(self):
        from integrate_wgs_rna.cpu_executor import WorkerBudget
        gib=1024**3;path=self.root/'leases/own.json'
        row=dict(cpus=4,total_cpus=4,memory_bytes=4*gib,manager_pid=os.getpid(),sample='B')
        path.write_text(json.dumps(row));waiting=self.pending('A',cpus=3)
        (self.root/'leases/other.json').write_text(json.dumps(waiting))
        with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path),'MFM_RAM_LIMIT_BYTES':str(200*gib)}),patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=180*gib):
            budget=WorkerBudget(4,4,True)
            self.assertEqual(4,budget.capacity(4))
            self.assertEqual(1,budget.capacity(4,running=0))

    def test_sequence_vcf_gets_two_chromosome_workers_without_exceeding_small_fixed_slot(self):
        with patch.dict(os.environ,{'MFM_WGS_PAYLOAD':'sequence_only','MFM_CPU_POOL_PROFILE':'throughput'},clear=True):
            self.assertEqual((6,68),POOL.stage_cap('WGS','03',68))
            self.assertEqual((4,4),POOL.stage_cap('WGS','03',4))
            self.assertEqual((2,2),POOL.stage_cap('WGS','03',2))

    def test_breadth_first_matches_nested_and_cyclic_processes(self):
        graph={1:[2,3],2:[4],4:[1],8:[9]}
        self.assertEqual({1,2,3,4},descendants(graph,[1,2]))
        self.assertEqual({8,9},descendants(graph,[8]))
    def test_timestamp_is_actual_collection_completion_and_cache_is_not_relabelled(self):
        o=PreprocessingCpuObserver();setting={'enabled':True,'reservation_reuse':{'enabled':True}}
        with patch('integrate_wgs_rna.pool_cpu_observer.time.monotonic',side_effect=[10,10.25]):
            first=o.observe(8,setting,reader=lambda roots:dict(scope=[(1,1)],usage_usec=100))
        self.assertEqual(10.25,first['sampled_at']);self.assertEqual(.25,first['collection_seconds'])
        with patch('integrate_wgs_rna.pool_cpu_observer.time.monotonic',return_value=11):
            cached=o.observe(8,setting,reader=lambda roots:self.fail('다시 읽지 않음'))
        self.assertEqual(10.25,cached['sampled_at'])

if __name__=='__main__':unittest.main()
