"""CPU 상한·실행 인수·장부 이관·저장 오류 재시도를 검증함."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import copy,errno,fcntl,hashlib,json,os,subprocess,tempfile,unittest

from handoff import cli,resources,planner,local_queue
from handoff.state_import import import_states
from handoff.state_bindings import StateBindings
from handoff.ready_queue import Scheduler
from handoff.queue_recovery import enqueue


class CpuPolicyTests(unittest.TestCase):
    def setUp(self):
        probe=patch.object(resources,'cgroup_limits',return_value=(None,None,None))
        probe.start();self.addCleanup(probe.stop)
        env=patch.dict(os.environ,{},clear=True);env.start();self.addCleanup(env.stop)

    def config(self, jobs):
        return dict(cpu_pool_limit=68,cpu_utilization='maximum',local_jobs=jobs,
                    sample_cpus='auto',download_jobs='auto',resource_profile='throughput',
                    slurm={'cpus':18,'concurrency':4})

    def test_auto와maximum은각각초기값과공유상한임(self):
        for jobs,expected in [(1,68),(2,34),(3,22),(4,17),(6,17)]:
            c=resources.resolve(self.config(jobs),cpu_count=72,
                available_memory_bytes=1024*resources.GIB,free_storage_bytes=20*resources.TIB)
            self.assertEqual(c['sample_cpus'],expected)
            self.assertEqual(c['resource_resolution']['effective_total_cpu_budget'],68)
        c=resources.resolve(self.config(6),cpu_count=72,
            available_memory_bytes=1024*resources.GIB,free_storage_bytes=20*resources.TIB)
        self.assertEqual(c['resource_resolution']['requested_sample_cpus'],'auto')

    def test_cgroup과Slurm이더작으면그범위를우선함(self):
        with patch.object(resources,'cgroup_limits',return_value=(None,None,24)):
            self.assertEqual(resources.cpu_budget(self.config(1),affinity=72),24)
        with patch.dict(os.environ,{'SLURM_JOB_ID':'test','SLURM_CPUS_PER_TASK':'12'}):
            self.assertEqual(resources.cpu_budget(self.config(1),affinity=72),12)
            c=resources.resolve(self.config(1),cpu_count=72,
                available_memory_bytes=1024*resources.GIB,free_storage_bytes=20*resources.TIB)
            self.assertEqual(c['sample_cpus'],12)
        for invalid in [0,-1,True]:
            with self.assertRaises(ValueError):
                resources.cpu_budget({**self.config(1),'cpu_pool_limit':invalid},affinity=72)


class StateTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.source=self.root/'old';self.source.mkdir()
        self.oldstate=self.source/'state_001'
        self.row=dict(sample_id='S',project_id='P',platform='hifi',data_group='wgs',
                      read_type='long_read',rna_path='',known_input_total_bytes='1000')
        extra={**self.row,'sample_id':'short','platform':'illumina','read_type':'short_read'}
        self.oldrows=[extra,self.row]
        self.config=dict(data_root=str(self.root/'data'),output_root=str(self.root/'output'),
             state_root=str(self.root/'state'),package_root=str(cli.ROOT),run_id='new',
             sample_cpus=17,cpu_pool_limit=68,local_jobs=6,download_jobs=3,resources_resolved=True,
             memory_limit_bytes=10**12,model_target_bytes=10**13,pipeline_env={},
             selection_fingerprint=hashlib.sha256(b'S').hexdigest(),
             resource_resolution={'requested_local_jobs':'auto','requested_download_jobs':'auto'})
        self.old={**self.config,'run_id':'old','state_root':str(self.oldstate),
                  'selection_fingerprint':hashlib.sha256(b'short\nS').hexdigest()}
        self.api=SimpleNamespace(ROOT=cli.ROOT,ORIGINAL_LAYOUT=True,atomic_json=cli.atomic_json,
             sample_root=lambda c,r:Path(c['output_root'])/'wgs/P/S',
             processing_complete=lambda *args:False,prepare=lambda *args:None)
        cli.atomic_json(self.source/'묶음001.json',self.old)
        planner.write_tsv(self.source/'묶음001.tsv',self.oldrows,fields=list(self.row))
        self.record={'sample':'S','status':'running','started':1,'model_bytes':0}
        cli.atomic_json(self.oldstate/'samples/S.json',self.record)
        cli.atomic_json(self.source/'작업재개요청.json',{'samples':{
            'S':{'stage':'12','requested_at':7},'short':{'stage':'01','requested_at':4}}})
        cli.atomic_json(self.source/'재시도상태.json',{'schema':'mfm_local_retry_v1',
            'workers':{},'downloads':{},'applied_recoveries':{'S':7,'short':4}})

    def imported(self):
        return import_states(self.config,[self.row],self.source,self.api)

    def test_장부이관은완료와재개표시보존하고외부코드를실행하지않음(self):
        result=self.imported();self.assertFalse(result['model_data_moved'])
        self.assertFalse(result['outside_python_imported'])
        root=Path(self.config['state_root'])
        self.assertEqual(json.loads((root/'samples/S.json').read_text()),self.record)
        ledger=json.loads((root/'queue/재시도상태.json').read_text())
        self.assertEqual(ledger['applied_recoveries'],{'S':7})
        self.assertEqual(self.imported()['samples'],1)
        binding=StateBindings(self.config,[self.row])
        config,selection,index=binding.launch('S','worker',{'stage':'12','overwrite_prepared':True},
                                             self.api,root/'queue/launch_configs')
        launched=json.loads(config.read_text())
        self.assertEqual(launched['run_id'],'old');self.assertEqual(launched['state_root'],str(self.oldstate))
        self.assertEqual(index,1);self.assertEqual(selection,self.source/'묶음001.tsv')
        self.assertEqual(launched['sample_cpus'],17)
        self.assertEqual(launched['cpu_allocation_mode'],'dynamic')
        binding.config['cpu_allocation_mode']='fixed_sample'
        again,_,_=binding.launch('S','worker',None,self.api,root/'queue/launch_configs')
        self.assertEqual(json.loads(again.read_text())['cpu_allocation_mode'],'fixed_sample')
        self.assertEqual(launched['sample_pipeline_env']['S']['WGS_START_AT'],'12')
        self.assertEqual(launched['sample_pipeline_env']['S']['STAGE04_OVERWRITE'],'1')

    def test_실행중worker와다운로드lock을인수하고중복투입하지않음(self):
        self.imported();a=local_queue.Adapter(self.config,[self.row],self.api)
        download=a.download_lock_state(self.row)
        paths=[self.oldstate/'samples/S.lock',download/'samples/S.lock']
        handles=[]
        try:
            for path in paths:
                path.parent.mkdir(parents=True,exist_ok=True);handle=path.open('a')
                fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB);handles.append(handle)
            s=Scheduler(a,[(1,[self.row])],[self.row]);s.adopt_workers()
            self.assertIn('S',s.active);self.assertIn('S',s.downloads)
            self.assertIsNone(s.active['S'][1].poll());self.assertIsNone(s.downloads['S'][1].poll())
            with patch.object(a,'start') as start:s.start_workers();s.start_downloads();start.assert_not_called()
        finally:
            for handle in handles:handle.close()

    def test_이관중인worker기록은변경된경우만집계본에반영함(self):
        self.imported();binding=StateBindings(self.config,[self.row])
        with patch.object(self.api,'atomic_json',wraps=cli.atomic_json) as save:
            binding.mirror(self.api);binding.mirror(self.api);self.assertEqual(save.call_count,1)
            cli.atomic_json(self.oldstate/'samples/S.json',{**self.record,'status':'success','model_bytes':123})
            binding.mirror(self.api);self.assertEqual(save.call_count,2)

    def test_이전제어큐가살아있으면이관을거부함(self):
        cli.atomic_json(self.source/'연속제어상태.json',{'pid':os.getpid(),'status':'running'})
        with self.assertRaisesRegex(ValueError,'draining'):self.imported()
        self.assertFalse((Path(self.config['state_root'])/'imported_samples.json').exists())

    def test_선택파일변경과출력경로변경을거부함(self):
        self.imported()
        with self.assertRaisesRegex(ValueError,'경로'):
            StateBindings({**self.config,'output_root':str(self.root/'other')},[self.row])
        with (self.source/'묶음001.tsv').open('a') as f:f.write('\n')
        with self.assertRaisesRegex(ValueError,'선택 파일'):StateBindings(self.config,[self.row])

    def test_일부샘플누락이면부분장부를남기지않음(self):
        missing={**self.row,'sample_id':'missing'}
        with self.assertRaisesRegex(ValueError,'찾지 못한'):
            import_states(self.config,[self.row,missing],self.source,self.api)
        self.assertFalse((Path(self.config['state_root'])/'imported_samples.json').exists())
        self.assertFalse((Path(self.config['state_root'])/'queue/재시도상태.json').exists())

    def test_resume는장부만갱신하며사용자가JSON을만들필요없음(self):
        self.imported();state=Path(self.config['state_root'])
        cli.atomic_json(state/'effective_config.json',self.config)
        planner.write_tsv(state/'selected_samples.tsv',[self.row],fields=list(self.row))
        result=enqueue(self.config,[self.row],['S'],'1',True,True,cli.atomic_json)
        self.assertEqual(result['from_stage'],'01')
        request=json.loads((state/'queue/작업재개요청.json').read_text())['samples']['S']
        self.assertTrue(request['rebuild_haplotag_partial']);self.assertNotEqual(request['requested_at'],7)
        with self.assertRaises(ValueError):enqueue(self.config,[self.row],['S'],'99',False,False,cli.atomic_json)

    def test_이관후이전제어기를다시켜도공유잠금으로차단함(self):
        path=self.source/'연속운영.lock';path.touch();self.imported()
        binding=StateBindings(self.config,[self.row])
        with binding.control_lock(),path.open('a') as handle:
            with self.assertRaises(BlockingIOError):fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)


class EntrypointTests(unittest.TestCase):
    def test_외부PYTHONPATH가없어도쉘단독사전검사성공함(self):
        root=Path(cli.ROOT)
        env=dict(os.environ,PYTHONPATH='',PYTHONDONTWRITEBYTECODE='1')
        result=subprocess.run(['bash',str(root/'run_pipeline.sh'),'dry-run',
            '--lab','medai','--batch','001','--config',str(root/'handoff_config.example.json')],
            env=env,text=True,capture_output=True,timeout=30)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertGreater(json.loads(result.stdout)['selected'],0)


class StorageTests(unittest.TestCase):
    def test_quota오류시실행객체를보존하며동일큐가재시도함(self):
        scheduler=object.__new__(Scheduler)
        events=[];scheduler.c=SimpleNamespace(event=lambda *a,**k:events.append(a[0]))
        scheduler.clock=lambda:100;scheduler.active={'S':object()}
        with patch.object(scheduler,'adopt_workers'),patch.object(scheduler,'step',side_effect=[
             OSError(errno.EDQUOT,'quota'),OSError(errno.ENOSPC,'full'),0]), \
             patch('handoff.ready_queue.time.sleep') as sleep:
            self.assertEqual(scheduler.run(),0);self.assertEqual(sleep.call_count,2)
        self.assertIn('S',scheduler.active);self.assertIn('storage_write_recovered',events)
        with patch.object(scheduler,'adopt_workers'),patch.object(scheduler,'step',side_effect=OSError(errno.EIO,'io')):
            with self.assertRaises(OSError):scheduler.run()


if __name__=='__main__':unittest.main()
