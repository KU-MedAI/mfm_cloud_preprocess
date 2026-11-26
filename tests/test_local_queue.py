"""공용 adapter의 자원·실제 상태 파일·동시 흐름을 네트워크 없이 검증함."""
from pathlib import Path
from types import SimpleNamespace
import copy
import fcntl
import hashlib
import json
import os
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src')]
from handoff import cli, local_queue as module, queue_resources as memory
from handoff.ready_queue import Scheduler as RealScheduler


class LocalTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        watcher=patch.object(module.Adapter,'ensure_watchdog')
        watcher.start();self.addCleanup(watcher.stop)  # 가짜 worker 검증은 실제 daemon을 띄우지 않음.
        self.base=Path(self.temp.name);self.done=set();self.events=[];self.clock=1000
        self.config={'data_root':str(self.base/'data'),'output_root':str(self.base/'output'),
                     'state_root':str(self.base/'state'),'package_root':str(ROOT),
                     'run_id':'test','selection_fingerprint':'a'*64,'sample_cpus':8,
                     'local_jobs':2,'download_jobs':3,'memory_limit_bytes':10**12,
                     'model_target_bytes':10**12,'pipeline_env':{
                         'MFM_SHORT_READ_OUTPUT_MODE':'full','MFM_LONG_READ_OUTPUT_MODE':'full'},
                     'resource_resolution':{'requested_local_jobs':'auto','requested_download_jobs':'auto'}}
        self.rows=[dict(sample_id=s,platform='illumina',known_input_total_bytes='1000',
                        data_group='wgs',project_id='P',rna_path='') for s in ['slow','fast','bad']]
        self.config['selection_fingerprint']=hashlib.sha256('\n'.join(r['sample_id'] for r in self.rows).encode()).hexdigest()
        self.finish_calls=[]
        self.api=SimpleNamespace(ROOT=ROOT,ORIGINAL_LAYOUT=True,atomic_json=cli.atomic_json,
            sample_root=lambda c,r:Path(c['output_root'])/'wgs'/'P'/r['sample_id'],
            processing_complete=lambda root,sample,mode:sample in self.done,
            model_bytes=lambda root:100,prepare=lambda *a:None,
            finish=lambda c,r:self.finish_calls.append(r) or {'success':len(self.done)})

    def test_Stage01완료후사전검증재개요청을수용(self):
        adapter=module.Adapter(self.config,self.rows,self.api)
        self.assertIn('02',adapter.RECOVERY_STAGES)
        self.api.atomic_json(adapter.BASE/'작업재개요청.json',{'samples':{'slow':{'stage':'02','requested_at':1,'reason':'누락index보완'}}})
        scheduler=RealScheduler(adapter,[(1,self.rows)],self.rows)
        scheduler.refresh_recoveries()
        self.assertEqual(scheduler.recoveries['slow']['stage'],'02')

    def test_실제상태파일로_빠른샘플우선_실패격리_버퍼진행(self):
        test=self
        class ClockScheduler(RealScheduler):
            def __init__(self,q,work,rows):super().__init__(q,work,rows,clock=lambda:test.clock)
        class Process:
            pid=100
            def __init__(self,kind,adapter,row):
                self.kind,self.adapter,self.row=kind,adapter,row
                self.started=test.clock;self.result=None
            def poll(self):
                if self.result is not None:return self.result
                delay=180 if self.row['sample_id']=='slow' and self.kind=='download-worker' else 30
                if test.clock-self.started<delay:return None
                sample=self.row['sample_id']
                if self.kind=='download-worker':
                    state=Path(self.adapter.download_config['state_root'])/'samples'/(sample+'.json')
                    if sample=='bad':
                        self.adapter.api.atomic_json(state,{'status':'failed','error':'의도한 HTTP 실패'})
                        self.result=1
                    else:
                        path=Path(self.adapter.config['data_root'])/(sample+'.bam');path.write_bytes(b'input')
                        self.adapter.api.atomic_json(state,{'status':'success','finished':test.clock,
                            'files':[{'status':'verified','local_path':str(path),'actual_bytes':5}]})
                        self.result=0
                        test.events.append(('download_done',sample))
                else:
                    test.done.add(sample)
                    self.adapter.api.atomic_json(self.adapter.state/'samples'/(sample+'.json'),
                                                 {'status':'success','model_bytes':100})
                    self.result=0
                return self.result
            def wait(self):return self.poll()
        def start(adapter,kind,tag,selection,config,**kwargs):
            row=test.rows[kwargs['task_index']]
            test.events.append((kind,row['sample_id']))
            return Process(kind,adapter,row)
        def sleep(seconds):test.clock+=60
        with patch.object(module,'Scheduler',ClockScheduler),patch.object(module.Adapter,'start',start), \
             patch('handoff.ready_queue.time.sleep',side_effect=sleep), \
             patch.object(memory,'admit',return_value=(True,'')), \
             patch.object(module,'admit',return_value=(True,'')):
            result=module.run(copy.deepcopy(self.config),self.rows,self.api)
        self.assertEqual(result,1)
        self.assertEqual(self.done,{'slow','fast'})
        self.assertLess(self.events.index(('worker','fast')),self.events.index(('download_done','slow')))
        self.assertEqual(self.events.count(('download-worker','bad')),3)
        self.assertEqual(len(self.finish_calls),1)
        state=json.loads((self.base/'state/queue/status.json').read_text())
        self.assertEqual(state['status'],'queue_finished')
        ledger=json.loads((self.base/'state/queue/재시도상태.json').read_text())
        self.assertTrue(ledger['downloads']['bad']['deferred'])

    def test_선택범위다르면기존설정저장전거부(self):
        a=module.Adapter(self.config,self.rows,self.api)
        old=dict(a.config,selection_fingerprint='b'*64)
        cli.atomic_json(a.config_path,old)
        with self.assertRaisesRegex(ValueError,'선택 범위'):
            module.run(self.config,self.rows,self.api)
        self.assertEqual(json.loads(a.config_path.read_text()),old)

    def test_작은cgroupRAM을실행과계획모두우선(self):
        with patch.object(memory,'cgroup_limits',return_value=(0,400_000_000_000,68)):
            a=module.Adapter(self.config,self.rows,self.api)
            self.assertEqual(a.config['memory_limit_bytes'],400_000_000_000)
            self.assertEqual(module.plan(self.config)['ram_reservation_bytes'],360_000_000_000)

    def test_준비CPU기본4와사용자명시값을보존(self):
        a=module.Adapter(self.config,self.rows,self.api)
        self.assertEqual(a.config['sample_pipeline_env']['slow']['MFM_REALIGN_PREPARE_MAX_CPUS'],'4')
        c=copy.deepcopy(self.config);c['pipeline_env']['MFM_REALIGN_PREPARE_MAX_CPUS']='8'
        a=module.Adapter(c,self.rows,self.api)
        self.assertNotIn('sample_pipeline_env',a.config)

    def test_완료메타데이터가있어도_worker잠금이끝나기전에는완료아님(self):
        self.done.add('slow');self.config['model_target_bytes']=50
        a=module.Adapter(self.config,self.rows,self.api)
        path=a.state/'samples/slow.lock';path.parent.mkdir(parents=True)
        cli.atomic_json(a.state/'samples/slow.json',{'status':'running','model_bytes':0})
        with path.open('a') as handle:
            fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
            self.assertEqual(a.completed_ids(self.rows),set())
            self.assertFalse(a.target_reached())
        self.assertEqual(a.completed_ids(self.rows),{'slow'})
        self.assertTrue(a.target_reached())

    def test_기존완료로용량목표에도달했으면새작업없음(self):
        self.done.add('slow');self.config['model_target_bytes']=50
        with patch.object(module.Adapter,'start') as start:
            result=module.run(self.config,self.rows,self.api)
        self.assertEqual(result,0);start.assert_not_called()
        state=json.loads((self.base/'state/queue/status.json').read_text())
        self.assertTrue(state['target_reached'])
        self.assertEqual(state['success'],1)

    def test_완료본정책이다르면덮어쓰기전중단(self):
        self.done.add('slow');self.config['pipeline_env']['WINDOW_STORAGE_POLICY']='annotated_acgtn_v1'
        with self.assertRaisesRegex(ValueError,'window 정책'):
            module.run(self.config,self.rows,self.api)
        self.assertFalse((self.base/'state/effective_config.json').exists())

    def test_선택목록순서가다르면상태생성전거부(self):
        with self.assertRaisesRegex(ValueError,'선택 목록/순서'):
            module.run(self.config,list(reversed(self.rows)),self.api)
        self.assertFalse((self.base/'state').exists())

    def test_중복제어기는준비이전에거부(self):
        adapter=module.Adapter(self.config,self.rows,self.api)
        with (adapter.BASE/'control.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            with self.assertRaises(BlockingIOError):module.run(self.config,self.rows,self.api)
        self.assertFalse(adapter.config_path.exists())

    def test_완료결과는원본이없어도완료로인식(self):
        self.done.add('slow');adapter=module.Adapter(self.config,self.rows,self.api)
        self.assertEqual(adapter.completed_ids(self.rows),{'slow'})
        self.assertFalse((self.base/'data').exists())

    def test_auto계산도확정된동시수상한을사용(self):
        a=module.Adapter(self.config,self.rows,self.api);self.assertEqual(self.config['local_jobs'],a.MAX_ACTIVE_SAMPLES)
        c=copy.deepcopy(self.config);c['resource_resolution']['requested_local_jobs']=2
        c['resource_resolution']['requested_download_jobs']=2;c['download_jobs']=2
        a=module.Adapter(c,self.rows,self.api)
        self.assertEqual(a.MAX_ACTIVE_SAMPLES,2);self.assertEqual(a.DOWNLOAD_TARGET,2)

    def test_요청보다작게확정된동시수와dryrun이일치함(self):
        c=copy.deepcopy(self.config);c['local_jobs']=2
        c['resource_resolution']['requested_local_jobs']=4
        a=module.Adapter(c,self.rows,self.api)
        self.assertEqual(a.MAX_ACTIVE_SAMPLES,2)
        self.assertEqual(module.plan(c)['compute_slots'],2)

    def test_legacy복구단계는08원본은10(self):
        a=module.Adapter(self.config,self.rows,self.api);self.assertIn('12',a.RECOVERY_STAGES)
        self.api.ORIGINAL_LAYOUT=False
        a=module.Adapter(self.config,self.rows,self.api);self.assertIn('08',a.RECOVERY_STAGES)
        self.assertNotIn('12',a.RECOVERY_STAGES)

    def test_dry계획은폴더를생성하지않음(self):
        result=module.plan(self.config)
        self.assertEqual(result['mode'],'ready_first_local')
        self.assertEqual(result['ram_reservation_bytes'],900_000_000_000)
        self.assertFalse((self.base/'state').exists())

    def test_cgroupCPU가affinity보다작으면작은상한(self):
        with patch.object(memory,'cpu_budget',return_value=72), \
             patch.object(memory,'cgroup_limits',return_value=(0,10**12,68)):
            self.assertEqual(memory.cpu_limit(self.config),68)

    def test_진행중단계실예약과_전환중전체예약구분(self):
        pool=self.base/'pool';(pool/'leases').mkdir(parents=True)
        row=self.rows[0]
        lease={'sample':'slow','manager_pid':os.getpid(),'child_pid':None,'cpus':1,
               'branch':'WGS','stage':'10','memory_bytes':110_000_000_000}
        (pool/'leases/a.json').write_text(json.dumps(lease))
        legacy=memory.active_ram(pool,[row],final_stage='10')
        full=memory.active_ram(pool,[row],final_stage='12')
        self.assertEqual(legacy,full)
        self.assertEqual(full,16*memory.GIB+110_000_000_000)
        self.assertGreater(memory.active_ram(self.base/"empty",[row]),full)
        pair=dict(row,data_group='wgs_rna_pair')
        self.assertEqual(memory.active_ram(pool,[pair],final_stage='10'),full)

if __name__=='__main__':unittest.main()
