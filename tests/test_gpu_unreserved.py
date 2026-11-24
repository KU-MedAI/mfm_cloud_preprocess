"""GPU 비교의 예약 예외가 명시 범위에만 적용되는지 검사함."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import json,os,signal,subprocess,sys,tempfile,threading,time,unittest
from integrate_wgs_rna import gpu_compare_resources as m
from integrate_wgs_rna.gpu_stage_resources import StageBudget
from integrate_wgs_rna.realign.phases import RealignPhases

ENV={m.FLAG:'1','RUN_ID':'gpu_compare_test','MFM_WGS_BACKEND':'parabricks'}
class UnreservedTests(unittest.TestCase):
    def test_명시GPU비교외에는예약예외를적용하지않음(self):
        self.assertFalse(m.enabled({}))
        with self.assertRaises(ValueError):m.enabled({m.FLAG:'1','RUN_ID':'production','MFM_WGS_BACKEND':'parabricks'})
        with self.assertRaises(ValueError):m.enabled({m.FLAG:'1','RUN_ID':'gpu_compare_test'})
        self.assertTrue(m.enabled(ENV))

    def test_예약예외도단계별실제병렬설정을전달함(self):
        from integrate_wgs_rna.stage_cpu_config import apply
        env={};apply(env,'WGS','06',68);self.assertEqual(env['MFM_HAPLOTAG_CPU_BUDGET'],'68')
        apply(env,'WGS','03',68);self.assertEqual(env['VCF_CHROM_JOBS'],'22');self.assertEqual(env['VCF_THREADS_PER_CHROM'],'2')
        plain={};apply(plain,'RNA','06',68);self.assertEqual(plain,{})

    def test_명시CPU대조군만같은예약예외를허용함(self):
        env={m.FLAG:'1','RUN_ID':'gpu_compare_cpu_control','MFM_CPU_COMPARISON_CONTROL':'1','MFM_WGS_BACKEND':'cpu','MFM_RNA_BACKEND':'cpu','MFM_BAM_BACKEND':'cpu'}
        self.assertTrue(m.enabled(env))
        env['RUN_ID']='production'
        with self.assertRaises(ValueError):m.enabled(env)
        env['RUN_ID']='gpu_compare_cpu_control';env.pop('MFM_CPU_COMPARISON_CONTROL')
        with self.assertRaises(ValueError):m.enabled(env)

    def test_기존예약파일이없어도단계와GPU3장확보는동작함(self):
        with patch.dict(os.environ,dict(ENV,MFM_CPU_POOL_LEASE_FILE='/missing/lease.json',MFM_GPU_HAP_ID='1'),clear=True):
            phase=RealignPhases();self.assertIsNone(phase.path)
            budget=StageBudget();self.assertIsNone(budget.path)
            lease=SimpleNamespace(close=lambda:None)
            with patch('integrate_wgs_rna.gpu_pool.acquire',return_value=lease) as acquire:
                with budget.gpu(1024,purpose='test',count=3):pass
            self.assertEqual(acquire.call_args.kwargs['count'],3);self.assertTrue(acquire.call_args.kwargs['exclusive'])
            self.assertTrue(all(e['reservation_applied'] is False for e in budget.events))

    def test_큐RAM예약입장도같은비교에만생략함(self):
        from handoff.queue_resources import admit,download_admit
        config={'run_id':'gpu_compare_test','pipeline_env':ENV}
        with patch('handoff.queue_resources.cgroup_limits',side_effect=AssertionError('예약검사 불필요')):
            self.assertTrue(admit(config,'/missing',[])[0]);self.assertTrue(download_admit(config,'/missing')[0])

    def test_실제CPUwrapper는대기장부없이실행하고예약0을기록함(self):
        root=Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as raw:
            p=Path(raw);audit=p/'audit.jsonl';pool=p/'pool'
            env=dict(os.environ,**ENV,MFM_GPU_HOST_MEMORY_BYTES='1')
            command=[sys.executable,str(root/'scripts/cpu_pool.py'),'--pool-dir',str(pool),'--total-cpus','68','--branch','WGS','--stage','10','--sample','S','--audit',str(audit),'--',sys.executable,'-c','import os;assert "MFM_CPU_POOL_LEASE_FILE" not in os.environ;assert int(os.environ["MFM_CPU_POOL_GRANTED"])<=68']
            result=subprocess.run(command,env=env,capture_output=True,text=True,timeout=8)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertFalse((pool/'leases').exists());self.assertFalse((pool/'waiters').exists())
            self.assertEqual(list((pool/'memory_observers').glob('*.json')),[])
            row=json.loads(audit.read_text());self.assertEqual(row['reserved_cpus'],0);self.assertEqual(row['reserved_memory_bytes'],0)
            self.assertEqual(row['queue_wait_seconds'],0)
            self.assertEqual(row['rss_limit_scope'],'preprocessing_pool_total')
            self.assertEqual(row['tool_host_memory_bytes'],1)

    def test_전체샘플은동시진입하고GPU단계잠금으로제어함(self):
        with tempfile.TemporaryDirectory() as raw:
            p=Path(raw);env=dict(ENV,MFM_GPU_POOL_DIR=str(p/'gpu'))
            config={'run_id':'gpu_compare_test','pipeline_env':env,'state_root':str(p/'state')}
            entered=threading.Event();finished=threading.Event()
            def other():
                with m.sample_slot(config,'B'):entered.set()
                finished.set()
            with m.sample_slot(config,'A'):
                thread=threading.Thread(target=other);thread.start()
                self.assertTrue(entered.wait(2))
            thread.join(3);self.assertTrue(finished.is_set())

    def test_Runtime종료신호는별도GPU도구세션도종료함(self):
        with tempfile.TemporaryDirectory() as raw:
            marker=Path(raw)/'pid'
            child="import os,pathlib,time;pathlib.Path("+repr(str(marker))+").write_text(str(os.getpid()));time.sleep(30)"
            driver="from integrate_wgs_rna.parabricks_runtime import Runtime;Runtime().run(['-c',"+repr(child)+"],None,host_limit=1024**3,cpu_limit=1)"
            env=dict(os.environ,MFM_GPU_RUNTIME='native',MFM_PBRUN=sys.executable)
            process=subprocess.Popen([sys.executable,'-c',driver],env=env)
            try:
                deadline=time.monotonic()+5
                while not marker.exists() and time.monotonic()<deadline:time.sleep(.02)
                self.assertTrue(marker.exists());pid=int(marker.read_text());time.sleep(.1)
                process.send_signal(signal.SIGTERM);self.assertEqual(process.wait(timeout=5),143)
                self.assertFalse(Path('/proc/'+str(pid)).exists())
            finally:
                if process.poll() is None:process.kill();process.wait()

if __name__=='__main__':unittest.main()
