"""다중 GPU 원자 예약과 실행환경 격리를 검증함."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import json,os,sys,tempfile,unittest
from integrate_wgs_rna import gpu_pool,gpu_wgs
from integrate_wgs_rna.parabricks_runtime import Runtime

class MultiGpuTests(unittest.TestCase):
    def test_종료된계산뒤재정렬대기는CPU를예약하지않음(self):
        from integrate_wgs_rna.pool_leases import regular_cpu_claim
        self.assertEqual(regular_cpu_claim({'cpus':1,'realign_phase':'align_restore_sort_qc:waiting'}),0)
        self.assertEqual(regular_cpu_claim({'cpus':8,'realign_phase':'align_restore_sort_qc'}),8)
        self.assertEqual(regular_cpu_claim({'cpus':1}),1)

    def test_GPU대여는명시정책과실측유휴조건이필요함(self):
        import time
        from integrate_wgs_rna.realign import idle_resources as m
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{'MFM_WGS_BACKEND':'parabricks','MFM_CPU_UTILIZATION':'maximum'},clear=True):
            root=Path(raw);policy={'cpu_backfill':{'enabled':True,'realign_phases':True,'nice':0}}
            (root/'policy.json').write_text(json.dumps(policy))
            (root/'watchdog.json').write_text(json.dumps({'updated':time.time(),'guard_action':'continue','cpu_backfill':{'ready':True,'quota':68,'observation_seconds':10,'idle_cpus':12,'used_cpus':40,'peak_used_cpus':40}}))
            sample={'branch':'WGS','stage':'10'}
            with patch.object(m,'cgroup_sample',return_value={'quota':68,'cpu_cgroup':'test'}),patch.object(m,'in_cpu_group',return_value=True),patch.object(m.os,'getpriority',return_value=0):
                self.assertIsNone(m.loan_capacity(root,68,sample,[],128*gpu_pool.GIB))
                policy['cpu_backfill']['gpu_realign_phases']=True;(root/'policy.json').write_text(json.dumps(policy))
                result=m.loan_capacity(root,68,sample,[],128*gpu_pool.GIB)
                self.assertEqual(result[2]['kind'],'gpu_realign_phase')
                self.assertGreaterEqual(result[0],12)

    def test_GPU단계반환후최초대여상한안에서만재확보함(self):
        from integrate_wgs_rna.gpu_stage_resources import StageBudget
        from integrate_wgs_rna.realign.phases import _atomic
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir();own=root/'leases/own.json'
            _atomic(own,{'cpus':12,'memory_bytes':4*gpu_pool.GIB,'total_cpus':68,'effective_alignment':{'hap_jobs':1},'cpu_backfill':{'cpus':12}})
            _atomic(root/'leases/other.json',{'cpus':68,'memory_bytes':4*gpu_pool.GIB})
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(own),'MFM_GPU_HAP_ID':'1'},clear=True):
                budget=StageBudget();budget.wait()
                with patch('integrate_wgs_rna.realign.idle_resources.loan_capacity',return_value=(12,300*gpu_pool.GIB,{})):
                    self.assertTrue(budget.try_phase('gpu_compute',12,150*gpu_pool.GIB))
                budget.wait()
                with patch('integrate_wgs_rna.realign.idle_resources.loan_capacity',return_value=(11,300*gpu_pool.GIB,{})):
                    self.assertFalse(budget.try_phase('gpu_compute',12,150*gpu_pool.GIB))
                self.assertEqual(json.loads(own.read_text())['cpus'],1)

    def test_GPU의pairedBAM준비도로컬에생성하고CPU기본경로는유지함(self):
        from integrate_wgs_rna.personal_realign import prepare_illumina_paired_input
        with tempfile.TemporaryDirectory() as raw:
            work=Path(raw)/'work';work.mkdir()
            for backend in ['cpu','parabricks']:
                with patch.dict(os.environ,{'MFM_WGS_BACKEND':backend,'MFM_GPU_SCRATCH_ROOT':str(Path(raw)/'scratch')},clear=True),patch('integrate_wgs_rna.realign.prepared_cache.prepare',return_value='ok') as prepare:
                    result=prepare_illumina_paired_input(SimpleNamespace(),'S',1,Path('/input.bam'),['chr1'],work,None,'all_to_both')
                self.assertEqual(result,'ok')
                target=Path(prepare.call_args.args[5])
                if backend=='cpu':self.assertEqual(target,work)
                else:self.assertIn(Path(raw)/'scratch',target.parents)
                self.assertEqual(Path(prepare.call_args.args[8]).parent,target)

    def test_공유잠금의일시EAGAIN만재시도하고다른오류는보고함(self):
        import importlib.util
        script=Path(__file__).resolve().parents[1]/'scripts/cpu_pool.py'
        spec=importlib.util.spec_from_file_location('gpu_lock_test',script);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
        with patch.object(m.fcntl,'flock',side_effect=[BlockingIOError(11,'일시 잠금 경합'),None]) as lock,patch.object(m.time,'sleep') as sleep:
            m.lock_pool(object());self.assertEqual(lock.call_count,2);sleep.assert_called_once_with(m.POLL_SECONDS)
        with patch.object(m.fcntl,'flock',side_effect=PermissionError(13,'권한 오류')),self.assertRaises(PermissionError):m.lock_pool(object())

    def devices(self,n=4):
        return [dict(index=i,uuid='GPU-'+str(i),pci=str(i),total=64*gpu_pool.GIB,used=0,utilization=0) for i in range(n)]

    def test_세GPU를모두예약하고다른정렬의부분예약을남기지않음(self):
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{'MFM_GPU_POOL_DIR':raw,'MFM_GPU_DEVICES':'auto'},clear=True),patch.object(gpu_pool,'inventory',return_value=self.devices()):
            one=gpu_pool.acquire(gpu_pool.GIB,purpose='three',exclusive=True,count=3,timeout=0)
            self.assertEqual([d['index'] for d in one.devices],[0,1,2])
            self.assertEqual(one.device['uuid'],'GPU-0,GPU-1,GPU-2')
            with self.assertRaises(TimeoutError):gpu_pool.acquire(gpu_pool.GIB,purpose='another',exclusive=True,count=3,timeout=0)
            self.assertEqual(len(list((Path(raw)/'leases').glob('*.json'))),3)
            one.close();self.assertEqual(list((Path(raw)/'leases').glob('*.json')),[])

    def test_예약파일저장실패도부분예약을모두회수함(self):
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{'MFM_GPU_POOL_DIR':raw,'MFM_GPU_DEVICES':'auto'},clear=True),patch.object(gpu_pool,'inventory',return_value=self.devices()):
            original=gpu_pool._save;calls=[]
            def save(path,row):
                calls.append(path)
                if len(calls)==3:raise OSError('저장 실패')
                original(path,row)
            with patch.object(gpu_pool,'_save',side_effect=save),self.assertRaises(OSError):gpu_pool.acquire(gpu_pool.GIB,purpose='fail',exclusive=True,count=3,timeout=0)
            self.assertEqual(list((Path(raw)/'leases').glob('*.json')),[])

    def test_외부작업과CUDA허용범위를보존함(self):
        devices=self.devices();devices[0].update(used=gpu_pool.GIB,utilization=100)
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{'MFM_GPU_POOL_DIR':raw,'MFM_GPU_DEVICES':'auto','CUDA_VISIBLE_DEVICES':'0,1,2'},clear=True),patch.object(gpu_pool,'inventory',return_value=devices):
            with self.assertRaises(TimeoutError):gpu_pool.acquire(gpu_pool.GIB,purpose='wait',exclusive=True,count=3,timeout=0)
            self.assertEqual(list((Path(raw)/'leases').glob('*.json')),[])

    def test_Docker에예약된세장만공개하고CLI개수도일치함(self):
        with patch.dict(os.environ,{'MFM_GPU_RUNTIME':'docker'}),patch('shutil.which',return_value='/bin/docker'):
            command=Runtime().command(['minimap2'],'GPU-a,GPU-b,GPU-c',cpu_limit=8,host_limit=16*1024**3)
        self.assertIn('CUDA_VISIBLE_DEVICES=0,1,2',command)
        import csv
        self.assertEqual(next(csv.reader([command[command.index('--gpus')+1]])),['device=GPU-a,GPU-b,GPU-c'])
        args=gpu_wgs.commands('/ref',['/reads'],'/output','hifi','map-hifi',8,'S',3)
        self.assertEqual(args[args.index('--num-gpus')+1],'3')

    def test_native자식의장치와메타데이터도세장에일치함(self):
        gpu=SimpleNamespace(device={'uuid':'GPU-a,GPU-b,GPU-c'},reserved=3*gpu_pool.GIB)
        with patch.dict(os.environ,{'MFM_GPU_RUNTIME':'native','MFM_PBRUN':sys.executable}):
            result=Runtime().run(['-c','import os;assert os.environ["CUDA_VISIBLE_DEVICES"]=="GPU-a,GPU-b,GPU-c"'],gpu,host_limit=gpu_pool.GIB,cpu_limit=1)
        self.assertEqual(result['gpu_count'],3);self.assertEqual(result['gpu_uuids'],['GPU-a','GPU-b','GPU-c'])

if __name__=='__main__':unittest.main()
