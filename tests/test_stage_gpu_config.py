"""GPU 공용 풀 설정 재반영과 장치 예약 경합을 검증함."""
from pathlib import Path
from unittest.mock import patch
import copy,json,os,sys,tempfile,unittest
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(os.environ.get('MFM_TEST_PACKAGE','/NHNHOME/korea/mfm_data/process_code/original'))
sys.path[:0]=[str(ROOT/'src'),str(ROOT)]
from integrate_wgs_rna.stage_gpu_config import apply_stage_environment
from integrate_wgs_rna import gpu_pool

class StageGpuConfigTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='stage_gpu_')
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.state=self.root/'state'
        (self.state/'sample_manifests/S').mkdir(parents=True)
        self.path=self.state/'effective_config.json'
        self.config={'run_id':'test','data_root':str(self.root/'raw'),'output_root':str(self.root/'out'),
                     'acceleration':{'gpu_devices':'0,1,2','gpu_pool_dir':str(self.root/'pool')},
                     'pipeline_env':{'CUDA_VISIBLE_DEVICES':'0,1,2','MFM_GPU_HAP_JOBS':'2','WGS_HAP_JOBS':'auto'}}
        self.env={'MFM_WGS_BACKEND':'parabricks','CUDA_VISIBLE_DEVICES':'1','MFM_GPU_DEVICES':'1',
                  'MFM_WGS_HAP_REQUESTED':'1','MFM_GPU_HAP_JOBS':'1','WGS_HAP_JOBS':'1',
                  'RUN_ID':'test','INPUT_ROOT':self.config['data_root'],'OUTPUT_ROOT':self.config['output_root'],
                  'PIPELINE_MANIFEST':str(self.state/'sample_manifests/S/preprocess_manifest.tsv')}
        self.save()
    def save(self):self.path.write_text(json.dumps(self.config))
    def test_CPU_only_환경_유지(self):
        env=dict(self.env,MFM_WGS_BACKEND='cpu');before=dict(env)
        self.assertIsNone(apply_stage_environment(env));self.assertEqual(env,before)
    def test_고정1번에서_공용3장으로_변경(self):
        info=apply_stage_environment(self.env)
        self.assertEqual(self.env['MFM_GPU_DEVICES'],'0,1,2')
        self.assertEqual(self.env['CUDA_VISIBLE_DEVICES'],'0,1,2')
        self.assertEqual(self.env['MFM_WGS_HAP_REQUESTED'],'auto')
        self.assertEqual(self.env['MFM_GPU_POOL_DIR'],str(self.root/'pool'))
        self.assertEqual(self.env['MFM_WGS_BACKEND'],'parabricks')
        self.assertEqual(info['source'],str(self.path))
    def test_다른_run_거부(self):
        self.env['RUN_ID']='other';before=dict(self.env)
        self.assertIsNone(apply_stage_environment(self.env));self.assertEqual(self.env,before)
    def test_다른_output_거부(self):
        self.env['OUTPUT_ROOT']=str(self.root/'other');before=dict(self.env)
        self.assertIsNone(apply_stage_environment(self.env));self.assertEqual(self.env,before)
    def test_다른_input_거부(self):
        self.env['INPUT_ROOT']=str(self.root/'other');self.assertIsNone(apply_stage_environment(self.env))
    def test_Slurm_CUDA범위와_컨테이너_범위_유지(self):
        self.env.update(SLURM_JOB_ID='123',NVIDIA_VISIBLE_DEVICES='1')
        info=apply_stage_environment(self.env)
        self.assertTrue(info['slurm_cuda_visibility_preserved'])
        self.assertEqual(self.env['CUDA_VISIBLE_DEVICES'],'1')
        self.assertEqual(self.env['NVIDIA_VISIBLE_DEVICES'],'1')
        devices=[{'index':i,'uuid':'GPU-'+str(i)} for i in range(3)]
        self.assertEqual([r['index'] for r in gpu_pool.allowed_devices(devices,environment=self.env)],[1])
    def test_새단계마다_설정_다시읽음(self):
        apply_stage_environment(self.env)
        self.config['acceleration']['gpu_devices']='0,2';self.config['pipeline_env']['CUDA_VISIBLE_DEVICES']='0,2';self.save()
        apply_stage_environment(self.env);self.assertEqual(self.env['MFM_GPU_DEVICES'],'0,2')
    def test_백엔드_변경은_반영하지않음(self):
        self.config['acceleration']['wgs_backend']='cpu';self.save()
        apply_stage_environment(self.env);self.assertEqual(self.env['MFM_WGS_BACKEND'],'parabricks')
    def test_잘못된_HAP상한_차단(self):
        self.config['pipeline_env']['MFM_GPU_HAP_JOBS']='3';self.save()
        with self.assertRaises(ValueError):apply_stage_environment(self.env)
    def test_RNA통합분리선택도다음단계설정으로반영함(self):
        for value in ['integrated','separate']:
            self.config['pipeline_env']['RNA_GPU_MARKDUP_MODE']=value;self.save()
            apply_stage_environment(self.env)
            self.assertEqual(self.env['RNA_GPU_MARKDUP_MODE'],value)
        self.config['pipeline_env']['RNA_GPU_MARKDUP_MODE']='typo';self.save()
        with self.assertRaises(ValueError):apply_stage_environment(self.env)
    def test_공유3장_동시예약과_반환재할당(self):
        devices=[dict(index=i,uuid='GPU-'+str(i),pci=str(i),total=160*gpu_pool.GIB,used=0,utilization=0) for i in range(3)]
        env={'MFM_GPU_POOL_DIR':str(self.root/'leases_pool'),'CUDA_VISIBLE_DEVICES':'0,1,2','MFM_GPU_DEVICES':'0,1,2','NVIDIA_VISIBLE_DEVICES':'all'}
        leases=[]
        with patch.dict(os.environ,env),patch.object(gpu_pool,'inventory',return_value=devices):
            try:
                with ThreadPoolExecutor(max_workers=3) as executor:
                    futures=[executor.submit(gpu_pool.acquire,gpu_pool.GIB,purpose='sample_'+str(i),exclusive=True,timeout=2) for i in range(3)]
                    leases=[f.result() for f in futures]
                self.assertEqual({v.device['index'] for v in leases},{0,1,2})
                with self.assertRaises(TimeoutError):gpu_pool.acquire(gpu_pool.GIB,purpose='fourth',exclusive=True,timeout=0)
                released=leases.pop();index=released.device['index'];released.close()
                reused=gpu_pool.acquire(gpu_pool.GIB,purpose='next_hap',exclusive=True,timeout=0);leases.append(reused)
                self.assertEqual(reused.device['index'],index)
            finally:
                for lease in leases:lease.close()
        self.assertEqual(list((self.root/'leases_pool/leases').glob('*.json')),[])

if __name__=='__main__':unittest.main()
