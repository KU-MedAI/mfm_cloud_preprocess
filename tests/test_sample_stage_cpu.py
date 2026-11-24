"""샘플 CPU 상한 및 병렬 검증의 동일성·실패 전파를 확인함."""
import hashlib,importlib.util,json,os,tempfile,unittest,subprocess,sys
from pathlib import Path
from unittest.mock import patch
import numpy as np
import zarr
from integrate_wgs_rna.validation_parallel import cpu_budget,sequence_checks,depth_checks,hash_workers
from integrate_wgs_rna.stage_cpu_config import apply
from integrate_wgs_rna.sample_cpu import sample_slot
from tools.personal_assets import verify
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('stage_caps_test',ROOT/'scripts/cpu_pool.py')
pool=importlib.util.module_from_spec(spec);spec.loader.exec_module(pool)

class SampleStageCpuTests(unittest.TestCase):
    def test_fixed_16_has_no_small_stage_caps(self):
        for layout in ['original','legacy']:
            with patch.dict(os.environ,{'MFM_CPU_ALLOCATION_MODE':'fixed_sample',
                'MFM_FIXED_SAMPLE_CPUS':'16','MFM_CPU_POOL_LAYOUT':layout,
                'MFM_CPU_UTILIZATION':'maximum','MFM_REALIGN_PREPARE_MAX_CPUS':'4'},clear=True):
                for branch,stages in [('WGS',list(range(1,13))+[23,24]),('RNA',range(11 if layout=='legacy' else 13,25 if layout=='legacy' else 27))]:
                    for stage in stages:self.assertEqual(pool.stage_cap(branch,f'{stage:02d}',68)[1],16,(layout,branch,stage))
                self.assertEqual(pool.stage_cap('WGS','23',8)[1],8)
                self.assertEqual(pool.stage_cap('IO','cram_crc',68)[1],16)
                with patch.dict(os.environ,MFM_CPU_POOL_STAGE_MAX='8'):
                    self.assertEqual(pool.stage_cap('RNA','19',68)[1],8)
                    self.assertEqual(pool.stage_cap('IO','cram_crc',68)[1],8)

    def test_rna_defaults_receive_stage_grant(self):
        for layout,index,align,fastq in [('original','15','16','13'),('legacy','13','14','11')]:
            for stage,key in [(index,'STAR_INDEX_THREADS'),(align,'STAR_ALIGN_THREADS'),(fastq,'RNA_COLLATE_THREADS')]:
                env={'MFM_CPU_ALLOCATION_MODE':'fixed_sample','MFM_CPU_POOL_LAYOUT':layout,key:'4'}
                apply(env,'RNA',stage,16);self.assertEqual(env[key],'16')
                apply(env,'RNA',stage,8);self.assertEqual(env[key],'8')

    def test_qc_budget_respects_affinity_and_shared_dag(self):
        with patch.dict(os.environ,{'MFM_SAMPLE_CPU_BUDGET':'16','MFM_CPU_POOL_GRANTED':'16'},clear=True),patch('os.sched_getaffinity',return_value=set(range(16))):
            self.assertEqual(cpu_budget(),16)
            with patch.dict(os.environ,MFM_CPU_POOL_STAGE_MAX='4'):self.assertEqual(cpu_budget(),4)
        with patch.dict(os.environ,{},clear=True):self.assertEqual(cpu_budget(),1)

    def test_depth_equal_and_failure_preserved_with_live_ram_lease(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);store=root/'depth.zarr';group=zarr.open_group(str(store),mode='w');expected=[];contigs=[]
            for i in range(4):
                values=(np.arange(10000,dtype=np.uint32)*(i+1))%37;contig=f'chr{i+1}';contigs.append(contig)
                group.create_dataset(contig,data=values,chunks=(1024,))
                expected.append((int(values.sum()),int(np.count_nonzero(values)),hashlib.sha256(values.astype('<u4').tobytes()).hexdigest()))
            config=dict(cpu_allocation_mode='fixed_sample',sample_cpus=4,cpu_pool_limit=4,memory_limit_bytes=256*(1<<30),
                        state_root=str(root),run_id='test',pipeline_env={'MFM_CPU_POOL_DIR':str(root/'pool')})
            with patch.dict(os.environ,{'MFM_RAM_LIMIT_BYTES':str(256*(1<<30)),'MFM_CPU_UTILIZATION':'maximum'},clear=True),sample_slot(config,'S') as env:
                lease=root/'pool/leases/stage.json'
                lease.write_text(json.dumps(dict(manager_pid=os.getpid(),sample='S',branch='RNA',stage='25',total_cpus=4,
                     cpus=4,memory_bytes=4*(1<<30),fixed_sample_token=env['MFM_FIXED_SAMPLE_TOKEN'])))
                with patch.dict(os.environ,{**env,'MFM_CPU_POOL_LEASE_FILE':str(lease),'MFM_CPU_POOL_GRANTED':'4'}):
                    self.assertEqual(depth_checks(store,contigs),expected)
                    self.assertLessEqual(json.loads(lease.read_text())['peak_cpus'],4)
                    self.assertEqual(json.loads(lease.read_text())['cpus'],1)
                    with hash_workers() as workers:
                        self.assertEqual(workers,3)
                        self.assertEqual(json.loads(lease.read_text())['cpus'],4)
                        self.assertGreater(json.loads(lease.read_text())['memory_bytes'],4*(1<<30))
                    group.create_dataset('wrong_dtype',data=np.ones(10,dtype=np.float32))
                    with self.assertRaisesRegex(ValueError,'dtype/shape'):depth_checks(store,['chr1','wrong_dtype'])
                    self.assertTrue(all(value is None and error is not None for value,error in sequence_checks([root/'missing1',root/'missing2'])))
                lease.unlink()
            self.assertEqual(list((root/'pool/leases').glob('*.json')),[])

    def test_qc_command_receives_all_16_and_sample_owner_stays_held(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            config=dict(cpu_allocation_mode='fixed_sample',sample_cpus=16,cpu_pool_limit=16,memory_limit_bytes=256*(1<<30),
                        state_root=str(root),run_id='test',pipeline_env={'MFM_CPU_POOL_DIR':str(root/'pool')})
            base={'PATH':os.environ['PATH'],'PYTHONPATH':str(ROOT)+':'+str(ROOT/'src'),
                  'MFM_RAM_LIMIT_BYTES':str(256*(1<<30)),'MFM_CPU_UTILIZATION':'maximum'}
            with patch.dict(os.environ,base,clear=True),sample_slot(config,'S') as env:
                code='import os;assert os.environ["MFM_CPU_POOL_GRANTED"]=="16";assert len(os.sched_getaffinity(0))==16'
                result=subprocess.run([sys.executable,str(ROOT/'scripts/cpu_pool.py'),'--pool-dir',str(root/'pool'),
                    '--total-cpus','16','--branch','WGS','--stage','23','--sample','S','--',sys.executable,'-c',code],
                    env={**base,**env},capture_output=True,text=True,timeout=15)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                leases=[json.loads(p.read_text()) for p in (root/'pool/leases').glob('*.json')]
                self.assertEqual(len(leases),1);self.assertTrue(leases[0]['sample_cpu_owner']);self.assertEqual(leases[0]['cpus'],16)

    def test_hash_parallel_preserves_tamper_detection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);files={}
            for i in range(10):
                data=(str(i)*10000).encode();p=root/str(i);p.write_bytes(data)
                files[p.name]={'size_bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
            (root/'permanent_assets.json').write_text(json.dumps({'status':'pass','files':files}))
            self.assertEqual(verify(root,workers=1),verify(root,workers=16))
            (root/'3').write_bytes(b'bad');(root/'8').unlink()
            serial=verify(root,workers=1);parallel=verify(root,workers=16)
            self.assertEqual(serial,parallel);self.assertEqual(parallel['status'],'fail');self.assertEqual(len(parallel['problems']),2)

if __name__=='__main__':unittest.main()
