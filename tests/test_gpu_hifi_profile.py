"""GPU HiFi 버퍼 산정·명시 설정 보존·recovery 진단을 검사함."""
from pathlib import Path
from unittest.mock import patch
import os,tempfile,unittest
from integrate_wgs_rna.gpu_wgs import profile_hifi_input,commands
from integrate_wgs_rna.parabricks_runtime import minimap_gpu_diagnostics

class HifiProfileTests(unittest.TestCase):
    def test_실제길이와최대길이로산정하고입력은보존함(self):
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{},clear=True):
            p=Path(raw)/'r.fq';p.write_text(''.join('@r'+str(i)+' HP:i:1 PS:i:10\n'+'A'*n+'\n+\n'+'I'*n+'\n' for i,n in enumerate([10000,20000,40000])))
            before=p.read_bytes();result=profile_hifi_input(p)
            self.assertEqual(result['p95_read_length'],40000);self.assertEqual(result['sampled_records'],3)
            self.assertEqual(result['effective_parameters']['MFM_GPU_MINIMAP_ESTIMATED_ANCHOR_PER_READ'],80000)
            self.assertEqual(p.read_bytes(),before)
            args=commands('/ref',[p],'/out','hifi','map-hifi',16,'S')
            self.assertEqual(args[args.index('--estimated-read-length')+1],'40000')
            self.assertNotIn('--stage2-large-chaining-size',args)
    def test_명시설정우선과비활성선택(self):
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{'MFM_GPU_MINIMAP_ESTIMATED_READ_LENGTH':'50000'},clear=True):
            p=Path(raw)/'r.fq';p.write_text('@r\nAAAA\n+\nIIII\n');r=profile_hifi_input(p)
            self.assertEqual(r['effective_parameters']['MFM_GPU_MINIMAP_ESTIMATED_READ_LENGTH'],50000)
            os.environ['MFM_GPU_HIFI_AUTO_PROFILE']='0';self.assertIsNone(profile_hifi_input('/missing'))
    def test_잘못된수치와Illumina비적용(self):
        with patch.dict(os.environ,{'MFM_GPU_MINIMAP_ESTIMATED_READ_LENGTH':'0'},clear=True):
            with self.assertRaises(ValueError):commands('/ref',['/fq'],'/out','hifi','map-hifi',16,'S')
            args=commands('/ref',['/fq1','/fq2'],'/out','illumina','map-hifi',16,'S')
            self.assertNotIn('--estimated-read-length',args)
    def test_로그에서복구와GPU커널시간기록함(self):
        with tempfile.TemporaryDirectory() as raw:
            p=Path(raw)/'p.log';p.write_text('GPU seed kernel time: 0.01 s\nGPU chain kernel time: 0.48 s\nGPU seeding/chaining recovery mode triggered: 8810 time(s)\n')
            self.assertEqual(minimap_gpu_diagnostics(p),{'gpu_seed_kernel_seconds':.01,'gpu_chain_kernel_seconds':.48,'gpu_seeding_chaining_recovery_count':8810})
            self.assertEqual(minimap_gpu_diagnostics(None),{})
    def test_같은실행설정의프로파일과host한도를새단계에반영함(self):
        import json
        from integrate_wgs_rna.stage_gpu_config import apply_stage_environment
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'sample_manifests/S').mkdir(parents=True)
            (root/'effective_config.json').write_text(json.dumps({'run_id':'test','data_root':'/raw','output_root':'/out','pipeline_env':{'MFM_GPU_HOST_MEMORY_BYTES':'549755813888','MFM_GPU_HIFI_AUTO_PROFILE':'1'}}))
            env={'RUN_ID':'test','INPUT_ROOT':'/raw','OUTPUT_ROOT':'/out','MFM_WGS_BACKEND':'parabricks','PIPELINE_MANIFEST':str(root/'sample_manifests/S/m.tsv')}
            apply_stage_environment(env);self.assertEqual(env['MFM_GPU_HOST_MEMORY_BYTES'],'549755813888');self.assertEqual(env['MFM_GPU_HIFI_AUTO_PROFILE'],'1')

if __name__=='__main__':unittest.main()
