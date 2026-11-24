"""Illumina 선택 옵션의 기본 보존·4GPU 범위·실행 설정 전달을 검증함."""
from pathlib import Path
from unittest.mock import patch
import json,os,tempfile,unittest
from integrate_wgs_rna.gpu_wgs import commands
from integrate_wgs_rna.stage_gpu_config import apply_stage_environment

class BwaTuningTests(unittest.TestCase):
 def command(self,platform='illumina',gpus=4):
  return commands('/ref',['/r1','/r2'],'/out',platform,'map-hifi',67,'S',gpus)
 def test_검증된4GPU기본은두옵션을함께사용함(self):
  with patch.dict(os.environ,{},clear=True):c=self.command()
  self.assertIn('--cigar-on-gpu',c);self.assertIn('--run-partition',c)
  for count in [1,2,3,8]:
   with patch.dict(os.environ,{},clear=True):other=self.command(gpus=count)
   self.assertNotIn('--cigar-on-gpu',other);self.assertNotIn('--run-partition',other)
 def test_CIGAR만선택가능함(self):
  with patch.dict(os.environ,{'MFM_GPU_CIGAR_ON_GPU':'1'},clear=True):c=self.command()
  self.assertIn('--cigar-on-gpu',c);self.assertNotIn('--run-partition',c)
 def test_4GPU를2GPU씩나누어두옵션을선택함(self):
  with patch.dict(os.environ,{'MFM_GPU_CIGAR_ON_GPU':'1','MFM_GPU_BWA_PARTITION':'1'},clear=True):c=self.command()
  for key in ['--cigar-on-gpu','--run-partition']:self.assertIn(key,c)
  self.assertEqual(c[c.index('--num-gpus')+1],'4')
  self.assertEqual(c[c.index('--bwa-gpu-num-per-partition')+1],'2')
 def test_잘못된값과홀수장치분할을차단함(self):
  for key in ['MFM_GPU_CIGAR_ON_GPU','MFM_GPU_BWA_PARTITION']:
   with patch.dict(os.environ,{key:'yes'},clear=True),self.assertRaises(ValueError):self.command()
  for count in [1,3]:
   with patch.dict(os.environ,{'MFM_GPU_BWA_PARTITION':'1'},clear=True),self.assertRaises(ValueError):self.command(gpus=count)
 def test_HiFi명령에는Illumina옵션을전달하지않음(self):
  with patch.dict(os.environ,{'MFM_GPU_CIGAR_ON_GPU':'1','MFM_GPU_BWA_PARTITION':'1'},clear=True):c=self.command('hifi')
  self.assertNotIn('--cigar-on-gpu',c);self.assertNotIn('--run-partition',c)
 def test_새단계설정으로두옵션을전달함(self):
  with tempfile.TemporaryDirectory() as t:
   p=Path(t);m=p/'sample_manifests/S/input.tsv';m.parent.mkdir(parents=True)
   opts={'MFM_GPU_CIGAR_ON_GPU':'1','MFM_GPU_BWA_PARTITION':'1'}
   (p/'effective_config.json').write_text(json.dumps({'run_id':'R','data_root':str(p/'raw'),'output_root':str(p/'out'),'pipeline_env':opts}))
   env={'RUN_ID':'R','INPUT_ROOT':str(p/'raw'),'OUTPUT_ROOT':str(p/'out'),'PIPELINE_MANIFEST':str(m),'MFM_WGS_BACKEND':'parabricks'}
   apply_stage_environment(env)
   for k,v in opts.items():self.assertEqual(env[k],v)

if __name__=='__main__':unittest.main()
