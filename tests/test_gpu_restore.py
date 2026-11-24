"""네이티브 태그 복원의 모든 필드·순서·오류 검출을 CPU 구현과 대조함."""
from pathlib import Path
from unittest.mock import patch
import os,tempfile,unittest,pysam
from integrate_wgs_rna.gpu_wgs import restore_bam
class RestoreTests(unittest.TestCase):
 def test_두구현의전체SAM필드와순서가동일함(self):
  with tempfile.TemporaryDirectory() as t:
   root=Path(t);src=root/'in.bam';header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000}],'RG':[{'ID':'old','SM':'old'}]}
   with pysam.AlignmentFile(str(src),'wb',header=header) as f:
    for i in range(20):
     r=pysam.AlignedSegment(f.header);r.query_name=f'r{i}';r.query_sequence='ACGT'*10;r.query_qualities=pysam.qualitystring_to_array('I'*40);r.flag=1024 if i%3==0 else 0;r.reference_id=0;r.reference_start=i;r.mapping_quality=60;r.cigarstring='40M'
     for k,v in [('RG','old'),('HP',i%2+1),('PS',17),('OD',i%2),('XU',1)]:r.set_tag(k,v)
     f.write(r)
   for policy in [None,'none','preserve_source']:
    records=[];headers=[]
    for backend in ['python','native']:
     out=root/f'{policy}_{backend}.bam'
     with patch.dict(os.environ,{'MFM_GPU_TAG_RESTORE':backend,'MFM_GPU_RESTORE_DECODE_THREADS':'2'}):restore_bam(src,'@RG\\tID:S\\tSM:S','S',policy,out)
     with pysam.AlignmentFile(str(out)) as f:headers.append(str(f.header));records.append([r.to_string() for r in f])
    self.assertEqual(headers[0],headers[1]);self.assertEqual(records[0],records[1])
 def test_최대비교의출력thread도68CPU안에포함함(self):
  from integrate_wgs_rna.realign.config import PersonalRealignConfig
  from integrate_wgs_rna.gpu_wgs import write_script
  with tempfile.TemporaryDirectory() as t,patch.dict(os.environ,{'RUN_ID':'gpu_compare_output','MFM_GPU_COMPARE_UNRESERVED':'1','MFM_WGS_BACKEND':'parabricks','MFM_CPU_BUDGET':'68'},clear=True):
   p=Path(t);target=p/'target.fa';target.write_text('>chr1\nACGT\n');pysam.faidx(str(target));bam=p/'in.bam';bam.touch()
   config=PersonalRealignConfig(cohort='test',output_root=p,sample_list=None,samples=['S'],chroms=['chr1'],sequence_variant_policy='snv_indel',dataset_role='pre_training',platform='hifi',threads=32,sort_threads=4,sort_mem='128M',hap_jobs=1,input_decode_threads=2)
   script=write_script(config,'S',1,bam,target,p/'out.cram',['chr1'],p,p,None).read_text()
   self.assertIn('export MFM_GPU_OUTPUT_CPUS=68',script);self.assertIn('--out-cram',script)
   from integrate_wgs_rna.gpu_cram import plan
   self.assertLessEqual(plan(68,44)['native_cpu_limit'],68)
if __name__=='__main__':unittest.main()
