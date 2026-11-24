"""완료가드·샘플경계·공용입력·열린파일·재실행을검사함."""
from pathlib import Path
from unittest.mock import patch
import tempfile,json,unittest,hashlib
from handoff.completion_cleanup import cleanup_external,plan
class Tests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.out=Path(self.tmp.name)/'out';self.data=self.out/'raw';self.data.mkdir(parents=True);self.write(self.data/'.mfm_handoff_store.json',{'schema':'mfm_handoff_store_v1'});self.root=self.out/'wgs/P/S1';(self.root/'metadata').mkdir(parents=True);self.write(self.root/'metadata/sample.json',{'sample_id':'S1','source_sample_id':'S1'});self.write(self.root/'metadata/download_inputs.json',[{'path':'public_vcf/P/S1.shared.vcf.gz','shared':True}]);self.write(self.root/'qc/final_package_validation.json',{'sample':'S1','mode':'wgs','status':'pass','failed_checks':[]});self.write(self.root/'qc/delivery_final_validation.json',{'status':'pass'});self.write(self.root/'metadata/cleanup_completed.json',{'schema':'mfm_delivery_cleanup_v1','sample':'S1','mode':'wgs','dry_run':False,'status':'complete','post_cleanup_validation':{'schema':'mfm_delivery_final_validation_v1','sample':'S1','mode':'wgs','status':'pass','errors':[]}})
  layout=self.root/'compact/layout';contig=layout/'chr1';contig.mkdir(parents=True)
  for name in ['hap_sequence.zarr','readset.HAP1.zarr','readset.HAP2.zarr','depth.zarr']:(contig/name).mkdir()
  (contig/'windows.tsv').write_text('model');self.write(layout/'manifest.json',{'status':'complete','validation_level':'full','records':{'chr1':{'contig':'chr1','status':'pass'}}});self.trace=self.root/'qc/stage12_trace/all/trace_manifest.json';self.write(self.trace,{'status':'pass','contigs':{'chr1':'pass'},'compact_root':str(layout)})
  self.raw=self.data/'wgs/P/S1/S1.bam';self.touch(self.raw);records=self.data/'.download_records';key=hashlib.sha256(str(self.raw).encode()).hexdigest();self.write(records/(key+'.json'),{'owned':True,'path':str(self.raw)});self.vcf=self.out/'_work/run/S1/persample_vcf/S1.seq.vcf.gz';self.touch(self.vcf);self.other=self.out/'_work/run/S2/persample_vcf/S2.seq.vcf.gz';self.touch(self.other);self.shared=self.data/'public_vcf/P/S1.shared.vcf.gz';self.touch(self.shared);self.private=self.data/'public_vcf/P/S1.private.vcf.gz';self.touch(self.private);key=hashlib.sha256(str(self.private).encode()).hexdigest();self.write(self.data/'.download_records'/(key+'.json'),{'owned':True,'path':str(self.private)})
 def write(self,p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d))
 def touch(self,p):p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b'input')
 def run_cleanup(self,**kw):return cleanup_external(self.root,'S1','wgs',self.data,self.out,**kw)
 def test_전용입력만삭제_공용모델다른샘플보존(self):
  with patch('handoff.completion_cleanup._opened',return_value=[]):r=self.run_cleanup()
  self.assertEqual(r['status'],'complete');self.assertFalse(self.raw.exists());self.assertFalse(self.vcf.exists());self.assertFalse(self.private.exists());self.assertTrue(self.other.exists());self.assertTrue(self.shared.exists());self.assertEqual((self.root/'compact/layout/chr1/windows.tsv').read_text(),'model')
 def test_재실행안전(self):
  with patch('handoff.completion_cleanup._opened',return_value=[]):self.run_cleanup();r=self.run_cleanup()
  self.assertEqual(r['logical_bytes'],0);self.assertEqual(r['removed'],[])
 def test_소유미등록원본은명시옵션필요(self):
  for p in (self.data/'.download_records').glob('*.json'):p.unlink()
  r=plan(self.root,'S1','wgs',self.data,self.out);self.assertNotIn(str(self.raw),[f['path'] for f in r['files']])
  with patch('handoff.completion_cleanup._opened',return_value=[]):self.run_cleanup(allow_unowned=True)
  self.assertFalse(self.raw.exists())
 def test_최종QC실패면삭제거부(self):
  self.write(self.root/'qc/final_package_validation.json',{'sample':'S1','mode':'wgs','status':'fail'})
  with self.assertRaises(ValueError):self.run_cleanup()
  self.assertTrue(self.raw.exists())
 def test_trace실패면삭제거부(self):
  self.write(self.trace,{'status':'fail'})
  with self.assertRaises(ValueError):self.run_cleanup()
  self.assertTrue(self.raw.exists())
 def test_열린원본삭제거부(self):
  with patch('handoff.completion_cleanup._opened',return_value=[{'pid':123,'path':str(self.raw)}]):
   with self.assertRaises(ValueError):self.run_cleanup()
  self.assertTrue(self.raw.exists())
 def test_샘플raw_symlink거부(self):
  self.raw.unlink();self.raw.symlink_to(self.other)
  with self.assertRaises(ValueError):self.run_cleanup(allow_unowned=True)
  self.assertTrue(self.other.exists())
 def test_work_symlink거부(self):
  self.vcf.unlink();self.vcf.symlink_to(self.other)
  with self.assertRaises(ValueError):self.run_cleanup()
  self.assertTrue(self.other.exists())
 def test_dryrun은파일보존(self):
  r=self.run_cleanup(dry_run=True);self.assertEqual(r['status'],'planned');self.assertTrue(self.raw.exists());self.assertTrue(self.vcf.exists())
 def test_모델필수파일누락거부(self):
  (self.root/'compact/layout/chr1/windows.tsv').unlink()
  with self.assertRaises(ValueError):self.run_cleanup()
  self.assertTrue(self.raw.exists())

 def test_샘플디렉터리내공유원본도보존(self):
  raw=self.data/'wgs/P/S1/shared.bam';self.touch(raw);d=json.loads((self.root/'metadata/download_inputs.json').read_text());d.append({'path':'wgs/P/S1/shared.bam','shared':True});self.write(self.root/'metadata/download_inputs.json',d)
  with patch('handoff.completion_cleanup._opened',return_value=[]):self.run_cleanup(allow_unowned=True)
  self.assertTrue(raw.exists())
 def test_계획후파일변경시삭제거부(self):
  def changed(paths):self.raw.write_bytes(b'changed-input');return []
  with patch('handoff.completion_cleanup._opened',side_effect=changed):
   with self.assertRaises(ValueError):self.run_cleanup()
  self.assertTrue(self.raw.exists());self.assertTrue(self.other.exists())
 def test_실제열린fd탐지(self):
  from handoff.completion_cleanup import _opened
  import os
  with self.raw.open('rb') as f:found=_opened({str(self.raw)})
  self.assertTrue(any(x['pid']==os.getpid() for x in found))

 def test_외부정리실패는완료아님_모델재생성없이재시도(self):
  from handoff.cleanup import cleanup_complete
  with patch('handoff.completion_cleanup._opened',return_value=[{'pid':1}]):
   with self.assertRaises(ValueError):self.run_cleanup()
  self.assertFalse(cleanup_complete(self.root,'S1','wgs'))
  with patch('handoff.completion_cleanup._opened',return_value=[]):self.run_cleanup()
  self.assertTrue(cleanup_complete(self.root,'S1','wgs'))
