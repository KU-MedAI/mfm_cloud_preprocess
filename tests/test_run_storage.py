"""중앙경로선택·명시경로보존·실행루트·이전배정선택의계약을검사함."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import hashlib,json,tempfile,unittest,csv
from handoff.storage_paths import storage_defaults
from handoff import cli
class Tests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name);self.package=self.root/'process_code/original';self.package.mkdir(parents=True);(self.package/'09_personal_fasta.sh').touch();p=self.root/'metadata/storage';p.mkdir(parents=True);self.run=self.root/'preprocessed_data/final/new_run';(p/'paths.json').write_text(json.dumps({'paths':[{'분류':'현재 MedAI 최종 입력','실제경로':str(self.run)},{'분류':'현재 MedAI 원본·캐시','실제경로':str(self.run/'raw')}]}))
 def test_중앙기본값과사용자명시경로(self):
  with patch.dict('os.environ',{'MFM_DATA_ROOT':str(self.root)}):
   c=storage_defaults({'data_root':'./data','output_root':'./output'},self.package);self.assertEqual(c['data_root'],str(self.run/'raw'));self.assertEqual(c['output_root'],str(self.run));self.assertEqual(c['output_layout'],'run_root')
   x=storage_defaults({'data_root':'/explicit_data','output_root':'/explicit_output'},self.package);self.assertEqual(x['data_root'],'/explicit_data');self.assertEqual(x['output_root'],'/explicit_output');self.assertNotIn('output_layout',x)
 def test_run_root는출력공유_state는배치분리(self):
  c={'output_root':str(self.run),'state_root':str(self.root/'state'),'output_layout':'run_root'};r=[{'assigned_lab':'MedAI','batch_id':'001','sample_id':'fixture'}];d=cli.batch_config(c,r);self.assertEqual(d['output_root'],str(self.run));self.assertEqual(d['state_root'],str(self.root/'state/medai_001'));self.assertEqual(c['state_root'],str(self.root/'state'))
 def test_기존portable배치출력유지(self):
  c={'output_root':str(self.run),'state_root':str(self.root/'state')};d=cli.batch_config(c,[{'assigned_lab':'MedAI','batch_id':'001','sample_id':'fixture'}]);self.assertEqual(d['output_root'],str(self.run/'medai_001'))
 def test_snapshot선택도URL은동결master만사용(self):
  folder=self.package/'manifests/분담';folder.mkdir(parents=True);row={'sample_id':'fixture','assigned_lab':'MedAI','batch_id':'002','project_id':'p','platform':'hifi','wgs_url':'trusted'}
  def write(p,r):
   with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(r),delimiter='\t');w.writeheader();w.writerow(r)
  write(folder/'전체_분담배치.tsv',row);frozen=dict(row,batch_id='001');write(folder/'old.tsv',frozen);selection=self.root/'selected.tsv';write(selection,dict(frozen,wgs_url='untrusted'));config=self.root/'config.json';config.write_text(json.dumps({'assignment_snapshot':{'relative_path':'manifests/분담/old.tsv','sha256':hashlib.sha256((folder/'old.tsv').read_bytes()).hexdigest()}}))
  with patch.object(cli,'ROOT',self.package):
   with self.assertRaises(ValueError):cli.selection(SimpleNamespace(selection=str(selection)))
   self.assertEqual(cli.selection(SimpleNamespace(selection=str(selection),effective_config=str(config)))[0]['wgs_url'],'trusted')
