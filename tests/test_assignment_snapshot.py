"""동결배정의SHA검증과신규배정경계를확인함."""
from pathlib import Path
from types import SimpleNamespace
import tempfile,unittest,json,hashlib
from handoff.assignment_snapshot import master_path
class Tests(unittest.TestCase):
 def setUp(self):
  self.t=tempfile.TemporaryDirectory();self.addCleanup(self.t.cleanup);self.r=Path(self.t.name);self.p=self.r/'manifests/분담/old.tsv';self.p.parent.mkdir(parents=True);self.p.write_text('frozen');self.c=self.r/'config.json';self.spec={'relative_path':'manifests/분담/old.tsv','sha256':hashlib.sha256(self.p.read_bytes()).hexdigest()};self.c.write_text(json.dumps({'assignment_snapshot':self.spec}))
 def test_명시effective_config만동결배정사용(self):
  self.assertEqual(master_path(self.r,SimpleNamespace(effective_config=str(self.c))),self.p)
  self.assertEqual(master_path(self.r,SimpleNamespace(config=str(self.c))),self.r/'manifests/분담/전체_분담배치.tsv')
 def test_변경된snapshot차단(self):
  self.p.write_text('tampered')
  with self.assertRaises(ValueError):master_path(self.r,SimpleNamespace(effective_config=str(self.c)))
 def test_분담밖경로차단(self):
  self.spec['relative_path']='../../outside';self.c.write_text(json.dumps({'assignment_snapshot':self.spec}))
  with self.assertRaises(ValueError):master_path(self.r,SimpleNamespace(effective_config=str(self.c)))
