"""배치 재배정의샘플보존·기관경계·상한·우선순위계약을검사함."""
import copy,unittest
from handoff.batch_priority import rebatch
class Tests(unittest.TestCase):
 def rows(self):
  specs=[('s1','001','short_read','wgs'),('s2','001','short_read','wgs'),('long_wgs','001','long_read','wgs'),('long_pair','002','long_read','wgs_rna_pair'),('s3','002','short_read','wgs'),('s4','002','short_read','wgs')]
  rows=[dict(sample_id=n,batch_id=b,read_type=t,data_group=g,assigned_lab='MedAI',estimated_model_bytes='10',known_input_total_bytes='100',batch_max_samples='3',batch_target_bytes='100') for n,b,t,g in specs]
  rows.append(dict(rows[0],sample_id='tgil',assigned_lab='TGIL'))
  return rows
 def test_기관배정원본인원보존과pair우선(self):
  rows=self.rows();saved=copy.deepcopy(rows);result=rebatch(rows)
  self.assertEqual(rows,saved);self.assertEqual({r['sample_id'] for r in rows},{r['sample_id'] for r in result});self.assertEqual(result[-1],saved[-1])
  first=[r for r in result if r['assigned_lab']=='MedAI' and r['batch_id']=='001'];self.assertEqual(len(first),3);self.assertEqual([r['sample_id'] for r in first[:2]],['long_pair','long_wgs']);self.assertTrue(all(r['batch_id']=='001' for r in result if r['read_type']=='long_read'))
 def test_재실행결정성(self):
  a=rebatch(self.rows());self.assertEqual(a,rebatch(a))
 def test_용량상한거부(self):
  rows=self.rows()
  for r in rows:r['batch_target_bytes']='20'
  with self.assertRaises(ValueError):rebatch(rows)
 def test_첫배치long수용불가거부(self):
  rows=self.rows()
  for r in rows[:4]:r['read_type']='long_read'
  with self.assertRaises(ValueError):rebatch(rows)
 def test_중복거부(self):
  rows=self.rows();rows.append(dict(rows[0]))
  with self.assertRaises(ValueError):rebatch(rows)
