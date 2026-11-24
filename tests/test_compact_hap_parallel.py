"""HAP 병렬 작업·압축·해독을 실제 CPU 예산 안에서 계산하는지 확인함."""
from unittest.mock import patch
import os,unittest
from integrate_wgs_rna.compact.resources import plan
class HapBudgetTests(unittest.TestCase):
 def test_68CPU에서17염색체34HAP를동시에예약함(self):
  with patch.dict(os.environ,{'MFM_CPU_BUDGET':'68','MFM_BLOSC_THREADS':'2'},clear=True),patch('os.sched_getaffinity',return_value=set(range(72))),patch('integrate_wgs_rna.compact.resources.cgroup_cpu_limit',return_value=68):
   p=plan(68,22,hap_split=True)
   self.assertEqual((p['workers'],p['hap_jobs'],p['alignment_threads'],p['cpu_per_worker']),(17,2,2,4));self.assertEqual(p['workers']*p['cpu_per_worker'],68)
 def test_작은grant와명시순차선택도지킴(self):
  for cpus in [1,2,3,4,5,16,68]:
   with self.subTest(cpus=cpus),patch.dict(os.environ,{'MFM_CPU_POOL_GRANTED':str(cpus)},clear=True):
    p=plan(68,22,hap_split=True);self.assertLessEqual(p['workers']*p['cpu_per_worker'],cpus)
  with patch.dict(os.environ,{'MFM_CPU_BUDGET':'68','MFM_COMPACT_HAP_JOBS':'1'},clear=True):self.assertEqual(plan(22,22,hap_split=True)['hap_jobs'],1)
if __name__=='__main__':unittest.main()
