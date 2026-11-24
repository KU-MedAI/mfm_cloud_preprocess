"""실제 CPU time과 예약 코어를 분리해 기록하는지 검증함."""
from pathlib import Path
import json,os,subprocess,sys,tempfile,unittest
ROOT=Path(__file__).resolve().parents[1]


class UsageTests(unittest.TestCase):
    def test_nested_child_cpu_and_sampled_rss(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);audit=root/'usage.jsonl'
            code='import time; data=bytearray(16000000); end=time.monotonic()+0.8\nwhile time.monotonic()<end: pass\ntime.sleep(0.3)'
            command=[sys.executable,str(ROOT/'scripts/cpu_pool.py'),'--pool-dir',str(root/'pool'),'--total-cpus','2','--sample','test','--branch','WGS','--stage','11','--audit',str(audit),'--',sys.executable,'-c',code]
            env=dict(os.environ,MFM_RAM_LIMIT_BYTES='1000000000000',MFM_CPU_UTILIZATION='maximum')
            subprocess.run(command,env=env,capture_output=True,text=True,check=True)
            row=json.loads(audit.read_text());self.assertEqual(row['granted_cpus'],1)
            self.assertGreater(row['actual_cpu_seconds'],0.5);self.assertLess(row['actual_cpu_seconds'],2)
            self.assertGreater(row['average_used_cores'],0.3);self.assertLess(row['average_used_cores'],1.5)
            self.assertGreater(row['peak_tree_rss_observed_bytes'],16000000)
            self.assertGreater(row['peak_single_process_rss_bytes'],16000000)


if __name__=='__main__':unittest.main()
