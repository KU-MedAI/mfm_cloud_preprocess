"""출력QC의 실제grant·동시HAP·RAM·직렬fallback 계획을 검증함."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import os,sys,tempfile,unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from integrate_wgs_rna.realign.output_scan import output_scan_plan


class ResourceTest(unittest.TestCase):
    def test_grant_hap_division_small_file_and_memory(self):
        with tempfile.TemporaryDirectory() as raw:
            p=Path(raw)/'input.cram';p.write_bytes(b'fixture')
            config=SimpleNamespace(jobs=1,hap_jobs=2,sort_threads=8)
            with patch.dict(os.environ,{'MFM_CPU_BUDGET':'6'},clear=True):
                plan=output_scan_plan(config,p,[('chr1',100),('chr2',100)],True,72)
                self.assertEqual(plan['workers_used'],3);self.assertEqual(plan['cpu_budget'],3)
                self.assertEqual(plan['decode_threads'],0)
                self.assertEqual(output_scan_plan(config,p,[('chr1',100)],True,None)['workers_used'],1)
                self.assertEqual(output_scan_plan(config,p,[('chr1',100)],False,72)['workers_used'],1)
            with patch.dict(os.environ,{'MFM_RAM_LIMIT_BYTES':str(128*1024**2)},clear=True):
                with self.assertRaisesRegex(ValueError,'RAM'):output_scan_plan(config,p,[('chr1',100)],True,2)


if __name__=='__main__':unittest.main()
