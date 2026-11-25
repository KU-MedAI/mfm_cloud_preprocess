"""실제 Stage01의 RNA 보류와 필수 입력 검증을 검사함."""
from pathlib import Path
import os
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]


class Stage01RnaDeferTests(unittest.TestCase):
    def test_deferred_modes_skip_missing_rna_but_required_modes_still_fail(self):
        with tempfile.TemporaryDirectory() as raw:
            base=Path(raw);missing=base/'not_downloaded_rna.bam'
            env=dict(os.environ,PIPELINE_MANIFEST='',RNA_PATH=str(missing),
                DEFAULT_PLATFORM='hifi',DEFAULT_READ_TYPE='long_read',OUTPUT_ROOT=str(base/'output'),
                OUTPUT_ROOT_OVERRIDE='1',INPUT_ROOT=str(base/'input'),
                PYTHON=sys.executable,PYTHON_BIN=sys.executable,
                MFM_CORE_ENV=str(Path(sys.executable).parent.parent),
                MFM_CPU_BUDGET='2',MFM_SAMPLE_CPU_BUDGET='2',
                MFM_WGS_BACKEND='cpu',MFM_RNA_BACKEND='cpu',MFM_RNA_OUTPUT_MODE='full')
            for mode,deferred,success in [('pileup_feature','1',True),('sequence_only','1',True),
                                          ('pileup_feature','0',False),('full','1',False)]:
                with self.subTest(mode=mode,deferred=deferred):
                    result=subprocess.run(['bash',str(ROOT/'01_prepare_analysis_bam.sh'),'S1'],cwd=ROOT,
                        env=dict(env,MFM_LONG_READ_OUTPUT_MODE=mode,MFM_SEQUENCE_DEFER_RNA=deferred),
                        capture_output=True,text=True,timeout=20)
                    self.assertEqual(result.returncode==0,success,result.stdout+result.stderr)
                    if success:self.assertIn('legacy 실행',result.stdout+result.stderr)
                    else:self.assertIn(str(missing),result.stdout+result.stderr)
            self.assertFalse(missing.exists())


if __name__=='__main__':unittest.main()
