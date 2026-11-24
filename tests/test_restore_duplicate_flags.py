"""OD 태그 복원 결과가 기존 pysam 경로와 같음을 확인한다."""
import json
import os
import subprocess
import sys
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "src/integrate_wgs_rna/restore_duplicate_flags.py"
SAMTOOLS = Path(sys.executable).parent / "samtools"
SAM = (
    b"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n"
    b"a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n"
    b"b\t16\tchr1\t11\t60\t4M\t*\t0\t0\tACGT\tIIII\tOD:i:1\tNM:i:1\n"
    b"c\t1024\tchr1\t12\t60\t4M\t*\t0\t0\tACGT\tIIII\tOD:i:0\n"
    b"d\t99\tchr1\t13\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tOD:i:1\n"
)


class RestoreDuplicateFlagsTests(unittest.TestCase):
    def run_engine(self, mode, engine):
        return subprocess.run(
            [sys.executable, "-B", str(SCRIPT), "--mode", mode], input=SAM,
            capture_output=True, check=True,
            env={**os.environ, "MFM_DUP_RESTORE_ENGINE": engine})

    @staticmethod
    def records(stream):
        result = subprocess.run([str(SAMTOOLS), "view", "-h", "-"], input=stream,
                                capture_output=True, check=True)
        return [line for line in result.stdout.splitlines() if not line.startswith(b"@")]

    def test_두_모드의_레코드와_집계_일치(self):
        for mode in ("preserve_source", "none"):
            with self.subTest(mode=mode):
                old = self.run_engine(mode, "pysam")
                fast = self.run_engine(mode, "sam")
                self.assertEqual(self.records(old.stdout), self.records(fast.stdout))
                self.assertEqual(json.loads(old.stderr), json.loads(fast.stderr))
                self.assertEqual(json.loads(fast.stderr)["od_records"], 2)

    def test_원시_SAM을_samtools_sort가_읽음(self):
        stream = self.run_engine("preserve_source", "sam").stdout
        result = subprocess.run([str(SAMTOOLS), "sort", "-O", "BAM", "-o", "/dev/null", "-"],
                                input=stream, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr.decode())


if __name__ == "__main__":
    unittest.main()
