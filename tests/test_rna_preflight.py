from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "rna_preflight.py"
SAMTOOLS = "samtools"


class RnaPreflightTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.bam = self.root / "rna.bam"
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 1000}],
            "RG": [{"ID": "rg1", "SM": "S1"}],
        }
        with pysam.AlignmentFile(self.bam, "wb", header=header) as out:
            for name, flag, start in (("pair", 99, 100), ("pair", 147, 180)):
                read = pysam.AlignedSegment()
                read.query_name = name
                read.query_sequence = "A" * 50
                read.query_qualities = pysam.qualitystring_to_array("I" * 50)
                read.flag = flag
                read.reference_id = 0
                read.reference_start = start
                read.mapping_quality = 60
                read.cigar = ((0, 50),)
                read.next_reference_id = 0
                read.next_reference_start = 180 if flag == 99 else 100
                read.template_length = 130 if flag == 99 else -130
                read.set_tag("RG", "rg1")
                out.write(read)
        pysam.index(str(self.bam))
        self.gtf = self.root / "genes.gtf"
        self.gtf.write_text('chr1\ttest\texon\t101\t150\t.\t+\t.\tgene_id "G1";\n')
        self.haps = self.root / "haps"
        for hap in ("HAP1", "HAP2"):
            directory = self.haps / hap
            directory.mkdir(parents=True)
            (directory / "chr1.fa").write_text(">chr1\n" + "A" * 1000 + "\n")
            (directory / "chr1.fa.fai").write_text("chr1\t1000\t6\t1000\t1001\n")
            (directory / "chr1.chain").write_text(
                "chain 0 chr1 1000 + 0 1000 chr1 1000 + 0 1000 1\n1000\n"
            )
        self.output = self.root / "rna_preflight.json"

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def command(
        self, sample: str = "S1", bam: Path | None = None,
        ratio_policy: str = "record_only",
    ) -> list[str]:
        return [
            sys.executable, str(SCRIPT), "--sample", sample,
            "--rna", str(bam or self.bam), "--gtf", str(self.gtf),
            "--hap-root", str(self.haps), "--contigs", "chr1",
            "--samtools", SAMTOOLS, "--hap-read-mode", "all_to_both",
            "--min-paired-fraction", "0.90", "--ratio-qc-policy", ratio_policy,
            "--output", str(self.output),
        ]

    def low_paired_bam(self) -> Path:
        path = self.root / "low_paired.bam"
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 1000}],
            "RG": [{"ID": "rg1", "SM": "S1"}],
        }
        with pysam.AlignmentFile(path, "wb", header=header) as out:
            for index in range(8):
                read = pysam.AlignedSegment()
                read.query_name = f"single{index}"
                read.query_sequence = "A" * 50
                read.query_qualities = pysam.qualitystring_to_array("I" * 50)
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 10 + index * 60
                read.mapping_quality = 60
                read.cigar = ((0, 50),)
                read.set_tag("RG", "rg1")
                out.write(read)
            for name, flag, start in (("pair", 99, 600), ("pair", 147, 680)):
                read = pysam.AlignedSegment()
                read.query_name = name
                read.query_sequence = "A" * 50
                read.query_qualities = pysam.qualitystring_to_array("I" * 50)
                read.flag = flag
                read.reference_id = 0
                read.reference_start = start
                read.mapping_quality = 60
                read.cigar = ((0, 50),)
                read.next_reference_id = 0
                read.next_reference_start = 680 if flag == 99 else 600
                read.template_length = 130 if flag == 99 else -130
                read.set_tag("RG", "rg1")
                out.write(read)
        pysam.index(str(path))
        return path

    def test_contract_changes_and_sm_mismatch_fails(self) -> None:
        subprocess.run(self.command(), check=True, stdout=subprocess.PIPE, text=True)
        first = json.loads(self.output.read_text())
        self.assertEqual(first["status"], "pass")
        self.assertEqual(first["checks"]["paired_fraction"], 1.0)
        self.gtf.write_text(self.gtf.read_text() + "# changed\n")
        subprocess.run(self.command(), check=True, stdout=subprocess.PIPE, text=True)
        second = json.loads(self.output.read_text())
        self.assertNotEqual(first["contract_sha256"], second["contract_sha256"])
        wrong = subprocess.run(self.command(sample="OTHER"), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertNotEqual(wrong.returncode, 0)
        self.assertIn("SM mismatch", wrong.stderr)

    def test_replaced_input_invalidates_quickcheck_cache_and_is_rechecked(self) -> None:
        cache = self.root / 'old_quickcheck.json'
        cache.write_text(json.dumps({'schema':'mfm_stage01_quickcheck_v1','sample':'S1',
                                    'status':'pass','input':{'path':'old_input.bam'}}))
        result = subprocess.run(self.command()+['--quickcheck-json',str(cache)],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('현재 입력 재검사',result.stderr)
        self.assertEqual(json.loads(self.output.read_text())['status'],'pass')
        mismatch = subprocess.run(self.command(sample='OTHER')+['--quickcheck-json',str(cache)],capture_output=True,text=True)
        self.assertNotEqual(mismatch.returncode,0)
        self.assertIn('SM mismatch',mismatch.stderr)
    def test_missing_cache_rechecks_and_rejects_corrupt_input(self) -> None:
        cache = self.root/'missing_quickcheck.json'
        result = subprocess.run(self.command()+['--quickcheck-json',str(cache)],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        broken=self.root/'broken.bam';broken.write_bytes(self.bam.read_bytes()[:30])
        shutil.copy2(Path(str(self.bam)+'.bai'),Path(str(broken)+'.bai'))
        result = subprocess.run(self.command(bam=broken)+['--quickcheck-json',str(cache)],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)

    def test_truncated_bam_is_rejected(self) -> None:
        broken = self.root / "broken.bam"
        data = self.bam.read_bytes()
        broken.write_bytes(data[: max(1, len(data) // 2)])
        shutil.copy2(Path(str(self.bam) + ".bai"), Path(str(broken) + ".bai"))
        result = subprocess.run(self.command(bam=broken), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertNotEqual(result.returncode, 0)

    def test_low_paired_fraction_is_recorded_but_not_filtered(self) -> None:
        low = self.low_paired_bam()
        subprocess.run(self.command(bam=low), check=True, stdout=subprocess.PIPE, text=True)
        checks = json.loads(self.output.read_text())["checks"]
        self.assertEqual(checks["paired_fraction"], 0.2)
        self.assertTrue(checks["paired_fraction_below_reference"])
        self.assertEqual(checks["paired_fraction_qc_action"], "record_only")
        strict = subprocess.run(
            self.command(bam=low, ratio_policy="strict"),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
        self.assertNotEqual(strict.returncode, 0)
        self.assertIn("paired fraction below threshold", strict.stderr)


if __name__ == "__main__":
    unittest.main()
