#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pysam

WORK = Path(__file__).resolve().parent
CODE = WORK if (WORK / "21_allelic_counts.py").is_file() else WORK.parent


def write_bam(path: Path, alignments: list[tuple[str, int]]) -> None:
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for query_name, score in alignments:
            read = pysam.AlignedSegment()
            read.query_name = query_name
            read.query_sequence = "A" * 50
            read.flag = 64
            read.reference_id = 0
            read.reference_start = 100
            read.mapping_quality = 255
            read.cigar = ((0, 50),)
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            read.set_tag("AS", score)
            out.write(read)



HV_COLS = ["contig", "pos", "ref", "alt", "gt", "kind", "status", "wgs_depth", "wgs_alt", "wgs_alt_frac",
           "hp1_alt", "hp1_depth", "hp2_alt", "hp2_depth", "h1_start", "h1_end", "h1_allele", "h2_start", "h2_end", "h2_allele", "genes"]


def write_hap_variants(path: Path, rows: list[dict]) -> None:
    """20_hap_variants.py 출력 형식. 기본값 = chr1:120 의 validated SNV (두 HAP 같은 좌표)."""
    base = dict(contig="chr1", pos=121, ref="A", alt="C", gt="0|1", kind="SNV", status="validated", wgs_depth=30,
                wgs_alt=15, wgs_alt_frac=0.5, hp1_alt=0, hp1_depth=15, hp2_alt=15, hp2_depth=15,
                h1_start=120, h1_end=121, h1_allele="A", h2_start=120, h2_end=121, h2_allele="C", genes="G1")
    with path.open("w") as fh:
        fh.write("\t".join(HV_COLS) + "\n")
        for r in rows:
            d = {**base, **r}
            fh.write("\t".join(str(d[k]) for k in HV_COLS) + "\n")

class TestAllelicOutputs(unittest.TestCase):
    def test_sidecar_counts_and_provenance(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            bam1, bam2 = root / "h1.bam", root / "h2.bam"
            write_bam(bam1, [("frag_h1", 50), ("frag_h2", 40),
                             ("frag_tie", 45), ("frag_only1", 45)])
            write_bam(bam2, [("frag_h1", 40), ("frag_h2", 50),
                             ("frag_tie", 45)])
            for name in ("e1.tsv", "e2.tsv"):
                (root / name).write_text("G1\tchr1\t0\t1000\n")
            out = root / "allelic_counts.tsv"
            write_hap_variants(root / "hv.tsv", [{}])
            subprocess.run([
                sys.executable, str(CODE / "21_allelic_counts.py"), "--hap-variants", str(root / "hv.tsv"),
                "--bam1", str(bam1), "--bam2", str(bam2),
                "--exons1", str(root / "e1.tsv"), "--exons2", str(root / "e2.tsv"),
                "--out", str(out),
            ], check=True)
            with out.open() as fh:
                main = next(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual((main["aCount"], main["bCount"]), ("2", "1"))
            with (root / "allelic_counts.both_haps.tsv").open() as fh:
                both = next(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual((both["aCount"], both["bCount"]), ("1", "1"))
            with (root / "allelic_counts.split_half.tsv").open() as fh:
                split = next(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(int(split["split0_totalCount"]) + int(split["split1_totalCount"]), 3)
            meta = json.loads((root / "allelic_counts.meta.json").read_text())
            self.assertEqual(meta["schema"], "mfm_allelic_counts_meta_v5")
            self.assertEqual(meta["fragments_assigned"], 3)
            self.assertEqual(meta["fragments_tie"], 1)
            self.assertEqual(meta["fragments_only_hap1"], 1)
            expected = {
                "counts": out,
                "split_half": root / "allelic_counts.split_half.tsv",
                "both_haps": root / "allelic_counts.both_haps.tsv",
                "mmap_adjusted": root / "allelic_counts.mmap_adjusted.tsv",
            }
            self.assertEqual(
                meta["output_sha256"],
                {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in expected.items()},
            )


if __name__ == "__main__":
    unittest.main()
