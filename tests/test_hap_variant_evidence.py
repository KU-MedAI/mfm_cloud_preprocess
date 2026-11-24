#!/usr/bin/env python3
"""20 HAP variant 검증 + 21 근거 규칙 (2026-10-02).

- 21: 검증된 HAP variant 를 덮는 fragment 만 allele 증거. RPS15 형(한쪽만 unique, variant 없음)·PKN1 형(unsupported
      indel 위) fragment 는 버리고, ZDBF2 형(validated SNV 위)은 센다.
- 20: WGS read 로 SNV 는 validated, read 에 없는 결실은 unsupported 로 판정하고 HAP 좌표 구간을 chain 으로 낸다.
"""
from __future__ import annotations

import csv, json, random, subprocess, sys, tempfile, unittest
from pathlib import Path

import pysam

WORK = Path(__file__).resolve().parent
CODE = WORK if (WORK / "21_allelic_counts.py").is_file() else WORK.parent
sys.path.insert(0, str(WORK))
from test_allelic_outputs import write_hap_variants  # noqa: E402


def write_bam(path: Path, reads: list[tuple[str, int, int]]) -> None:
    """(name, start, AS) — 50 bp read1, MAPQ 255."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    unsorted = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for name, start, score in reads:
            r = pysam.AlignedSegment()
            r.query_name = name; r.query_sequence = "A" * 50; r.flag = 64
            r.reference_id = 0; r.reference_start = start; r.mapping_quality = 255
            r.cigar = ((0, 50),); r.query_qualities = pysam.qualitystring_to_array("I" * 50)
            r.set_tag("AS", score); out.write(r)
    pysam.sort("-o", str(path), str(unsorted)); unsorted.unlink()


class EvidenceRuleTest(unittest.TestCase):
    def test_only_validated_variant_fragments_count(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            # frag_good    : 두 HAP, HAP1 점수 높음, validated SNV(120) 위        → aCount
            # frag_only_var: HAP2 에서만 unique, validated SNV(120) 위           → bCount (HAP 한쪽에만 있는 서열 같은 경우)
            # frag_only_nov: HAP2 에서만 unique, variant 없는 곳(300)  [RPS15 형] → 버림
            # frag_bad     : 두 HAP, HAP2 점수 높음, unsupported 결실(620) 위 [PKN1 형] → 버림
            # frag_novar   : 두 HAP, 점수 다름, variant 없는 곳(800)             → 버림
            write_bam(root / "h1.bam", [("frag_good", 100, 50), ("frag_bad", 600, 40), ("frag_novar", 780, 50)])
            write_bam(root / "h2.bam", [("frag_good", 100, 48), ("frag_only_var", 100, 50), ("frag_only_nov", 300, 50),
                                        ("frag_bad", 600, 50), ("frag_novar", 780, 46)])
            for name in ("e1.tsv", "e2.tsv"):
                (root / name).write_text("G1\tchr1\t0\t1000\n")
            write_hap_variants(root / "hv.tsv", [{}, dict(pos=621, ref="TC", alt="T", gt="1|0", kind="INDEL", status="unsupported",
                                                         wgs_alt=0, wgs_alt_frac=0.0, h1_start=620, h1_end=621, h1_allele="T",
                                                         h2_start=620, h2_end=622, h2_allele="TC")])
            out = root / "allelic_counts.tsv"
            subprocess.run([sys.executable, str(CODE / "21_allelic_counts.py"), "--hap-variants", str(root / "hv.tsv"),
                            "--bam1", str(root / "h1.bam"), "--bam2", str(root / "h2.bam"),
                            "--exons1", str(root / "e1.tsv"), "--exons2", str(root / "e2.tsv"), "--out", str(out)],
                           check=True, capture_output=True)
            row = next(csv.DictReader(out.open(), delimiter="\t"))
            self.assertEqual((row["aCount"], row["bCount"]), ("1", "1"))
            meta = json.loads((root / "allelic_counts.meta.json").read_text())
            self.assertEqual(meta["schema"], "mfm_allelic_counts_meta_v5")
            self.assertEqual(meta["fragments_dropped_unvalidated_variant"], 1)
            self.assertEqual(meta["fragments_dropped_no_validated_variant"], 2)
            self.assertEqual(meta["hap_variants_status_counts"], {"validated": 1, "unsupported": 1})


class HapVariantStageTest(unittest.TestCase):
    def test_snv_validated_and_absent_deletion_unsupported(self):
        random.seed(1)
        ref = "".join(random.choice("ACGT") for _ in range(400))
        ref = ref[:100] + "A" + ref[101:200] + "TC" + ref[202:]          # SNV pos 101 A>G, 결실 pos 201 TC>T
        hap1 = ref[:201] + ref[202:]                                        # HAP1 = 결실 (GT 1|0)
        hap2 = ref[:100] + "G" + ref[101:]                                  # HAP2 = SNV (GT 0|1)
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            (root / "ref.fa").write_text(f">chr1\n{ref}\n"); pysam.faidx(str(root / "ref.fa"))
            for h, seq, chain in ((1, hap1, "201 1 0\n198\n"), (2, hap2, "400\n")):
                d = root / f"HAP{h}"; d.mkdir()
                (d / "chr1.fa").write_text(f">chr1\n{seq}\n"); pysam.faidx(str(d / "chr1.fa"))
                (d / "chr1.chain").write_text(f"chain 1 chr1 400 + 0 400 chr1 {len(seq)} + 0 {len(seq)} 1\n{chain}\n")
            vcf = root / "v.vcf"
            vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=400>\n"
                           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                           "chr1\t101\t.\tA\tG\t.\t.\t.\tGT\t0|1\n"
                           "chr1\t201\t.\tTC\tT\t.\t.\t.\tGT\t1|0\n")
            pysam.tabix_index(str(vcf), preset="vcf", force=True)
            header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 400}]}
            bam = root / "wgs.bam"
            with pysam.AlignmentFile(bam, "wb", header=header) as out:
                for i in range(20):                                         # 결실 없는 read 20 개, 그중 10 개가 SNV alt
                    seq = hap2 if i % 2 else ref
                    r = pysam.AlignedSegment(); r.query_name = f"r{i}"; r.flag = 0
                    r.reference_id = 0; r.reference_start = 50; r.mapping_quality = 60
                    r.query_sequence = seq[50:300]; r.cigar = ((0, 250),)
                    r.query_qualities = pysam.qualitystring_to_array("I" * 250); r.set_tag("HP", 2 if i % 2 else 1)
                    out.write(r)
            pysam.index(str(bam))
            (root / "g.gtf").write_text('chr1\tT\texon\t1\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
            res = subprocess.run([sys.executable, str(CODE / "20_hap_variants.py"), "--vcf", str(vcf) + ".gz",
                                  "--wgs-bam", str(bam), "--hap-dir1", str(root / "HAP1"), "--hap-dir2", str(root / "HAP2"),
                                  "--gtf", str(root / "g.gtf"), "--contigs", "chr1", "--threads", "1",
                                  "--out", str(root / "hv.tsv"), "--qc-json", str(root / "hv.json")],
                                 capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, res.stderr)
            rows = {int(r["pos"]): r for r in csv.DictReader((root / "hv.tsv").open(), delimiter="\t")}
            self.assertEqual(rows[101]["status"], "validated")
            self.assertEqual(rows[101]["wgs_alt"], "10")
            self.assertEqual(rows[201]["status"], "unsupported")
            # HAP 좌표: SNV 는 두 HAP 모두 100, 결실은 HAP1 에서 anchor 1 bp, HAP2 에서 ref 2 bp
            self.assertEqual((rows[101]["h1_start"], rows[101]["h2_start"]), ("100", "100"))
            self.assertEqual((rows[201]["h1_start"], rows[201]["h1_end"], rows[201]["h2_end"]), ("200", "201", "202"))
            qc = json.loads((root / "hv.json").read_text())
            self.assertEqual(qc["variants"], 2)


if __name__ == "__main__":
    unittest.main()
