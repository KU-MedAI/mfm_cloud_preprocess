#!/usr/bin/env python3
from __future__ import annotations

import errno
import hashlib
import importlib.util
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from integrate_wgs_rna import personal_realign
from integrate_wgs_rna.sample_vcf import (
    PerSampleVcfConfig,
    summarize_sequence_vcf,
    summarize_tag_genotypes,
    validate_seq_vcf,
    validate_tag_vcf,
)

spec = importlib.util.spec_from_file_location("v2_stage17_mark_duplicates", ROOT / "17_mark_duplicates.py")
assert spec and spec.loader
markdup = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = markdup
spec.loader.exec_module(markdup)

SAMTOOLS = shutil.which("samtools") or "/home/aibio_1/miniconda3/envs/jj_mfm/bin/samtools"
BCFTOOLS = shutil.which("bcftools") or "/home/aibio_1/miniconda3/envs/jj_mfm/bin/bcftools"


def make_read(header, name, flag=0, tid=0, pos=10, mtid=-1, mpos=-1, hp=None, mapq=30):
    read = pysam.AlignedSegment(header)
    read.query_name = name
    read.query_sequence = "A" * 10
    read.query_qualities = pysam.qualitystring_to_array("I" * 10)
    read.flag = flag
    read.mapping_quality = mapq
    if flag & 4:
        read.reference_id = -1
        read.reference_start = -1
    else:
        read.reference_id = tid
        read.reference_start = pos
        read.cigar = [(0, 10)]
    if flag & 1:
        read.next_reference_id = mtid
        read.next_reference_start = mpos
    if hp is not None:
        read.set_tag("HP", hp)
    return read


def write_bam(path: Path, reads, *, reference: Path | None = None, index: bool = True) -> None:
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 1000}]}
    mode = "wc" if path.suffix == ".cram" else "wb"
    kwargs = {"reference_filename": str(reference)} if reference else {}
    with pysam.AlignmentFile(str(path), mode, header=header, **kwargs) as out:
        for read in reads(out.header):
            out.write(read)
    if index:
        pysam.index(str(path))


def make_reference(root: Path) -> Path:
    reference = root / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 1000 + "\n>chr2\n" + "A" * 1000 + "\n")
    pysam.faidx(str(reference))
    return reference


def legacy_samtools_count(alignment: Path, reference: Path, flag: str, expr: str | None, regions: list[str]) -> int:
    command = [SAMTOOLS, "view", "-c", "-F", flag]
    if expr:
        command.extend(["-e", expr])
    if alignment.suffix == ".cram":
        command.extend(["-T", str(reference)])
    command.append(str(alignment))
    command.extend(regions)
    return int(subprocess.run(command, check=True, text=True, capture_output=True).stdout.strip())


def make_vcf(root: Path, name: str, records: list[tuple[int, str, str, str]]) -> Path:
    plain = root / f"{name}.vcf"
    with plain.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n")
        handle.write('##ALT=<ID=DEL,Description="Deletion">\n')
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        for pos, ref, alt, gt in records:
            handle.write(f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")
    compressed = root / f"{name}.vcf.gz"
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    return compressed


def legacy_sequence_metrics(vcf: Path) -> dict[str, int]:
    def lines(command):
        return subprocess.run(command, check=True, text=True, capture_output=True).stdout.splitlines()

    indels = len(lines([BCFTOOLS, "view", "-H", "-v", "indels", str(vcf)]))
    hom_refs = sum(value.strip() in {"0|0", "0/0"} for value in lines([BCFTOOLS, "query", "-f", "[%GT]\n", str(vcf)]))
    symbolic = len(lines([BCFTOOLS, "view", "-H", "-i", 'ALT~"<.*>"', str(vcf)]))
    chars = subprocess.run(
        [BCFTOOLS, "query", "-f", "%REF\t%ALT\n", str(vcf)],
        check=True,
        text=True,
        capture_output=True,
    ).stdout
    allowed = set("ACGTN\t\n")
    bad = sum(character not in allowed for character in chars)
    return {"indel_count": indels, "hom_ref_count": hom_refs, "symbolic_alt_count": symbolic, "bad_base_char_count": bad}


def legacy_tag_metrics(vcf: Path) -> dict[str, int]:
    rows = subprocess.run([BCFTOOLS, "query", "-f", "[%GT]\n", str(vcf)], check=True, text=True, capture_output=True).stdout.splitlines()
    counts = {"records": 0, "complete_phased_het_count": 0, "unphased_het_count": 0,
              "missing_gt_count": 0, "other_gt_count": 0, "hom_count": 0}
    for value in rows:
        counts["records"] += 1
        if value in {"0|1", "1|0"}:
            counts["complete_phased_het_count"] += 1
        elif "." in value or not value:
            counts["missing_gt_count"] += 1
        elif value in {"0/1", "1/0"}:
            counts["unphased_het_count"] += 1
        else:
            counts["other_gt_count"] += 1
            if value in {"0|0", "0/0", "1|1", "1/1"}:
                counts["hom_count"] += 1
    counts["indel_count"] = len(subprocess.run([BCFTOOLS, "view", "-H", "-v", "indels", str(vcf)], check=True, text=True, capture_output=True).stdout.splitlines())
    return counts


class StageInternalDedupTest(unittest.TestCase):
    def test_realign_qc_is_single_source_and_input_counts_match_samtools(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            reference = make_reference(root)

            def reads(header):
                values = [
                    make_read(header, "hp1", hp=1, pos=10),
                    make_read(header, "hp2", hp=2, pos=20),
                    make_read(header, "hp0", hp=0, pos=30),
                    make_read(header, "hp3", hp=3, pos=40),
                    make_read(header, "missing", pos=50),
                    make_read(header, "other_chr", tid=1, pos=60, hp=2),
                    make_read(header, "secondary", flag=256, pos=70, hp=1),
                    make_read(header, "supplementary", flag=2048, pos=80, hp=2),
                    make_read(header, "unmapped", flag=4, hp=None),
                    make_read(header, "unmapped_secondary", flag=260),
                    make_read(header, "unmapped_supplementary", flag=2052),
                ]
                return sorted(values, key=lambda read: (read.reference_id < 0, read.reference_id, read.reference_start))

            for suffix in (".bam", ".cram"):
                alignment = root / f"input{suffix}"
                write_bam(alignment, reads, reference=reference if suffix == ".cram" else None)
                for platform, expected_untagged in (("ont", 2), ("illumina", 3)):
                    config = personal_realign.PersonalRealignConfig(
                        cohort="test", output_root=root, sample_list=None, samples=["S"], chroms=["chr1"],
                        sequence_variant_policy="snv_indel", dataset_role="test", platform=platform,
                        reference_fasta=reference, samtools_bin=SAMTOOLS, sort_threads=2,
                    )
                    result = personal_realign.count_input_reads(config, alignment, ["chr1"])
                    self.assertEqual(result["hp1"], legacy_samtools_count(alignment, reference, "0x904", "[HP]==1", ["chr1"]))
                    self.assertEqual(result["hp2"], legacy_samtools_count(alignment, reference, "0x904", "[HP]==2", ["chr1"]))
                    self.assertEqual(result["untagged"], expected_untagged)
                    self.assertEqual(result["all_primary_mapped"], legacy_samtools_count(alignment, reference, "0x904", None, ["chr1"]))
                    self.assertEqual(result["unmapped"], legacy_samtools_count(alignment, reference, "0x900", None, ["*"]))
            self.assertEqual(personal_realign.scan_realigned_output.__module__, "integrate_wgs_rna.realign.qc")

    def test_flagstat_text_matches_samtools_from_same_scan(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            bam = root / "reads.bam"
            rng = random.Random(261006)
            bit_values = (4, 256, 2048, 1024, 512, 1, 2, 64, 128, 8, 16, 32)
            flags = [0, 1024, 99, 147, 73, 133, 256, 2048, 512, 97, 97 | 256, 97 | 2048]
            flags.extend(sum(bit for bit in bit_values if rng.randrange(2)) for _ in range(120))

            def reads(header):
                result = []
                for i, flag in enumerate(flags):
                    unmapped = bool(flag & 4)
                    mate_unmapped = bool(flag & 8)
                    result.append(make_read(
                        header, f"r{i}", flag=flag, tid=i % 2, pos=10 + i,
                        mtid=(-1 if mate_unmapped else (i + 1) % 2), mpos=(20 + i if not mate_unmapped else -1),
                        mapq=rng.choice((0, 4, 5, 20, 60)),
                    ))
                return result

            write_bam(bam, reads, index=False)
            config = personal_realign.PersonalRealignConfig(
                cohort="test", output_root=root, sample_list=None, samples=["S"], chroms=["chr1"],
                sequence_variant_policy="snv_indel", dataset_role="test", output_format="bam", samtools_bin=SAMTOOLS,
            )
            scan = personal_realign.scan_realigned_output(config, bam, root / "unused.fa", 1)
            old = subprocess.run([SAMTOOLS, "flagstat", str(bam)], check=True, text=True, capture_output=True).stdout
            self.assertEqual(scan["flagstat_text"], old)

    def test_markdup_single_pass_counts_and_atomic_link_copy_fallback(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            bam = root / "marked.bam"

            def reads(header):
                return [
                    make_read(header, "primary", pos=10),
                    make_read(header, "duplicate", flag=1024, pos=20),
                    make_read(header, "secondary", flag=256, pos=30),
                    make_read(header, "supplementary", flag=2048, pos=40),
                    make_read(header, "unmapped", flag=4),
                    make_read(header, "qcfail", flag=512, pos=50),
                ]

            write_bam(bam, reads, index=False)
            observed = markdup.primary_record_counts(bam)
            expected_primary = int(subprocess.run([SAMTOOLS, "view", "-c", "-F", "2308", str(bam)], check=True, text=True, capture_output=True).stdout)
            expected_dups = int(subprocess.run([SAMTOOLS, "view", "-c", "-F", "2308", "-f", "1024", str(bam)], check=True, text=True, capture_output=True).stdout)
            self.assertEqual(observed, (expected_dups, expected_primary))

            linked = root / "linked.bam"
            was_linked = markdup.link_or_copy_atomic(bam, linked)
            self.assertEqual(hashlib.sha256(bam.read_bytes()).digest(), hashlib.sha256(linked.read_bytes()).digest())
            if was_linked:
                self.assertEqual(bam.stat().st_ino, linked.stat().st_ino)

            copied = root / "copied.bam"
            with patch.object(markdup.os, "link", side_effect=OSError(errno.EXDEV, "cross-device link")):
                self.assertFalse(markdup.link_or_copy_atomic(bam, copied))
            self.assertEqual(hashlib.sha256(bam.read_bytes()).digest(), hashlib.sha256(copied.read_bytes()).digest())
            self.assertFalse(Path(str(copied) + ".part").exists())

    def test_seq_and_tag_vcf_metrics_match_previous_bcftools_queries(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            seq = make_vcf(root, "seq", [
                (1, "A", "C", "0|1"),
                (2, "AT", "A", "0|0"),
                (3, "A", "AT", "1/1"),
                (4, "A", "<DEL>", "0|1"),
                (5, "A", "C,G", "1|2"),
                (6, "R", "C", "0/1"),
                (7, "AT", "A,ATT", "0|1"),
            ])
            self.assertEqual(summarize_sequence_vcf(BCFTOOLS, seq), legacy_sequence_metrics(seq))
            with self.assertRaisesRegex(RuntimeError, "hom-ref records remain"):
                validate_seq_vcf(
                    PerSampleVcfConfig("test", root, root, ["chr1"], "snv_indel", bcftools_bin=BCFTOOLS),
                    seq,
                    "snv_indel",
                )
            valid_seq = make_vcf(root, "valid_seq", [
                (1, "A", "C", "0|1"),
                (2, "AT", "A", "0|1"),
                (3, "A", "AT", "1|0"),
            ])
            seq_metrics = validate_seq_vcf(
                PerSampleVcfConfig("test", root, root, ["chr1"], "snv_indel", bcftools_bin=BCFTOOLS),
                valid_seq,
                "snv_indel",
            )
            self.assertEqual(seq_metrics["variant_count"], int(subprocess.run([BCFTOOLS, "index", "-n", str(valid_seq)], check=True, text=True, capture_output=True).stdout))

            tag = make_vcf(root, "tag", [
                (1, "A", "C", "0|1"),
                (2, "AT", "A", "1|0"),
                (3, "A", "AT", "0|1"),
                (4, "AT", "A,ATT", "0|1"),
            ])
            observed = summarize_tag_genotypes(BCFTOOLS, tag)
            self.assertEqual(observed, legacy_tag_metrics(tag))
            validated = validate_tag_vcf(
                PerSampleVcfConfig("test", root, root, ["chr1"], "snv_indel", bcftools_bin=BCFTOOLS), tag
            )
            self.assertEqual(validated["indel_count"], 3)
            self.assertEqual(summarize_tag_genotypes(BCFTOOLS, tag, include_indels=False)["records"], 4)


if __name__ == "__main__":
    unittest.main()
