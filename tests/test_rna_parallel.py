"""RNA 최적화 회귀: 독립 depth/bedcov 기준선, BAM/CRAM, 실패와 입력 변조."""
from __future__ import annotations

import importlib.util
import csv
import os
import statistics
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pysam
import zarr

ROOT = Path(os.environ.get("RNA_TEST_ROOT", Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT))
def stage_path(name):
    if (ROOT/'09_personal_fasta.sh').exists():return ROOT/name
    return ROOT/f'{int(name[:2])-2:02d}{name[2:]}'

from rna_resources import worker_count  # noqa: E402


class RnaParallelTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="rna_parallel_test_")
        self.base = Path(self.temp.name)
        self.fasta = self.base / "target.fa"
        self.fasta.write_text("".join(f">chr{i}\n{'A' * 600}\n" for i in range(1, 4)))
        pysam.faidx(str(self.fasta))
        self.fai = Path(str(self.fasta) + ".fai")
        self.bam = self.base / "reads.bam"
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": f"chr{i}", "LN": 600} for i in range(1, 4)]}
        with pysam.AlignmentFile(self.bam, "wb", header=header) as output:
            records = [(0, 0, 255, [(0, 100)]), (0, 1024, 255, [(0, 100)]),
                       (20, 256, 255, [(0, 100)]), (20, 2048, 255, [(0, 100)]),
                       (20, 512, 255, [(0, 100)]), (20, 0, 254, [(0, 100)]),
                       (120, 0, 255, [(0, 10), (2, 10), (0, 10)]),
                       (160, 0, 255, [(0, 10), (3, 20), (0, 10)]),
                       (220, 99, 255, [(0, 100)]), (250, 147, 255, [(0, 100)])]
            for refid in (0, 1):
                for n, (start, flag, mapq, cigar) in enumerate(records):
                    read = pysam.AlignedSegment()
                    read.query_name = "pair" if flag in (99, 147) else f"r{refid}_{n}"
                    qlen = sum(length for op, length in cigar if op in (0, 1, 4, 7, 8))
                    read.query_sequence = "A" * qlen
                    read.query_qualities = pysam.qualitystring_to_array("I" * qlen)
                    read.flag, read.reference_id, read.reference_start = flag, refid, start
                    read.mapping_quality, read.cigartuples = mapq, cigar
                    if flag in (99, 147):
                        read.next_reference_id = refid
                        read.next_reference_start = 250 if flag == 99 else 220
                        read.template_length = 130 if flag == 99 else -130
                    output.write(read)
        pysam.index(str(self.bam))
        self.cram = self.base / "reads.cram"
        with pysam.AlignmentFile(self.bam, "rb") as source:
            with pysam.AlignmentFile(self.cram, "wc", template=source,
                                     reference_filename=str(self.fasta)) as output:
                for read in source:
                    output.write(read)
        pysam.index(str(self.cram))
        self.exons = self.base / "exons.tsv"
        self.exons.write_text(
            "g1\tchr1\t0\t80\ng1\tchr1\t70\t100\ng1\tchr1\t100\t130\n"
            "g2\tchr1\t20\t120\ng3\tchr2\t220\t350\n"
            "cross\tchr1\t0\t100\ncross\tchr2\t0\t100\n"
            "short\tchr2\t160\t180\nzero\tchr3\t0\t110\n")
        self.windows = self.base / "windows.tsv"
        self.windows.write_text(
            "contig\thaplotype\twindow_index\tpersonal_start\tpersonal_end\tkept\n"
            "chr2\tHAP1\t8\t160\t360\t1\nchr1\tHAP1\t2\t0\t150\t1\n"
            "chr3\tHAP1\t12\t0\t600\t1\nchr1\tHAP1\t4\t160\t360\t1\n"
            "chr1\tHAP1\t5\t400\t600\t0\n")

    def tearDown(self):
        self.temp.cleanup()

    def run_script(self, name, args, check=True):
        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join([str(ROOT), str(ROOT / "src"),
            env.get("PYTHONPATH", "")])
        env["MFM_SAMPLE_CPU_BUDGET"] = "14"
        result = subprocess.run([sys.executable, str(stage_path(name)), *map(str, args)],
                                env=env, capture_output=True, text=True)
        if check:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def cnv_args(self, alignment, out, jobs):
        args = []
        for _ in range(2):
            args += ["--exons", self.exons, "--bam", alignment, "--reference", self.fasta]
        return [*args, "--out", out, "--jobs", jobs]

    def test_cnv_serial_parallel_and_original(self):
        outputs = []
        for alignment in (self.bam, self.cram):
            for jobs in (1, 2, 13):
                out = self.base / f"cnv_{alignment.suffix}_{jobs}.tsv"
                result = self.run_script("22_cnv_filter.py", self.cnv_args(alignment, out, jobs))
                self.assertIn("[cnv:progress] 6/6", result.stdout)
                outputs.append(out.read_bytes())
        self.assertTrue(all(data == outputs[0] for data in outputs))
        # 수정 전 파일이 제공되면 byte 단위 기준선도 대조한다(배포 전 필수 실행).
        baseline = os.environ.get("RNA_TEST_BASELINE")
        if baseline:
            original = Path(baseline) / "22_cnv_filter.py"
            out = self.base / "original.tsv"
            args = self.cnv_args(self.cram, out, 1)[:-2]
            subprocess.run([sys.executable, str(original), *map(str, args)], check=True,
                           capture_output=True, text=True)
            # 승인된 결측 정책 변경만 허용한다. rd와 관측 가능한 유전자 판정은 불변.
            with out.open() as handle:
                old = {r["gene_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
            new = {r["gene_id"]: r for r in csv.DictReader(outputs[0].decode().splitlines(), delimiter="\t")}
            for gene, row in old.items():
                for key, value in row.items():
                    if key == "pass_cnv" and row["rd"] == "NA":
                        self.assertEqual(new[gene][key], "0")
                    else:
                        self.assertEqual(new[gene][key], value)
        text = outputs[0].decode()
        self.assertIn("short\t20\tNA\tNA\tNA\t0", text)
        self.assertIn("zero\t110\t0.0000\t0.0000\t0.0000\t0", text)

    def test_cnv_independent_bedcov_reference(self):
        spec = importlib.util.spec_from_file_location("cnv_stage", stage_path("22_cnv_filter.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        genes = module.load_exons(self.exons)
        rows, lengths = [], {}
        for gene, chromosomes in genes.items():
            for contig, intervals in chromosomes.items():
                for start, end in module.merge(intervals):
                    rows.append((contig, start, end, gene))
                    lengths[gene] = lengths.get(gene, 0) + end - start
        rows.sort(key=lambda r: (r[0], r[1], r[2]))
        bed = self.base / "baseline.bed"
        bed.write_text("".join(f"{c}\t{s}\t{e}\t{g}\n" for c, s, e, g in rows))
        proc = subprocess.run(["samtools", "bedcov", str(bed), str(self.bam)],
                              check=True, capture_output=True, text=True)
        sums = {}
        for line in proc.stdout.splitlines():
            fields = line.split("\t")
            sums[fields[3]] = sums.get(fields[3], 0.0) + float(fields[4])
        means = {gene: sums[gene] / length for gene, length in lengths.items() if length >= 100}
        median = statistics.median(means.values())
        out = self.base / "ref.tsv"
        self.run_script("22_cnv_filter.py", self.cnv_args(self.bam, out, 2))
        for line in out.read_text().splitlines()[1:]:
            gene, _, rd1, _, _, passed = line.split("\t")
            if gene in means:
                rd = means[gene] / median
                self.assertEqual(rd1, f"{rd:.4f}")
                self.assertEqual(passed, str(int(0.5 <= rd <= 1.5)))

    def test_cnv_failure_does_not_publish(self):
        out = self.base / "failed.tsv"
        result = self.run_script("22_cnv_filter.py",
            [*self.cnv_args(self.bam, out, 2), "--samtools", "/bin/false"], check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(out.exists())
        self.exons.write_text("zero\tchr3\t0\t110\n")
        result = self.run_script("22_cnv_filter.py", self.cnv_args(self.bam, out, 1), check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("중앙값이 0", result.stderr)
        self.assertFalse(out.exists())

    def genome(self, jobs):
        out = self.base / f"genome_{jobs}.zarr"
        self.run_script("18_coverage_zarr.py", ["--bam", self.bam, "--fai", self.fai,
            "--out", out, "--contigs", "chr1 chr2 chr3", "--sample", "S1", "--hap", "HAP1",
            "--jobs", jobs, "--chunk", 64])
        return out

    def test_genome_serial_parallel_and_raw_filter(self):
        first, second = self.genome(1), self.genome(2)
        a, b = (zarr.open_group(str(p), mode="r") for p in (first, second))
        self.assertEqual(dict(a.attrs), dict(b.attrs))
        proc = subprocess.run(["samtools", "depth", "-a", "-J", "-Q", "255", "-G", "2820", "-g", "1024",
                               str(self.bam)], check=True, capture_output=True, text=True)
        expected = {f"chr{i}": np.zeros(600, dtype=np.uint32) for i in range(1, 4)}
        for line in proc.stdout.splitlines():
            c, pos, depth = line.split("\t")
            expected[c][int(pos) - 1] = int(depth)
        for c in expected:
            np.testing.assert_array_equal(a[c][:], expected[c])
            np.testing.assert_array_equal(b[c][:], expected[c])
        # GTEx target과 같이 -g 1024로 DUP를 포함한다.
        self.assertEqual(a["chr1"][10], 2)
        self.assertEqual(a["chr1"][135], 1)  # D 포함.
        self.assertEqual(a["chr1"][175], 0)  # N 제외.
        self.assertEqual(a["chr1"][270], 2)  # mate overlap 두 번.
        self.assertEqual(a["chr3"][:].sum(), 0)
        baseline = os.environ.get("RNA_TEST_BASELINE")
        if baseline:
            old = self.base / "old.zarr"
            subprocess.run([sys.executable, str(Path(baseline) / "18_coverage_zarr.py"),
                "--bam", str(self.bam), "--fai", str(self.fai), "--out", str(old),
                "--contigs", "chr1 chr2 chr3", "--sample", "S1", "--hap", "HAP1", "--chunk", "64"],
                env={**os.environ, "PYTHONPATH": os.pathsep.join([str(ROOT), os.environ.get("PYTHONPATH", "")])},
                check=True, capture_output=True, text=True)
            original = zarr.open_group(str(old), mode="r")
            for c in expected:
                np.testing.assert_array_equal(original[c][:], a[c][:])
            for key, value in original.attrs.items():
                if key not in {"read_filter", "coverage_policy", "duplicate_policy",
                               "effective_duplicate_policy", "legacy_filter_metadata_mismatch"}:
                    self.assertEqual(a.attrs[key], value, key)
            self.assertEqual(a.attrs["duplicate_policy"], "included")
            self.assertFalse(a.attrs["legacy_filter_metadata_mismatch"])

    def test_coverage_change_does_not_invalidate_star_cache(self):
        import rna_profile
        gtf = self.base / "g.gtf"
        gtf.write_text('chr1\tx\texon\t1\t10\t.\t+\t.\tgene_id "G";\n')
        full = rna_profile.fingerprint(str(gtf), str(self.fai))["sha256"]
        align = rna_profile.fingerprint(str(gtf), str(self.fai), align_only=True)["sha256"]
        old = list(rna_profile.COVERAGE_DEPTH_ARGS)
        rna_profile.COVERAGE_DEPTH_ARGS[:] = old[:-2]
        try:
            self.assertNotEqual(full, rna_profile.fingerprint(str(gtf), str(self.fai))["sha256"])
            self.assertEqual(align, rna_profile.fingerprint(str(gtf), str(self.fai), align_only=True)["sha256"])
        finally:
            rna_profile.COVERAGE_DEPTH_ARGS[:] = old

    def test_genome_failure_has_no_completion_metadata(self):
        out = self.base / "bad_genome.zarr"
        result = self.run_script("18_coverage_zarr.py", ["--bam", self.bam, "--fai", self.fai,
            "--out", out, "--contigs", "chr1 chr2", "--jobs", "2", "--samtools", "/bin/false"], check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("total_depth", zarr.open_group(str(out), mode="r").attrs)

    def test_worker_budget(self):
        with patch.dict(os.environ, {"MFM_SAMPLE_CPU_BUDGET": "14", "MFM_CPU_BUDGET": "72"}), \
                patch("rna_resources.os.sched_getaffinity", return_value=set(range(72))):
            self.assertEqual(worker_count("auto", 44), 13)
            self.assertEqual(worker_count("auto", 22, 2), 6)
            self.assertEqual(worker_count("100", 44), 13)
            self.assertEqual(worker_count("2", 44), 2)
            with self.assertRaises(ValueError):
                worker_count("0", 44)


if __name__ == "__main__":
    unittest.main()
