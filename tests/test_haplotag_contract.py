"""파일럿 복구의 donor/전체 reference/read 보존 회귀 시험."""
import csv
import os
import shlex
import shutil
import tempfile
import subprocess
import unittest
from pathlib import Path

import pysam

from integrate_wgs_rna.haplotag import HaplotagConfig, run_whatshap
from integrate_wgs_rna.haplotag_contract import validate_haplotag_contract


class HaplotagContractTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.ref = self.root / "ref.fa"
        self.ref.write_text(">chr1\n" + "A" * 100 + "\n>chrX\n" + "A" * 100 + "\n")
        pysam.faidx(str(self.ref))
        self.bam = self.root / "input.bam"
        self.write_bam(["HG001"])
        self.vcf = self.root / "tag.vcf.gz"
        self.write_vcf(["NA12878"])

    def tearDown(self):
        self.tmp.cleanup()

    def write_bam(self, samples):
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": "chr1", "LN": 100}, {"SN": "chrX", "LN": 100}],
                  "RG": [{"ID": f"rg{i}", "SM": s} for i, s in enumerate(samples)]}
        with pysam.AlignmentFile(str(self.bam), "wb", header=header) as out:
            for name, tid in [("chr1_read", 0), ("chrX_read", 1), ("unmapped_read", -1)]:
                read = pysam.AlignedSegment(out.header)
                read.query_name = name
                read.query_sequence = "AAAAAGAAAA"
                read.query_qualities = pysam.qualitystring_to_array("IIIIIIIIII")
                read.reference_id = tid
                read.reference_start = 0 if tid >= 0 else -1
                read.mapping_quality = 60 if tid >= 0 else 0
                read.flag = 0 if tid >= 0 else 4
                if tid >= 0:
                    read.cigarstring = "10M"
                read.set_tag("RG", "rg0")
                out.write(read)
        pysam.index(str(self.bam))

    def write_vcf(self, samples):
        header = pysam.VariantHeader()
        header.contigs.add("chr1", length=100)
        header.contigs.add("chrX", length=100)
        header.formats.add("GT", 1, "String", "Genotype")
        header.formats.add("PS", 1, "Integer", "Phase set")
        for sample in samples:
            header.add_sample(sample)
        with pysam.VariantFile(str(self.vcf), "wz", header=header) as out:
            rec = out.new_record(contig="chr1", start=5, stop=6, alleles=("A", "G"))
            for sample in samples:
                rec.samples[sample]["GT"] = (0, 1)
                rec.samples[sample].phased = True
                rec.samples[sample]["PS"] = 1
            out.write(rec)
        pysam.tabix_index(str(self.vcf), preset="vcf", force=True)

    def contract(self, alias="HG001"):
        return validate_haplotag_contract(self.bam, self.vcf, self.ref, "NA12878", alias)

    def test_explicit_alias_is_allowed(self):
        self.assertTrue(self.contract()["ignore_read_groups"])

    def test_canonical_sample_keeps_read_group_check(self):
        self.write_bam(["NA12878"])
        self.assertFalse(self.contract()["ignore_read_groups"])

    def test_alias_not_in_manifest_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "단일 donor"):
            self.contract(alias=None)

    def test_multiple_bam_donors_are_rejected(self):
        self.write_bam(["HG001", "NA12878"])
        with self.assertRaisesRegex(ValueError, "단일 donor"):
            self.contract()

    def test_missing_rg_sm_is_rejected(self):
        self.write_bam([""])
        with self.assertRaisesRegex(ValueError, "SM"):
            self.contract()

    def test_multiple_vcf_donors_are_rejected(self):
        self.write_vcf(["NA12878", "OTHER"])
        with self.assertRaisesRegex(ValueError, "단일 donor VCF"):
            self.contract()

    def test_autosome_only_reference_is_rejected(self):
        ref = self.root / "only.fa"
        ref.write_text(">chr1\n" + "A" * 100 + "\n")
        pysam.faidx(str(ref))
        with self.assertRaisesRegex(ValueError, "missing=.*chrX"):
            validate_haplotag_contract(self.bam, self.vcf, ref, "NA12878", "HG001")

    def test_length_mismatch_is_rejected(self):
        ref = self.root / "short.fa"
        ref.write_text(">chr1\n" + "A" * 99 + "\n>chrX\n" + "A" * 100 + "\n")
        pysam.faidx(str(ref))
        with self.assertRaisesRegex(ValueError, "length_mismatch=.*chr1"):
            validate_haplotag_contract(self.bam, self.vcf, ref, "NA12878", "HG001")

    def test_real_whatshap_preserves_all_reads_and_original_sm(self):
        config = HaplotagConfig(cohort="test", raw_root=self.root, tag_vcf_dir=self.root,
                                output_root=self.root, reference_fasta=self.ref,
                                platform="hifi", source_sample_id="HG001",
                                allow_longread_vcf_haplotag=True, output_threads=1,
                                tag_supplementary=True)
        output = self.root / "tagged.bam"
        run_whatshap(config, "NA12878", self.bam, self.vcf, self.ref,
                    output, self.root / "taglist.tsv.gz")
        with pysam.AlignmentFile(str(self.bam)) as src, pysam.AlignmentFile(str(output)) as dst:
            before = [(r.query_name, r.query_sequence, list(r.query_qualities), r.flag) for r in src]
            after_reads = list(dst)
            after = [(r.query_name, r.query_sequence, list(r.query_qualities), r.flag) for r in after_reads]
            self.assertEqual(before, after)
            self.assertEqual(dst.header.to_dict()["RG"][0]["SM"], "HG001")
            self.assertEqual(after_reads[0].get_tag("HP"), 2)

    def test_partial_backup_includes_taglist_qc_and_part(self):
        script = Path(__file__).resolve().parents[1] / "06_haplotag.sh"
        text = script.read_text()
        function = text[text.index("HAPLOTAG_APPEND_OVERWRITE=0"):text.index("\nif is_long_read_sample")]
        hap = self.root / "haplotagged"
        hap.mkdir()
        names = ["S.public.bam.part", "S.taglist.tsv.gz", "S.haplotag_qc.json",
                 "S.haplotag.summary.tsv", "S.haplotag.metadata.json"]
        for name in names:
            (hap / name).write_text("보존할 이전 결과")
        env_setup = f'S=S; HAP_DIR="{hap}"; BAM="$HAP_DIR/S.public.bam"; SUMMARY_TSV="$HAP_DIR/S.haplotag.summary.tsv"; METADATA_JSON="$HAP_DIR/S.haplotag.metadata.json"; HAPLOTAG_REBUILD_PARTIAL=1; have_index() {{ return 1; }}; log() {{ :; }}; die() {{ exit 1; }};\n'
        subprocess.run(["bash", "-c", "set -euo pipefail\n" + env_setup + function + "\nprepare_haplotag_output"], check=True)
        backup = list((hap / "partial_backup").iterdir())
        self.assertEqual(len(backup), 1)
        for name in names:
            self.assertFalse((hap / name).exists())
            self.assertEqual((backup[0] / name).read_text(), "보존할 이전 결과")



class HaplotagEntrypointTest(unittest.TestCase):
    """실제 Stage06과 공용 config의 입력 경로 처리·누락 중단을 확인한다."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="mfm-stage06-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.pipeline = Path(__file__).resolve().parents[1]
        self.sample = "ILM_FIXTURE"
        self.bam = self.root / "prepared input" / "ILM_FIXTURE.analysis.bam"
        self.bam.parent.mkdir()
        self.bam.write_text("표준 BAM 존재 검사 fixture")
        self.vcf = self.root / "persample" / "ILM_FIXTURE.tag.vcf.gz"
        self.vcf.parent.mkdir()
        self.vcf.write_text("VCF 존재 검사 fixture")
        self.capture = self.root / "cli.args"
        self.runner = self.root / "capture_cli.sh"
        self.runner.write_text('#!/usr/bin/env bash\nset -euo pipefail\nprintf "%s\\0" "$@" > "$CLI_CAPTURE"\n')
        self.stage = self.root / "06_haplotag.sh"
        shutil.copy2(self.pipeline / self.stage.name, self.stage)
        # 경로·manifest·파일검사는 실제 config를 사용하고 무거운 CLI 실행만 대체한다.
        (self.root / "config_v2.sh").write_text(
            f"source {shlex.quote(str(self.pipeline / 'config_v2.sh'))}\n"
            'run_cli() { bash "$CAPTURE_RUNNER" "$@"; }\n'
        )
        self.manifest = self.root / "samples.tsv"
        self.env = dict(os.environ)
        self.env.pop("RAW_ROOT", None)
        self.env.update({"PIPELINE_MANIFEST": str(self.manifest),
                         "INPUT_ROOT": str(self.root), "OUTPUT_ROOT": str(self.root / "output"),
                         "CLI_CAPTURE": str(self.capture), "CAPTURE_RUNNER": str(self.runner)})

    def run_stage(self, raw_root="", env_raw_root=None):
        # 앞 14개 열은 공용 manifest 계약의 순서를 유지한다.
        columns = ["sample_id", "platform", "read_type", "dataset_role", "project_id", "cohort",
                   "alignment_path", "raw_root", "rna_path", "final_root", "work_root",
                   "common_vcf_dir", "persample_vcf_dir", "rna_bigwigs",
                   "prepared_alignment_path", "source_decode_reference_fasta", "vcf_source"]
        row = dict.fromkeys(columns, "")
        row.update(sample_id=self.sample, platform="illumina", read_type="short_read",
                   dataset_role="pre_training", project_id="fixture_project", cohort="fixture",
                   raw_root=raw_root, final_root=str(self.root / "output"),
                   work_root=str(self.root / "work"), persample_vcf_dir=str(self.vcf.parent),
                   prepared_alignment_path=str(self.bam),
                   source_decode_reference_fasta=str(self.root / "reference.fa"),
                   vcf_source="illumina_panel")
        with self.manifest.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(row)
        env = dict(self.env)
        if env_raw_root is not None:
            env["RAW_ROOT"] = str(env_raw_root)
        return subprocess.run(["bash", str(self.stage), self.sample], env=env,
                              text=True, capture_output=True, timeout=30)

    def assert_cli_paths(self, result, raw_root):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        args = self.capture.read_bytes().decode().split("\0")[:-1]
        self.assertEqual(args[0], "create-haplotagged-alignments")
        self.assertEqual(args[args.index("--raw-root") + 1], str(raw_root))
        self.assertEqual(args[args.index("--alignment") + 1], str(self.bam))
        self.assertEqual(args[args.index("--sample") + 1], self.sample)

    def test_missing_raw_root_uses_prepared_bam_parent_with_spaces(self):
        self.assert_cli_paths(self.run_stage(), self.bam.parent)

    def test_manifest_raw_root_overrides_environment(self):
        explicit = self.root / "manifest raw"
        self.assert_cli_paths(self.run_stage(raw_root="manifest raw",
                                            env_raw_root=self.root / "environment raw"), explicit)

    def test_environment_raw_root_is_preserved(self):
        explicit = self.root / "environment raw"
        self.assert_cli_paths(self.run_stage(env_raw_root=explicit), explicit)

    def test_missing_bam_stops_before_cli(self):
        self.bam.unlink()
        result = self.run_stage()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(f"missing file: {self.bam}", result.stderr)
        self.assertFalse(self.capture.exists())

    def test_missing_vcf_stops_before_cli(self):
        self.vcf.unlink()
        result = self.run_stage()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(f"missing file: {self.vcf}", result.stderr)
        self.assertFalse(self.capture.exists())


if __name__ == "__main__":
    unittest.main()
