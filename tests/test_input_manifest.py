from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from integrate_wgs_rna.input_manifest import FILE_COLUMNS, InputManifestError, validate_wgs_input_manifest
from integrate_wgs_rna.library_provenance import PROVENANCE_COLUMNS, project_manifest_fields


SAMPLE_COLUMNS = [
    "sample_id", "platform", "read_type", "dataset_role", "project_id", "cohort", "alignment_path",
    "read_files_manifest", "prepared_alignment_path", "source_decode_reference_fasta",
    "analysis_reference_fasta", "vcf_source", "public_vcf_path", "public_vcf_sample_id",
    "vcf_reference_fasta", "vcf_callset", "phase_source",
    "rna_path", "data_group",
    *PROVENANCE_COLUMNS,
]


class InputManifestTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        for name in ("reads.bam", "r1.fq.gz", "r2.fq.gz", "sample.cram", "ref.fa", "sample.vcf.gz"):
            (self.root / name).touch()

    def tearDown(self) -> None:
        self.temp.cleanup()

    def write_tsv(self, name: str, columns: list[str] | tuple[str, ...], rows: list[dict[str, str]]) -> Path:
        path = self.root / name
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        return path

    def sample_row(self, **changes: str) -> dict[str, str]:
        row = {column: "-" for column in SAMPLE_COLUMNS}
        row.update({
            "sample_id": "S1", "platform": "hifi", "read_type": "long_read",
            "dataset_role": "pre_training", "project_id": "test", "cohort": "test",
            "alignment_path": "reads.bam", "vcf_source": "clair3",
        })
        row.update(changes)
        return row

    def file_row(self, **changes: str) -> dict[str, str]:
        row = {column: "-" for column in FILE_COLUMNS}
        row.update({
            "sample_id": "S1", "run_id": "run1", "file_role": "alignment",
            "format": "ubam", "alignment_state": "unaligned", "path": "reads.bam",
        })
        row.update(changes)
        return row

    def test_legacy_alignment_is_backward_compatible(self) -> None:
        manifest = self.write_tsv("samples.tsv", SAMPLE_COLUMNS, [self.sample_row()])
        contract = validate_wgs_input_manifest(manifest, "S1")
        self.assertEqual(contract["input_mode"], "legacy_alignment")
        self.assertIsNone(contract["prepared_alignment"])
        self.assertEqual(contract["library_provenance"]["wgs"]["amplification"], "unspecified")

    def test_known_project_provenance_is_required_and_returned(self) -> None:
        missing = self.sample_row(
            project_id="1KGP_ILLUMINA", platform="illumina", read_type="short_read",
            vcf_source="illumina_panel",
        )
        manifest = self.write_tsv("missing.tsv", SAMPLE_COLUMNS, [missing])
        with self.assertRaisesRegex(InputManifestError, "provenance 불일치"):
            validate_wgs_input_manifest(manifest, "S1")
        valid = dict(missing)
        valid.update(project_manifest_fields("1KGP_ILLUMINA", False))
        manifest = self.write_tsv("valid.tsv", SAMPLE_COLUMNS, [valid])
        contract = validate_wgs_input_manifest(manifest, "S1")
        self.assertEqual(contract["schema"], "mfm_wgs_input_contract_v2")
        self.assertEqual(contract["library_provenance"]["wgs"]["amplification"],
                         "pcr_free_confirmed")

    def test_known_rna_pair_requires_confirmed_ten_cycle_provenance(self) -> None:
        sample = self.sample_row(
            project_id="HPRC_R2", rna_path="rna.bam", data_group="wgs_rna_pair",
        )
        sample.update(project_manifest_fields("HPRC_R2", True))
        manifest = self.write_tsv("pair.tsv", SAMPLE_COLUMNS, [sample])
        contract = validate_wgs_input_manifest(manifest, "S1")
        self.assertEqual(contract["library_provenance"]["wgs"]["amplification"],
                         "pcr_free_inferred")
        self.assertEqual(contract["library_provenance"]["rna"]["pcr_cycles"], 10)
        sample["rna_library_pcr_cycles"] = "8"
        manifest = self.write_tsv("pair_bad.tsv", SAMPLE_COLUMNS, [sample])
        with self.assertRaisesRegex(InputManifestError, "provenance 불일치"):
            validate_wgs_input_manifest(manifest, "S1")

    def test_legacy_cram_requires_decode_reference(self) -> None:
        manifest = self.write_tsv(
            "samples.tsv", SAMPLE_COLUMNS,
            [self.sample_row(platform="illumina", read_type="short_read",
                             alignment_path="sample.cram", vcf_source="illumina_panel")],
        )
        with self.assertRaisesRegex(InputManifestError, "source_decode_reference_fasta"):
            validate_wgs_input_manifest(manifest, "S1")

    def test_legacy_cram_with_decode_reference_is_valid(self) -> None:
        manifest = self.write_tsv(
            "samples.tsv", SAMPLE_COLUMNS,
            [self.sample_row(platform="illumina", read_type="short_read",
                             alignment_path="sample.cram", source_decode_reference_fasta="ref.fa",
                             analysis_reference_fasta="ref.fa", vcf_source="illumina_panel")],
        )
        self.assertEqual(validate_wgs_input_manifest(manifest, "S1")["input_mode"], "legacy_alignment")

    def test_long_ubam_with_public_vcf_is_valid(self) -> None:
        files = self.write_tsv("files.tsv", FILE_COLUMNS, [self.file_row()])
        sample = self.sample_row(
            alignment_path="-", read_files_manifest=files.name, analysis_reference_fasta="ref.fa",
            vcf_source="public", public_vcf_path="sample.vcf.gz", vcf_reference_fasta="ref.fa",
            public_vcf_sample_id="S1", vcf_callset="test_callset", phase_source="provided",
        )
        manifest = self.write_tsv("samples.tsv", SAMPLE_COLUMNS, [sample])
        contract = validate_wgs_input_manifest(manifest, "S1")
        self.assertEqual(contract["source"]["formats"], ["ubam"])
        self.assertIsNone(contract["prepared_alignment"])

    def test_clr_fastq_with_public_vcf_is_valid(self) -> None:
        files = self.write_tsv("clr_files.tsv", FILE_COLUMNS, [self.file_row(
            file_role="long_reads", format="fastq_gz", path="r1.fq.gz",
        )])
        sample = self.sample_row(
            platform="pacbio_clr", read_type="long_read", alignment_path="-",
            read_files_manifest=files.name, analysis_reference_fasta="ref.fa",
            vcf_source="public", public_vcf_path="sample.vcf.gz",
            public_vcf_sample_id="S1", vcf_reference_fasta="ref.fa",
            vcf_callset="test_callset", phase_source="provided",
        )
        manifest = self.write_tsv("samples.tsv", SAMPLE_COLUMNS, [sample])
        contract = validate_wgs_input_manifest(manifest, "S1")
        self.assertEqual(contract["platform"], "clr")
        self.assertEqual(contract["source"]["formats"], ["fastq_gz"])

    def test_long_read_run_must_have_one_raw_file(self) -> None:
        rows = [
            self.file_row(path="reads.bam"),
            self.file_row(path="r1.fq.gz", file_role="long_reads", format="fastq_gz"),
        ]
        files = self.write_tsv("files.tsv", FILE_COLUMNS, rows)
        sample = self.sample_row(
            alignment_path="-", read_files_manifest=files.name,
            analysis_reference_fasta="ref.fa",
        )
        manifest = self.write_tsv("samples.tsv", SAMPLE_COLUMNS, [sample])
        with self.assertRaisesRegex(InputManifestError, "파일 하나"):
            validate_wgs_input_manifest(manifest, "S1")

    def test_paired_fastq_missing_mate_fails(self) -> None:
        files = self.write_tsv("files.tsv", FILE_COLUMNS, [self.file_row(
            file_role="read1", format="fastq_gz", path="r1.fq.gz",
        )])
        sample = self.sample_row(
            platform="illumina", read_type="short_read", alignment_path="-",
            read_files_manifest=files.name, analysis_reference_fasta="ref.fa", vcf_source="illumina_panel",
        )
        manifest = self.write_tsv("samples.tsv", SAMPLE_COLUMNS, [sample])
        with self.assertRaisesRegex(InputManifestError, "read1/read2"):
            validate_wgs_input_manifest(manifest, "S1")

    def test_paired_fastq_is_valid(self) -> None:
        rows = [
            self.file_row(file_role="read1", format="fastq_gz", path="r1.fq.gz"),
            self.file_row(file_role="read2", format="fastq_gz", path="r2.fq.gz"),
        ]
        files = self.write_tsv("files.tsv", FILE_COLUMNS, rows)
        sample = self.sample_row(
            platform="illumina", read_type="short_read", alignment_path="-",
            read_files_manifest=files.name, analysis_reference_fasta="ref.fa", vcf_source="illumina_panel",
        )
        manifest = self.write_tsv("samples.tsv", SAMPLE_COLUMNS, [sample])
        contract = validate_wgs_input_manifest(manifest, "S1")
        self.assertEqual(contract["source"]["file_count"], 2)

    def test_illumina_ubam_is_valid(self) -> None:
        files = self.write_tsv("files.tsv", FILE_COLUMNS, [self.file_row()])
        sample = self.sample_row(
            platform="illumina", read_type="short_read", alignment_path="-",
            read_files_manifest=files.name, analysis_reference_fasta="ref.fa", vcf_source="illumina_panel",
        )
        manifest = self.write_tsv("samples.tsv", SAMPLE_COLUMNS, [sample])
        contract = validate_wgs_input_manifest(manifest, "S1")
        self.assertEqual(contract["source"]["formats"], ["ubam"])

    def test_prepared_cram_is_rejected(self) -> None:
        manifest = self.write_tsv(
            "samples.tsv", SAMPLE_COLUMNS,
            [self.sample_row(prepared_alignment_path="sample.cram")],
        )
        with self.assertRaisesRegex(InputManifestError, "canonical BAM"):
            validate_wgs_input_manifest(manifest, "S1")

    def test_single_aligned_cram_is_reused(self) -> None:

        files = self.write_tsv("files.tsv", FILE_COLUMNS, [self.file_row(
            file_role="alignment", format="cram", alignment_state="aligned",
            path="sample.cram", source_reference_fasta="ref.fa",
        )])
        sample = self.sample_row(
            platform="illumina", read_type="short_read", alignment_path="-",
            read_files_manifest=files.name, analysis_reference_fasta="ref.fa",
            vcf_source="illumina_panel",
        )
        manifest = self.write_tsv("samples.tsv", SAMPLE_COLUMNS, [sample])
        contract = validate_wgs_input_manifest(manifest, "S1")
        self.assertIsNone(contract["prepared_alignment"])

    def test_legacy_and_file_manifest_are_mutually_exclusive(self) -> None:
        files = self.write_tsv("files.tsv", FILE_COLUMNS, [self.file_row()])
        manifest = self.write_tsv(
            "samples.tsv", SAMPLE_COLUMNS, [self.sample_row(read_files_manifest=files.name)]
        )
        with self.assertRaisesRegex(InputManifestError, "정확히 하나"):
            validate_wgs_input_manifest(manifest, "S1")


if __name__ == "__main__":
    unittest.main()
