"""재개 산출물의 최종 검증·정리 전후 완료 집계 계약을 확인한다."""
from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import zarr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
from tools.compact_manifest import completed_contigs
from scripts.validate_sample_package import Audit, validate_compact
from handoff.cleanup import cleanup
from handoff.final_check import validate as validate_delivery


class CompactManifestTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name) / "S1"
        self.layout = self.root / "compact" / "layout"
        self.directory = self.layout / "chr1"
        self.directory.mkdir(parents=True)
        self.metadata = {
            "sample_id": "S1", "contig": "chr1", "validation_level": "full",
            "validation": {"level": "full", "by_haplotype": {
                hap: {"level": "full", "status": "pass"} for hap in ("HAP1", "HAP2")}},
            "input_provenance": {"contig": "chr1", "validation_level": "full",
                                 "fingerprint_sha256": "a" * 64},
        }
        self.manifest = {"status": "complete", "validation_level": "full", "records": {
            "S1/chr1": {"sample": "S1", "contig": "chr1", "status": "skip_existing",
                         "zarr": str(self.directory)},
        }}
        self.write_metadata(self.metadata)
        (self.layout / "manifest.json").write_text(json.dumps(self.manifest))

    def tearDown(self):
        self.temp.cleanup()

    def write_metadata(self, metadata):
        (self.directory / "contig_metadata.json").write_text(json.dumps(metadata))

    def make_package(self):
        windows = self.root / "windows" / "windows.tsv"
        windows.parent.mkdir()
        windows.write_text(
            "contig\thaplotype\twindow_index\traw_start\tactual_length\tpersonal_start\tpersonal_end\tref_start\tref_end\tkept\n"
            "chr1\tHAP1\t0\t0\t5\t0\t5\t0\t5\t1\n"
            "chr1\tHAP2\t0\t0\t5\t0\t5\t0\t5\t1\n")
        (self.directory / "windows.tsv").write_text(
            "contig\thaplotype\twindow_index\traw_start\tactual_length\n"
            "chr1\tHAP1\t0\t0\t5\nchr1\tHAP2\t0\t0\t5\n")
        (self.directory / "filter_qc.json").write_text("{}\n")
        for name in ("hap_sequence.zarr", "depth.zarr", "readset.HAP1.zarr", "readset.HAP2.zarr"):
            group = zarr.open_group(str(self.directory / name), mode="w")
            group.create_dataset("values", data=np.arange(5, dtype="u4"))
        trace = self.root / "qc" / "stage12_trace" / "all_to_both" / "trace_manifest.json"
        trace.parent.mkdir(parents=True)
        trace.write_text(json.dumps({"status": "pass", "contigs": {"chr1": "pass"},
                                     "compact_root": str(self.layout)}))

    def test_mixed_new_and_valid_reused_contigs(self):
        self.manifest["records"]["S1/chr2"] = {"contig": "chr2", "status": "pass"}
        self.assertEqual(completed_contigs(self.manifest, self.layout, "S1"), {"chr1", "chr2"})
        self.manifest["records"]["S1/chr2"]["status"] = "failed"
        self.assertEqual(completed_contigs(self.manifest, self.layout, "S1"), {"chr1"})

    def test_reused_metadata_must_prove_full_validation_and_identity(self):
        cases = []
        for key, value in (("sample_id", "other"), ("contig", "chr2"),
                           ("validation_level", "fast"), ("input_provenance", {})):
            metadata = copy.deepcopy(self.metadata); metadata[key] = value; cases.append(metadata)
        for key, value in (("status", "fail"), ("level", "standard")):
            metadata = copy.deepcopy(self.metadata)
            metadata["validation"]["by_haplotype"]["HAP2"][key] = value; cases.append(metadata)
        for metadata in cases:
            with self.subTest(metadata=metadata):
                self.write_metadata(metadata)
                with self.assertRaises(ValueError):completed_contigs(self.manifest, self.layout, "S1")

    def test_reused_record_must_point_to_current_sample(self):
        for key, value in (("sample", "other"), ("zarr", str(self.directory / "other"))):
            manifest = copy.deepcopy(self.manifest); manifest["records"]["S1/chr1"][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                completed_contigs(manifest, self.layout, "S1")
        (self.directory / "contig_metadata.json").unlink()
        with self.assertRaises(FileNotFoundError):completed_contigs(self.manifest, self.layout, "S1")

    def test_final_compact_validation_accepts_resume_but_rejects_missing_files(self):
        self.make_package()
        args = argparse.Namespace(root=self.root, layout="layout", sample="S1")
        audit = Audit();validate_compact(args, audit, {"chr1"});self.assertTrue(audit.passed, audit.checks)
        (self.directory / "windows.tsv").unlink()
        audit = Audit();validate_compact(args, audit, {"chr1"});self.assertFalse(audit.passed)

    def test_cleanup_dry_run_and_post_cleanup_include_reused_contigs(self):
        self.make_package()
        with patch("handoff.cleanup.validate_predelete", return_value={"status": "pass"}):
            result = cleanup(self.root, "S1", "wgs", dry_run=True)
        self.assertEqual(result["status"], "planned")
        self.assertTrue((self.root / "windows/windows.tsv").exists())
        # 실제 삭제 대신 보존 경로만 검사하도록 임시 fixture의 원본 window 폴더를 제거한다.
        (self.root / "windows/windows.tsv").unlink();(self.root / "windows").rmdir()
        self.assertEqual(validate_delivery(self.root, "S1", "wgs")["status"], "pass")
        (self.directory / "filter_qc.json").unlink()
        with self.assertRaisesRegex(ValueError, "필수 모델 입력 누락"):
            validate_delivery(self.root, "S1", "wgs")

    def test_catalog_includes_verified_reused_contigs_and_rejects_invalid_metadata(self):
        name = "27_make_catalog.py" if (ROOT / "09_personal_fasta.sh").is_file() else "25_make_catalog.py"
        spec = importlib.util.spec_from_file_location("resume_catalog", ROOT / name)
        catalog = importlib.util.module_from_spec(spec);spec.loader.exec_module(catalog)
        self.make_package()
        for contig in catalog.CHROMS[1:]:
            self.manifest["records"]["S1/" + contig] = {"contig": contig, "status": "pass"}
        (self.layout / "manifest.json").write_text(json.dumps(self.manifest))
        (self.root / "metadata").mkdir()
        (self.root / "metadata/sample.json").write_text(json.dumps({
            "sequence_variant_policy": "snv_indel", "dataset_role": "pre_training", "layout_name": "layout"}))
        (self.root / "qc/final_package_validation.json").write_text(json.dumps({"status": "pass", "failed_checks": []}))
        (self.root / "qc/stage12_trace/all_to_both/trace_manifest.json").write_text(json.dumps({
            "status": "pass", "compact_root": str(self.layout),
            "contigs": {contig: "pass" for contig in catalog.CHROMS}}))
        row = catalog.scan_sample(self.root, "test", {}, {})
        self.assertEqual(row["compact_contigs"], 22)
        self.assertEqual(row["processing_status"], "complete")
        self.assertEqual(row["has_compact"], "True")
        self.metadata["validation"]["by_haplotype"]["HAP1"]["status"] = "fail"
        self.write_metadata(self.metadata)
        row = catalog.scan_sample(self.root, "test", {}, {})
        self.assertEqual(row["has_compact"], "False")
        self.assertEqual(row["processing_status"], "incomplete")


if __name__ == "__main__":
    unittest.main()
