#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path


CODE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE))
from scripts.cleanup_guard import validate_predelete  # noqa: E402


class CleanupGuardTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name) / "S1"
        (self.root / "qc").mkdir(parents=True)
        assets = self.root / "individual_fasta"
        assets.mkdir()
        self.asset_manifest = assets / "permanent_assets.json"
        self.asset_manifest.write_text('{"schema":"mfm_personal_permanent_assets_v1","status":"pass"}\n')

    def tearDown(self) -> None:
        self.temp.cleanup()

    def write_report(self, mode: str) -> Path:
        path = self.root / "qc/final_package_validation.json"
        path.write_text(json.dumps({
            "schema": "mfm_final_package_validation_v1",
            "contract_revision": "2026-10-03_integrity_v2",
            "sample": "S1", "mode": mode, "status": "pass", "failed_checks": [],
            "personal_assets_manifest_sha256": hashlib.sha256(
                self.asset_manifest.read_bytes()
            ).hexdigest(),
        }))
        return path

    def build_rna(self, *, ase_hash: str | None = None) -> tuple[Path, Path]:
        rna = self.root / "rna_personal/panel"
        rna.mkdir(parents=True)
        hap = rna / "hap_variants.tsv"
        hap.write_text("contig\tstatus\tgenes\nchr1\tvalidated\tG1\n")
        digest = hashlib.sha256(hap.read_bytes()).hexdigest()
        (rna / "hap_variants_qc.json").write_text(json.dumps({
            "schema": "mfm_hap_variants_qc_v1", "output_sha256": digest,
        }))
        for name in (
            "allelic_counts.tsv", "allelic_counts.split_half.tsv",
            "allelic_counts.both_haps.tsv", "allelic_counts.mmap_adjusted.tsv",
            "ase_training_eligibility.tsv", "rna_realign_summary.json",
        ):
            (rna / name).write_text("x\n")
        count_paths = {
            "counts": rna / "allelic_counts.tsv",
            "split_half": rna / "allelic_counts.split_half.tsv",
            "both_haps": rna / "allelic_counts.both_haps.tsv",
            "mmap_adjusted": rna / "allelic_counts.mmap_adjusted.tsv",
        }
        count_hashes = {name: hashlib.sha256(path.read_bytes()).hexdigest()
                        for name, path in count_paths.items()}
        count_meta = rna / "allelic_counts.meta.json"
        count_meta.write_text(json.dumps({
            "schema": "mfm_allelic_counts_meta_v5",
            "script_revision": "2026-10-02_validated_hap_variant_evidence",
            "hap_variants_sha256": digest,
            "output_sha256": count_hashes,
        }))
        (rna / "ase_qc.json").write_text(json.dumps({
            "schema": "mfm_ase_qc_v3",
            "script_revision": "2026-10-02_validated_hap_variant_evidence_cnv_hap_consistency",
            "inputs": {
                "hap_variants_sha256": ase_hash if ase_hash is not None else digest,
                "counts_sha256": count_hashes["counts"],
                "split_half_sha256": count_hashes["split_half"],
                "both_haps_sha256": count_hashes["both_haps"],
                "mmap_adjusted_sha256": count_hashes["mmap_adjusted"],
                "meta_sha256": hashlib.sha256(count_meta.read_bytes()).hexdigest(),
            },
            "sample_qc": {"provenance_ok": True, "duplicates_excluded": True},
        }))
        for hap_name in ("HAP1", "HAP2"):
            (rna / hap_name / "coverage.zarr").mkdir(parents=True)
            (rna / hap_name / "coverage_windows.zarr").mkdir()
        report = self.write_report("rna")
        return rna, report

    def test_valid_rna_contract_passes(self) -> None:
        _, _ = self.build_rna()
        result = validate_predelete(self.root, "S1", "rna")
        self.assertEqual(result["schema"], "mfm_cleanup_predelete_guard_v1")
        self.assertTrue(result["validation_report_fresh"])

    def test_old_pass_does_not_hide_missing_hap_variants(self) -> None:
        rna, _ = self.build_rna()
        sentinel = self.root / "would_be_deleted.txt"
        sentinel.write_text("keep\n")
        (rna / "hap_variants.tsv").unlink()
        with self.assertRaisesRegex(ValueError, "필수 파일"):
            validate_predelete(self.root, "S1", "rna")
        self.assertTrue(sentinel.exists())

    def test_ase_hap_hash_mismatch_is_rejected(self) -> None:
        self.build_rna(ase_hash="0" * 64)
        with self.assertRaisesRegex(ValueError, "ASE schema/revision/input SHA-256"):
            validate_predelete(self.root, "S1", "rna")

    def test_count_content_tampering_is_rejected(self) -> None:
        rna, _ = self.build_rna()
        (rna / "allelic_counts.tsv").write_text("tampered\n")
        with self.assertRaisesRegex(ValueError, "output SHA-256"):
            validate_predelete(self.root, "S1", "rna")

    def test_stale_validation_report_is_rejected(self) -> None:
        rna, report = self.build_rna()
        future = report.stat().st_mtime_ns + 1_000_000_000
        os.utime(rna / "allelic_counts.tsv", ns=(future, future))
        with self.assertRaisesRegex(ValueError, "보고서가 RNA evidence보다 오래됨"):
            validate_predelete(self.root, "S1", "rna")

    def test_wgs_contract_checks_pass_report(self) -> None:
        self.write_report("wgs")
        self.assertEqual(validate_predelete(self.root, "S1", "wgs")["mode"], "wgs")

    def test_personal_asset_manifest_tampering_is_rejected(self) -> None:
        self.write_report("wgs")
        self.asset_manifest.write_text("tampered\n")
        with self.assertRaisesRegex(ValueError, "영구 HAP 자산 manifest SHA-256"):
            validate_predelete(self.root, "S1", "wgs")

    def test_v2_backend_calls_guard_before_first_remove(self) -> None:
        text = (CODE / "20_cleanup_intermediates.sh").read_text()
        self.assertLess(
            text.index("assert_cleanup_predelete_ready || die"),
            text.index('remove_target "$ROOT/$REALIGNED_SUBDIR"'),
        )
        self.assertIn("personal_assets.py verify", text)
        self.assertNotIn("personal_assets.py manifest --individual-fasta-dir", text)

    def test_stage09_anchors_personal_assets_before_cleanup(self) -> None:
        text = (CODE / "09_personal_fasta.sh").read_text()
        self.assertIn("personal_assets.py manifest", text)
        self.assertLess(text.index("personal_assets.py manifest"), text.index('log "09 완료: $S"'))


if __name__ == "__main__":
    unittest.main()
