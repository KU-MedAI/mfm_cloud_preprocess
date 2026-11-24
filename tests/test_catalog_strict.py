from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("catalog_v2", ROOT / "21_make_catalog.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class StrictCatalogTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.sample = Path(self.tmp.name) / "S1"
        (self.sample / "metadata").mkdir(parents=True)
        self.layout = "full_shared_hap_4bit_acgtnp_depth__all_to_both"
        self.meta = {
            "sequence_variant_policy": "snv_indel",
            "dataset_role": "pre_training",
            "read_type": "short_read",
            "project_id": "test",
            "platform": "illumina",
            "reference_build": "hg38",
            "sample_status": "healthy",
            "layout_name": self.layout,
            "hap_read_mode": "all_to_both",
        }
        self.write_meta()
        compact = self.sample / "compact" / self.layout
        compact.mkdir(parents=True)
        records = {f"S1/chr{i}": {"status": "pass", "contig": f"chr{i}"} for i in range(1, 23)}
        (compact / "manifest.json").write_text(json.dumps({
            "status": "complete", "validation_level": "full", "records": records,
        }))
        trace_dir = self.sample / "qc" / "stage12_trace" / "all_to_both"
        trace_dir.mkdir(parents=True)
        (trace_dir / "trace_manifest.json").write_text(json.dumps({
            "status": "pass", "compact_root": str(compact.resolve()),
            "contigs": {f"chr{i}": "pass" for i in range(1, 23)},
        }))

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def write_meta(self) -> None:
        (self.sample / "metadata" / "sample.json").write_text(json.dumps(self.meta))

    def scan(self) -> dict:
        return MODULE.scan_sample(self.sample, "cohort", {}, {})

    def test_full_22_and_trace_are_required(self) -> None:
        row = self.scan()
        self.assertEqual(row["training_eligible"], "True")
        self.assertEqual(row["compact_contigs"], 22)
        manifest = self.sample / "compact" / self.layout / "manifest.json"
        payload = json.loads(manifest.read_text())
        payload["validation_level"] = "standard"
        manifest.write_text(json.dumps(payload))
        row = self.scan()
        self.assertEqual(row["training_eligible"], "False")
        self.assertIn("compact_full_22", row["missing"])

    def test_rna_qc_fail_is_not_technical_failure(self) -> None:
        phase = self.sample / "rna_personal" / "panel"
        for hap in ("HAP1", "HAP2"):
            (phase / hap / "coverage.zarr").mkdir(parents=True)
        (phase / "allelic_counts.tsv").write_text("gene_id\tcount\nG1\t1\n")
        (phase / "allelic_counts.meta.json").write_text(json.dumps({"counting_unit": "fragment"}))
        (phase / "ase_qc.json").write_text(json.dumps({"sample_training_ready": False}))
        (phase / "ase_training_eligibility.tsv").write_text("sample_training_ready\nFalse\n")
        self.meta["dataset_role"] = "post_training"
        self.meta["rna_personal"] = {
            "rna_personal_root": str(phase),
            "phase_mode": "panel",
            "coverage_zarr": {
                "HAP1": str(phase / "HAP1" / "coverage.zarr"),
                "HAP2": str(phase / "HAP2" / "coverage.zarr"),
            },
            "allelic_counts": str(phase / "allelic_counts.tsv"),
            "ase_qc_json": str(phase / "ase_qc.json"),
            "ase_training_eligibility": str(phase / "ase_training_eligibility.tsv"),
        }
        self.write_meta()
        row = self.scan()
        self.assertEqual(row["qc_all_paths_exist"], "True")
        self.assertEqual(row["processing_status"], "processed_qc_fail")
        self.assertEqual(row["training_eligible"], "False")


if __name__ == "__main__":
    unittest.main()
