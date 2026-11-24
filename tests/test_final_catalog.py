from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def write_sample(root: Path, sample: str, status: str) -> None:
    sample_root = root / "wgs" / "P" / sample
    (sample_root / "metadata").mkdir(parents=True)
    (sample_root / "qc").mkdir()
    (sample_root / "windows").mkdir()
    (sample_root / "metadata" / "sample.json").write_text(json.dumps({
        "sample_id": sample, "project_id": "P", "cohort": "C",
        "dataset_role": "pre_training", "read_type": "short_read", "platform": "illumina",
        "layout_name": "layout", "sequence_variant_policy": "snv_indel",
    }))
    (sample_root / "windows" / "windows.tsv").write_text("header\nrow\n")
    compact = sample_root / "compact" / "layout"
    compact.mkdir(parents=True)
    chroms = [f"chr{i}" for i in range(1, 23)]
    (compact / "manifest.json").write_text(json.dumps({
        "status": "complete", "validation_level": "full",
        "records": {c: {"contig": c, "status": "pass"} for c in chroms},
    }))
    trace = sample_root / "qc" / "stage12_trace" / "all_to_both"
    trace.mkdir(parents=True)
    (trace / "trace_manifest.json").write_text(json.dumps({
        "status": "pass", "compact_root": str(compact),
        "contigs": {c: "pass" for c in chroms},
    }))
    (sample_root / "qc" / "final_package_validation.json").write_text(json.dumps({
        "schema": "mfm_final_package_validation_v1", "sample": sample, "mode": "wgs",
        "status": status, "failed_checks": [] if status == "pass" else ["tampered"],
        "layout": "layout", "hap_read_mode": "all_to_both",
    }))


def test_pass_only_catalog() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        write_sample(root, "PASS", "pass")
        write_sample(root, "FAIL", "fail")
        result = subprocess.run(
            [sys.executable, str(ROOT / "27_make_catalog.py"), str(root)],
            capture_output=True, text=True,
        )
        assert result.returncode != 0
        with (root / "catalog" / "sample_catalog.tsv").open(newline="") as handle:
            catalog = list(csv.DictReader(handle, delimiter="\t"))
        with (root / "catalog" / "incomplete_samples.tsv").open(newline="") as handle:
            incomplete = list(csv.DictReader(handle, delimiter="\t"))
        assert [row["sample_id"] for row in catalog] == ["PASS"]
        assert [row["sample_id"] for row in incomplete] == ["FAIL"]
        assert incomplete[0]["failed_checks"] == "tampered"


if __name__ == "__main__":
    test_pass_only_catalog()
