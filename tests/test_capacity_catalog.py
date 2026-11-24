from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_capacity_catalog_keeps_shared_raw_separate_and_uses_standard_compact(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    columns = [
        "sample_id", "platform", "read_type", "project_id", "data_group",
        "wgs_url", "wgs_input_bytes", "wgs_index_url", "wgs_index_bytes",
        "rna_url", "rna_input_bytes", "rna_index_url", "rna_index_bytes", "reference_id",
    ]
    rows = []
    for index in range(5):
        sample = f"S{index}"
        rows.append({
            "sample_id": sample, "platform": "illumina", "read_type": "short_read",
            "project_id": "TEST", "data_group": "wgs", "wgs_url": f"https://example.invalid/{sample}.cram",
            "wgs_input_bytes": "1000", "wgs_index_url": "", "wgs_index_bytes": "0",
            "rna_url": "", "rna_input_bytes": "0", "rna_index_url": "", "rna_index_bytes": "0", "reference_id": "TEST_REF",
        })
        layout = tmp_path / "output" / "pre_training" / "short_read" / "test" / sample / "compact" / "full_shared_hap_4bit_acgtnp_depth"
        layout.mkdir(parents=True)
        (layout / "chunk.bin").write_bytes(b"x" * 2048)
        (layout / "manifest.json").write_text(json.dumps({"status": "complete", "validation_level": "standard"}))
    write_tsv(manifest, columns, rows)
    for name, header in {
        "read.tsv": ["sample_id", "source_url", "path", "size_bytes", "md5"],
        "vcf.tsv": ["project", "sample_id", "file_url", "size_bytes", "file_type", "checksum_type", "checksum_value"],
        "ref.tsv": ["reference_id", "file_url", "size_bytes", "resource_type", "checksum_type", "checksum_value"],
        "sources.tsv": ["project", "file_url_or_template", "index_url_or_template"],
    }.items():
        write_tsv(tmp_path / name, header, [])
    output = tmp_path / "catalog"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/build_capacity_catalog.py"),
        "--manifest", str(manifest), "--read-files", str(tmp_path / "read.tsv"),
        "--vcf-files", str(tmp_path / "vcf.tsv"), "--reference-files", str(tmp_path / "ref.tsv"),
        "--vcf-sources", str(tmp_path / "sources.tsv"), "--output-root", str(tmp_path / "output"),
        "--out-dir", str(output),
    ], check=True, capture_output=True, text=True)
    summary = json.loads((output / "capacity_summary.json").read_text())
    calibration = read_tsv(output / "zarr_calibration.tsv")
    capacity = read_tsv(output / "sample_capacity_manifest.tsv")
    assert summary["samples"] == 5
    assert summary["observed_compact_count"] == 5
    illumina = next(row for row in calibration if row["platform"] == "illumina")
    assert illumina["eligible"] == "1"
    assert all(row["raw_size_complete"] == "1" for row in capacity)
    assert all(row["wgs_estimate_method"] == "observed_illumina_allocated_compact_ratio" for row in capacity)

