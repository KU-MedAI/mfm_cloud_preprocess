from __future__ import annotations

import csv
import json
import subprocess
import tempfile
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_portable_manifest_generation(tmp_path: Path) -> None:
    subprocess.run(
        [sys.executable, str(ROOT / "scripts/build_portable_manifest.py"), "--out-dir", str(tmp_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    manifest = tmp_path / "manifests/preprocess_manifest.tsv"
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_sample = {row["sample_id"]: row for row in rows}
    summary = json.loads((tmp_path / "manifests/portable_summary.json").read_text())

    assert len(rows) == 3212
    assert len(by_sample) == 3212
    assert summary["schema"] == "mfm_portable_manifest_v3"
    assert summary["wgs_only"] == 2763
    assert summary["wgs_rna_pair"] == 449
    assert summary["wgs_primary_bytes"] == 86_820_574_655_817
    assert summary["rna_primary_bytes"] == 1_930_541_687_593
    assert summary["primary_input_bytes"] == 88_751_116_343_410
    assert summary["unique_known_resource_bytes"] == 88_798_744_565_917

    required = {
        "data_group", "wgs_input_bytes", "wgs_input_gib", "rna_input_bytes",
        "rna_input_gib", "wgs_index_bytes", "rna_index_bytes", "vcf_input_bytes",
        "reference_input_bytes", "known_input_total_bytes", "known_input_total_gib",
        "size_completeness", "size_basis",
    }
    assert required <= set(rows[0])
    assert {row["data_group"] for row in rows} == {"wgs", "wgs_rna_pair"}

    pair_first = (tmp_path / "sample_lists/all_pair_first_3212.txt").read_text().splitlines()
    pair = (tmp_path / "sample_lists/wgs_rna_pair_449.txt").read_text().splitlines()
    wgs = (tmp_path / "sample_lists/wgs_only_2763.txt").read_text().splitlines()
    assert len(pair_first) == 3212
    assert pair_first == pair + wgs
    assert all(by_sample[sample]["data_group"] == "wgs_rna_pair" for sample in pair_first[:449])
    assert all(by_sample[sample]["data_group"] == "wgs" for sample in pair_first[449:])

    manifest_text = manifest.read_text()
    assert "/NHNHOME/" not in manifest_text
    assert "/home/aibio_1/" not in manifest_text
    assert all(not Path(row["alignment_path"]).is_absolute() for row in rows if row["alignment_path"])
    assert all(not Path(row["rna_path"]).is_absolute() for row in rows if row["rna_path"])

    with (tmp_path / "manifests/read_files.tsv").open(newline="") as handle:
        read_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(read_rows) == 507
    assert all(not Path(row["path"]).is_absolute() for row in read_rows)


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as directory:
        test_portable_manifest_generation(Path(directory))
