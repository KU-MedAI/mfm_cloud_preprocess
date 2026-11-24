from __future__ import annotations

import csv
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_real_manifest_stratified_shards(tmp_path: Path) -> None:
    output = tmp_path / "delivery"
    subprocess.run(
        [
            sys.executable, str(ROOT / "scripts/build_stratified_shards.py"),
            "--out-dir", str(output),
        ],
        check=True, capture_output=True, text=True,
    )
    summary = json.loads((output / "catalog/shard_summary.json").read_text())
    assignments = read_tsv(output / "catalog/sample_assignment.tsv")
    index = read_tsv(output / "catalog/shard_index.tsv")
    qc = read_tsv(output / "catalog/stratification_qc.tsv")
    eta = read_tsv(output / "catalog/runtime_estimate.tsv")
    shared = read_tsv(output / "shared/download_plan.tsv")

    assert summary["samples"] == 3212
    assert summary["wgs_rna_pair"] == 449
    assert summary["wgs_only"] == 2763
    assert summary["shards"] == 37
    assert len(assignments) == len({row["sample_id"] for row in assignments}) == 3212
    assert len(index) == 37
    assert all(row["hard_cap_pass"] == "1" for row in index)
    assert max(float(row["estimated_output_tib"]) for row in index) <= 10.0
    assert len(shared) == len({row["source_url"] for row in shared})
    assert sum(row["resource_kind"] == "illumina_panel_vcf" for row in shared) == 22
    assert sum(row["resource_kind"] == "illumina_panel_index" for row in shared) == 22

    by_stratum: dict[tuple[str, str], int] = defaultdict(int)
    for row in qc:
        if row["dimension"] in {"data_group", "platform", "project_id"}:
            by_stratum[(row["dimension"], row["value"])] += int(row["sample_count"])
    assert by_stratum[("data_group", "wgs_rna_pair")] == 449
    assert by_stratum[("data_group", "wgs")] == 2763
    assert by_stratum[("platform", "illumina")] == 2918
    assert by_stratum[("platform", "hifi")] == 294
    assert by_stratum.get(("platform", "clr"), 0) == 0

    group_counts: dict[str, list[int]] = defaultdict(list)
    for shard in index:
        group_counts["wgs_rna_pair"].append(int(shard["wgs_rna_pair"]))
        group_counts["wgs"].append(int(shard["wgs_only"]))
    assert max(group_counts["wgs_rna_pair"]) - min(group_counts["wgs_rna_pair"]) <= 1
    assert max(group_counts["wgs"]) - min(group_counts["wgs"]) <= 1

    all_eta = next(row for row in eta if row["scope_type"] == "all")
    assert all_eta["sample_count"] == "3212"
    assert all_eta["days_central"] == "222.15"

    first = output / "shards/shard_001"
    ordered = (first / "sample_list.txt").read_text().splitlines()
    pair = (first / "wgs_rna_pair.txt").read_text().splitlines()
    wgs = (first / "wgs_only.txt").read_text().splitlines()
    assert ordered == pair + wgs
    assert json.loads((first / "SHARD_SUMMARY.json").read_text())["hard_cap_pass"]


if __name__ == "__main__":
    import tempfile
    with tempfile.TemporaryDirectory() as directory:
        test_real_manifest_stratified_shards(Path(directory))
