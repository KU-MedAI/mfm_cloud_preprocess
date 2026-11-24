#!/usr/bin/env python3
"""대표 5개 sample의 portable manifest/read-files와 실행 전 용량 보고서를 만든다."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_rows(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(f".{path.name}.part.{os.getpid()}")
    with part.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(part, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("manifests/preprocess_manifest.tsv"))
    parser.add_argument("--read-files", type=Path, default=Path("manifests/read_files.tsv"))
    parser.add_argument("--samples", type=Path, default=Path("manifests/pilot5_samples.txt"))
    parser.add_argument("--out-dir", type=Path, default=Path("manifests/pilot5"))
    args = parser.parse_args()
    selected = [line.strip() for line in args.samples.read_text().splitlines() if line.strip() and not line.startswith("#")]
    selected_set = set(selected)
    if len(selected) != 5 or len(selected_set) != 5:
        raise SystemExit("pilot sample은 중복 없는 5개여야 합니다")

    fields, rows = read_rows(args.manifest)
    by_sample = {row["sample_id"]: row for row in rows}
    missing = [sample for sample in selected if sample not in by_sample]
    if missing:
        raise SystemExit(f"manifest sample 누락: {missing}")
    pilot = [by_sample[sample] for sample in selected]
    read_fields, read_rows_all = read_rows(args.read_files)
    pilot_reads = [row for row in read_rows_all if row.get("sample_id") in selected_set]
    write_rows(args.out_dir / "preprocess_manifest.tsv", fields, pilot)
    write_rows(args.out_dir / "read_files.tsv", read_fields, pilot_reads)
    (args.out_dir / "sample_list.txt").write_text("\n".join(selected) + "\n")
    summary = {
        "schema": "mfm_pilot5_manifest_v1", "sample_count": 5,
        "samples": selected,
        "wgs_rna_pair": sum(row.get("rna_pair_status") == "wgs_rna_pair" for row in pilot),
        "wgs_only": sum(row.get("rna_pair_status") == "wgs_only" for row in pilot),
        "wgs_input_bytes": sum(int(row.get("wgs_input_bytes") or 0) for row in pilot),
        "rna_input_bytes": sum(int(row.get("rna_input_bytes") or 0) for row in pilot),
        "formats": {row["sample_id"]: row.get("wgs_input_format", "") for row in pilot},
        "execution_order": selected,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
