#!/usr/bin/env python3
"""main/read-files manifest에서 Stage 01 입력 inventory를 생성하고 main manifest를 보강한다."""
from __future__ import annotations

import argparse
import csv
import json
import os
from collections import Counter, defaultdict
from pathlib import Path


EXTRA_COLUMNS = [
    "wgs_input_format", "wgs_input_file_count", "rna_pair_status",
    "rna_input_format", "expected_pipeline_branch", "manifest_schema_version",
]


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv_atomic(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    tmp = path.with_name(f".{path.name}.part.{os.getpid()}")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def infer_extension(path: str) -> str:
    lower = path.lower()
    for suffix, label in (
        (".fastq.gz", "fastq_gz"), (".fq.gz", "fastq_gz"),
        (".cram", "cram"), (".ubam", "ubam"), (".bam", "bam"),
    ):
        if lower.endswith(suffix):
            return label
    return "unknown"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--read-files", required=True, type=Path)
    ap.add_argument("--output-tsv", required=True, type=Path)
    ap.add_argument("--summary-json", required=True, type=Path)
    ap.add_argument("--augment-manifest", action="store_true")
    ap.add_argument("--expected-samples", type=int, default=3212)
    a = ap.parse_args()

    fields, rows = read_tsv(a.manifest)
    _, file_rows = read_tsv(a.read_files)
    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in file_rows:
        by_sample[row.get("sample_id", "")].append(row)

    seen: set[str] = set()
    inventory: list[dict[str, object]] = []
    format_counts: Counter[str] = Counter()
    platform_counts: Counter[str] = Counter()
    branch_counts: Counter[str] = Counter()
    project_counts: Counter[str] = Counter()
    for row in rows:
        sample = row.get("sample_id", "")
        if not sample or sample in seen:
            raise SystemExit(f"duplicate/empty sample_id: {sample}")
        seen.add(sample)
        files = by_sample.get(sample, []) if row.get("read_files_manifest") else []
        if files:
            formats = sorted({f.get("format", "").lower() for f in files if f.get("format")})
            source_paths = [f.get("path", "") for f in files if f.get("path")]
            file_count = len(files)
        else:
            fmt = infer_extension(row.get("alignment_path", ""))
            formats = [fmt]
            source_paths = [row.get("alignment_path", "")] if row.get("alignment_path") else []
            file_count = 1 if source_paths else 0
        format_label = "+".join(formats) if formats else "unknown"
        has_rna = bool(row.get("rna_path") and row.get("rna_path") != "-")
        rna_format = infer_extension(row.get("rna_path", "")) if has_rna else "not_applicable"
        branch = "wgs_rna_pair" if has_rna else "wgs_only"
        row.update({
            "wgs_input_format": format_label,
            "wgs_input_file_count": str(file_count),
            "rna_pair_status": branch,
            "rna_input_format": rna_format,
            "expected_pipeline_branch": branch,
            "manifest_schema_version": "mfm_preprocess_manifest_v3",
        })
        inventory.append({
            "sample_id": sample, "project_id": row.get("project_id", ""),
            "cohort": row.get("cohort", ""), "platform": row.get("platform", ""),
            "read_type": row.get("read_type", ""), "wgs_input_format": format_label,
            "wgs_input_file_count": file_count,
            "wgs_input_bytes": row.get("wgs_input_bytes", ""),
            "rna_pair_status": branch, "rna_input_format": rna_format,
            "source_paths": ";".join(source_paths),
        })
        format_counts[format_label] += 1
        platform_counts[row.get("platform", "")] += 1
        branch_counts[branch] += 1
        project_counts[row.get("project_id", "")] += 1

    if len(rows) != a.expected_samples:
        raise SystemExit(f"sample count mismatch: {len(rows)} != {a.expected_samples}")
    a.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    write_tsv_atomic(a.output_tsv, list(inventory[0]), inventory)
    summary = {
        "schema": "mfm_input_inventory_v1", "sample_count": len(rows),
        "format_counts": dict(sorted(format_counts.items())),
        "platform_counts": dict(sorted(platform_counts.items())),
        "branch_counts": dict(sorted(branch_counts.items())),
        "project_counts": dict(sorted(project_counts.items())),
    }
    a.summary_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
    if a.augment_manifest:
        write_tsv_atomic(a.manifest, fields + [x for x in EXTRA_COLUMNS if x not in fields], rows)
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
