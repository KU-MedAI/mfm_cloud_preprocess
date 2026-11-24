from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

FULL_LAYOUT_NAME = "full_shared_hap_4bit_acgtnp_depth"

def upsert_sample_catalog(
    catalog_path: Path,
    manifest_rows: list[dict[str, str]],
    results: list[dict[str, Any]],
) -> None:
    row_by_sample = {row["sample_id"]: row for row in manifest_rows}
    catalog_rows, columns = read_catalog(catalog_path)
    existing_index = {row.get("sample_id", ""): index for index, row in enumerate(catalog_rows)}
    for result in results:
        sample_id = str(result.get("sample_id", ""))
        if not sample_id or sample_id not in row_by_sample:
            continue
        manifest_row = row_by_sample[sample_id]
        catalog_row = catalog_row_from_result(manifest_row, result)
        for column in catalog_row:
            if column not in columns:
                columns.append(column)
        if sample_id in existing_index:
            catalog_rows[existing_index[sample_id]].update(catalog_row)
        else:
            existing_index[sample_id] = len(catalog_rows)
            catalog_rows.append(catalog_row)
    write_tsv(catalog_path, catalog_rows, columns)


def read_catalog(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    if not path.exists():
        columns = [
            "sample_id",
            "dataset_role",
            "read_type",
            "project_id",
            "status",
            "data_state",
            "health_context",
            "has_wgs",
            "has_rna",
            "has_bigwig",
            "source_sample_id",
            "source_wgs",
            "source_vcfs",
            "rna_bigwigs",
            "final_sample_root",
            "size_human_from_contigs",
            "read_count_sum",
            "hap_invalid_policy",
            "hap_invalid_base_count",
            "notes",
        ]
        return [], columns
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
        columns = list(reader.fieldnames or [])
    return rows, columns


def catalog_row_from_result(manifest_row: dict[str, str], result: dict[str, Any]) -> dict[str, str]:
    sample_root = Path(str(result.get("output_sample_root") or manifest_row.get("output_sample_root") or ""))
    summary = read_optional_json(sample_root / "metadata" / "preprocess.json")
    compact_summary = summary.get("compact_summary", {}) if isinstance(summary, dict) else {}
    contigs = compact_summary.get("contigs", []) if isinstance(compact_summary, dict) else []
    read_count_sum = sum(int(item.get("read_count", 0) or 0) for item in contigs if isinstance(item, dict))
    hap_invalid_base_count = sum(int(item.get("hap_invalid_base_count", 0) or 0) for item in contigs if isinstance(item, dict))
    size_human = summary.get("size_human_after_cleanup") or result.get("size_human", "")
    source_vcfs = manifest_row.get("source_vcfs", "")
    if not source_vcfs and manifest_row.get("common_vcf_dir"):
        source_vcfs = manifest_row.get("common_vcf_dir", "")
    notes = (
        "ONT WGS-only pre-training; full shared-HAP compact 4bit output complete; "
        "common phased VCF split into sequence nonref VCF for HAP FASTA and het SNV+indel VCF for whatshap haplotagging; "
        "personal FASTA uses sequence_reference_fasta, CRAM/whatshap uses reference_fasta; "
        "PS is preserved as source-provided and missing/zero remains unknown unless synthetic PS is explicitly enabled; "
        "raw_intermediate_removed"
    )
    return {
        "sample_id": manifest_row.get("sample_id", ""),
        "dataset_role": manifest_row.get("dataset_role", "pre_training"),
        "read_type": manifest_row.get("read_type", "long_read"),
        "project_id": manifest_row.get("project_id", ""),
        "status": "complete",
        "data_state": FULL_LAYOUT_NAME,
        "health_context": "healthy_individual" if manifest_row.get("sample_status", "healthy") == "healthy" else manifest_row.get("sample_status", ""),
        "has_wgs": "True",
        # has_rna 는 bigwig 또는 rna_personal(개인좌표 STAR) 중 하나라도 있으면 True.
        # has_bigwig 은 bigwig 전용이므로 그대로 둔다 — 둘을 같은 식으로 두면 구분이 사라진다.
        "has_rna": "True" if (manifest_row.get("rna_bigwigs") or detect_rna_personal(sample_root)) else "False",
        "has_bigwig": "True" if manifest_row.get("rna_bigwigs") else "False",
        "source_sample_id": manifest_row.get("source_sample_id") or manifest_row.get("sample_id", ""),
        "source_wgs": manifest_row.get("source_wgs") or manifest_row.get("alignment_path", ""),
        "source_vcfs": source_vcfs,
        "rna_bigwigs": manifest_row.get("rna_bigwigs", ""),
        "final_sample_root": str(sample_root),
        "size_human_from_contigs": str(size_human),
        "read_count_sum": str(read_count_sum),
        "hap_invalid_policy": "to_n",
        "hap_invalid_base_count": str(hap_invalid_base_count),
        "notes": notes,
    }


def read_optional_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))
def detect_rna_personal(sample_root: Path) -> dict[str, Any]:
    """`rna_personal/<phase>/` 산출물을 찾아 경로를 돌려준다. 없으면 `{}`.

    ★ 왜 필요한가 (2026-09-04)
      `has_rna` 를 `rna_bigwigs` 유무로만 판정하던 것이 버그였다. STAR 개인좌표 재정렬 경로
      (snv_indel)는 bigwig 을 쓰지 않고 `rna_personal/<phase>/HAP{1,2}/coverage.zarr` 를 만든다.
      그래서 RNA 가 실제로 있는데도 `sample.json` 에 `has_rna: false` 가 박혔다 — 신규 10명 전부.
      이 필드로 필터링하면 RNA 가 통째로 누락된다.
      근거: `post-training/journal/2609/post_training_data_issue_260903.md` §5

    ⚠ RNA 파이프라인(`run_rna_sample.sh`)이 WGS 뒤에 도는 경우가 있어, 이 함수가 불리는 시점에
      아직 `rna_personal/` 이 없을 수 있다. 그 경우를 위해 RNA 쪽에서도 완료 시 `sample.json` 을
      갱신한다(v2 RNA Stage 24). 여기서는 "있으면 반영"만 한다.
    """
    root = sample_root / "rna_personal"
    if not root.is_dir():
        return {}
    for phase_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        haps = {
            hap: str(phase_dir / hap / "coverage.zarr")
            for hap in ("HAP1", "HAP2")
            if (phase_dir / hap / "coverage.zarr").is_dir()
        }
        if not haps:
            continue
        out: dict[str, Any] = {
            "rna_personal_root": str(phase_dir),
            "phase_mode": phase_dir.name,
            "coverage_zarr": haps,
        }
        allelic = phase_dir / "allelic_counts.tsv"
        if allelic.is_file():
            out["allelic_counts"] = str(allelic)
        return out
    return {}

def write_tsv(path: Path, rows: list[dict[str, Any]], cols: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=cols, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in cols})
