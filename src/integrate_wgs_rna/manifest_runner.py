from __future__ import annotations

import csv
import json
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .legacy_shared_hap_converter import (
    LegacyConvertConfig,
    convert_contig,
    human_bytes,
    prepare_sample_root,
    resolve_contigs,
    run_legacy_convert,
    write_metadata,
)


@dataclass
class ManifestRunConfig:
    manifest: Path
    status_tsv: Path
    overwrite: bool = False
    jobs: int = 1


@dataclass
class ManifestContigRunConfig:
    manifest: Path
    status_tsv: Path
    overwrite: bool = False
    jobs: int = 24


def run_manifest(config: ManifestRunConfig) -> list[dict[str, Any]]:
    rows = read_manifest(config.manifest)
    config.status_tsv.parent.mkdir(parents=True, exist_ok=True)
    if config.jobs <= 1:
        results = [run_manifest_row(row, config.overwrite) for row in rows]
    else:
        results = []
        with ProcessPoolExecutor(max_workers=config.jobs) as executor:
            future_to_sample = {
                executor.submit(run_manifest_row, row, config.overwrite): row.get("sample_id", "")
                for row in rows
            }
            for future in as_completed(future_to_sample):
                results.append(future.result())
    results.sort(key=lambda row: row.get("manifest_order", 0))
    write_sample_status(config.status_tsv, results)
    print(f"[MANIFEST] wrote status={config.status_tsv}")
    return results


def run_manifest_contigs(config: ManifestContigRunConfig) -> list[dict[str, Any]]:
    rows = read_manifest(config.manifest)
    config.status_tsv.parent.mkdir(parents=True, exist_ok=True)

    sample_configs: dict[str, LegacyConvertConfig] = {}
    sample_order: dict[str, int] = {}
    tasks: list[tuple[int, str, str, LegacyConvertConfig]] = []
    for row in rows:
        sample_config = build_legacy_config(row, config.overwrite)
        prepare_sample_root(sample_config)
        contigs = parse_contigs(row.get("contigs", "")) or resolve_contigs(sample_config.legacy_root, None)
        sample_configs[sample_config.sample_id] = sample_config
        sample_order[sample_config.sample_id] = int(row.get("_manifest_order", 0))
        for contig in contigs:
            tasks.append((sample_order[sample_config.sample_id], sample_config.sample_id, contig, sample_config))

    print(f"[MANIFEST] contig_tasks={len(tasks)} samples={len(sample_configs)} jobs={config.jobs}")
    results: list[dict[str, Any]] = []
    if config.jobs <= 1:
        for task in tasks:
            result = run_contig_task(task)
            print_contig_result(result)
            results.append(result)
    else:
        with ProcessPoolExecutor(max_workers=config.jobs) as executor:
            future_to_task = {executor.submit(run_contig_task, task): task for task in tasks}
            for future in as_completed(future_to_task):
                result = future.result()
                print_contig_result(result)
                results.append(result)

    results.sort(key=lambda row: (row.get("manifest_order", 0), contig_sort_key(str(row.get("contig", "")))))
    write_contig_status(config.status_tsv, results)
    write_sample_metadata_from_contigs(sample_configs, sample_order, results)
    print(f"[MANIFEST] wrote contig status={config.status_tsv}")
    return results


def run_contig_task(task: tuple[int, str, str, LegacyConvertConfig]) -> dict[str, Any]:
    manifest_order, sample_id, contig, sample_config = task
    start = time.monotonic()
    try:
        summary = convert_contig(sample_config, contig)
        elapsed = time.monotonic() - start
        return {
            "manifest_order": manifest_order,
            "sample_id": sample_id,
            "contig": contig,
            "status": summary.get("status", "ok"),
            "elapsed_seconds": f"{elapsed:.3f}",
            "elapsed_hms": format_hms(elapsed),
            "read_count": int(summary.get("read_count", 0)),
            "read_bases": int(summary.get("read_bases", 0)),
            "size_bytes": int(summary.get("size_bytes", 0)),
            "size_human": summary.get("size_human", human_bytes(int(summary.get("size_bytes", 0)))),
            "hap_invalid_base_count": int(summary.get("hap_invalid_base_count", 0)),
            "hap_invalid_base_symbols": json.dumps(summary.get("hap_invalid_base_symbols", {}), sort_keys=True),
            "source": summary.get("source", str(sample_config.legacy_root / f"{contig}.zarr")),
            "target": summary.get(
                "target",
                str(sample_config.output_sample_root / "compact" / "legacy_shared_hap_4bit_acgtnp" / contig),
            ),
            "error": "",
        }
    except Exception as exc:  # noqa: BLE001 - status file should preserve full failure context.
        elapsed = time.monotonic() - start
        return {
            "manifest_order": manifest_order,
            "sample_id": sample_id,
            "contig": contig,
            "status": "failed",
            "elapsed_seconds": f"{elapsed:.3f}",
            "elapsed_hms": format_hms(elapsed),
            "read_count": 0,
            "read_bases": 0,
            "size_bytes": 0,
            "size_human": "0.00 B",
            "source": str(sample_config.legacy_root / f"{contig}.zarr"),
            "target": str(sample_config.output_sample_root / "compact" / "legacy_shared_hap_4bit_acgtnp" / contig),
            "error": f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}",
        }


def print_contig_result(result: dict[str, Any]) -> None:
    print(
        "[CONTIG] "
        f"sample={result['sample_id']} contig={result['contig']} status={result['status']} "
        f"elapsed={result['elapsed_seconds']}s size={result['size_human']}"
    )


def write_sample_metadata_from_contigs(
    sample_configs: dict[str, LegacyConvertConfig],
    sample_order: dict[str, int],
    results: list[dict[str, Any]],
) -> None:
    grouped: dict[str, list[dict[str, Any]]] = {sample_id: [] for sample_id in sample_configs}
    for result in results:
        grouped.setdefault(str(result["sample_id"]), []).append(result)
    for sample_id, sample_config in sample_configs.items():
        contig_results = grouped.get(sample_id, [])
        timing = [
            {
                "name": f"convert_{result['contig']}",
                "elapsed_seconds": float(result.get("elapsed_seconds") or 0),
                "status": result.get("status", ""),
            }
            for result in contig_results
        ]
        summary = {
            "sample_id": sample_config.sample_id,
            "source_sample_id": sample_config.source_sample_id,
            "dataset_role": sample_config.dataset_role,
            "read_type": sample_config.read_type,
            "project_id": sample_config.project_id,
            "legacy_root": str(sample_config.legacy_root),
            "output_sample_root": str(sample_config.output_sample_root),
            "encoding": "4bit_acgtnp",
            "hap_invalid_policy": sample_config.hap_invalid_policy,
            "sequence_allowed_symbols": "A,C,G,T,N,P",
            "read_allowed_symbols": "A,C,G,T,N",
            "manifest_order": sample_order.get(sample_id, 0),
            "contigs": contig_results,
            "timing": timing,
            "status": "complete" if contig_results and all(row.get("status") in {"converted", "skipped_existing", "ok"} for row in contig_results) else "partial_or_failed",
        }
        write_metadata(sample_config, summary)


def build_legacy_config(row: dict[str, str], overwrite: bool) -> LegacyConvertConfig:
    return LegacyConvertConfig(
        sample_id=row["sample_id"],
        source_sample_id=row["source_sample_id"],
        legacy_root=Path(row["legacy_root"]),
        output_sample_root=Path(row["output_sample_root"]),
        dataset_role=row["dataset_role"],
        read_type=row["read_type"],
        project_id=row["project_id"],
        contigs=parse_contigs(row.get("contigs", "")),
        layout_name=row.get("layout_name") or "legacy_shared_hap_4bit_acgtnp",
        overwrite=overwrite or parse_bool(row.get("overwrite", "")),
        chunk_records=int(row.get("chunk_records") or 1024),
        hap_invalid_policy=row.get("hap_invalid_policy") or "fail",
        rna_bigwigs=parse_paths(row.get("rna_bigwigs", "")),
        source_wgs=parse_paths(row.get("source_wgs", "")),
        source_vcfs=parse_paths(row.get("source_vcfs", "")),
    )


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
    for index, row in enumerate(rows):
        row["_manifest_order"] = str(index)
    return rows


def run_manifest_row(row: dict[str, str], overwrite: bool) -> dict[str, Any]:
    start = time.monotonic()
    sample_id = row["sample_id"]
    try:
        config = build_legacy_config(row, overwrite)
        summary = run_legacy_convert(config)
        elapsed = time.monotonic() - start
        size_bytes = sum(int(item.get("size_bytes", 0)) for item in summary.get("contigs", []))
        return {
            "manifest_order": int(row.get("_manifest_order", 0)),
            "sample_id": sample_id,
            "status": "ok",
            "elapsed_seconds": f"{elapsed:.3f}",
            "elapsed_hms": format_hms(elapsed),
            "contig_count": len(summary.get("contigs", [])),
            "read_count": sum(int(item.get("read_count", 0)) for item in summary.get("contigs", [])),
            "read_bases": sum(int(item.get("read_bases", 0)) for item in summary.get("contigs", [])),
            "size_bytes": size_bytes,
            "size_human": human_bytes(size_bytes),
            "output_sample_root": row["output_sample_root"],
            "error": "",
        }
    except Exception as exc:  # noqa: BLE001 - status file should preserve full failure context.
        elapsed = time.monotonic() - start
        return {
            "manifest_order": int(row.get("_manifest_order", 0)),
            "sample_id": sample_id,
            "status": "failed",
            "elapsed_seconds": f"{elapsed:.3f}",
            "elapsed_hms": format_hms(elapsed),
            "contig_count": 0,
            "read_count": 0,
            "read_bases": 0,
            "size_bytes": 0,
            "size_human": "0.00 B",
            "output_sample_root": row.get("output_sample_root", ""),
            "error": f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}",
        }


def parse_paths(value: str) -> list[Path]:
    return [Path(item) for item in value.split(";") if item]


def parse_contigs(value: str) -> list[str] | None:
    values = [item.strip() for item in value.split(",") if item.strip()]
    return values or None


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def write_sample_status(path: Path, rows: list[dict[str, Any]]) -> None:
    cols = [
        "manifest_order",
        "sample_id",
        "status",
        "elapsed_seconds",
        "elapsed_hms",
        "contig_count",
        "read_count",
        "read_bases",
        "size_bytes",
        "size_human",
        "output_sample_root",
        "error",
    ]
    write_rows(path, rows, cols)


def write_contig_status(path: Path, rows: list[dict[str, Any]]) -> None:
    cols = [
        "manifest_order",
        "sample_id",
        "contig",
        "status",
        "elapsed_seconds",
        "elapsed_hms",
        "read_count",
        "read_bases",
        "size_bytes",
        "size_human",
        "hap_invalid_base_count",
        "hap_invalid_base_symbols",
        "source",
        "target",
        "error",
    ]
    write_rows(path, rows, cols)


def write_rows(path: Path, rows: list[dict[str, Any]], cols: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=cols, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in cols})


def contig_sort_key(contig: str) -> tuple[int, str]:
    if contig.startswith("chr") and contig[3:].isdigit():
        return (int(contig[3:]), "")
    return (10_000, contig)


def format_hms(seconds: float) -> str:
    total = int(round(seconds))
    hours, rem = divmod(total, 3600)
    minutes, sec = divmod(rem, 60)
    return f"{hours:02d}:{minutes:02d}:{sec:02d}"
