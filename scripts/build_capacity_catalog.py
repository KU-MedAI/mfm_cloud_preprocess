#!/usr/bin/env python3
"""원본 resource 원장과 Compact Zarr 용량 추정 catalog를 생성한다.

기존 preprocess_manifest.tsv는 실행 계약으로 보존한다. 이 도구는 그 manifest를
입력으로 받아 파일 단위 원장과 확장된 sample_capacity_manifest.tsv를 별도 생성한다.
공용 reference/VCF는 URL 단위로 한 번만 합산하고, sample별 Zarr 용량은 실제 complete
Compact layout 관측값이 충분한 플랫폼에서만 보정한다.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import urlparse

GIB = 1 << 30
TIB = 1 << 40
RESOURCE_COLUMNS = [
    "resource_id", "resource_kind", "scope", "project_id", "source_url",
    "relative_path", "size_bytes", "size_status", "size_source",
    "checksum_type", "checksum_value", "consumer_count",
]
LINK_COLUMNS = ["consumer_type", "consumer_id", "resource_id", "use_role"]
OBSERVATION_COLUMNS = [
    "sample_id", "layout", "output_root", "manifest_path", "pipeline_read_type",
    "manifest_status", "validation_level", "apparent_bytes", "allocated_bytes",
    "observation_status", "exclude_reason",
]
CALIBRATION_COLUMNS = [
    "platform", "observation_n", "ratio_low", "ratio_central", "ratio_high",
    "ratio_max", "method", "eligible",
]
CAPACITY_EXTRA_COLUMNS = [
    "raw_wgs_primary_bytes", "raw_wgs_index_bytes", "raw_rna_primary_bytes",
    "raw_rna_index_bytes", "raw_sample_vcf_bytes", "raw_sample_auxiliary_bytes",
    "sample_raw_total_bytes", "sample_raw_total_gib", "raw_size_complete",
    "raw_missing_resource_count", "shared_required_resource_count",
    "estimated_compact_wgs_bytes_low", "estimated_compact_wgs_bytes_central",
    "estimated_compact_wgs_bytes_high", "estimated_compact_rna_bytes_low",
    "estimated_compact_rna_bytes_central", "estimated_compact_rna_bytes_high",
    "estimated_compact_total_bytes_low", "estimated_compact_total_bytes_central",
    "estimated_compact_total_bytes_high", "wgs_estimate_method",
    "rna_estimate_method", "estimate_confidence", "estimate_version",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, columns: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    temporary.write_text(value, encoding="utf-8")
    os.replace(temporary, path)


def as_int(value: object) -> int | None:
    text = str(value or "").strip()
    return int(text) if text.isdigit() else None


def is_remote(value: str) -> bool:
    return urlparse(value or "").scheme.lower() in {"http", "https", "ftp", "s3"}


def resource_id(url: str) -> str:
    return "res_" + hashlib.sha256(url.encode("utf-8")).hexdigest()[:20]


def quantile(values: list[float], fraction: float) -> float:
    if not values:
        raise ValueError("empty quantile")
    if len(values) == 1:
        return values[0]
    position = (len(values) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return values[lower]
    return values[lower] + (values[upper] - values[lower]) * (position - lower)


def directory_bytes(path: Path) -> tuple[int, int]:
    """Zarr directory의 apparent/allocated bytes를 한 번의 metadata walk로 센다."""
    apparent = 0
    allocated = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            entry = Path(root) / name
            try:
                stat = entry.stat()
            except FileNotFoundError:
                continue
            apparent += stat.st_size
            allocated += getattr(stat, "st_blocks", 0) * 512
    return apparent, allocated


class ResourceLedger:
    def __init__(self, size_cache: dict[str, dict[str, str]]) -> None:
        self.size_cache = size_cache
        self.rows: dict[str, dict[str, Any]] = {}
        self.links: set[tuple[str, str, str, str]] = set()

    def size(self, url: str, reported: object, source: str) -> tuple[int | None, str, str]:
        cached = self.size_cache.get(url, {})
        value = as_int(cached.get("size_bytes"))
        if value is not None:
            return value, cached.get("size_status") or "remote_content_length", cached.get("size_source") or "remote_probe"
        value = as_int(reported)
        if value is not None:
            return value, "source_reported", source
        return None, "missing", "unresolved"

    def add(
        self,
        *,
        consumer_type: str,
        consumer_id: str,
        use_role: str,
        resource_kind: str,
        scope: str,
        project_id: str,
        source_url: str,
        relative_path: str = "",
        reported_size: object = "",
        size_source: str = "source_manifest",
        checksum_type: str = "",
        checksum_value: str = "",
    ) -> tuple[int | None, str]:
        if not is_remote(source_url):
            return None, "not_remote"
        identifier = resource_id(source_url)
        size, status, resolved_source = self.size(source_url, reported_size, size_source)
        prior = self.rows.get(identifier)
        if prior is None:
            self.rows[identifier] = {
                "resource_id": identifier,
                "resource_kind": resource_kind,
                "scope": scope,
                "project_id": project_id,
                "source_url": source_url,
                "relative_path": relative_path,
                "size_bytes": "" if size is None else size,
                "size_status": status,
                "size_source": resolved_source,
                "checksum_type": checksum_type,
                "checksum_value": checksum_value,
                "consumer_count": 0,
            }
        else:
            if prior["scope"] == "sample" and scope == "shared":
                prior["scope"] = "shared"
            if as_int(prior["size_bytes"]) is None and size is not None:
                prior.update({"size_bytes": size, "size_status": status, "size_source": resolved_source})
        self.links.add((consumer_type, consumer_id, identifier, use_role))
        return size, status

    def finalize(self) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
        consumers: Counter[str] = Counter(link[2] for link in self.links)
        rows = []
        for row in self.rows.values():
            result = dict(row)
            result["consumer_count"] = consumers[result["resource_id"]]
            if result["scope"] == "sample" and result["consumer_count"] > 1:
                result["scope"] = "shared"
            rows.append(result)
        links = [
            {"consumer_type": a, "consumer_id": b, "resource_id": c, "use_role": d}
            for a, b, c, d in sorted(self.links)
        ]
        return sorted(rows, key=lambda row: (row["scope"], row["project_id"], row["resource_kind"], row["source_url"])), links


def source_map(rows: list[dict[str, str]], key: str) -> dict[str, list[dict[str, str]]]:
    result: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        result[row.get(key, "")].append(row)
    return result


def expand_panel_template(template: str) -> list[str]:
    return [template.replace("{1..22}", str(chrom)) for chrom in range(1, 23)] if "{1..22}" in template else [template]


def scan_compact_outputs(roots: list[Path], manifest_by_sample: dict[str, dict[str, str]]) -> list[dict[str, Any]]:
    observed: dict[tuple[str, str], dict[str, Any]] = {}
    for root in roots:
        if not root.is_dir():
            continue
        for manifest_path in root.glob("**/compact/*/manifest.json"):
            layout_root = manifest_path.parent
            if not layout_root.name.startswith("full_shared_hap_4bit_acgtnp_depth"):
                continue
            sample = layout_root.parent.parent.name
            source = manifest_by_sample.get(sample)
            try:
                payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as error:
                payload = {"status": "invalid", "validation_level": "", "error": str(error)}
            status = str(payload.get("status", ""))
            validation = str(payload.get("validation_level", ""))
            pipeline_type = "long_read" if "/long_read/" in str(layout_root) else "short_read" if "/short_read/" in str(layout_root) else "unknown"
            reason = ""
            if source is None:
                reason = "sample_not_in_active_manifest"
            elif source.get("read_type") != pipeline_type:
                reason = f"platform_path_mismatch:{source.get('read_type')}!={pipeline_type}"
            elif status != "complete" or validation not in {"full", "standard"}:
                reason = f"not_full_complete:status={status},validation={validation}"
            if reason:
                apparent = allocated = 0
                observation_status = "excluded"
            else:
                apparent, allocated = directory_bytes(layout_root)
                observation_status = "eligible_full" if validation == "full" else "eligible_standard"
            row = {
                "sample_id": sample,
                "layout": layout_root.name,
                "output_root": str(layout_root),
                "manifest_path": str(manifest_path),
                "pipeline_read_type": pipeline_type,
                "manifest_status": status,
                "validation_level": validation,
                "apparent_bytes": apparent,
                "allocated_bytes": allocated,
                "observation_status": observation_status,
                "exclude_reason": reason,
            }
            key = (sample, str(layout_root))
            observed[key] = row
    return sorted(observed.values(), key=lambda row: (row["sample_id"], row["output_root"]))


def make_calibration(
    observations: list[dict[str, Any]], manifest_by_sample: dict[str, dict[str, str]]
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    ratios: dict[str, list[float]] = defaultdict(list)
    for row in observations:
        if not str(row["observation_status"]).startswith("eligible"):
            continue
        source = manifest_by_sample[row["sample_id"]]
        raw = as_int(source.get("wgs_input_bytes")) or 0
        allocated = as_int(row["allocated_bytes"]) or 0
        if raw > 0 and allocated > 0:
            ratios[source["platform"].lower()].append(allocated / raw)
    calibration: dict[str, dict[str, Any]] = {}
    output: list[dict[str, Any]] = []
    for platform in sorted(set(ratios) | {"illumina", "hifi"}):
        values = sorted(ratios.get(platform, []))
        eligible = len(values) >= 5
        result = {
            "platform": platform,
            "observation_n": len(values),
            "ratio_low": quantile(values, 0.10) if values else "",
            "ratio_central": statistics.median(values) if values else "",
            "ratio_high": quantile(values, 0.90) if values else "",
            "ratio_max": max(values) if values else "",
            "method": "allocated_compact_bytes/raw_wgs_primary_bytes; p10/median/p90" if eligible else "insufficient_complete_same_platform_observations",
            "eligible": int(eligible),
        }
        calibration[platform] = result
        output.append(result)
    return calibration, output


def legacy_wgs_estimate(row: dict[str, str]) -> int:
    raw = as_int(row.get("wgs_input_bytes")) or 0
    if row.get("platform", "").lower() == "illumina":
        return max(int(93.028 * GIB), int(raw * 5.5))
    return max(64 * GIB, int(raw * 0.9))


def legacy_rna_estimate(row: dict[str, str]) -> int:
    if row.get("data_group") != "wgs_rna_pair":
        return 0
    return max(8 * GIB, int((as_int(row.get("rna_input_bytes")) or 0) * 4.0))


def aggregate_rows(
    rows: list[dict[str, str]], ledger: ResourceLedger, read_by_sample: dict[str, list[dict[str, str]],],
    vcf_by_project_sample: dict[tuple[str, str], list[dict[str, str]]],
    ref_by_id: dict[str, list[dict[str, str]]], panel_sources: list[dict[str, str]],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, int]]]:
    sample_sizes: dict[str, dict[str, int]] = {}
    for source in rows:
        sample = source["sample_id"]
        project = source["project_id"]
        values = defaultdict(int)
        missing = 0
        def add_sample(role: str, kind: str, url: str, path: str, reported: object, checksum_type: str = "", checksum_value: str = "") -> None:
            nonlocal missing
            if not url:
                return
            size, status = ledger.add(
                consumer_type="sample", consumer_id=sample, use_role=role, resource_kind=kind,
                scope="sample", project_id=project, source_url=url, relative_path=path,
                reported_size=reported, checksum_type=checksum_type, checksum_value=checksum_value,
            )
            if size is None:
                missing += 1
            else:
                values[role] += size

        read_rows = read_by_sample.get(sample, [])
        if read_rows:
            for item in read_rows:
                add_sample("wgs_primary", "wgs_read", item.get("source_url", ""), item.get("path", ""), item.get("size_bytes", ""), "md5" if item.get("md5") else "", item.get("md5", ""))
        else:
            add_sample("wgs_primary", "wgs_alignment", source.get("wgs_url", ""), source.get("alignment_path", ""), source.get("wgs_input_bytes", ""))
        add_sample("wgs_index", "wgs_index", source.get("wgs_index_url", ""), source.get("alignment_path", ""), source.get("wgs_index_bytes", ""))
        if source.get("data_group") == "wgs_rna_pair":
            add_sample("rna_primary", "rna_alignment", source.get("rna_url", ""), source.get("rna_path", ""), source.get("rna_input_bytes", ""))
            add_sample("rna_index", "rna_index", source.get("rna_index_url", ""), source.get("rna_path", ""), source.get("rna_index_bytes", ""))
        for item in vcf_by_project_sample.get((project, sample), []):
            add_sample("sample_vcf", f"sample_{item.get('file_type', 'vcf')}", item.get("file_url", ""), "", item.get("size_bytes", ""), item.get("checksum_type", ""), item.get("checksum_value", ""))

        for item in vcf_by_project_sample.get((project, ""), []):
            ledger.add(
                consumer_type="project", consumer_id=project, use_role="shared_vcf", resource_kind=f"shared_{item.get('file_type', 'vcf')}",
                scope="shared", project_id=project, source_url=item.get("file_url", ""), reported_size=item.get("size_bytes", ""),
                checksum_type=item.get("checksum_type", ""), checksum_value=item.get("checksum_value", ""),
            )
        reference_id = source.get("reference_id", "")
        for item in ref_by_id.get(reference_id, []):
            ledger.add(
                consumer_type="project", consumer_id=project, use_role="reference", resource_kind=f"reference_{item.get('resource_type', 'file')}",
                scope="shared", project_id=project, source_url=item.get("file_url", ""), reported_size=item.get("size_bytes", ""),
                checksum_type=item.get("checksum_type", ""), checksum_value=item.get("checksum_value", ""),
            )
        if source.get("platform") == "illumina":
            for item in panel_sources:
                for key, kind in (("file_url_or_template", "illumina_panel_vcf"), ("index_url_or_template", "illumina_panel_index")):
                    for url in expand_panel_template(item.get(key, "")):
                        ledger.add(
                            consumer_type="project", consumer_id=project, use_role="shared_panel", resource_kind=kind,
                            scope="shared", project_id=project, source_url=url,
                        )
        sample_sizes[sample] = dict(values)
        sample_sizes[sample]["missing"] = missing

    return [dict(row) for row in rows], sample_sizes


def capacity_rows(
    rows: list[dict[str, str]], sample_sizes: dict[str, dict[str, int]], calibration: dict[str, dict[str, Any]], ledger_rows: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    project_shared: Counter[str] = Counter()
    project_shared_count: Counter[str] = Counter()
    for resource in ledger_rows:
        if resource["scope"] == "shared":
            value = as_int(resource["size_bytes"])
            if value is not None:
                project_shared[resource["project_id"]] += value
            project_shared_count[resource["project_id"]] += 1
    result: list[dict[str, Any]] = []
    for source in rows:
        sample = source["sample_id"]
        sizes = sample_sizes[sample]
        platform = source["platform"].lower()
        model = calibration.get(platform, {})
        raw_wgs = as_int(source.get("wgs_input_bytes")) or 0
        if int(model.get("eligible", 0)):
            ratios = [float(model[key]) for key in ("ratio_low", "ratio_central", "ratio_high")]
            wgs = [max(1, int(raw_wgs * ratio)) for ratio in ratios]
            wgs_method = f"observed_{platform}_allocated_compact_ratio"
            confidence = "measured_platform"
        else:
            fallback = legacy_wgs_estimate(source)
            wgs = [fallback, fallback, fallback]
            wgs_method = "legacy_conservative_no_same_platform_observation"
            confidence = "provisional"
        rna = [legacy_rna_estimate(source)] * 3
        rna_method = "legacy_conservative_no_final_rna_observation" if rna[0] else "not_applicable"
        total = [wgs[index] + rna[index] for index in range(3)]
        raw_total = sum(sizes.get(key, 0) for key in ("wgs_primary", "wgs_index", "rna_primary", "rna_index", "sample_vcf"))
        row: dict[str, Any] = dict(source)
        row.update({
            "raw_wgs_primary_bytes": sizes.get("wgs_primary", 0), "raw_wgs_index_bytes": sizes.get("wgs_index", 0),
            "raw_rna_primary_bytes": sizes.get("rna_primary", 0), "raw_rna_index_bytes": sizes.get("rna_index", 0),
            "raw_sample_vcf_bytes": sizes.get("sample_vcf", 0),
            "raw_sample_auxiliary_bytes": sizes.get("wgs_index", 0) + sizes.get("rna_index", 0) + sizes.get("sample_vcf", 0),
            "sample_raw_total_bytes": raw_total, "sample_raw_total_gib": f"{raw_total / GIB:.3f}",
            "raw_size_complete": int(sizes.get("missing", 0) == 0), "raw_missing_resource_count": sizes.get("missing", 0),
            "shared_required_resource_count": project_shared_count[source["project_id"]],
            "estimated_compact_wgs_bytes_low": wgs[0], "estimated_compact_wgs_bytes_central": wgs[1], "estimated_compact_wgs_bytes_high": wgs[2],
            "estimated_compact_rna_bytes_low": rna[0], "estimated_compact_rna_bytes_central": rna[1], "estimated_compact_rna_bytes_high": rna[2],
            "estimated_compact_total_bytes_low": total[0], "estimated_compact_total_bytes_central": total[1], "estimated_compact_total_bytes_high": total[2],
            "wgs_estimate_method": wgs_method, "rna_estimate_method": rna_method,
            "estimate_confidence": confidence, "estimate_version": "raw_resource_compact_v1",
        })
        result.append(row)
    return result


def aggregate_scope(rows: list[dict[str, Any]], dimensions: list[str]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[tuple(str(row.get(dimension, "")) for dimension in dimensions)].append(row)
    output = []
    for key, members in sorted(groups.items()):
        item = {dimension: value for dimension, value in zip(dimensions, key)}
        item.update({
            "sample_count": len(members),
            "sample_raw_total_bytes": sum(as_int(row["sample_raw_total_bytes"]) or 0 for row in members),
            "estimated_compact_total_bytes_low": sum(as_int(row["estimated_compact_total_bytes_low"]) or 0 for row in members),
            "estimated_compact_total_bytes_central": sum(as_int(row["estimated_compact_total_bytes_central"]) or 0 for row in members),
            "estimated_compact_total_bytes_high": sum(as_int(row["estimated_compact_total_bytes_high"]) or 0 for row in members),
        })
        output.append(item)
    return output


def main() -> int:
    here = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=here / "manifests/preprocess_manifest.tsv")
    parser.add_argument("--read-files", type=Path, default=here / "manifests/read_files.tsv")
    parser.add_argument("--vcf-files", type=Path, default=here / "share_package/manifests/public_vcf_file_manifest.tsv")
    parser.add_argument("--reference-files", type=Path, default=here / "share_package/manifests/reference_file_manifest.tsv")
    parser.add_argument("--vcf-sources", type=Path, default=here / "share_package/manifests/public_vcf_sources.tsv")
    parser.add_argument("--size-cache", type=Path, default=here / "manifests/resource_size_cache.tsv")
    parser.add_argument("--output-root", type=Path, action="append", default=[])
    parser.add_argument("--out-dir", type=Path, default=here / "manifests/capacity_catalog")
    args = parser.parse_args()

    rows = read_tsv(args.manifest)
    if not rows or len(rows) != len({row.get("sample_id") for row in rows}):
        raise SystemExit("manifest가 비었거나 sample_id가 중복됨")
    fields = list(rows[0])
    read_by_sample = source_map(read_tsv(args.read_files), "sample_id")
    vcf_by_project_sample: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(args.vcf_files):
        vcf_by_project_sample[(row.get("project", ""), row.get("sample_id", ""))].append(row)
    ref_by_id = source_map(read_tsv(args.reference_files), "reference_id")
    panel_sources = [row for row in read_tsv(args.vcf_sources) if row.get("project") == "1KGP_ILLUMINA"]
    cache = {row.get("source_url", ""): row for row in read_tsv(args.size_cache) if row.get("source_url")}
    ledger = ResourceLedger(cache)
    _raw_rows, sample_sizes = aggregate_rows(rows, ledger, read_by_sample, vcf_by_project_sample, ref_by_id, panel_sources)
    ledger_rows, links = ledger.finalize()
    manifest_by_sample = {row["sample_id"]: row for row in rows}
    roots = args.output_root or [Path(__file__).resolve().parents[1] / "output"]
    observations = scan_compact_outputs(roots, manifest_by_sample)
    calibration, calibration_rows = make_calibration(observations, manifest_by_sample)
    capacity = capacity_rows(rows, sample_sizes, calibration, ledger_rows)
    out = args.out_dir
    write_tsv(out / "raw_resource_inventory.tsv", RESOURCE_COLUMNS, ledger_rows)
    write_tsv(out / "raw_resource_links.tsv", LINK_COLUMNS, links)
    write_tsv(out / "zarr_observations.tsv", OBSERVATION_COLUMNS, observations)
    write_tsv(out / "zarr_calibration.tsv", CALIBRATION_COLUMNS, calibration_rows)
    write_tsv(out / "sample_capacity_manifest.tsv", fields + CAPACITY_EXTRA_COLUMNS, capacity)
    by_project = aggregate_scope(capacity, ["project_id", "platform", "data_group"])
    by_all = aggregate_scope(capacity, ["data_group"])
    write_tsv(out / "cohort_capacity.tsv", list(by_project[0]), by_project)
    write_tsv(out / "data_group_capacity.tsv", list(by_all[0]), by_all)
    known_raw = sum(as_int(row["size_bytes"]) or 0 for row in ledger_rows)
    shared_raw = sum(as_int(row["size_bytes"]) or 0 for row in ledger_rows if row["scope"] == "shared")
    summary = {
        "schema": "mfm_raw_resource_compact_capacity_v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "samples": len(capacity),
        "resource_count_unique_url": len(ledger_rows),
        "resource_missing_size_count": sum(as_int(row["size_bytes"]) is None for row in ledger_rows),
        "raw_unique_known_bytes": known_raw,
        "raw_shared_unique_known_bytes": shared_raw,
        "raw_sample_owned_known_bytes": known_raw - shared_raw,
        "observed_compact_count": sum(str(row["observation_status"]).startswith("eligible") for row in observations),
        "observed_compact_excluded_count": sum(not str(row["observation_status"]).startswith("eligible") for row in observations),
        "estimated_compact_total_bytes_low": sum(as_int(row["estimated_compact_total_bytes_low"]) or 0 for row in capacity),
        "estimated_compact_total_bytes_central": sum(as_int(row["estimated_compact_total_bytes_central"]) or 0 for row in capacity),
        "estimated_compact_total_bytes_high": sum(as_int(row["estimated_compact_total_bytes_high"]) or 0 for row in capacity),
        "estimate_note": "Illumina는 eligible complete Compact 관측 비율을 사용하며 validation level은 observation table에 명시; HiFi와 RNA는 동종 완료 관측이 충분할 때까지 legacy conservative를 유지.",
    }
    write_text(out / "capacity_summary.json", json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
