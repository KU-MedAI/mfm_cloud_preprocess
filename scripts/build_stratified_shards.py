#!/usr/bin/env python3
"""최종 산출물 크기를 기준으로 portable manifest를 계층화 shard로 분할한다."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlparse

GIB = 1 << 30
TIB = 1 << 40
EXTRA_COLUMNS = [
    "shard_id", "stratum", "estimated_wgs_output_bytes",
    "estimated_rna_output_bytes", "estimated_output_bytes", "estimated_output_gib",
]
DOWNLOAD_COLUMNS = [
    "sample_id", "resource_kind", "source_url", "relative_path",
    "size_bytes", "checksum_type", "checksum_value", "size_status",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_text(value, encoding="utf-8")
    os.replace(tmp, path)


def as_int(value: str | int | None) -> int:
    text = str(value or "").strip()
    return int(text) if text.isdigit() else 0


def is_remote(value: str) -> bool:
    return urlparse(value or "").scheme.lower() in {"http", "https", "ftp", "s3"}


def stable_tie(seed: int, sample: str) -> str:
    return hashlib.sha256(f"{seed}|{sample}".encode()).hexdigest()


def estimate_output(row: dict[str, str]) -> tuple[int, int, int]:
    """용량 catalog가 있으면 실측 기반 중앙 추정치를 우선하고, 없으면 legacy quota를 쓴다.

    기본 manifest에는 확장 열이 없으므로 기존 전달 shard의 분할 결과는 바뀌지 않는다.
    """
    measured_wgs = as_int(row.get("estimated_compact_wgs_bytes_central"))
    measured_rna = as_int(row.get("estimated_compact_rna_bytes_central"))
    if measured_wgs > 0 and measured_rna >= 0:
        return measured_wgs, measured_rna, measured_wgs + measured_rna
    wgs_input = as_int(row.get("wgs_input_bytes"))
    rna_input = as_int(row.get("rna_input_bytes"))
    platform = row.get("platform", "").lower()
    if platform == "illumina":
        wgs_output = max(int(93.028 * GIB), int(wgs_input * 5.5))
    else:
        # 현재 HiFi 전구간 실측 전까지 기존 long-read 관측치와 입력량의 0.9배 중 큰 값 사용.
        wgs_output = max(64 * GIB, int(wgs_input * 0.9))
    rna_output = 0
    if row.get("data_group") == "wgs_rna_pair":
        rna_output = max(8 * GIB, int(rna_input * 4.0))
    return wgs_output, rna_output, wgs_output + rna_output


def stratum(row: dict[str, str]) -> str:
    return "|".join((row.get("data_group", ""), row.get("platform", ""), row.get("project_id", "")))


def allocate(rows: list[dict[str, object]], shard_count: int, seed: int) -> list[list[dict[str, object]]]:
    shards: list[list[dict[str, object]]] = [[] for _ in range(shard_count)]
    loads = [0] * shard_count
    by_group: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_group[str(row["data_group"])].append(row)
    for group in ("wgs_rna_pair", "wgs"):
        group_rows = by_group.get(group, [])
        quotient, remainder = divmod(len(group_rows), shard_count)
        remaining = [quotient] * shard_count
        for index in sorted(range(shard_count), key=lambda i: (loads[i], i))[:remainder]:
            remaining[index] += 1
        strata: dict[str, list[dict[str, object]]] = defaultdict(list)
        for row in group_rows:
            strata[str(row["stratum"])].append(row)
        stratum_counts: dict[str, list[int]] = {key: [0] * shard_count for key in strata}
        for key in sorted(strata, key=lambda value: (-len(strata[value]), value)):
            members = sorted(
                strata[key],
                key=lambda row: (-int(row["estimated_output_bytes"]), stable_tie(seed, str(row["sample_id"]))),
            )
            for row in members:
                candidates = [i for i, quota in enumerate(remaining) if quota > 0]
                selected = min(candidates, key=lambda i: (stratum_counts[key][i], loads[i], len(shards[i]), i))
                shards[selected].append(row)
                loads[selected] += int(row["estimated_output_bytes"])
                remaining[selected] -= 1
                stratum_counts[key][selected] += 1
        if any(remaining):
            raise RuntimeError(f"unfilled {group} quotas: {remaining}")
    return shards


def remote_downloads(
    row: dict[str, object],
    read_files: dict[str, list[dict[str, str]]],
    sample_vcfs: dict[tuple[str, str], list[dict[str, str]]],
) -> list[dict[str, object]]:
    sample = str(row["sample_id"])
    project = str(row["project_id"])
    result: list[dict[str, object]] = []
    pacbio = read_files.get(sample, [])
    if pacbio:
        for item in pacbio:
            if is_remote(item.get("source_url", "")):
                result.append({
                    "sample_id": sample, "resource_kind": "wgs_read",
                    "source_url": item["source_url"], "relative_path": item["path"],
                    "size_bytes": item.get("size_bytes", ""),
                    "checksum_type": "md5" if item.get("md5") else "",
                    "checksum_value": item.get("md5", ""),
                    "size_status": "known" if item.get("size_bytes") else "unknown",
                })
    elif is_remote(str(row.get("wgs_url", ""))):
        result.append({
            "sample_id": sample, "resource_kind": "wgs_alignment",
            "source_url": row["wgs_url"], "relative_path": row["alignment_path"],
            "size_bytes": row.get("wgs_input_bytes", ""), "checksum_type": "",
            "checksum_value": "", "size_status": "known",
        })
    for kind, url_key, path_key, size_key in (
        ("wgs_index", "wgs_index_url", "alignment_path", "wgs_index_bytes"),
        ("rna_alignment", "rna_url", "rna_path", "rna_input_bytes"),
        ("rna_index", "rna_index_url", "rna_path", "rna_index_bytes"),
    ):
        url = str(row.get(url_key, ""))
        if not is_remote(url):
            continue
        base_path = str(row.get(path_key, ""))
        relative = base_path + (".crai" if kind == "wgs_index" and base_path.endswith(".cram") else ".bai" if kind.endswith("index") else "")
        result.append({
            "sample_id": sample, "resource_kind": kind, "source_url": url,
            "relative_path": relative, "size_bytes": row.get(size_key, ""),
            "checksum_type": "", "checksum_value": "",
            "size_status": "known" if as_int(row.get(size_key)) else "unknown",
        })
    for item in sample_vcfs.get((project, sample), []):
        result.append({
            "sample_id": sample, "resource_kind": "sample_vcf",
            "source_url": item["file_url"],
            "relative_path": f"public_vcf/{project}/{Path(urlparse(item['file_url']).path).name}",
            "size_bytes": item.get("size_bytes", ""),
            "checksum_type": item.get("checksum_type", ""),
            "checksum_value": item.get("checksum_value", ""),
            "size_status": "known" if item.get("size_bytes") else "unknown",
        })
    unique: dict[str, dict[str, object]] = {}
    for item in result:
        unique[str(item["source_url"])] = item
    return list(unique.values())

def shared_downloads(
    reference_rows: list[dict[str, str]],
    vcf_rows: list[dict[str, str]],
    source_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for item in reference_rows:
        url = item.get("file_url", "")
        if not is_remote(url):
            continue
        result.append({
            "sample_id": "", "resource_kind": f"reference_{item.get('resource_type', 'file')}",
            "source_url": url, "relative_path": f"references/source/{Path(urlparse(url).path).name}",
            "size_bytes": item.get("size_bytes", ""), "checksum_type": item.get("checksum_type", ""),
            "checksum_value": item.get("checksum_value", ""),
            "size_status": "known" if item.get("size_bytes") else "unknown",
        })
    for item in vcf_rows:
        if item.get("sample_id") or not is_remote(item.get("file_url", "")):
            continue
        url = item["file_url"]
        result.append({
            "sample_id": "", "resource_kind": f"shared_{item.get('file_type', 'vcf')}",
            "source_url": url,
            "relative_path": f"public_vcf/{item['project']}/{Path(urlparse(url).path).name}",
            "size_bytes": item.get("size_bytes", ""), "checksum_type": item.get("checksum_type", ""),
            "checksum_value": item.get("checksum_value", ""),
            "size_status": "known" if item.get("size_bytes") else "unknown",
        })
    for item in source_rows:
        if item.get("project") != "1KGP_ILLUMINA":
            continue
        for key, kind in (("file_url_or_template", "illumina_panel_vcf"), ("index_url_or_template", "illumina_panel_index")):
            template = item.get(key, "")
            if not is_remote(template):
                continue
            urls = [template.replace("{1..22}", str(chrom)) for chrom in range(1, 23)] if "{1..22}" in template else [template]
            for url in urls:
                result.append({
                    "sample_id": "", "resource_kind": kind, "source_url": url,
                    "relative_path": f"vcf/1KGP_ILLUMINA/{Path(urlparse(url).path).name}",
                    "size_bytes": "", "checksum_type": "", "checksum_value": "", "size_status": "unknown",
                })
    return list({str(item["source_url"]): item for item in result}.values())


def runtime_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """72코어 한 대의 cohort throughput 모델. 다운로드 시간은 제외한다."""
    scopes: list[tuple[str, str, list[dict[str, object]]]] = [("all", "all_3212", rows)]
    for group in sorted({str(row["data_group"]) for row in rows}):
        scopes.append(("data_group", group, [row for row in rows if row["data_group"] == group]))
    for project in sorted({str(row["project_id"]) for row in rows}):
        scopes.append(("project", project, [row for row in rows if row["project_id"] == project]))
    output = []
    for scope_type, value, members in scopes:
        short = sum(str(row["platform"]) == "illumina" for row in members)
        long = len(members) - short
        pair = sum(str(row["data_group"]) == "wgs_rna_pair" for row in members)
        low = short / 18.0 + long / 18.0 + pair / 288.0
        central = short / 15.0 + long / 12.0 + pair / 144.0
        high = short / 12.0 + long / 8.0 + pair / 72.0
        output.append({
            "scope_type": scope_type, "scope_value": value, "sample_count": len(members),
            "short_read_samples": short, "long_read_samples": long, "rna_pair_samples": pair,
            "days_low": f"{low:.2f}", "days_central": f"{central:.2f}", "days_high": f"{high:.2f}",
            "confidence": "medium" if long == 0 else "low",
            "basis": "72core; short WGS 12-18/day (central 15); long WGS 8-18/day (central 12); RNA add-on 72-288/day",
        })
    return output


def main() -> int:
    here = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=here / "manifests/preprocess_manifest.tsv")
    parser.add_argument("--read-files", type=Path, default=here / "manifests/read_files.tsv")
    parser.add_argument("--vcf-files", type=Path, default=here / "share_package/manifests/public_vcf_file_manifest.tsv")
    parser.add_argument("--reference-files", type=Path, default=here / "share_package/manifests/reference_file_manifest.tsv")
    parser.add_argument("--vcf-sources", type=Path, default=here / "share_package/manifests/public_vcf_sources.tsv")
    parser.add_argument("--out-dir", type=Path, default=here / "share_package/sharded_delivery")
    parser.add_argument("--target-output-tib", type=float, default=8.5)
    parser.add_argument("--hard-cap-tib", type=float, default=10.0)
    parser.add_argument("--shards", type=int, default=0, help="0이면 예상 출력/목표치로 자동 산정")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if args.target_output_tib <= 0 or args.hard_cap_tib <= 0 or args.target_output_tib > args.hard_cap_tib:
        raise SystemExit("0 < target-output-tib <= hard-cap-tib 이어야 함")
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise SystemExit(f"기존 비어 있지 않은 out-dir은 덮어쓰지 않음: {args.out_dir}")

    raw_rows = read_tsv(args.manifest)
    if len(raw_rows) != len({row["sample_id"] for row in raw_rows}):
        raise SystemExit("manifest sample_id 중복")
    fields = list(raw_rows[0])
    rows: list[dict[str, object]] = []
    for source in raw_rows:
        row: dict[str, object] = dict(source)
        wgs_out, rna_out, total = estimate_output(source)
        row.update({
            "stratum": stratum(source), "estimated_wgs_output_bytes": wgs_out,
            "estimated_rna_output_bytes": rna_out, "estimated_output_bytes": total,
            "estimated_output_gib": f"{total / GIB:.3f}",
        })
        rows.append(row)

    target = int(args.target_output_tib * TIB)
    cap = int(args.hard_cap_tib * TIB)
    predicted_total = sum(int(row["estimated_output_bytes"]) for row in rows)
    shard_count = args.shards or max(1, math.ceil(predicted_total / target))
    while True:
        shards = allocate(rows, shard_count, args.seed)
        if max(sum(int(row["estimated_output_bytes"]) for row in shard) for shard in shards) <= cap:
            break
        if args.shards:
            raise SystemExit("지정 shard 수로 hard cap을 만족하지 못함")
        shard_count += 1

    read_by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(args.read_files):
        read_by_sample[row["sample_id"]].append(row)
    vcf_file_rows = read_tsv(args.vcf_files)
    sample_vcfs: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in vcf_file_rows:
        if row.get("sample_id") and row.get("file_type") == "vcf":
            sample_vcfs[(row["project"], row["sample_id"])].append(row)

    assignments: list[dict[str, object]] = []
    shard_index: list[dict[str, object]] = []
    strat_qc: list[dict[str, object]] = []
    seen: set[str] = set()
    for number, shard in enumerate(shards, 1):
        shard_id = f"shard_{number:03d}"
        directory = args.out_dir / "shards" / shard_id
        for row in shard:
            row["shard_id"] = shard_id
            sample = str(row["sample_id"])
            if sample in seen:
                raise SystemExit(f"sample 중복 배정: {sample}")
            seen.add(sample)
            assignments.append({
                "sample_id": sample, "shard_id": shard_id, "stratum": row["stratum"],
                "estimated_output_bytes": row["estimated_output_bytes"],
            })
        ordered = sorted(shard, key=lambda row: (0 if row["data_group"] == "wgs_rna_pair" else 1, str(row["sample_id"])))
        write_tsv(directory / "shard_manifest.tsv", fields + EXTRA_COLUMNS, ordered)
        write_text(directory / "sample_list.txt", "".join(f"{row['sample_id']}\n" for row in ordered))
        write_text(directory / "wgs_rna_pair.txt", "".join(f"{row['sample_id']}\n" for row in ordered if row["data_group"] == "wgs_rna_pair"))
        write_text(directory / "wgs_only.txt", "".join(f"{row['sample_id']}\n" for row in ordered if row["data_group"] == "wgs"))
        downloads = []
        for row in ordered:
            downloads.extend(remote_downloads(row, read_by_sample, sample_vcfs))
        by_url = {str(item["source_url"]): item for item in downloads}
        write_tsv(directory / "download_plan.tsv", DOWNLOAD_COLUMNS, list(by_url.values()))
        (directory / "output").mkdir(parents=True, exist_ok=True)
        (directory / "qc").mkdir(parents=True, exist_ok=True)
        predicted = sum(int(row["estimated_output_bytes"]) for row in shard)
        known_input = sum(as_int(row.get("wgs_input_bytes")) + as_int(row.get("rna_input_bytes")) for row in shard)
        counts = Counter(str(row["data_group"]) for row in shard)
        summary = {
            "schema": "mfm_preprocess_shard_v1", "shard_id": shard_id,
            "sample_count": len(shard), "wgs_only": counts["wgs"],
            "wgs_rna_pair": counts["wgs_rna_pair"], "primary_input_bytes": known_input,
            "estimated_output_bytes": predicted, "estimated_output_tib": predicted / TIB,
            "hard_cap_bytes": cap, "hard_cap_pass": predicted <= cap,
            "execution_order": "wgs_rna_pair_first",
        }
        write_text(directory / "SHARD_SUMMARY.json", json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
        shard_index.append({
            "shard_id": shard_id, "sample_count": len(shard), "wgs_rna_pair": counts["wgs_rna_pair"],
            "wgs_only": counts["wgs"], "primary_input_bytes": known_input,
            "primary_input_tib": f"{known_input / TIB:.4f}", "estimated_output_bytes": predicted,
            "estimated_output_tib": f"{predicted / TIB:.4f}", "hard_cap_pass": int(predicted <= cap),
        })
        for dimension in ("data_group", "platform", "project_id"):
            for value in sorted({str(row[dimension]) for row in shard}):
                members = [row for row in shard if str(row[dimension]) == value]
                strat_qc.append({
                    "shard_id": shard_id, "dimension": dimension, "value": value,
                    "sample_count": len(members),
                    "estimated_output_bytes": sum(int(row["estimated_output_bytes"]) for row in members),
                })

    if len(seen) != len(rows):
        raise SystemExit("누락된 sample이 있음")
    catalog = args.out_dir / "catalog"
    write_tsv(catalog / "sample_assignment.tsv", ["sample_id", "shard_id", "stratum", "estimated_output_bytes"], sorted(assignments, key=lambda row: str(row["sample_id"])))
    write_tsv(catalog / "shard_index.tsv", list(shard_index[0]), shard_index)
    write_tsv(catalog / "capacity_summary.tsv", list(shard_index[0]), shard_index)
    write_tsv(catalog / "stratification_qc.tsv", ["shard_id", "dimension", "value", "sample_count", "estimated_output_bytes"], strat_qc)
    eta = runtime_rows(rows)
    write_tsv(catalog / "runtime_estimate.tsv", list(eta[0]), eta)
    shared = shared_downloads(
        read_tsv(args.reference_files), vcf_file_rows, read_tsv(args.vcf_sources)
    )
    write_tsv(args.out_dir / "shared/download_plan.tsv", DOWNLOAD_COLUMNS, shared)
    write_text(
        args.out_dir / "shared/README.md",
        "공통 reference와 population VCF는 모든 shard에서 재사용합니다. source 압축 해제·index 준비는 pipeline 문서를 따르십시오.\n",
    )
    summary = {
        "schema": "mfm_stratified_shards_v1", "samples": len(rows), "shards": shard_count,
        "target_output_tib": args.target_output_tib, "hard_cap_tib": args.hard_cap_tib,
        "predicted_output_bytes": predicted_total, "predicted_output_tib": predicted_total / TIB,
        "wgs_rna_pair": sum(row["data_group"] == "wgs_rna_pair" for row in rows),
        "wgs_only": sum(row["data_group"] == "wgs" for row in rows),
        "estimation_status": "provisional_until_full_illumina_pair_and_hifi_pilots",
        "seed": args.seed,
    }
    write_text(catalog / "shard_summary.json", json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
