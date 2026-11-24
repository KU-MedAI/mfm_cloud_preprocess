#!/usr/bin/env python3
"""PacBio 295명 read/VCF/reference freeze를 네트워크 없이 엄격히 검증한다."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


class FreezeError(RuntimeError):
    pass


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FreezeError(f"파일이 없거나 비었습니다: {path}")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FreezeError(message)


def validate_checksum(kind: str, value: str, label: str) -> None:
    patterns = {
        "md5": r"[0-9a-f]{32}",
        "md5_s3_etag": r"[0-9a-f]{32}",
        "s3_multipart_etag_64MiB": r"[0-9a-f]{32}-[1-9][0-9]*",
    }
    require(kind in patterns, f"지원하지 않는 checksum_type: {label} {kind}")
    require(bool(re.fullmatch(patterns[kind], value)), f"checksum 형식 오류: {label} {kind}={value}")


def validate_reads(root: Path, cohort: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = read_tsv(root / "manifests/pacbio_file_manifest.tsv")
    expected = {item["sample_id"] for item in cohort}
    observed = {item["sample_id"] for item in rows}
    require(expected == observed, f"PacBio donor set 불일치: missing={sorted(expected-observed)} extra={sorted(observed-expected)}")
    require(len(rows) == 511, f"PacBio file rows expected=511 actual={len(rows)}")
    require(len({(item["sample_id"], item["file_url"]) for item in rows}) == len(rows), "sample/file URL 중복")
    require(all(item["status"] == "ready" for item in rows), "ready 아닌 PacBio 파일 행 존재")
    for item in rows:
        label = f"{item['sample_id']}:{item['source_record_id']}"
        require(item["file_url"].startswith("https://"), f"HTTPS URL 아님: {label}")
        require(item["size_bytes"].isdigit() and int(item["size_bytes"]) > 0, f"size 오류: {label}")
        validate_checksum(item["checksum_type"], item["checksum_value"], label)
        if item["index_url"]:
            require(item["index_url"].startswith("https://"), f"index HTTPS URL 아님: {label}")
            require(item["index_size_bytes"].isdigit() and int(item["index_size_bytes"]) > 0,
                    f"index size 오류: {label}")
            validate_checksum(item["index_checksum_type"], item["index_checksum_value"], label + ":index")
    donors_by_project = {project: len({item["sample_id"] for item in rows if item["project"] == project})
                         for project in {item["project"] for item in rows}}
    require(donors_by_project == {"HPRC_R2": 232, "HGSVC3": 60, "HGSVC2": 2, "HPRC_R1": 1},
            f"project donor 수 불일치: {donors_by_project}")
    hprc2 = [item for item in rows if item["project"] == "HPRC_R2"]
    require(len(hprc2) == 232 and all(item["file_type"] == "aligned_bam" and item["index_url"]
                                      and item["reference"] == "GRCh38_no_alt" for item in hprc2),
            "HPRC R2 aligned BAM/BAI/reference 계약 오류")
    h3 = [item for item in rows if item["project"] == "HGSVC3"]
    require(len(h3) == 267 and all(item["read_technology"] == "PacBio_HiFi" for item in h3),
            "HGSVC3 HiFi run 계약 오류")
    require(len({item["run_accession"] for item in h3}) == len(h3), "HGSVC3 run accession 중복")
    h2 = [item for item in rows if item["project"] == "HGSVC2"]
    require(Counter(item["read_technology"] for item in h2) == Counter({"PacBio_CLR": 4, "PacBio_HiFi": 1}),
            "HGSVC2 CLR/HiFi 기술 구분 오류")
    require({item["sample_id"] for item in h2 if item["read_technology"] == "PacBio_CLR"} == {"HG03125"},
            "HG03125 CLR 표기 오류")
    return rows


def validate_active_reads(root: Path, active_cohort: list[dict[str, str]], historical_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Historical freeze와 분리된 active HiFi-only file-level manifest를 검증한다."""
    rows = read_tsv(root / "manifests/pacbio_active_file_manifest.tsv")
    expected = {item["sample_id"] for item in active_cohort}
    observed = {item["sample_id"] for item in rows}
    require(len(active_cohort) == 294, f"active PacBio donor expected=294 actual={len(active_cohort)}")
    require(expected == observed, f"active PacBio donor set 불일치: missing={sorted(expected-observed)} extra={sorted(observed-expected)}")
    require(len(rows) == 507, f"active PacBio file rows expected=507 actual={len(rows)}")
    require(all(item["read_technology"] == "PacBio_HiFi" for item in rows), "active manifest에 HiFi 외 기술이 있음")
    require("HG03125" not in observed, "active manifest에 historical CLR donor HG03125가 있음")
    require(all(item["status"] == "ready" for item in rows), "active manifest에 ready 아닌 행 존재")
    historical_expected = [item for item in historical_rows if item["sample_id"] != "HG03125"]
    key = lambda item: (item["sample_id"], item["source_record_id"], item["file_url"])
    require({key(item) for item in rows} == {key(item) for item in historical_expected},
            "active 507-file manifest가 historical freeze의 HiFi subset과 다름")
    return rows


def validate_vcf(root: Path, cohort: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = read_tsv(root / "manifests/public_vcf_file_manifest.tsv")
    require(len(rows) == 704, f"public VCF file rows expected=704 actual={len(rows)}")
    require(all(item["status"] == "ready" for item in rows), "ready 아닌 VCF resource 행 존재")
    require(len({(item["project"], item["sample_id"], item["file_url"]) for item in rows}) == len(rows),
            "VCF resource 중복")
    for item in rows:
        label = f"{item['project']}:{item['sample_id']}:{item['file_type']}"
        require(item["file_url"].startswith("https://"), f"VCF HTTPS URL 아님: {label}")
        require(item["size_bytes"].isdigit() and int(item["size_bytes"]) > 0, f"VCF size 오류: {label}")
        validate_checksum(item["checksum_type"], item["checksum_value"], label)
        require(item["reference_id"] == "GRCh38_no_alt_GCA_000001405.15", f"reference ID 오류: {label}")
    selected_hprc2 = {item["sample_id"] for item in cohort if item["selected_project"] == "HPRC_R2"}
    h2 = [item for item in rows if item["project"] == "HPRC_R2"]
    require(len(h2) == 696 and {item["sample_id"] for item in h2} == selected_hprc2,
            "HPRC R2 VCF donor/file coverage 오류")
    by_sample: dict[str, Counter[str]] = defaultdict(Counter)
    for item in h2:
        by_sample[item["sample_id"]][item["file_type"]] += 1
    require(all(counts == Counter({"vcf": 1, "tbi": 1, "tsv": 1}) for counts in by_sample.values()),
            "HPRC R2 sample당 VCF/TBI/phase TSV가 1개씩 아님")
    require(Counter(item["project"] for item in rows if item["project"] != "HPRC_R2") ==
            Counter({"HGSVC3": 4, "HGSVC2": 2, "HPRC_R1": 2}), "shared VCF file 수 오류")
    return rows


def validate_references(root: Path) -> list[dict[str, str]]:
    rows = read_tsv(root / "manifests/reference_file_manifest.tsv")
    require(len(rows) == 3, f"reference resource rows expected=3 actual={len(rows)}")
    require(all(item["status"] == "ready" for item in rows), "ready 아닌 reference 행 존재")
    require(any(item["resource_type"] == "fasta_gz" and item["checksum_type"] == "md5" for item in rows),
            "exact FASTA+MD5 행 없음")
    for item in rows:
        validate_checksum(item["checksum_type"], item["checksum_value"], item["resource_type"])
        require(item["size_bytes"].isdigit() and int(item["size_bytes"]) > 0, "reference size 오류")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-root", type=Path, required=True)
    args = parser.parse_args()
    cohort = read_tsv(args.package_root / "pacbio_295_samples.tsv")
    active_cohort = read_tsv(args.package_root / "pacbio_294_samples.tsv")
    require(len(cohort) == 295, f"PacBio cohort expected=295 actual={len(cohort)}")
    require(Counter(item["selected_platform"] for item in cohort) ==
            Counter({"PacBio_HiFi": 294, "PacBio_CLR": 1}), "PacBio HiFi/CLR donor 구분 오류")
    require({item["sample_id"] for item in cohort if item["selected_platform"] == "PacBio_CLR"} == {"HG03125"},
            "HG03125 CLR catalog 표기 오류")
    require(all(item["wgs_download_status"] == "ready_frozen_manifest" for item in cohort),
            "PacBio catalog에 download-ready가 아닌 donor 존재")
    reads = validate_reads(args.package_root, cohort)
    active_reads = validate_active_reads(args.package_root, active_cohort, reads)
    vcfs = validate_vcf(args.package_root, cohort)
    refs = validate_references(args.package_root)
    print(
        f"[OK] PacBio freeze: historical_donors=295 historical_read_files={len(reads)} "
        f"active_hifi_donors=294 active_hifi_read_files={len(active_reads)} "
        f"vcf_files={len(vcfs)} reference_files={len(refs)}"
    )


if __name__ == "__main__":
    main()
