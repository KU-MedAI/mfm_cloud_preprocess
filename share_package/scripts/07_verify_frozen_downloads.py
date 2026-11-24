#!/usr/bin/env python3
"""phase 다운로드 파일의 크기와 동결 checksum을 URL 기준으로 검증한다."""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
from urllib.parse import unquote, urlparse

CHUNK_64_MIB = 64 * 1024 * 1024


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def multipart_etag(path: Path, part_size: int = CHUNK_64_MIB) -> str:
    parts: list[bytes] = []
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(part_size), b""):
            parts.append(hashlib.md5(chunk).digest())
    if len(parts) <= 1:
        return hashlib.md5(path.read_bytes()).hexdigest()
    return f"{hashlib.md5(b''.join(parts)).hexdigest()}-{len(parts)}"


def add_expected(index: dict[str, dict[str, str]], url: str, size: str = "",
                 checksum_type: str = "", checksum: str = "") -> None:
    if not url:
        return
    value = {"size": size, "checksum_type": checksum_type, "checksum": checksum}
    previous = index.get(url)
    if previous is not None and previous != value:
        raise RuntimeError(f"같은 URL의 검증 정보가 충돌합니다: {url}")
    index[url] = value


def build_index(root: Path) -> dict[str, dict[str, str]]:
    index: dict[str, dict[str, str]] = {}
    for row in rows(root / "manifests/pacbio_file_manifest.tsv"):
        add_expected(index, row["file_url"], row["size_bytes"], row["checksum_type"], row["checksum_value"])
        add_expected(index, row["index_url"], row["index_size_bytes"], row["index_checksum_type"], row["index_checksum_value"])
    for manifest in ("public_vcf_file_manifest.tsv", "reference_file_manifest.tsv"):
        for row in rows(root / "manifests" / manifest):
            add_expected(index, row["file_url"], row["size_bytes"], row["checksum_type"], row["checksum_value"])
    for manifest in ("illumina_30x_files.tsv", "rna_geuvadis_files.tsv"):
        for row in rows(root / "manifests" / manifest):
            add_expected(index, row["file_url"], row.get("size_bytes", ""), "md5" if row.get("md5") else "", row.get("md5", ""))
            add_expected(index, row["index_url"])
    return index


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-root", required=True, type=Path)
    parser.add_argument("--download-root", required=True, type=Path)
    parser.add_argument("--url-list", required=True, type=Path)
    args = parser.parse_args()
    expected = build_index(args.package_root)
    listed = rows(args.url_list)
    errors: list[str] = []
    checked_size = checked_checksum = metadata_missing = 0
    seen: set[str] = set()
    for row in listed:
        url = row["url"]
        if url in seen:
            continue
        seen.add(url)
        local = args.download_root / Path(unquote(urlparse(url).path)).name
        if not local.is_file():
            errors.append(f"다운로드 파일 없음: {local} ({url})")
            continue
        metadata = expected.get(url)
        if metadata is None:
            metadata_missing += 1
            continue
        if metadata["size"]:
            checked_size += 1
            if local.stat().st_size != int(metadata["size"]):
                errors.append(f"크기 불일치: {local} expected={metadata['size']} actual={local.stat().st_size}")
                continue
        kind, value = metadata["checksum_type"], metadata["checksum"]
        if value:
            checked_checksum += 1
            if kind in {"md5", "md5_s3_etag"}:
                actual = md5_file(local)
            elif kind == "s3_multipart_etag_64MiB":
                actual = multipart_etag(local)
            else:
                errors.append(f"지원하지 않는 checksum: {kind} {url}")
                continue
            if actual != value:
                errors.append(f"checksum 불일치: {local} type={kind} expected={value} actual={actual}")
    for error in errors:
        print(f"[ERROR] {error}")
    print(f"[INFO] listed_urls={len(seen)} size_checked={checked_size} checksum_checked={checked_checksum} metadata_missing={metadata_missing}")
    if errors:
        return 1
    print("[OK] phase download size/checksum validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
