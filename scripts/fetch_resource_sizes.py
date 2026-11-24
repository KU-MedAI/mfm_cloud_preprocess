#!/usr/bin/env python3
"""원격 raw resource의 Content-Length를 cache로 동결한다.

다운로드하지 않고 HTTP/FTP HEAD를 먼저 요청하고, HEAD가 거부된 경우 Range 0-0 요청의
Content-Range만 이용한다. 실패는 size=0으로 치환하지 않고 cache에 상태와 오류를 남긴다.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

CACHE_COLUMNS = ["source_url", "size_bytes", "size_status", "size_source", "http_status", "message", "checked_at_utc"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CACHE_COLUMNS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def as_int(value: object) -> int | None:
    text = str(value or "").strip()
    return int(text) if text.isdigit() else None


def checked_row(url: str, size: int | None, status: str, source: str, http_status: str, message: str) -> dict[str, Any]:
    return {
        "source_url": url, "size_bytes": "" if size is None else size,
        "size_status": status, "size_source": source, "http_status": http_status,
        "message": message[:300], "checked_at_utc": datetime.now(timezone.utc).isoformat(),
    }


def probe(url: str, timeout: float) -> dict[str, Any]:
    try:
        request = Request(url, method="HEAD", headers={"User-Agent": "curl/8.5.0"})
        with urlopen(request, timeout=timeout) as response:
            length = as_int(response.headers.get("Content-Length"))
            if length is not None:
                return checked_row(url, length, "remote_content_length", "HEAD", str(response.status), "")
            status = str(response.status)
    except (HTTPError, URLError, TimeoutError, OSError) as error:
        status = str(getattr(error, "code", ""))
        head_error = str(error)
    else:
        head_error = "HEAD Content-Length missing"
    try:
        request = Request(url, headers={"Range": "bytes=0-0", "User-Agent": "curl/8.5.0"})
        with urlopen(request, timeout=timeout) as response:
            content_range = response.headers.get("Content-Range", "")
            matched = re.search(r"/(\d+)$", content_range)
            if matched:
                return checked_row(url, int(matched.group(1)), "remote_content_range", "Range 0-0", str(response.status), "")
            length = as_int(response.headers.get("Content-Length"))
            # 200/full-body 응답은 실수로 파일을 읽지 않기 위해 size로 사용하지 않는다.
            return checked_row(url, None, "unresolved", "Range 0-0", str(response.status), f"Content-Range missing; Content-Length={length}")
    except (HTTPError, URLError, TimeoutError, OSError) as error:
        return checked_row(url, None, "probe_failed", "HEAD+Range", status, f"{head_error}; {error}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=12)
    parser.add_argument("--timeout", type=float, default=20.0)
    parser.add_argument("--retry-failed", action="store_true", help="기존 failed/unresolved URL도 재시도")
    args = parser.parse_args()
    if args.jobs < 1 or args.jobs > 32:
        raise SystemExit("jobs는 1..32")
    inventory = read_tsv(args.inventory)
    existing = {row.get("source_url", ""): row for row in read_tsv(args.out) if row.get("source_url")}
    pending = []
    for item in inventory:
        url = item.get("source_url", "")
        if not url or as_int(item.get("size_bytes")) is not None:
            continue
        prior = existing.get(url, {})
        if as_int(prior.get("size_bytes")) is not None:
            continue
        if prior and not args.retry_failed:
            continue
        pending.append(url)
    print(f"size_probe_pending={len(pending)} existing_cache={len(existing)}", file=sys.stderr)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        future_to_url = {pool.submit(probe, url, args.timeout): url for url in pending}
        for index, future in enumerate(as_completed(future_to_url), 1):
            url = future_to_url[future]
            try:
                existing[url] = future.result()
            except Exception as error:  # pragma: no cover - final defensive record
                existing[url] = checked_row(url, None, "probe_failed", "internal", "", repr(error))
            if index % 100 == 0 or index == len(pending):
                write_tsv(args.out, sorted(existing.values(), key=lambda row: row["source_url"]))
                print(f"size_probe_complete={index}/{len(pending)}", file=sys.stderr)
    write_tsv(args.out, sorted(existing.values(), key=lambda row: row["source_url"]))
    resolved = sum(as_int(row.get("size_bytes")) is not None for row in existing.values())
    print(f"size_probe_cache={len(existing)} resolved={resolved} unresolved={len(existing)-resolved}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
