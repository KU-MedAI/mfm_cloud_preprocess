#!/usr/bin/env python3
"""shard download_plan.tsv를 상대경로 그대로 내려받고 크기/checksum을 검증한다."""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def digest(path: Path, algorithm: str) -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            value.update(chunk)
    return value.hexdigest()


def process_row(row: dict[str, str], args: argparse.Namespace) -> tuple[bool, list[str], str]:
    target = args.data_root / row["relative_path"]
    expected = int(row["size_bytes"]) if row.get("size_bytes", "").isdigit() else None
    checksum_type = row.get("checksum_type", "").lower()
    checksum_value = row.get("checksum_value", "").lower()
    messages: list[str] = []
    valid = target.is_file() and (expected is None or target.stat().st_size == expected)
    if valid and checksum_type in {"md5", "sha256"} and checksum_value:
        valid = digest(target, checksum_type) == checksum_value
    if valid:
        return True, [f"[OK] {target}"], str(target)
    if args.verify_only:
        return False, [f"[FAIL] missing-or-invalid {target}"], str(target)
    messages.append(f"[PLAN] {row['source_url']} -> {target}")
    if args.dry_run:
        return True, messages, str(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(target.name + ".part")
    command = [
        "curl", "--fail", "--location", "--retry", "5", "--retry-delay", "10",
        "--continue-at", "-", "--output", str(partial), row["source_url"],
    ]
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        return False, messages + [f"[FAIL] curl rc={completed.returncode} {target}"], str(target)
    if expected is not None and partial.stat().st_size != expected:
        return False, messages + [
            f"[FAIL] size expected={expected} observed={partial.stat().st_size} {partial}"
        ], str(target)
    if checksum_type in {"md5", "sha256"} and checksum_value:
        observed = digest(partial, checksum_type)
        if observed != checksum_value:
            return False, messages + [
                f"[FAIL] {checksum_type} expected={checksum_value} observed={observed} {partial}"
            ], str(target)
    os.replace(partial, target)
    return True, messages + [f"[DOWNLOADED] {target}"], str(target)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be >= 1")
    rows = list(csv.DictReader(args.plan.open(newline="", encoding="utf-8"), delimiter="\t"))
    failed = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        for ok, messages, target in executor.map(lambda row: process_row(row, args), rows):
            for message in messages:
                print(message)
            if not ok:
                failed.append(target)
    if failed:
        print(f"[SUMMARY] failed={len(failed)} total={len(rows)} jobs={args.jobs}")
        return 1
    print(f"[SUMMARY] pass={len(rows)} jobs={args.jobs}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
