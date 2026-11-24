#!/usr/bin/env python3
"""main manifest의 library provenance를 기존 shard manifest에 원자적으로 동기화한다."""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from integrate_wgs_rna.library_provenance import PROVENANCE_COLUMNS  # noqa: E402


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV header 없음: {path}")
        return list(reader.fieldnames), list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    temporary = path.with_suffix(path.suffix + ".part")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=ROOT / "manifests/preprocess_manifest.tsv")
    parser.add_argument("--shard-root", type=Path,
                        default=ROOT / "share_package/sharded_delivery/shards")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()

    _, main_rows = read_tsv(args.manifest)
    main = {row["sample_id"]: row for row in main_rows}
    if len(main) != len(main_rows):
        raise SystemExit("main manifest sample_id 중복")
    missing_columns = [column for column in PROVENANCE_COLUMNS if column not in main_rows[0]]
    if missing_columns:
        raise SystemExit(f"main manifest provenance 열 누락: {missing_columns}")

    prepared: list[tuple[Path, list[str], list[dict[str, str]]]] = []
    seen: set[str] = set()
    for path in sorted(args.shard_root.glob("shard_*/shard_manifest.tsv")):
        fields, rows = read_tsv(path)
        base_fields = [field for field in fields if field not in PROVENANCE_COLUMNS]
        insert_at = base_fields.index("data_group") + 1 if "data_group" in base_fields else len(base_fields)
        output_fields = base_fields[:insert_at] + list(PROVENANCE_COLUMNS) + base_fields[insert_at:]
        for row in rows:
            sample = row["sample_id"]
            if sample in seen:
                raise SystemExit(f"shard sample 중복: {sample}")
            if sample not in main:
                raise SystemExit(f"main manifest에 없는 shard sample: {sample}")
            seen.add(sample)
            for column in PROVENANCE_COLUMNS:
                row[column] = main[sample][column]
        prepared.append((path, output_fields, rows))
    if not prepared:
        raise SystemExit(f"shard manifest 없음: {args.shard_root}")
    if seen != set(main):
        missing = sorted(set(main) - seen)
        raise SystemExit(f"shard에 없는 main sample: count={len(missing)} first={missing[:5]}")

    if not args.check_only:
        for path, fields, rows in prepared:
            write_tsv(path, fields, rows)
    print(f"shards={len(prepared)} samples={len(seen)} mode={'check' if args.check_only else 'sync'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
