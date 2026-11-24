#!/usr/bin/env python3
"""기존 execution_provenance.jsonl에서 사람이 읽는 stage_runtime.tsv를 재생성한다."""
from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict
from datetime import datetime
from pathlib import Path


FIELDS = [
    "run_id", "sample_id", "branch", "stage", "stage_name", "status",
    "started_at", "finished_at", "elapsed_seconds", "elapsed_minutes",
    "exit_code", "reason",
]


def key_of(record: dict) -> tuple[str, ...]:
    return tuple(str(record.get(name, "")) for name in (
        "run_id", "sample_id", "branch", "stage", "stage_name",
    ))


def elapsed_seconds(started_at: str, finished_at: str) -> float:
    try:
        return max(0.0, (datetime.fromisoformat(finished_at) - datetime.fromisoformat(started_at)).total_seconds())
    except ValueError:
        return 0.0


def clean(value: object) -> str:
    return str(value if value is not None else "").replace("\t", " ").replace("\r", " ").replace("\n", " ")


def make_row(finish: dict, start: dict | None) -> dict[str, str]:
    started_at = str((start or finish).get("timestamp", ""))
    finished_at = str(finish.get("timestamp", started_at))
    seconds = finish.get("elapsed_seconds")
    if seconds is None:
        seconds = elapsed_seconds(started_at, finished_at)
    seconds = float(seconds)
    return {
        "run_id": clean(finish.get("run_id", "")),
        "sample_id": clean(finish.get("sample_id", "")),
        "branch": clean(finish.get("branch", "")),
        "stage": clean(finish.get("stage", "")),
        "stage_name": clean(finish.get("stage_name", "")),
        "status": clean(finish.get("status", "skipped" if finish.get("event") == "skip" else "")),
        "started_at": clean(started_at),
        "finished_at": clean(finished_at),
        "elapsed_seconds": f"{seconds:g}",
        "elapsed_minutes": f"{seconds / 60:.3f}",
        "exit_code": clean(finish.get("exit_code", 0 if finish.get("event") == "skip" else "")),
        "reason": clean(finish.get("reason", "")),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="execution_provenance.jsonl")
    parser.add_argument("--output", required=True, type=Path, help="재생성할 stage_runtime.tsv")
    args = parser.parse_args()

    starts: dict[tuple[str, ...], list[dict]] = defaultdict(list)
    rows: list[dict[str, str]] = []
    with args.input.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise SystemExit(f"{args.input}:{line_number}: JSON 오류: {exc}")
            event = record.get("event")
            key = key_of(record)
            if event == "start":
                starts[key].append(record)
            elif event == "finish":
                start = starts[key].pop() if starts[key] else None
                rows.append(make_row(record, start))
            elif event == "skip":
                rows.append(make_row(record, None))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f".{args.output.name}.part.{os.getpid()}")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
