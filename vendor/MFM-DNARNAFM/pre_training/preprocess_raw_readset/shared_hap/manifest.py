"""Manifest helpers for shared-HAP preprocessing runs.

The CLI updates manifest.json after every chromosome so failed or interrupted
runs remain auditable and resumable. Keeping the write/update logic here makes
the runner easier to read and avoids copy-pasted timestamp/status handling.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .schema import MANIFEST_SCHEMA, OUTPUT_LAYOUT, SCHEMA_VERSION


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    """Atomically write JSON so manifest/QC files are not left half-written."""
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    try:
        partial.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
        os.replace(partial, path)
    finally:
        partial.unlink(missing_ok=True)


def record_key(result: dict[str, Any]) -> str:
    return f"{result['sample']}/{result['contig']}"


def init_manifest(
    *,
    record_tsv: Path,
    output_root: Path,
    windows_per_chunk: int,
    jobs: int,
    schedule: str,
    max_windows: int | None,
    validation_level: str,
    config: dict[str, Any],
) -> dict[str, Any]:
    now = utc_now()
    return {
        "schema": MANIFEST_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "output_layout": OUTPUT_LAYOUT,
        "record_tsv": str(record_tsv),
        "output_root": str(output_root),
        "windows_per_chunk": int(windows_per_chunk),
        "jobs": int(jobs),
        "schedule": schedule,
        "max_windows": max_windows,
        "validation_level": validation_level,
        "config": config,
        "started_at": now,
        "updated_at": now,
        "status": "running",
        "records": {},
    }


def update_manifest_record(manifest: dict[str, Any], result: dict[str, Any]) -> None:
    manifest["records"][record_key(result)] = result
    manifest["updated_at"] = utc_now()


def mark_manifest_failed(manifest: dict[str, Any], exc: BaseException) -> None:
    manifest["status"] = "failed"
    manifest["error"] = repr(exc)
    manifest["updated_at"] = utc_now()


def mark_manifest_complete(manifest: dict[str, Any]) -> None:
    manifest["status"] = "complete"
    manifest["completed_at"] = utc_now()
    manifest["updated_at"] = manifest["completed_at"]
