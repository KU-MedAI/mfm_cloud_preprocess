#!/usr/bin/env python3
"""v2 catalog 생성기들이 공유하는 파일·집단 메타데이터 reader."""
from __future__ import annotations

import csv
import json
from pathlib import Path


def first_existing(directory: Path, patterns, **fmt) -> Path | None:
    for pattern in patterns:
        candidate = directory / pattern.format(**fmt)
        if candidate.exists():
            return candidate
    return None


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def load_population_meta(path: Path | None) -> dict[str, dict]:
    """집단/가족 메타 TSV를 sample ID 기준으로 색인한다."""
    if path is None or not path.is_file():
        return {}
    out: dict[str, dict] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row.get("sample") or row.get("sample_id")
            if key:
                out[key] = row
    return out


def load_splits(pp_metadata_dir: Path) -> dict[str, str]:
    split_of: dict[str, str] = {}
    split_dir = pp_metadata_dir / "splits"
    for name in ("train", "valid", "test", "ase_holdout"):
        path = split_dir / f"{name}.txt"
        if path.exists():
            for token in path.read_text().split():
                split_of[token.strip()] = name
    return split_of


def read_clair3_qc(sample_root: Path) -> dict:
    """구 catalog schema용 Clair3 QC reader. 파일이 없으면 빈 dict를 돌려준다."""
    path = sample_root / "qc" / "clair3_output_qc.tsv"
    if not path.is_file():
        return {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            return row
    return {}
