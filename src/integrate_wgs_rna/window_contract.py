"""WGS/RNA 공통 개인 HAP window 계약. 기존 window 선택·순서는 바꾸지 않는다."""
from __future__ import annotations

import csv
import math
import os
import shutil
import tempfile
import time
from collections import defaultdict
from pathlib import Path

N_ANNOTATION_SCHEMA = "mfm_window_n_annotations_v3"
LEGACY_N_FLAG_FIELDS = ("has_n", "n_gt_1pct", "training_n5_eligible")
N_VALUE_FIELDS = ("n_count", "n_nonacgt", "n_fraction", "n_max_run")
N_VALUE_DEFINITIONS = {
    "n_count": "number of N bases in actual window; padding P excluded",
    "n_nonacgt": "number of non-ACGT bases; equals n_count for validated ACGTN",
    "n_fraction": "n_nonacgt / actual_length; padding P excluded",
    "n_max_run": "longest consecutive N run in actual window",
}


def validate_n_values(row: dict) -> None:
    """고정 학습 임계값 없이 실제 길이·N 원값의 일관성만 검증함."""
    length, count, non_acgt, run = (int(row[key]) for key in
                                  ("actual_length", "n_count", "n_nonacgt", "n_max_run"))
    if (length <= 0 or not 0 <= count <= non_acgt <= length
            or not 0 <= run <= count or (count > 0 and run == 0)
            or not math.isclose(float(row["n_fraction"]), non_acgt / length, rel_tol=1e-9, abs_tol=1e-12)):
        raise ValueError("window N 원값/실제 길이 불일치")


def n_value_metadata(metadata: dict) -> dict:
    """파생 flag·고정 학습 cutoff 정의를 제거하고 원값 계약만 남김."""
    updated = dict(metadata)
    for field in ("n_flag_definitions", "training_max_n_fraction"):
        updated.pop(field, None)
    updated.update(n_annotation_schema=N_ANNOTATION_SCHEMA, n_value_definitions=N_VALUE_DEFINITIONS,
                   training_n_selection="consumer_defined_from_n_fraction")
    return updated


def ensure_window_n_annotations(path: Path) -> bool:
    """구캐시의 원값·좌표·선택을 보존하며 중복 N flag 열만 제거함."""
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    if not set(N_VALUE_FIELDS + ("actual_length",)) <= set(fields):
        raise ValueError("window N 원값 열이 누락됨")
    for row in rows:
        if None in row or any(value is None for value in row.values()):
            raise ValueError("window N 원값 TSV 열 개수 불일치")
        validate_n_values(row)
    if not set(LEGACY_N_FLAG_FIELDS).intersection(fields):
        return False
    fields = [field for field in fields if field not in LEGACY_N_FLAG_FIELDS]
    shutil.copy2(path, path.with_name(f"{path.name}.nvalues.{time.time_ns()}"))
    descriptor, name = tempfile.mkstemp(prefix=".windows.nvalues.", suffix=".part", dir=path.parent)
    part = Path(name)
    try:
        with os.fdopen(descriptor, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows({field: row[field] for field in fields} for row in rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(part, path)
    finally:
        part.unlink(missing_ok=True)
    return True


def read_windows(path: Path, *, compact: bool = False) -> list[dict[str, str]]:
    """구형 raw 좌표와 신형 personal 좌표를 검증하고 동일한 표현으로 읽는다.

    compact의 window_index는 kept 행을 재번호화하므로 전역 index와 비교하지 않는다.
    REF 좌표는 chain 없이 추정하지 않으며, 없는 값은 호출부의 -1 sentinel로 남긴다.
    """
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        required = {"contig", "haplotype", "window_index"}
        if not compact:
            required.add("kept")
        if required - fields:
            raise ValueError(f"windows.tsv 필수 열 누락: {sorted(required-fields)}")
        raw_coords = {"raw_start", "actual_length"} <= fields
        personal_coords = {"personal_start", "personal_end"} <= fields
        if not raw_coords and not personal_coords:
            raise ValueError("windows.tsv에 완전한 raw 또는 personal 좌표 쌍이 필요함")
        rows = []
        for line, raw in enumerate(reader, 2):
            if None in raw or any(value is None for value in raw.values()):
                raise ValueError(f"windows.tsv {line}행: 열 개수 불일치")
            row = dict(raw)
            start = int(row["raw_start"] if raw_coords else row["personal_start"])
            end = start + int(row["actual_length"]) if raw_coords else int(row["personal_end"])
            if personal_coords and (int(row["personal_start"]), int(row["personal_end"])) != (start, end):
                raise ValueError(f"windows.tsv {line}행: raw/personal 좌표 불일치")
            if start < 0 or end <= start or int(row["window_index"]) < 0:
                raise ValueError(f"windows.tsv {line}행: 잘못된 좌표 또는 index")
            if row["haplotype"] not in {"HAP1", "HAP2"} or not row["contig"]:
                raise ValueError(f"windows.tsv {line}행: 잘못된 contig/HAP")
            if not compact and row["kept"].strip().lower() not in {"0", "1", "false", "true"}:
                raise ValueError(f"windows.tsv {line}행: 잘못된 kept 값")
            row["personal_start"], row["personal_end"] = str(start), str(end)
            if compact:
                row["kept"] = "1"
            rows.append(row)
    if not rows:
        raise ValueError("windows.tsv 데이터 행이 없음")
    intervals = defaultdict(list)
    indices = set()
    for row in rows:
        if row["kept"].strip().lower() not in {"1", "true"}:
            continue
        identity = (row["contig"], row["haplotype"], int(row["window_index"]))
        if identity in indices:
            raise ValueError(f"kept window index 중복: {identity}")
        indices.add(identity)
        intervals[identity[:2]].append((int(row["personal_start"]), int(row["personal_end"])))
    for identity, spans in intervals.items():
        spans.sort()
        if any(current[0] < previous[1] for previous, current in zip(spans, spans[1:])):
            raise ValueError(f"kept window 중복 또는 overlap: {identity}")
    return rows


def ensure_window_contract(path: Path) -> bool:
    """재사용 캐시도 검사한다. 구형 열만 원본 백업 후 원자적으로 보충한다."""
    rows = read_windows(path)
    with path.open(newline="") as handle:
        fields = next(csv.reader(handle, delimiter="\t"))
    additions = [name for name in ("personal_start", "personal_end") if name not in fields]
    if not additions:
        return False
    backup = path.with_name(f"{path.name}.legacy.{time.time_ns()}")
    shutil.copy2(path, backup)
    descriptor, name = tempfile.mkstemp(prefix=".windows.", suffix=".part", dir=path.parent)
    part = Path(name)
    try:
        with os.fdopen(descriptor, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields + additions, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({field: row[field] for field in fields + additions})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(part, path)
    finally:
        part.unlink(missing_ok=True)
    return True
