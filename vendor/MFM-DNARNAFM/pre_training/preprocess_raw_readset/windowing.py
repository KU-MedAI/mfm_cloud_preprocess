from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass

@dataclass(frozen=True)
class TrimmedSequence:
    sequence: str
    source_start: int
    source_end: int


@dataclass(frozen=True)
class SequenceWindow:
    index: int
    raw_start: int
    raw_end: int
    actual_length: int
    sequence: str


def trim_terminal_n(sequence: str) -> TrimmedSequence:
    seq = sequence.upper()
    left = len(seq) - len(seq.lstrip("N"))
    right = len(seq.rstrip("N"))
    return TrimmedSequence(seq[left:right], left, right)


def deterministic_offset(seed: int, sample: str, contig: str, unit: int) -> int:
    if unit <= 0:
        raise ValueError("offset unit must be positive")
    digest = hashlib.sha256(f"{seed}:{sample}:{contig}".encode()).digest()
    return int.from_bytes(digest[:8], "big") % unit


def make_windows(
    sequence: str, *, source_start: int, offset: int, window_size: int,
    max_n_ratio: float = 1.0, pad: str = "P", max_windows: int | None = None,
) -> list[SequenceWindow]:
    windows: list[SequenceWindow] = []
    start = offset
    while start < len(sequence):
        actual = sequence[start : start + window_size]
        if actual and any(base in "ACGT" for base in actual) and actual.count("N") / len(actual) <= max_n_ratio:
            windows.append(SequenceWindow(
                index=len(windows), raw_start=source_start + start,
                raw_end=source_start + start + len(actual), actual_length=len(actual),
                sequence=actual + pad * (window_size - len(actual)),
            ))
            if max_windows is not None and len(windows) >= max_windows:
                break
        start += window_size
    return windows


def sequence_annotations(sequence: str) -> dict:
    """패딩을 뺀 실제 구간의 N 원값만 계산함."""
    seq = sequence.upper()
    n_count = seq.count("N")
    non_acgt = len(seq) - sum(seq.count(base) for base in "ACGT")
    return {
        "n_count": n_count, "n_nonacgt": non_acgt,
        "n_fraction": non_acgt / len(seq) if seq else 1.0,
        "n_max_run": max((m.end() - m.start() for m in re.finditer("N+", seq)), default=0) if n_count else 0,
        "storage_policy": "annotated_acgtn_v1",
    }
