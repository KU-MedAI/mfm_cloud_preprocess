#!/usr/bin/env python3
"""RNA depth/CNV 단계에서 공유하는 0-based half-open 구간 helper."""
from __future__ import annotations


def merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """겹치거나 맞닿은 구간을 병합해 중복 exon 길이의 이중 계수를 막는다."""
    out: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if out and start <= out[-1][1]:
            if end > out[-1][1]:
                out[-1] = (out[-1][0], end)
        else:
            out.append((start, end))
    return out
