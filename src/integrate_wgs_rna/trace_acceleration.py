"""독립 source trace 전용 가속: writer의 projection/depth 구현과 공유하지 않는다."""

from __future__ import annotations

from bisect import bisect_right
from typing import Any

import numpy as np


class WindowLookup:
    """후보만 줄이고, 실제 segment 좌표 결정은 기존 독립 CIGAR projector에 맡긴다."""

    def __init__(self, windows: list[tuple[int, int, int]]) -> None:
        self.starts = [start for start, _actual, _stored in windows]
        self.ends = [start + actual for start, actual, _stored in windows]
        # 손상된 window metadata의 기존 fail/audit 경로를 최적화가 가리지 않게 한다.
        self.searchable = (
            all(actual > 0 for _start, actual, _stored in windows)
            and self.starts == sorted(self.starts)
            and self.ends == sorted(self.ends)
        )

    def candidates(self, read: Any) -> range:
        # 기존 전수 projection은 CIGAR를 만나면 반드시 검증했다. 후보가 0개여도
        # unsupported/0-length op를 조용히 무시하지 않는다.
        for op, length in read.cigartuples or ():
            if int(length) <= 0:
                raise ValueError("invalid source CIGAR length")
            if int(op) not in (0, 1, 2, 3, 4, 5, 6, 7, 8):
                raise ValueError(f"unsupported source CIGAR op: {int(op)}")
        if not self.searchable:
            return range(len(self.starts))
        start = int(read.reference_start)
        end = int(read.reference_end or start)
        # trailing insertion의 anchor는 reference_end 자체다. half-open overlap만
        # 쓰면 정확히 다음 window 시작점에 있는 insertion을 누락할 수 있다.
        return range(bisect_right(self.ends, start), bisect_right(self.starts, end))


class DepthAccumulator:
    """구간 endpoints를 더한 뒤 window별 누적합. uint32 기존 wrap 의미를 유지한다."""

    def __init__(self, windows: list[tuple[int, int, int]]) -> None:
        self.buffers = [np.zeros(stored + 1, dtype=np.uint32) for _start, _actual, stored in windows]
        self.pending: dict[int, list[tuple[int, int]]] = {}
        self.pending_count = 0
        self.active: set[int] = set()

    def add(self, index: int, start: int, end: int) -> None:
        self.pending.setdefault(index, []).append((start, end))
        self.active.add(index)
        self.pending_count += 1
        if self.pending_count >= 8192:
            for window_index in list(self.pending):
                self._flush(window_index)

    def _flush(self, index: int) -> None:
        rows = self.pending.pop(index, [])
        if rows:
            positions = np.asarray(rows, dtype=np.intp)
            # repeated endpoints가 있어 fancy indexing += 는 잘못된 결과를 낸다.
            # at은 같은 위치의 모든 증감을 보존하고 uint32 wrap도 기존과 같다.
            np.add.at(self.buffers[index], positions[:, 0], np.uint32(1))
            np.subtract.at(self.buffers[index], positions[:, 1], np.uint32(1))
            self.pending_count -= len(rows)

    def finish_window(self, index: int) -> np.ndarray:
        # 한 window를 한 번만 finalize한다. 추가 chr-size 복사/메모리 증가는 없다.
        self._flush(index)
        values = self.buffers[index][:-1]
        # supplementary/sparse 입력에서 read가 없는 window를 쓰기 방문하면
        # zero-page도 물리 메모리에 올라가고 무의미한 누적합 비용이 생긴다.
        if index in self.active:
            np.cumsum(values, dtype=np.uint32, out=values)
        return values
