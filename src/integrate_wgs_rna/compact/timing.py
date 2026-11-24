from __future__ import annotations

import time
from collections import Counter
from typing import Any

import numpy as np


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(v) for v in value]
    if isinstance(value, Counter):
        return {str(k): int(v) for k, v in value.items()}
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


class DirectTiming:
    def __init__(self) -> None:
        self.steps: dict[str, float] = {}
        self.chunks: list[dict[str, Any]] = []
        self.counts: dict[str, int] = {}

    def add(self, name: str, elapsed: float) -> None:
        self.steps[name] = self.steps.get(name, 0.0) + float(elapsed)

    def add_count(self, name: str, value: int) -> None:
        self.counts[name] = self.counts.get(name, 0) + int(value)

    def step(self, name: str) -> "_DirectTimingStep":
        return _DirectTimingStep(self, name)

    def add_chunk(self, chunk: dict[str, Any]) -> None:
        self.chunks.append(json_ready(chunk))

    def to_dict(self) -> dict[str, Any]:
        return {
            "steps": {
                key: round(value, 6)
                for key, value in sorted(self.steps.items())
            },
            "counts": dict(sorted(self.counts.items())),
            "chunks": self.chunks,
        }


class _DirectTimingStep:
    def __init__(self, timing: DirectTiming, name: str) -> None:
        self.timing = timing
        self.name = name
        self.started = 0.0

    def __enter__(self) -> "_DirectTimingStep":
        self.started = time.perf_counter()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.timing.add(self.name, time.perf_counter() - self.started)


def timing_add(
    timing: dict[str, Any] | None, name: str, elapsed: float
) -> None:
    if timing is None:
        return
    steps = timing.setdefault("steps", {})
    steps[name] = float(steps.get(name, 0.0)) + float(elapsed)


def timing_set_count(
    timing: dict[str, Any] | None, name: str, value: int
) -> None:
    if timing is None:
        return
    counts = timing.setdefault("counts", {})
    counts[name] = int(value)


def timing_add_count(
    timing: dict[str, Any] | None, name: str, value: int
) -> None:
    if timing is None:
        return
    counts = timing.setdefault("counts", {})
    counts[name] = int(counts.get(name, 0)) + int(value)


def merge_timing_dict(target: DirectTiming, source: dict[str, Any]) -> None:
    for name, elapsed in source.get("steps", {}).items():
        target.add(str(name), float(elapsed))
    for name, value in source.get("counts", {}).items():
        target.add_count(str(name), int(value))
