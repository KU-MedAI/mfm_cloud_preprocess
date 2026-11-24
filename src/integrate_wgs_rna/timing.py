from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, field
from time import perf_counter
from typing import Iterator


def format_seconds(seconds: float) -> str:
    total = int(round(seconds))
    return f"{total // 3600:02d}:{(total % 3600) // 60:02d}:{total % 60:02d}"


@dataclass
class StepTimer:
    records: list[dict] = field(default_factory=list)

    @contextmanager
    def step(self, name: str) -> Iterator[None]:
        print(f"[TIME] START {name}", flush=True)
        start = perf_counter()
        status = "ok"
        try:
            yield
        except BaseException:
            status = "failed"
            raise
        finally:
            elapsed = perf_counter() - start
            self.records.append({
                "step": name,
                "elapsed_seconds": elapsed,
                "elapsed_hms": format_seconds(elapsed),
                "status": status,
            })
            print(f"[TIME] END {name} elapsed={elapsed:.3f}s ({format_seconds(elapsed)}) status={status}", flush=True)

