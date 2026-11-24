"""RNA worker 수를 실제 CPU affinity와 샘플 예산 안에서 결정한다."""
from __future__ import annotations

import os
from pathlib import Path


def worker_count(requested: str, tasks: int, cpu_cost: int = 1) -> int:
    """조정 프로세스 1코어를 남기며 중첩 병렬화의 CPU 비용도 계산한다."""
    if tasks < 1 or cpu_cost < 1:
        raise ValueError("tasks와 cpu_cost는 양수여야 합니다")
    available = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    from integrate_wgs_rna.gpu_compare_resources import cpu_limit
    available=cpu_limit(available)
    for name in ("MFM_SAMPLE_CPU_BUDGET", "MFM_CPU_BUDGET"):
        if os.environ.get(name):
            budget = int(os.environ[name])
            if budget < 1:
                raise ValueError(f"{name}은 양수여야 합니다")
            available = min(available, budget)
    cap = max(1, (available - 1) // cpu_cost)
    chosen = cap if requested == "auto" else int(requested)
    if chosen < 1:
        raise ValueError("worker 수는 auto 또는 양수여야 합니다")
    return min(chosen, cap, tasks)


def file_signature(path: str | Path) -> dict:
    """재사용 직전 파일 동일성 검사. 내용 해시가 아닌 실행 중 변경 감지용이다."""
    resolved = Path(path).resolve()
    st = resolved.stat()
    return {"path": str(resolved), "size": st.st_size,
            "mtime_ns": st.st_mtime_ns, "inode": st.st_ino}
