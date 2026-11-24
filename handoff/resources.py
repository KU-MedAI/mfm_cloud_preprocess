"""로컬/Slurm 실행 자원을 과구독 없이 결정한다."""
from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Any
from .host_limits import cgroup_limits


GIB = 1024**3
TIB = 1024**4
PROFILES = {
    "auto": {"base_max_jobs": 4, "elastic_max_jobs": 6,
             "min_cpus_per_sample": 16, "memory_per_sample_gib": 128},
    "throughput": {"base_max_jobs": 6, "elastic_max_jobs": 6,
                   "min_cpus_per_sample": 12, "memory_per_sample_gib": 128},
    "balanced": {"base_max_jobs": 4, "elastic_max_jobs": 4,
                 "min_cpus_per_sample": 16, "memory_per_sample_gib": 128},
    "conservative": {"base_max_jobs": 1, "elastic_max_jobs": 1,
                     "min_cpus_per_sample": 8, "memory_per_sample_gib": 128},
}
ELASTIC_MIN_CPUS = 64
ELASTIC_MIN_MEMORY_GIB = 768
ELASTIC_MIN_FREE_STORAGE_BYTES = 12 * TIB
DEFAULT_MEMORY_LIMIT_BYTES = 1_800_000_000_000


def memory_budget(config: dict, executor: str = 'local') -> int:
    """전체 전처리 RAM 운영 상한. Slurm에서는 실제 할당량을 넘지 않는다."""
    raw = config.get('memory_limit_bytes', DEFAULT_MEMORY_LIMIT_BYTES)
    if isinstance(raw, bool) or not isinstance(raw, int) or raw < 1:
        raise ValueError('memory_limit_bytes는 양의 정수여야 합니다.')
    limit = min(raw, DEFAULT_MEMORY_LIMIT_BYTES)
    if executor == 'slurm' or os.environ.get('SLURM_JOB_ID'):
        if os.environ.get('SLURM_MEM_PER_NODE'):
            limit = min(limit, int(os.environ['SLURM_MEM_PER_NODE']) * 1024 ** 2)
        elif os.environ.get('SLURM_MEM_PER_CPU'):
            limit = min(limit, int(os.environ['SLURM_MEM_PER_CPU']) * 1024 ** 2 * cpu_budget(config, 'slurm'))
    # Cloud VM의 실제 RAM/cgroup 할당은 요청한 운영 예산보다 작을 수 있음.
    from integrate_wgs_rna.memory_policy import allocation_memory_limit
    allocated = allocation_memory_limit()
    return min(limit, allocated) if allocated is not None else limit


def cpu_ceiling(config, affinity=None):
    available = max(1, int(affinity if affinity is not None else detected_cpus()))
    limit = _value(config.get('cpu_pool_limit', 'auto'), 'cpu_pool_limit')
    if limit != 'auto':
        available = min(available, limit)
    _, _, quota = cgroup_limits()
    return min(available, quota) if quota else available


def cpu_budget(config: dict, executor: str = "local", *, affinity: int | None = None) -> int:
    """동시 샘플 수와 별개로 실제 허용 CPU 전체를 공유한다. Slurm 할당은 존중한다."""
    available = cpu_ceiling(config, affinity)
    mode = config.get("cpu_utilization", "maximum")
    if mode not in {"maximum", "bounded"}:
        raise ValueError("cpu_utilization은 maximum|bounded 중 하나여야 합니다.")
    if executor == "slurm" or os.environ.get("SLURM_JOB_ID"):
        raw = os.environ.get("SLURM_CPUS_PER_TASK") or config.get("slurm", {}).get("cpus", available)
        allocated = available if raw == "auto" else int(raw)
        return max(1, min(available, allocated))
    if mode == "maximum" or config.get('cpu_allocation_mode') == 'fixed_sample':
        return available
    return max(1, min(available, int(config["sample_cpus"]) * int(config["local_jobs"])))


def detected_cpus() -> int:
    return max(1, len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity")
               else (os.cpu_count() or 1))


def detected_available_memory_bytes() -> int:
    try:
        with open("/proc/meminfo", encoding="utf-8") as stream:
            values = {parts[0].rstrip(":"): int(parts[1]) * 1024
                      for line in stream if len(parts := line.split()) >= 2}
        return max(1, values["MemAvailable"])
    except (OSError, KeyError, ValueError):
        # 메모리 정보를 읽지 못하면 CPU만으로 계산하되 최소 한 sample은 허용한다.
        return 128 * GIB


def _existing_ancestor(path: Path) -> Path:
    path = path.expanduser().resolve()
    while not path.exists() and path != path.parent:
        path = path.parent
    return path


def detected_storage(config: dict) -> tuple[int, dict[str, dict[str, int | str]]]:
    """data/output/state가 놓일 파일시스템별 최소 가용 공간을 반환한다."""
    report: dict[str, dict[str, int | str]] = {}
    by_device: dict[int, int] = {}
    for name in ("data_root", "output_root", "state_root"):
        anchor = _existing_ancestor(Path(config[name]))
        usage = shutil.disk_usage(anchor)
        device = anchor.stat().st_dev
        by_device[device] = usage.free
        report[name] = {"requested_path": str(config[name]), "probe_path": str(anchor),
                        "device": device, "free_bytes": usage.free}
    return min(by_device.values()), report


def _value(value: Any, name: str) -> int | str:
    if value == "auto" or value is None:
        return "auto"
    if isinstance(value, bool):
        raise ValueError(f"{name}은 auto 또는 양의 정수여야 합니다.")
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name}은 auto 또는 양의 정수여야 합니다.") from exc
    if number < 1:
        raise ValueError(f"{name}은 auto 또는 양의 정수여야 합니다.")
    return number


def resolve(config: dict, overrides: dict | None = None, *, cpu_count: int | None = None,
            available_memory_bytes: int | None = None,
            free_storage_bytes: int | None = None) -> dict:
    """설정과 CLI override를 합쳐 재현 가능한 effective resource config를 만든다."""
    c = dict(config)
    overrides = {k: v for k, v in (overrides or {}).items() if v is not None}
    mode = c['cpu_allocation_mode'] = overrides.get('cpu_allocation_mode', c.get('cpu_allocation_mode','dynamic'))
    if mode not in {'dynamic', 'fixed_sample'}:
        raise ValueError('cpu_allocation_mode은 dynamic|fixed_sample이어야 합니다.')
    declared = c.get('pipeline_env', {}).get('MFM_CPU_ALLOCATION_MODE', mode)
    if declared != mode:
        raise ValueError('cpu_allocation_mode과 pipeline_env의 모드 설정이 충돌합니다.')
    c.setdefault("cpu_utilization", "maximum")
    c.setdefault('memory_limit_bytes', DEFAULT_MEMORY_LIMIT_BYTES)
    c["memory_limit_bytes"] = memory_budget(c)
    if c["cpu_utilization"] not in {"maximum", "bounded"}:
        raise ValueError("cpu_utilization은 maximum|bounded 중 하나여야 합니다.")
    # worker/finish가 effective_config를 다시 읽을 때 자원을 재계산하지 않는다.
    if c.get("resources_resolved") and not overrides:
        return c

    profile = overrides.get("resource_profile", c.get("resource_profile", "auto"))
    if profile not in PROFILES:
        raise ValueError("resource_profile은 auto|throughput|balanced|conservative 중 하나여야 합니다.")
    policy = PROFILES[profile]
    c['cpu_pool_limit'] = _value(overrides.get('cpu_pool_limit', c.get('cpu_pool_limit', 'auto')), 'cpu_pool_limit')
    detected = max(1, int(cpu_count if cpu_count is not None else detected_cpus()))
    cpus = cpu_ceiling(c, detected)
    if os.environ.get('SLURM_JOB_ID'):
        cpus = cpu_budget(c, 'slurm', affinity=cpus)
    memory = max(1, int(available_memory_bytes if available_memory_bytes is not None
                        else detected_available_memory_bytes()))
    memory = min(memory, memory_budget(c))
    storage_report: dict[str, dict[str, int | str]] = {}
    if free_storage_bytes is None:
        storage, storage_report = detected_storage(c)
    else:
        storage = max(0, int(free_storage_bytes))
    memory_cap = max(1, memory // (policy["memory_per_sample_gib"] * GIB))
    cpu_job_cap = max(1, cpus // policy["min_cpus_per_sample"])
    elastic_requested = _value(
        overrides.get("elastic_max_jobs", c.get("elastic_max_jobs", "auto")),
        "elastic_max_jobs")
    elastic_cap = policy["elastic_max_jobs"] if elastic_requested == "auto" else elastic_requested
    elastic_cap = max(1, min(elastic_cap, memory_cap, cpus))
    elastic_eligible = (cpus >= ELASTIC_MIN_CPUS and
                        memory >= ELASTIC_MIN_MEMORY_GIB * GIB and
                        storage >= ELASTIC_MIN_FREE_STORAGE_BYTES)
    if profile == "auto":
        profile_cap = elastic_cap if elastic_eligible else policy["base_max_jobs"]
    else:
        profile_cap = min(elastic_cap, policy["base_max_jobs"])
    elastic_cpu_cap = elastic_cap if (profile == "throughput" or
                                      (profile == "auto" and elastic_eligible)) else cpu_job_cap
    auto_job_cap = max(1, min(profile_cap, memory_cap, elastic_cpu_cap))

    requested_jobs = _value(overrides.get("local_jobs", c.get("local_jobs", "auto")), "local_jobs")
    requested_cpus = _value(overrides.get("sample_cpus", c.get("sample_cpus", "auto")), "sample_cpus")
    if mode == 'fixed_sample' and (requested_cpus == 'auto' or requested_cpus > cpus):
        raise ValueError('fixed_sample은 실제 CPU 한도 이내의 sample_cpus를 명시해야 합니다.')
    adjustments: list[str] = []

    if requested_jobs == "auto":
        if requested_cpus == "auto":
            jobs = auto_job_cap
        else:
            # worker는 CPU 예약이 아니라 stage 대기열이다. 실제 thread 합은 lease pool이 제한한다.
            jobs = max(1, min(profile_cap, memory_cap))
    else:
        jobs = min(requested_jobs, elastic_cap, memory_cap, cpus)
        if jobs != requested_jobs:
            adjustments.append(f"local_jobs {requested_jobs}→{jobs}: 탄력 worker/RAM 상한")

    if requested_cpus == "auto":
        # 6개 worker가 모두 CPU-heavy일 때는 lease가 약 12개씩 배정되지만,
        # 일부가 I/O stage일 때 CPU-heavy stage가 최대 18개까지 확장할 여지를 둔다.
        sample_cpus = max(1, cpus // min(jobs, 4))
    else:
        sample_cpus = min(requested_cpus, cpus)
        if sample_cpus != requested_cpus:
            adjustments.append(f"sample_cpus {requested_cpus}→{sample_cpus}: CPU affinity 상한")

    if mode == 'fixed_sample':
        fixed_jobs = max(1, cpus // sample_cpus)
        if requested_jobs != 'auto' and requested_jobs > fixed_jobs:
            raise ValueError('fixed_sample의 local_jobs × sample_cpus가 CPU 풀을 초과합니다.')
        jobs = min(jobs, fixed_jobs)

    requested_downloads = _value(
        overrides.get("download_jobs", c.get("download_jobs", "auto")), "download_jobs")
    download_jobs = min(4, jobs) if requested_downloads == "auto" else requested_downloads
    forced_elastic = jobs > 4 and not elastic_eligible
    if forced_elastic:
        adjustments.append("탄력 worker를 명시적으로 사용함: CPU/RAM/저장공간 preflight를 확인하세요.")

    slurm = dict(c.get("slurm", {}))
    if "slurm_cpus" in overrides:
        slurm["cpus"] = _value(overrides["slurm_cpus"], "slurm_cpus")
    if "slurm_concurrency" in overrides:
        slurm["concurrency"] = _value(overrides["slurm_concurrency"], "slurm_concurrency")
    for name in ("cpus", "concurrency"):
        slurm[name] = _value(slurm.get(name), f"slurm.{name}")
        if slurm[name] == "auto":
            raise ValueError(f"slurm.{name}은 양의 정수여야 합니다.")

    c.update(resource_profile=profile, local_jobs=jobs, sample_cpus=sample_cpus,
             download_jobs=download_jobs, elastic_max_jobs=elastic_cap,
             slurm=slurm, resources_resolved=True,
             resource_resolution={
                 "detected_cpus": detected,
                 "effective_cpu_ceiling": cpus,
                 "cpu_pool_limit": c["cpu_pool_limit"],
                 "effective_memory_limit_bytes": memory_budget(c),
                 "detected_available_memory_gib": round(memory / GIB, 2),
                 "detected_min_free_storage_bytes": storage,
                 "detected_min_free_storage_tib": round(storage / TIB, 3),
                 "storage_probes": storage_report,
                 "memory_per_sample_gib": policy["memory_per_sample_gib"],
                 "base_max_jobs": policy["base_max_jobs"],
                 "elastic_max_jobs": elastic_cap,
                 "elastic_eligible": elastic_eligible,
                 "elastic_enabled": jobs > 4,
                 "elastic_reason": ("CPU/RAM/storage 기준 충족" if jobs > 4 and elastic_eligible else
                                    "throughput/명시 override" if jobs > 4 else
                                    "안정 기본 또는 CPU/RAM/storage 기준 미충족"),
                 "requested_local_jobs": requested_jobs,
                 "requested_sample_cpus": requested_cpus,
                 "requested_download_jobs": requested_downloads,
                 "requested_elastic_max_jobs": elastic_requested,
                 "effective_local_jobs": jobs,
                 "effective_sample_cpus": sample_cpus,
                 "effective_total_cpu_budget": cpu_budget({**c, "local_jobs": jobs, "sample_cpus": sample_cpus, "slurm": slurm}, affinity=cpus),
                 "cpu_utilization": c["cpu_utilization"],
                 "worker_cpu_overcommit_ratio": round((jobs * sample_cpus) / cpus, 3),
                 "effective_download_jobs": download_jobs,
                 "adjustments": adjustments,
             })
    from integrate_wgs_rna.processing_buffers import apply
    return apply(c)
