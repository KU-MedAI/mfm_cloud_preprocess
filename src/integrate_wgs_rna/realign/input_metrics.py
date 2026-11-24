"""QNAME 경계를 보존하는 WGS 입력 통계 병렬 집계. SAM은 파일로 저장하지 않음."""
from __future__ import annotations

from collections import deque
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import io
import json
import math
import multiprocessing
import os
from pathlib import Path
import resource
import subprocess
import threading
import time

from ..memory_policy import limit_worker_address_space, memory_limit

MIB = 1024 ** 2
GIB = 1024 ** 3
DEFAULT_CHUNK_BYTES = 8 * MIB
MAX_GROUP_BYTES = 64 * MIB
METRIC_KEYS = (
    "paired_fragments", "alignment_input_primary", "alignment_input_routed_primary",
    "alignment_input_rescue_primary", "alignment_input_hp1", "alignment_input_hp2",
    "alignment_input_untagged", "alignment_input_duplicate",
    "singleton_discard_records", "read_other_discard_records",
)


def metric_plan(path: Path, decode_threads: int, *, workers: int | None = None,
                chunk_bytes: int = DEFAULT_CHUNK_BYTES, concurrency: int = 1) -> dict:
    """실제 grant와 RAM 예약에서 최대 worker 수·bounded buffer 피크를 계산함."""
    if chunk_bytes < 1 or concurrency < 1 or (workers is not None and workers < 1):
        raise ValueError("집계 chunk/workers/concurrency는 양수여야 함")
    affinity = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count() or 1
    lease = {}
    if os.environ.get("MFM_CPU_POOL_LEASE_FILE"):
        lease = json.loads(Path(os.environ["MFM_CPU_POOL_LEASE_FILE"]).read_text())
    cpu_caps = [affinity]
    cpu_keys = ("MFM_CPU_BUDGET", "SLURM_CPUS_PER_TASK")
    if lease.get("cpus"):
        cpu_caps.append(int(lease["cpus"]))
    else:
        cpu_keys += ("MFM_CPU_POOL_GRANTED", "MFM_SAMPLE_CPU_BUDGET")
    for key in cpu_keys:
        if os.environ.get(key):
            cpu_caps.append(int(os.environ[key]))
    cpu = max(1, min(cpu_caps) // concurrency)
    ram = memory_limit()
    if int(os.environ.get("SLURM_MEM_PER_NODE", "0")) > 0:
        ram = min(ram, int(os.environ["SLURM_MEM_PER_NODE"]) * MIB)
    if int(os.environ.get("SLURM_MEM_PER_CPU", "0")) > 0:
        ram = min(ram, int(os.environ["SLURM_MEM_PER_CPU"]) * MIB * min(cpu_caps))
    ram = ram * 9 // 10
    if lease.get("memory_bytes"):
        ram = min(ram, int(lease["memory_bytes"]))
    ram //= concurrency
    aux = max(0, min(int(decode_threads), cpu - 4))
    cpu_workers = max(1, cpu - aux - 2)
    max_payload = chunk_bytes + MAX_GROUP_BYTES
    fixed = 512 * MIB + 3 * max_payload
    per_worker = 512 * MIB + 6 * max_payload  # worker 해석4배와 대기2건을 예약함.
    ram_workers = max(0, (ram * 10 // 11 - fixed) // per_worker)
    input_bytes = path.stat().st_size
    tasks = max(1, math.ceil(input_bytes / chunk_bytes))
    chosen = min(workers if workers is not None else tasks, cpu_workers, ram_workers)
    if chosen < 2 or cpu < 3:
        chosen, aux = 1, max(0, min(int(decode_threads), cpu - 1))
    predicted = (fixed + chosen * per_worker) * 11 // 10
    if chosen == 1:
        predicted = 512 * MIB
    if predicted > ram:
        raise ValueError("입력 집계의 예측 RAM이 예약 예산보다 큼")
    return {
        "mode": "parallel_sam_stream" if chosen > 1 else "serial_pysam",
        "input_bytes": input_bytes, "cpu_budget": cpu, "workers_requested": chosen,
        "samtools_aux_threads": aux, "chunk_bytes": chunk_bytes,
        "max_group_bytes": MAX_GROUP_BYTES, "max_payload_bytes": max_payload,
        "max_inflight": 2 * chosen, "ram_budget_bytes": ram,
        "estimated_peak_rss_bytes": predicted,
        "worker_address_space_bytes": min(ram, max(4 * GIB, max_payload * 8 + GIB)),
    }


def complete_qname_batches(stream, chunk_bytes: int, max_group_bytes: int = MAX_GROUP_BYTES):
    """마지막 QNAME 묶음을 보류하여 primary pair 경계를 보존함."""
    pending = b""
    while True:
        data = stream.read(chunk_bytes)
        if not data:
            if pending:
                yield pending
            return
        pending += data
        newline = pending.rfind(b"\n")
        if newline >= 0:
            end = newline + 1
            start = pending.rfind(b"\n", 0, newline) + 1
            tab = pending.find(b"\t", start, end)
            if tab < 0:
                raise ValueError("SAM QNAME 열이 없음")
            name = pending[start:tab]
            while start:
                previous = pending.rfind(b"\n", 0, start - 1) + 1
                previous_tab = pending.find(b"\t", previous, start)
                if previous_tab < 0 or pending[previous:previous_tab] != name:
                    break
                start = previous
            if start:
                yield pending[:start]
                pending = pending[start:]
        if len(pending) > max_group_bytes:
            raise ValueError("단일 QNAME/미완성 SAM record가 집계 buffer 상한을 초과함")


def _integer_tag(tags: dict[bytes, tuple[bytes, bytes]], name: bytes) -> int:
    kind, value = tags[name]
    return int(float(value)) if kind in {b"f", b"d"} else int(value)


def count_sam_batch(data: bytes) -> tuple[dict, int, float, int]:
    """기존 첫 R1/첫 R2 선택, HP/XU/OD·singleton 의미를 그대로 집계함."""
    started = time.process_time()
    metrics = dict.fromkeys(METRIC_KEYS, 0)
    current = None
    first1 = first2 = None
    n1 = n2 = other = 0

    def finish():
        if first1 is not None and first2 is not None:
            metrics["paired_fragments"] += 1
            metrics["alignment_input_primary"] += 2
            for selected in (first1, first2):
                tags = {}
                for token in selected.split(b"\t", 11)[11:]:
                    for tag in token.rstrip(b"\r\n").split(b"\t"):
                        parts = tag.split(b":", 2)
                        if len(parts) == 3 and parts[0] in {b"HP", b"XU", b"OD"}:
                            tags.setdefault(parts[0], (parts[1], parts[2]))
                if b"XU" in tags and _integer_tag(tags, b"XU") == 1:
                    metrics["alignment_input_rescue_primary"] += 1
                else:
                    metrics["alignment_input_routed_primary"] += 1
                    try:
                        hp = _integer_tag(tags, b"HP") if b"HP" in tags else 0
                    except (ValueError, TypeError):
                        hp = -1
                    key = f"alignment_input_hp{hp}" if hp in {1, 2} else "alignment_input_untagged"
                    metrics[key] += 1
                metrics["alignment_input_duplicate"] += int(b"OD" in tags and _integer_tag(tags, b"OD") == 1)
            metrics["singleton_discard_records"] += n1 + n2 - 2
        else:
            metrics["singleton_discard_records"] += n1 + n2
        metrics["read_other_discard_records"] += other

    for line in io.BytesIO(data):
        fields = line.split(b"\t", 2)
        if len(fields) < 3 or not fields[0] or fields[0] == b"*":
            raise ValueError("입력 집계 SAM record/QNAME가 잘못됨")
        flag = int(fields[1])
        if flag & (0x100 | 0x800):
            continue
        if current is not None and fields[0] != current:
            finish()
            first1 = first2 = None
            n1 = n2 = other = 0
        current = fields[0]
        if flag & 0x40 and not flag & 0x80:
            n1 += 1
            if first1 is None:
                first1 = line
        elif flag & 0x80 and not flag & 0x40:
            n2 += 1
            if first2 is None:
                first2 = line
        else:
            other += 1
    if current is not None:
        finish()
    return metrics, os.getpid(), time.process_time() - started, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024


def parallel_input_metrics(path: Path, samtools_bin: str, plan: dict) -> tuple[dict, dict]:
    """bounded queue를 이용해 합산함. 모든 worker와 samtools 성공 후에만 반환함."""
    command = [samtools_bin, "view", "-F", "2304", "-@", str(plan["samtools_aux_threads"]), str(path)]
    stderr_tail = deque(maxlen=32)
    started = time.monotonic()
    totals = dict.fromkeys(METRIC_KEYS, 0)
    worker_rss = {}
    worker_cpu = 0.0
    batches = 0
    peak_inflight = 0
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def drain_stderr():
        for line in iter(lambda: process.stderr.read(2048), b""):
            stderr_tail.append(line)

    reader = threading.Thread(target=drain_stderr, daemon=True)
    reader.start()
    pending = set()
    executor = ProcessPoolExecutor(
        max_workers=plan["workers_requested"], mp_context=multiprocessing.get_context("spawn"),
        initializer=limit_worker_address_space, initargs=(plan["worker_address_space_bytes"],),
    )

    def collect(completed):
        nonlocal worker_cpu, batches
        for future in completed:
            metrics, pid, cpu, rss = future.result()
            for key in METRIC_KEYS:
                totals[key] += metrics[key]
            worker_rss[pid] = max(worker_rss.get(pid, 0), rss)
            worker_cpu += cpu
            batches += 1

    try:
        for batch in complete_qname_batches(process.stdout, plan["chunk_bytes"], plan["max_group_bytes"]):
            pending.add(executor.submit(count_sam_batch, batch))
            peak_inflight = max(peak_inflight, len(pending))
            if len(pending) >= plan["max_inflight"]:
                completed, pending = wait(pending, return_when=FIRST_COMPLETED)
                collect(completed)
        while pending:
            completed, pending = wait(pending, return_when=FIRST_COMPLETED)
            collect(completed)
        code = process.wait()
        reader.join()
        if code:
            raise RuntimeError(f"samtools 입력 집계 실패({code}): " + b"".join(stderr_tail).decode(errors="replace"))
    finally:
        process.stdout.close()
        if process.poll() is None:
            process.terminate()
        process.wait()
        reader.join()
        process.stderr.close()
        for future in pending:
            future.cancel()
        executor.shutdown(wait=True, cancel_futures=True)
    execution = {
        **plan, "workers_used": len(worker_rss), "batches": batches,
        "peak_inflight": peak_inflight, "wall_seconds": round(time.monotonic() - started, 3),
        "worker_cpu_seconds": round(worker_cpu, 3),
        "worker_peak_rss_sum_bytes": sum(worker_rss.values()),
    }
    return totals, execution
