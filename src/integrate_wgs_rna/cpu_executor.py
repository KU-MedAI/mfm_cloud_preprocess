"""독립 작업 수에 맞춰 공유 CPU 임대를 반환하고 여유 CPU를 재사용한다."""
from __future__ import annotations

import fcntl
import json
import os
import uuid
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from contextlib import ExitStack
from multiprocessing import get_context
from pathlib import Path
from .pool_leases import prune_leases, fitting_waiters, regular_cpu_claim, pool_cpu_limit, read_pool_records, live_waiters
from .memory_policy import GIB, limit_worker_address_space, memory_active, usable_memory, memory_limit, available_memory
from .cpu_loan_policy import available_loan_cpus, available_regular_cpus
from .resource_fairness import fair_cpu_capacity


def elastic_cpu_enabled(request_name: str, *, thread_cost_accounted: bool = False) -> bool:
    thread_name = ("MFM_TRACE_BLOSC_THREADS" if request_name == "MFM_TRACE_JOBS_REQUESTED"
                   else "MFM_BLOSC_THREADS")
    if request_name=="MFM_COMPACT_JOBS_REQUESTED" and any(
            os.environ.get('MFM_'+stage+'_BACKEND')=='cupy' for stage in ['ZARR','QC','COVERAGE']):
        thread_cost_accounted=False
    return (os.environ.get("MFM_CPU_UTILIZATION", "maximum") == "maximum"
            and os.environ.get(request_name, "auto") == "auto"
            and (thread_cost_accounted or os.environ.get(thread_name, "1") == "1")
            and bool(os.environ.get("MFM_CPU_POOL_LEASE_FILE")))


class WorkerBudget:
    """한 CPU를 사용하는 작업들의 실제 동시 실행 상한을 관리한다."""

    def __init__(self, workers: int, limit: int, elastic: bool, memory_per_worker: int = 0,
                 cpus_per_worker: int = 1, parent_cpus: int = 0, *, wait_memory_bytes: int | None = None):
        self.initial = max(1, workers)
        self.limit = max(1, limit)
        value = os.environ.get("MFM_CPU_POOL_LEASE_FILE")
        root = Path(value).parent.parent if value else None
        self.memory_enabled = memory_active(root)
        self.path = Path(value) if value and (elastic or self.memory_enabled) else None
        self.elastic = elastic
        self.memory_per_worker = memory_per_worker if self.memory_enabled else 0
        self.base_memory = int(os.environ.get('MFM_MEMORY_STAGE_BASE_BYTES', 4 * GIB))
        self.cpu_cost = max(1, cpus_per_worker)
        self.parent_cpus = max(0, int(parent_cpus))
        # 순수 예측 대기에서만 명시함. 보존 중인 버퍼는 호출자가 이 최소량에 포함함.
        if wait_memory_bytes is not None and wait_memory_bytes < 0:
            raise ValueError('대기 RAM 최소량은 음수가 될 수 없음')
        self.wait_memory_bytes = wait_memory_bytes
        self.resident_workers = 0
        self.last_write = 0

    def capacity(self, remaining: int, grow: bool = True, *, running: int | None = None) -> int:
        idle_count_known = running is not None
        running = max(0,int(running or 0))
        target = max(1, min(self.limit, remaining))
        if not self.elastic:
            target = min(target, self.initial)
        if self.path is None:
            if self.memory_per_worker:
                cap = (memory_limit() * 9 // 10 - self.base_memory) // self.memory_per_worker
                if cap < 1:
                    raise ValueError("계산 작업 한 개의 RAM 예상량이 예산을 초과합니다.")
                target = min(target, cap)
            return min(self.initial, target)
        root = self.path.parent.parent
        with (root / "pool.lock").open("a+") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            prune_leases(root)
            row = json.loads(self.path.read_text())
            before_claim=(row.get('cpus'),row.get('memory_bytes'))
            previous_allocation=row.get('worker_allocation',{})
            total = pool_cpu_limit(root,int(row['total_cpus']))
            policy_path=root/'policy.json'
            setting=json.loads(policy_path.read_text()).get('cpu_backfill',{}) if policy_path.exists() else {}
            capacity = min(max(0,int(row["cpus"])-self.parent_cpus) // self.cpu_cost, target)
            other_leases = [value for path, value in read_pool_records(root / 'leases') if path != self.path]
            other_memory = sum(int(value.get('memory_bytes', 0)) for value in other_leases)
            available_ram = available_memory(root,other_leases)
            memory_cap = target
            base_fits = not self.memory_enabled or self.base_memory <= available_ram
            if not base_fits:
                target = capacity = 0
            if self.memory_per_worker:
                memory_cap = max(0, (available_ram - self.base_memory) // self.memory_per_worker)
                target = min(target, memory_cap)
                capacity = min(capacity, target)
            if grow:
                leases = other_leases
                free = available_regular_cpus(total,setting,[*leases,row])
                free_memory = max(0, available_ram - int(row.get('memory_bytes', 0)))
                waiters = sorted(live_waiters(root / "waiters"), key=lambda item: str(item[0]))
                regular_room = fair_cpu_capacity(root,total,[*leases,row],
                    available_regular_cpus(total,setting,leases),request=row,waiters=waiters)
                if idle_count_known and not row.get('cpu_backfill'):
                    # 새 묶음에만 적용하며 이미 실행 중인 worker는 아래 running 하한으로 보존함.
                    capacity = min(capacity,max(0,regular_room-self.parent_cpus)//self.cpu_cost)
                reacquire = int(row.get('cpu_backfill', {}).get('phase_reacquire_cpus', 0))
                if reacquire and row.get('branch')=='WGS' and row.get('stage') in {'02','06'}:
                    # 검증 묶음 사이에 반환한 대여만 재확보한다. 정상 풀과 RAM을 침범하지 않는다.
                    room = available_loan_cpus(total,setting,leases,requested=reacquire)
                    phase_flag='haplotag_phases' if row.get('stage')=='06' else 'validation_phases'
                    if (setting.get(phase_flag) is not True
                            or os.environ.get('MFM_CPU_UTILIZATION','maximum')!='maximum'):room=0
                    capacity = min(target, max(capacity, max(0,room-self.parent_cpus)//self.cpu_cost))
                # 현재 RAM/최소 CPU로 입장할 수 없는 요청은 실제 검증 작업을 막지 않음.
                elif not row.get('cpu_backfill') and not fitting_waiters(waiters, free, free_memory):
                    capacity = min(target, max(capacity, max(0, regular_room-self.parent_cpus) // self.cpu_cost))
            if remaining == 0 and self.memory_enabled:
                capacity = 0
            # 이미 계산 중인 worker의 CPU/RAM은 join/완료 전에 다른 작업에 반환하지 않음.
            capacity=max(capacity,min(self.limit,max(0,int(running))))
            row["cpus"] = 1 if remaining == 0 else (self.parent_cpus + capacity * self.cpu_cost if capacity else 0)
            if self.memory_enabled:
                row['memory_bytes'] = (self.base_memory + capacity * self.memory_per_worker if base_fits
                                       else min(self.base_memory, int(row.get('memory_bytes', self.base_memory))))
            if self.memory_enabled and self.wait_memory_bytes is not None and capacity == 0 and not running:
                rss = int(Path('/proc/self/statm').read_text().split()[1]) * os.sysconf('SC_PAGE_SIZE')
                row['memory_bytes'] = max(self.wait_memory_bytes, (rss * 11 + 9) // 10)
            if self.memory_enabled and (running or self.resident_workers):
                row['memory_bytes']=max(row['memory_bytes'],self.base_memory+max(running,self.resident_workers)*self.memory_per_worker)
            available_cpu = max(0, available_regular_cpus(total,setting,other_leases) - self.parent_cpus)
            limits = []
            if self.memory_enabled and (not base_fits or memory_cap < min(self.limit, remaining)):
                limits.append('RAM 예약')
            if available_cpu // self.cpu_cost < min(self.limit, remaining):
                limits.append('CPU 예약')
            if not limits and capacity < min(self.limit, remaining):
                limits.append('명시 worker 상한 또는 대기 작업 배분')
            previous_wait = row.get('worker_allocation', {}).get('queued_at')
            queued_at = (previous_wait if previous_wait is not None else time.time()) if remaining and capacity==0 else None
            row['worker_allocation'] = {
                'queued_at': queued_at, 'minimum_cpus': self.parent_cpus+self.cpu_cost,
                'minimum_memory_bytes': self.base_memory+self.memory_per_worker,
                'requested_workers': min(self.limit, remaining), 'granted_workers': capacity,
                'running_workers': running, 'cpus_per_worker': self.cpu_cost,
                'memory_per_worker_bytes': self.memory_per_worker,
                'other_memory_bytes': other_memory, 'available_cpus': available_cpu,
                'available_memory_bytes': max(0, available_ram - self.base_memory),
                'limiting_factors': limits,
            }
            row["peak_cpus"] = max(int(row.get("peak_cpus", 0)), row["cpus"])
            important=('requested_workers','granted_workers','running_workers','queued_at','minimum_cpus','minimum_memory_bytes')
            changed=before_claim!=(row.get('cpus'),row.get('memory_bytes')) or any(
                previous_allocation.get(key)!=row['worker_allocation'].get(key) for key in important)
            if changed or time.monotonic()-self.last_write>=1:
                candidate = self.path.with_name(self.path.name + "." + uuid.uuid4().hex + ".tmp")
                candidate.write_text(json.dumps(row) + "\n")
                os.replace(candidate, self.path)
                self.last_write=time.monotonic()
            return capacity


def cpu_results(function, tasks: list[tuple], workers: int, *, elastic: bool = False,
                mp_context=None, max_workers: int | None = None, memory_per_worker_bytes: int = 0,
                cpus_per_worker: int = 1, max_tasks_per_child: int | None = None,
                retain_worker_memory: bool = False, parent_cpus: int = 0,
                address_space_per_worker_bytes: int | None = None,
                base_memory_bytes: int | None = None, wait_memory_bytes: int | None = None):
    """결과 또는 예외를 원래 작업 번호와 반환한다. 미실행 작업은 대기열에 둔다."""
    if memory_per_worker_bytes < 0 or (address_space_per_worker_bytes is not None and
            address_space_per_worker_bytes < memory_per_worker_bytes):
        raise ValueError('worker 주소공간 한도는 RAM 예약 이상이어야 합니다.')
    if not tasks:
        return
    address_space = memory_per_worker_bytes if address_space_per_worker_bytes is None else address_space_per_worker_bytes
    limit = min(len(tasks), max_workers or (len(tasks) if elastic else max(1, workers)))
    budget = WorkerBudget(workers, limit, elastic, 0 if retain_worker_memory else memory_per_worker_bytes, cpus_per_worker, parent_cpus, wait_memory_bytes=wait_memory_bytes)
    if base_memory_bytes is not None:
        if base_memory_bytes < 0:raise ValueError('부모 RAM 예약은 음수일 수 없습니다.')
        budget.base_memory=max(budget.base_memory,base_memory_bytes)
    base_memory=budget.base_memory
    if retain_worker_memory and budget.memory_enabled:
        # fork된 유휴 worker도 공유 context/allocator RSS를 보유하므로 join까지 RAM을 유지함.
        budget.base_memory+=limit*memory_per_worker_bytes
    if limit == 1 and not elastic and not budget.memory_enabled and max_tasks_per_child is None and address_space_per_worker_bytes is None:
        for index, args in enumerate(tasks):
            try:
                yield index, function(*args), None
            except Exception as exc:
                yield index, None, exc
        return
    pending = iter(enumerate(tasks))
    submitted = completed = 0
    futures = {}
    try:
        with ExitStack() as stack:
            def create_executor():
                # spawn은 실제 제출 수만큼 자식을 생성하며 fork의 전체 pool 선생성을 피함.
                return stack.enter_context(ProcessPoolExecutor(max_workers=limit,
                    mp_context=mp_context or get_context('spawn'), initializer=limit_worker_address_space,
                    initargs=(address_space if budget.memory_enabled or address_space_per_worker_bytes is not None else 0,),
                    max_tasks_per_child=max_tasks_per_child))
            executor=create_executor()
            while completed < len(tasks):
                budget.resident_workers=0 if retain_worker_memory else len(getattr(executor,'_processes',None) or {})
                capacity = budget.capacity(len(tasks) - completed,running=len(futures))
                while submitted < len(tasks) and len(futures) < capacity:
                    index, args = next(pending)
                    futures[executor.submit(function, *args)] = index
                    submitted += 1
                if not futures:
                    if budget.resident_workers:
                        # 모든 작업 완료가 확인된 유휴 자식만 join한 뒤 RAM을 반환함.
                        executor.shutdown(wait=True)
                        budget.resident_workers=0
                        budget.capacity(len(tasks)-completed,grow=False)
                        executor=create_executor()
                    time.sleep(0.5)
                    continue
                done, _ = wait(futures, timeout=0.2, return_when=FIRST_COMPLETED)
                for future in done:
                    index = futures.pop(future)
                    completed += 1
                    budget.resident_workers=0 if retain_worker_memory else len(getattr(executor,'_processes',None) or {})
                    budget.capacity(len(tasks) - completed, grow=False,running=len(futures))
                    try:
                        value, error = future.result(), None
                    except Exception as exc:
                        value, error = None, exc
                    yield index, value, error
    finally:
        # 최종 metadata 기록·다음 단계 전환에 필요한 한 CPU만 유지한다.
        budget.base_memory=base_memory
        budget.resident_workers=0
        budget.capacity(0, grow=False)
