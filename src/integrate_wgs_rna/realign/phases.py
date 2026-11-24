"""재정렬 하위 단계가 사용하는 CPU/RAM만 동일 lease에 예약함."""
from __future__ import annotations

import fcntl
import json
import os
import tempfile
import time
from dataclasses import replace
from pathlib import Path

from ..memory_policy import GIB, memory_limit, usable_memory, available_memory
from ..pool_leases import prune_leases, fitting_waiters, pool_cpu_limit, regular_cpu_claim, read_pool_records
from ..cpu_loan_policy import available_regular_cpus
from ..resource_fairness import fair_cpu_capacity
from .align_resources import alignment_plan, parse_memory, alignment_cpu_overhead
from .align_memory import alignment_memory
from .idle_resources import loan_capacity, commit_loan


def _atomic(path, row):
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as handle:
            json.dump(row, handle, ensure_ascii=False)
            handle.write('\n')
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def _living_waiters(root, cpus, memory_bytes):
    rows = []
    for path in sorted((root / 'waiters').glob('*.json')):
        try:
            row = json.loads(path.read_text())
            os.kill(int(row['pid']), 0)
            rows.append((path, row))
        except (FileNotFoundError, ProcessLookupError):
            pass
    policy_path=root/'policy.json'
    policy=json.loads(policy_path.read_text()) if policy_path.exists() else {}
    active=[row for _,row in read_pool_records(root/'leases')]
    return len(fitting_waiters(rows, max(0, cpus), max(0, memory_bytes),
                              stage_limits=policy.get('stage_concurrency'),active=active,
                              prepare_read_only=policy.get('prepare_phased') is True))


def reference_bases(hap_dirs, contigs):
    """FASTA 바이트 크기 대신 FAI의 실제 sequence 길이로 피크를 계산함."""
    lengths = {}
    for hap, folder in hap_dirs.items():
        total = 0
        for contig in contigs:
            source = folder / (contig + '.fa')
            fai = Path(str(source) + '.fai')
            if fai.exists():
                total += sum(int(line.split('\t')[1]) for line in fai.read_text().splitlines())
            else:
                # FAI가 없는 직접 입력도 과소 추정하지 않음.
                total += source.stat().st_size
        lengths[hap] = total
    return lengths


class RealignPhases:
    """완료된 하위 단계 사이에서만 예약을 바꿈. 실행 중 도구는 선점하지 않음."""

    def __init__(self):
        value = os.environ.get('MFM_CPU_POOL_LEASE_FILE')
        self.path = Path(value) if value else None
        from ..gpu_compare_resources import enabled
        self.unreserved=enabled()
        if self.unreserved:self.path=None
        affinity = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else os.cpu_count() or 1
        caps = [affinity]
        for key in ['MFM_CPU_BUDGET', 'SLURM_CPUS_PER_TASK']:
            if os.environ.get(key):
                caps.append(int(os.environ[key]))
        if self.path:
            caps.append(int(json.loads(self.path.read_text())['total_cpus']))
        self.total = max(1, min(caps))
        if self.path:
            self.total = pool_cpu_limit(self.path.parent.parent, self.total)
        self.ram_limit = memory_limit(self.path.parent.parent if self.path else None)
        for key, factor in [('SLURM_MEM_PER_NODE', 1), ('SLURM_MEM_PER_CPU', self.total)]:
            if int(os.environ.get(key, '0')) > 0:
                self.ram_limit = min(self.ram_limit, int(os.environ[key]) * 1024 ** 2 * factor)
        self.events = []
        self.started = None

    @property
    def resource_limits(self):
        event=self.events[-1]
        return {'cpus':event['reserved_cpus'],'memory_bytes':event['reserved_memory_bytes']}

    def ceiling(self):
        self.ram_limit = memory_limit(self.path.parent.parent if self.path else None)
        for key, factor in [('SLURM_MEM_PER_NODE', 1), ('SLURM_MEM_PER_CPU', self.total)]:
            if int(os.environ.get(key, '0')) > 0:
                self.ram_limit = min(self.ram_limit, int(os.environ[key]) * 1024 ** 2 * factor)
        return min(self.ram_limit * 9 // 10,
                   usable_memory(self.path.parent.parent) if self.path else self.ram_limit * 9 // 10)

    def alignment_plan(self, config, lengths, cpus, ram):
        baseline = alignment_plan(config, lengths, cpus, ram)
        if os.environ.get('MFM_WGS_BACKEND', 'cpu') == 'cpu':
            policy_path = self.path.parent.parent / 'policy.json' if self.path else None
            policy = json.loads(policy_path.read_text()) if policy_path and policy_path.exists() else {}
            buffer = int(policy.get('realign_sort_memory_bytes', 0))
            # CPU/HAP 병렬도를 유지할 RAM이 있을 때만 기본 sort 버퍼를 확대함.
            if (not buffer and os.environ.get('MFM_WGS_SORT_MEM_REQUESTED') == 'auto'
                    and os.environ.get('WGS_AUTO_THREADS', '1') == '1'
                    and os.environ.get('MFM_MEMORY_PROFILE_RESOLVED',
                        os.environ.get('MFM_MEMORY_PROFILE', 'balanced')) == 'balanced'):
                buffer = 8 * GIB
            while buffer > parse_memory(config.sort_mem):
                larger = replace(config, sort_mem=str(buffer))
                plan = alignment_plan(larger, lengths, cpus, ram)
                if plan is not None and (baseline is None or (plan['hap_jobs'] >= baseline['hap_jobs'] and plan['reserved_cpus'] >= baseline['reserved_cpus'] and plan['threads'] >= baseline['threads'] and plan['sort_threads'] >= baseline['sort_threads'])):
                    return dict(plan, sort_mem=larger.sort_mem)
                buffer //= 2
        return dict(baseline, sort_mem=config.sort_mem) if baseline is not None else None

    def reserve(self, phase, *, max_units, min_units=1, cpu_per_unit=1,
                memory_per_unit=0, base_memory=4 * GIB, base_cpus=0):
        """부모/자식 lease를 중첩하지 않고 가능한 동시 작업 수를 원자 예약함."""
        if not (1 <= min_units <= max_units) or cpu_per_unit < 1 or base_cpus < 0:
            raise ValueError('재정렬 단계의 최소/최대 자원이 올바르지 않음')
        required = base_memory + min_units * memory_per_unit
        ceiling = self.ceiling()
        if base_cpus + min_units * cpu_per_unit > self.total or required > ceiling:
            raise ValueError('재정렬 하위 단계의 최소 CPU/RAM이 전체 예산보다 큼')
        self.finish()
        queued = time.time()
        while True:
            if self.path is None:
                units = min(max_units, (self.total-base_cpus) // cpu_per_unit,
                            (ceiling - base_memory) // memory_per_unit if memory_per_unit else max_units)
                row = None
            else:
                root = self.path.parent.parent
                with (root / 'pool.lock').open('a+') as lock:
                    fcntl.flock(lock, fcntl.LOCK_EX)
                    prune_leases(root)
                    row = json.loads(self.path.read_text())
                    others = [other for path, other in read_pool_records(root / 'leases') if path != self.path]
                    policy_path=root/'policy.json'
                    setting=json.loads(policy_path.read_text()).get('cpu_backfill',{}) if policy_path.exists() else {}
                    available = available_regular_cpus(self.total,setting,others)
                    available = fair_cpu_capacity(root,self.total,[*others,row],available,request=row)
                    available_ram = available_memory(root,others,ceiling_bytes=self.ceiling())
                    ram = max(0, available_ram - base_memory)
                    loan = None
                    units = min(max_units, max(0,available-base_cpus) // cpu_per_unit,
                                ram // memory_per_unit if memory_per_unit else max_units)
                    if available_ram < base_memory:
                        units = 0
                    if units < max_units:
                        capacity = loan_capacity(root, self.total, row, others, available_ram)
                        if capacity is not None:
                            loan_cpu, loan_ram, candidate = capacity
                            loan_units = min(max_units, max(0,loan_cpu-base_cpus) // cpu_per_unit,
                                             max(0, loan_ram - base_memory) // memory_per_unit if memory_per_unit else max_units)
                            if loan_units >= min_units and loan_units > units and loan_ram >= base_memory:
                                units, loan = loan_units, candidate
                    waiting = _living_waiters(root, available - base_cpus - min_units * cpu_per_unit,
                                              ram - min_units * memory_per_unit)
                    if waiting and units >= min_units and loan is None:
                        units = min(units, max(min_units, max(0,available-base_cpus) // (waiting + 1) // cpu_per_unit))
                    if units >= min_units:
                        row.update(cpus=base_cpus + units * cpu_per_unit,
                                   memory_bytes=base_memory + units * memory_per_unit,
                                   realign_phase=phase)
                        row.pop('realign_wait', None)
                        commit_loan(root, row, loan)
                    else:
                        # 모든 이전 하위 작업이 종료한 뒤이므로 계산 없이 부모만 대기함.
                        row.update(cpus=min(1, available), memory_bytes=4 * GIB, realign_phase=phase + ':waiting')
                        row.pop('cpu_backfill', None)
                        row['realign_wait'] = {'queued_at': queued, 'minimum_cpus': base_cpus + min_units * cpu_per_unit,
                            'memory_bytes': required, 'available_regular_cpus': available,
                            'available_memory_bytes': available_ram,
                            'reason': '정상 CPU/RAM 또는 실측 유휴 대여 조건 대기'}
                    row['peak_cpus'] = max(int(row.get('peak_cpus', 0)), int(row['cpus']))
                    _atomic(self.path, row)
            if units >= min_units:
                event = {'phase': phase, 'queued_at': queued, 'started_at': time.time(),
                         'queue_wait_seconds': time.time() - queued, 'units': int(units),
                         'reserved_cpus': int(base_cpus + units * cpu_per_unit),
                         'reserved_memory_bytes': int(base_memory + units * memory_per_unit),
                         'ram_limit_bytes': self.ram_limit,
                         'cpu_pool_total': self.total}
                if self.unreserved:event.update(reservation_applied=False,cpu_reservation=0,memory_reservation_bytes=0)
                self.events.append(event)
                self.started = event
                print('[REALIGN-PHASE] ' + json.dumps(event, ensure_ascii=False), flush=True)
                return int(units)
            time.sleep(.5)

    def finish(self):
        if self.started is not None:
            self.started['finished_at'] = time.time()
            self.started['elapsed_seconds'] = self.started['finished_at'] - self.started['started_at']
            self.started = None

    def release(self):
        """전체 하위 작업 join 뒤 호출함. 실패해도 큰 CPU/RAM 예약을 남기지 않음."""
        self.reserve('finalize', max_units=1)
        self.finish()

    def index_jobs(self, lengths, requested):
        if requested not in ['auto', '1', '2']:
            raise ValueError('MFM_BWA_INDEX_HAP_JOBS는 auto|1|2여야 함')
        maximum = len(lengths) if requested == 'auto' else min(len(lengths), int(requested))
        # 직렬 index 생성은 BWA의 -t 옵션이 없으며 한 CPU를 예약함.
        per_index = (max(lengths.values()) * 40 + 2 * GIB) * 11 // 10
        self.index_address_space = per_index
        return self.reserve('bwa_index', max_units=maximum, memory_per_unit=per_index)

    def hap_completed(self, config, remaining_lengths):
        """한 HAP의 정렬+index+QC가 모두 끝난 뒤 그 몫만 반환함."""
        if not remaining_lengths or self.path is None:
            return
        haps = len(remaining_lengths)
        model = alignment_memory(config, remaining_lengths, haps, config.threads, config.sort_threads)
        cpus = haps * (config.threads + config.sort_threads + alignment_cpu_overhead(config) + getattr(config,'input_decode_threads',0))
        root = self.path.parent.parent
        with (root / 'pool.lock').open('a+') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            row = json.loads(self.path.read_text())
            if row.get('realign_phase') != 'align_restore_sort_qc':
                raise RuntimeError('HAP 예약 반환 중 현재 단계가 다름')
            row['cpus'] = min(int(row['cpus']), cpus)
            row['memory_bytes'] = min(int(row['memory_bytes']), model['reserved_memory_bytes'])
            row['effective_alignment'] = dict(row.get('effective_alignment', {}),
                hap_jobs=haps, reserved_cpus=row['cpus'], reserved_memory_bytes=row['memory_bytes'],
                memory_model=model)
            _atomic(self.path, row)
        event = {'phase': 'completed_hap_release', 'remaining_haps': haps,
                 'reserved_cpus': row['cpus'], 'reserved_memory_bytes': row['memory_bytes'],
                 'time': time.time()}
        self.events.append(event)
        print('[REALIGN-PHASE] ' + json.dumps(event, ensure_ascii=False), flush=True)

    def align_config(self, config, lengths):
        """실행 직전 auto 요청을 재산정하고 계획과 예약을 한 lock에서 확정함."""
        if self.alignment_plan(config, lengths, self.total, self.ceiling()) is None:
            raise ValueError('재정렬 최소 CPU/RAM 또는 명시 thread 요청이 전체 예산보다 큼')
        minimum_cpus = next(cpus for cpus in range(1, self.total+1)
                            if self.alignment_plan(config, lengths, cpus, self.ceiling()) is not None)
        self.finish()
        queued = time.time()
        while True:
            if self.path is None:
                plan = self.alignment_plan(config, lengths, self.total, self.ceiling())
            else:
                root = self.path.parent.parent
                with (root / 'pool.lock').open('a+') as lock:
                    fcntl.flock(lock, fcntl.LOCK_EX)
                    prune_leases(root)
                    row = json.loads(self.path.read_text())
                    others = [value for path, value in read_pool_records(root / 'leases')
                              if path != self.path]
                    policy_path=root/'policy.json'
                    setting=json.loads(policy_path.read_text()).get('cpu_backfill',{}) if policy_path.exists() else {}
                    available = available_regular_cpus(self.total,setting,others)
                    available = fair_cpu_capacity(root,self.total,[*others,row],available,request=row)
                    ram = available_memory(root,others,ceiling_bytes=self.ceiling())
                    minimum_plan = self.alignment_plan(config, lengths, minimum_cpus, ram)
                    waiting = _living_waiters(root, available - minimum_cpus,
                                              ram - minimum_plan['reserved_memory_bytes']) if minimum_plan else 0
                    budget = min(available, max(minimum_cpus, available // (waiting + 1))) if waiting else available
                    plan = self.alignment_plan(config, lengths, budget, ram)
                    loan = None
                    if (len(lengths)>1 and os.environ.get('MFM_WGS_HAP_REQUESTED')=='auto'
                            and (plan is None or plan['hap_jobs']==1)):
                        # 공정 분할만으로 HAP를 직렬화하지 않고 두 HAP 최소 몫을 확인한다.
                        for candidate in range(budget+1, available+1):
                            parallel=self.alignment_plan(config,lengths,candidate,ram)
                            if parallel is not None and parallel['hap_jobs']==2:
                                plan=parallel
                                break
                    if plan is None or plan['reserved_cpus'] < self.total:
                        capacity = loan_capacity(root, self.total, row, others, ram)
                        if capacity is not None:
                            loan_cpu, loan_ram, candidate = capacity
                            larger = self.alignment_plan(config, lengths, loan_cpu, loan_ram)
                            if larger is not None and (plan is None or
                                    larger['reserved_cpus'] > plan['reserved_cpus'] and larger['hap_jobs'] >= plan['hap_jobs']):
                                plan, loan = larger, candidate
                    if plan is None:
                        row.update(cpus=min(1, available), memory_bytes=4 * GIB,
                                   realign_phase='align_restore_sort_qc:waiting')
                        row.pop('cpu_backfill', None)
                        row['realign_wait'] = {'queued_at': queued, 'minimum_cpus': minimum_cpus,
                            'memory_bytes': self.alignment_plan(config, lengths, minimum_cpus, self.ceiling())['reserved_memory_bytes'],
                            'available_regular_cpus': available, 'available_memory_bytes': ram,
                            'reason': '정상 CPU/RAM 또는 실측 유휴 대여 조건 대기'}
                    else:
                        row.pop('hap_resources',None)
                        row.update(cpus=plan['reserved_cpus'], memory_bytes=plan['reserved_memory_bytes'],
                                   realign_phase='align_restore_sort_qc', effective_alignment=plan)
                        row.pop('realign_wait', None)
                        commit_loan(root, row, loan)
                    row['peak_cpus'] = max(int(row.get('peak_cpus', 0)), int(row['cpus']))
                    _atomic(self.path, row)
            if plan is not None:
                event = {'phase': 'align_restore_sort_qc', 'queued_at': queued,
                         'started_at': time.time(), 'queue_wait_seconds': time.time() - queued,
                         'ram_limit_bytes': self.ram_limit, 'cpu_pool_total': self.total,
                         'requested_hap_jobs': config.hap_jobs, 'requested_threads': config.threads,
                         'requested_sort_threads': config.sort_threads, **plan}
                if self.unreserved:event.update(reservation_applied=False,cpu_reservation=0,memory_reservation_bytes=0)
                self.events.append(event)
                self.started = event
                print('[REALIGN-PHASE] ' + json.dumps(event, ensure_ascii=False), flush=True)
                return replace(config, **{key: plan[key] for key in
                               ['hap_jobs', 'threads', 'sort_threads', 'fastq_threads', 'input_decode_threads', 'sort_mem']})
            time.sleep(.5)
