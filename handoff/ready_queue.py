"""공용 local 실행의 선다운로드·준비 우선·실패 격리 스케줄러임."""
from pathlib import Path
import json
import errno
import os
import time
from .queue_priority import PriorityPolicy

MIN_DOWNLOAD_LIMIT = 3
MAX_DOWNLOAD_ATTEMPTS = 3
RETRY_DELAYS = (30, 120, 480)
POLL_SECONDS = 1.0


class Scheduler:
    def __init__(self, queue, work, all_rows, clock=time.time):
        self.q = queue
        self.c = queue.common
        self.work = work
        self.rows = all_rows
        self.clock = clock
        self.contexts = {}
        self.completed = self.c.completed_ids(all_rows)
        self.active = {}
        self.downloads = {}
        self.finished = set()
        self.record_cache = {}
        self.recoveries = {}
        self.last_slot_released = None
        self.ledger_path = queue.BASE / '재시도상태.json'
        self.ledger = json.loads(self.ledger_path.read_text()) if self.ledger_path.exists() else {
            'schema': 'mfm_local_retry_v1', 'downloads': {}, 'workers': {}}
        self.entries = [(position, index, row) for position, (_, wave) in enumerate(work)
                        for index, row in enumerate(wave)]
        self.priority = PriorityPolicy(self)

    def save(self):
        self.c.save(self.ledger_path, self.ledger)

    def refresh_recoveries(self):
        path = self.q.BASE / '작업재개요청.json'
        if not path.exists():
            return
        requests = json.loads(path.read_text()).get('samples', {})
        known = {e[2]['sample_id'] for e in self.entries}
        applied = self.ledger.setdefault('applied_recoveries', {})
        for sample, request in requests.items():
            if sample not in known or request.get('stage') not in getattr(self.q, 'RECOVERY_STAGES', {'06', '10', '12'}):
                raise ValueError('고정 배치 또는 승인된 복구 단계와 재개 요청이 다름')
            self.recoveries[sample] = request
            stamp = request['requested_at']
            if applied.get(sample) == stamp or sample in self.active:
                continue
            if request.get('retry_download'):
                retry = self.ledger['downloads'].get(sample)
                if retry:
                    previous = {k: v for k, v in retry.items() if k != 'recovery_history'}
                    retry.setdefault('recovery_history', []).append({
                        'time': self.clock(), 'previous': previous,
                        'reason': '명시적 다운로드 재개 요청. 실행 중 수신과 기존 부분 파일 유지.'})
                    retry.update(attempts=1 if sample in self.downloads else 0,
                                 deferred=False, next_retry_at=0, last_error='')
            self.ledger['workers'].pop(sample, None)
            applied[sample] = stamp
            self.c.event('recovery_requested', sample=sample, stage=request['stage'])
        self.save()

    def context(self, position):
        if position not in self.contexts:
            number, wave = self.work[position]
            tag, selection, config = self.c.files(number, wave)
            self.contexts[position] = dict(index=number, wave=wave, tag=tag, selection=selection,
                                           config=config, state=self.q.state_root(config), prepared=False)
        return self.contexts[position]

    def record(self, entry):
        position, _, row = entry
        sample = row['sample_id']
        if sample in self.record_cache:
            return self.record_cache[sample]
        ctx = self.context(position)
        value = self.q.download_record(ctx['state'], ctx['wave'], sample)
        self.record_cache[sample] = value
        return value

    def ready(self, entry):
        info = self.record(entry)
        if info.get('status') not in {'success', 'success_with_missing_indexes'}:
            return False
        # 정리된 완료 샘플은 이 함수 호출 전에 제외함. 남은 성공 영수증의 실파일도 확인함.
        return all(Path(f['local_path']).is_file() and
                   Path(f['local_path']).stat().st_size == f.get('actual_bytes')
                   for f in info.get('files', []) if f.get('status') == 'verified')

    def attempts(self, entry):
        sample = entry[2]['sample_id']
        if sample not in self.ledger['downloads']:
            info = self.record(entry)
            old_failure = info.get('status') == 'failed'
            self.ledger['downloads'][sample] = {
                'attempts': 1 if old_failure else 0,
                'next_retry_at': float(info.get('finished', 0)) + RETRY_DELAYS[0] if old_failure else 0,
                'last_error': info.get('error', '') if old_failure else '',
                'deferred': False}
        return self.ledger['downloads'][sample]

    def deferred(self, entry, *, worker=True):
        sample = entry[2]['sample_id']
        return (self.priority.paused(sample) or bool(self.priority.blocked(sample, worker=worker)) or sample in self.ledger['workers'] or
                bool(self.ledger['downloads'].get(sample, {}).get('deferred')))

    def detached_rows(self):
        running = getattr(self.q, 'recovery_stage_running', lambda request: False)
        return [e[2] for e in self.entries
                if e[2]['sample_id'] not in self.completed
                and e[2]['sample_id'] not in self.active
                and (running(self.recoveries.get(e[2]['sample_id'], {})) or
                     getattr(self.q, 'detached_running', lambda row: False)(e[2]))]

    def active_rows(self):
        return [entry[2] for entry, _ in self.active.values()] + self.detached_rows()

    def prefetch_rows(self):
        consumed = set(self.active) | {r['sample_id'] for r in self.detached_rows()}
        candidates = [entry[2] for entry, _ in self.downloads.values()]
        candidates += [e[2] for e in self.entries if e[2]['sample_id'] not in self.completed
                       and e[2]['sample_id'] not in consumed and not self.deferred(e, worker=False) and self.ready(e)]
        return list({row['sample_id']: row for row in candidates}.values())

    def download_target(self):
        configured = getattr(self.q, 'DOWNLOAD_TARGET', None)
        return configured if configured is not None else max(MIN_DOWNLOAD_LIMIT, len(self.active_rows()))

    def prefetch_limit(self):
        return getattr(self.q, 'PREFETCH_LIMIT', None) or 2 * self.download_target()

    def adopt_workers(self):
        for entry in self.entries:
            sample = entry[2]['sample_id']
            if sample in self.completed:
                continue
            ctx = self.context(entry[0])
            if self.q.sample_locked(ctx['state'], sample):
                self.active[sample] = (entry, self.q.ExternalWorker(entry[2], ctx['state']))
                self.c.event('adopted_existing_worker', sample=sample)
            download_state = self.download_lock_state(entry)
            if self.q.sample_locked(download_state, sample):
                self.downloads[sample] = (entry, ExistingDownload(self, entry))
                self.c.event('adopted_existing_download', sample=sample)

    def download_lock_state(self, entry):
        hook = getattr(self.q, 'download_lock_state', None)
        ctx = self.context(entry[0])
        return hook(entry[2]) if hook else self.q.download_state(ctx['state'], ctx['wave'])

    def collect_workers(self):
        for sample, (entry, process) in list(self.active.items()):
            result = process.poll()
            if result is None:
                continue
            del self.active[sample]
            self.last_slot_released = self.clock()
            if result == 0 and self.c.verified(entry[2], self.context(entry[0])['state']):
                self.completed.add(sample)
                self.c.event('sample_complete', sample=sample, active_workers=list(self.active))
            else:
                self.ledger['workers'][sample] = {'returncode': result, 'time': self.clock(),
                                                  'reason': 'worker 실패 또는 QC/정리 검증 미통과'}
                self.c.event('worker_deferred', sample=sample, returncode=result)
            self.save()

    def collect_downloads(self):
        for sample, (entry, process) in list(self.downloads.items()):
            result = process.poll()
            if result is None:
                continue
            del self.downloads[sample]
            retry = self.attempts(entry)
            if result == 0 and self.ready(entry):
                retry.update(next_retry_at=0, last_error='', deferred=False)
                self.c.event('download_complete', sample=sample, attempt=retry['attempts'])
            else:
                info = self.record(entry)
                retry.update(last_error=info.get('error', f'반환값 {result}, 다운로드 성공/실파일 검증 실패'),
                             next_retry_at=self.clock() + RETRY_DELAYS[min(retry['attempts'] - 1, 2)],
                             deferred=retry['attempts'] >= MAX_DOWNLOAD_ATTEMPTS)
                self.c.event('download_deferred' if retry['deferred'] else 'download_retry_wait',
                             sample=sample, attempt=retry['attempts'], next_retry_at=retry['next_retry_at'],
                             error=retry['last_error'])
            self.save()

    def start_workers(self):
        reasons = []
        for entry in self.priority.ordered(self.entries, worker=True):
            limit = self.q.MAX_ACTIVE_SAMPLES
            if limit is not None and len(self.active) + len(self.detached_rows()) >= limit:
                break
            sample = entry[2]['sample_id']
            if sample in self.completed or sample in self.active or self.deferred(entry):
                continue
            if (getattr(self.q, 'recovery_stage_running', lambda request: False)(self.recoveries.get(sample, {}))
                    or getattr(self.q, 'detached_running', lambda row: False)(entry[2])):
                reasons.append({'sample': sample, 'reason': '분리 실행 중인 기존 단계 완료 대기. 중복 계산하지 않음.'})
                continue
            priority_reason = self.priority.reason(entry[2], worker=True)
            if priority_reason:
                reasons.append({'sample': sample, 'reason': priority_reason})
                continue
            ctx = self.context(entry[0])
            if self.c.verified(entry[2], ctx['state']):
                self.completed.add(sample)
                continue
            if not self.ready(entry):
                continue
            sample_state = getattr(self.q, 'sample_state', lambda s: ctx['state'])(sample)
            record = sample_state / 'samples' / (sample + '.json')
            if record.exists() and json.loads(record.read_text()).get('status') == 'failed' and sample not in self.recoveries:
                self.ledger['workers'][sample] = {'time': self.clock(), 'reason': '기존 worker 실패 기록 보존'}
                self.c.event('worker_deferred', sample=sample, reason='기존 실패 기록 보존')
                self.save()
                continue
            if self.q.sample_locked(ctx['state'], sample):
                continue
            allowed, reason = self.q.resource_check(self.active_rows(), candidate=entry[2],
                                                     prefetch=self.prefetch_rows())
            if not allowed:
                reasons.append({'sample': sample, 'reason': reason})
                continue
            ready_since = self.ledger.setdefault('ready_since', {}).setdefault(
                sample, self.record(entry).get('finished') or self.clock())
            try:
                if not ctx['prepared']:
                    self.q.prepare(ctx['tag'], ctx['selection'], ctx['config'])
                    ctx['prepared'] = True
                config = json.loads(ctx['config'].read_text())
                budget = min(config['sample_cpus'], self.q.cpu_quota())
                process = self.q.start('worker', ctx['tag'], ctx['selection'], ctx['config'],
                                       task_index=entry[1], sample_cpus=budget,
                                       worker_config_path=self.q.worker_config(ctx, sample, budget))
            except Exception as error:
                if isinstance(error, OSError) and error.errno in {errno.ENOSPC, errno.EDQUOT}:
                    raise
                self.ledger['workers'][sample] = {'time': self.clock(), 'reason': '작업 시작 실패',
                                                  'error': repr(error)}
                self.save()
                self.c.event('worker_start_deferred', sample=sample, error=repr(error))
                reasons.append({'sample': sample, 'reason': '작업 시작 실패: ' + repr(error)})
                continue
            self.active[sample] = (entry, process)
            self.c.event('processing_ready_sample', sample=sample, worker_pid=process.pid,
                         active_workers=list(self.active), shared_cpu_pool_limit=self.q.cpu_quota(),
                         ram_peak_bytes=getattr(self.q, 'admission_ram', self.q.reserved_ram)(self.active_rows()),
                         ready_wait_seconds=max(0, self.clock() - float(ready_since)),
                         slot_released_to_start_seconds=(None if self.last_slot_released is None else
                                                        max(0, self.clock() - self.last_slot_released)))
        return reasons

    def initialize_download(self, ctx):
        state = self.q.download_state(ctx['state'], ctx['wave'])
        if not (state / 'effective_config.json').exists():
            # 묶음은 상태 식별용으로만 유지하며 실제 전송은 샘플별 worker로 실행함.
            from handoff.cli import load_config
            from handoff.download_only import configure, initialize
            initialize(configure(load_config(ctx['config'])), ctx['wave'])

    def start_downloads(self):
        reasons = []
        for entry in self.priority.ordered(self.entries, worker=False):
            if len(self.downloads) >= self.download_target():
                break
            sample = entry[2]['sample_id']
            if sample in self.completed or sample in self.active or sample in self.downloads or self.deferred(entry, worker=False):
                continue
            if self.priority.prefetch_count(entry[2]) >= self.prefetch_limit():
                continue
            priority_reason = self.priority.reason(entry[2], worker=False)
            if priority_reason:
                continue
            if self.ready(entry):
                continue
            retry = self.attempts(entry)
            if retry['attempts'] >= MAX_DOWNLOAD_ATTEMPTS:
                retry['deferred'] = True
                self.c.event('download_deferred', sample=sample, attempt=retry['attempts'],
                             error=retry['last_error'])
                self.save()
                continue
            if self.clock() < retry['next_retry_at']:
                continue
            download_check = getattr(self.q, 'download_resource_check', None)
            if download_check is None:
                allowed, reason = self.q.resource_check(self.active_rows(),
                                                         prefetch=[*self.prefetch_rows(), entry[2]])
            else:
                allowed, reason = download_check(self.active_rows(), candidate=entry[2])
            if not allowed:
                reasons.append({'sample': sample, 'reason': reason})
                continue
            ctx = self.context(entry[0])
            retry['attempts'] += 1
            retry['last_started'] = self.clock()
            self.save()  # 프로세스 생성 전에 기록하여 재시작으로 횟수가 초기화되지 않도록 함.
            try:
                self.initialize_download(ctx)
                process = self.q.start('download-worker', ctx['tag'], ctx['selection'], ctx['config'],
                                       task_index=entry[1], wave=ctx['wave'])
            except Exception as error:
                if isinstance(error, OSError) and error.errno in {errno.ENOSPC, errno.EDQUOT}:
                    raise
                retry.update(last_error=repr(error),
                             next_retry_at=self.clock() + RETRY_DELAYS[min(retry['attempts'] - 1, 2)],
                             deferred=retry['attempts'] >= MAX_DOWNLOAD_ATTEMPTS)
                self.save()
                self.c.event('download_start_deferred' if retry['deferred'] else 'download_start_retry_wait',
                             sample=sample, error=repr(error), attempt=retry['attempts'])
                reasons.append({'sample': sample, 'reason': '다운로드 시작 실패: ' + repr(error)})
                continue
            self.downloads[sample] = (entry, process)
            self.c.event('prefetching_sample', sample=sample, download_pid=process.pid,
                         attempt=retry['attempts'], active_downloads=list(self.downloads))
        return reasons

    def finalize(self):
        if getattr(self.q, 'DEFER_FINALIZE', False):
            return
        if self.active or self.downloads:
            return
        for position, (_, wave) in enumerate(self.work):
            if position in self.finished or not all(r['sample_id'] in self.completed for r in wave):
                continue
            ctx = self.context(position)
            process = self.q.start('finish', ctx['tag'], ctx['selection'], ctx['config'])
            if process.wait() != 0:
                raise RuntimeError(f'완료 묶음 집계 실패: {ctx["tag"]}')
            self.finished.add(position)

    def step(self):
        getattr(self.q, 'maintenance', lambda: None)()
        self.priority.refresh()
        self.record_cache.clear()
        self.refresh_recoveries()
        self.collect_workers()
        self.collect_downloads()
        limited = bool(getattr(self.q, 'target_reached', lambda: False)())
        worker_wait = [] if limited else self.start_workers()  # 준비 완료 샘플이 다운로드 재시도보다 우선함.
        download_wait = [] if limited else self.start_downloads()
        pending = [e for e in self.entries if e[2]['sample_id'] not in self.completed]
        deferred = [e[2]['sample_id'] for e in pending
                    if self.deferred(e) and not self.priority.paused(e[2]['sample_id'])]
        future = [e for e in pending if not self.deferred(e) or not self.deferred(e, worker=False)]
        self.save()
        running = {r['sample_id'] for r in self.active_rows()}
        ready_waiting = [e[2]['sample_id'] for e in pending
                         if e[2]['sample_id'] not in running and not self.deferred(e) and self.ready(e)]
        used_slots = len(running)
        limit = self.q.MAX_ACTIVE_SAMPLES
        if ready_waiting and limit is not None and used_slots >= limit:
            worker_wait.append({'sample': ','.join(ready_waiting[:3]),
                                'reason': f'실행 슬롯 사용 중: {used_slots}/{limit}'})
        self.c.event('sample_ready_queue', completed=len(self.completed), total=len(self.rows),
                     remaining=len(pending), active_workers=list(self.active),
                     active_downloads=list(self.downloads), deferred_samples=deferred,
                     ready_waiting=ready_waiting, used_compute_slots=used_slots,
                     max_compute_slots=limit,
                     compute_admission='RAM·cgroup·RSS 검사' if limit is None else '고정 슬롯+자원 검사',
                     prefetch_buffer_samples=[r['sample_id'] for r in self.prefetch_rows()],
                     prefetch_limit=self.prefetch_limit(), download_target=self.download_target(),
                     poll_seconds=POLL_SECONDS, priority_policy=self.priority.status(),
                     detached_stage_workers=[r['sample_id'] for r in self.detached_rows()],
                     recovery_waiting=[name for name in self.recoveries if name not in self.completed and name not in self.active],
                     retry_waiting=[e[2]['sample_id'] for e in future
                                    if self.ledger['downloads'].get(e[2]['sample_id'], {}).get('next_retry_at', 0) > self.clock()],
                     waiting=(worker_wait + download_wait)[:8])
        if not self.active and not self.downloads and limited:
            self.c.event('target_reached', completed=len(self.completed), total=len(self.rows),
                         remaining_samples=[e[2]['sample_id'] for e in pending])
            return 0
        if not self.active and not self.downloads and not future:
            self.finalize()
            self.c.event('complete' if not pending else 'complete_with_pending',
                         completed=len(self.completed), total=len(self.rows),
                         remaining_samples=[e[2]['sample_id'] for e in pending], deferred_samples=deferred)
            return 0 if not pending else 1
        return None

    def run(self):
        adopted = False
        storage_error_since = None
        while True:
            try:
                if not adopted:
                    self.adopt_workers()
                    adopted = True
                result = self.step()
                if storage_error_since is not None:
                    self.c.event('storage_write_recovered', since=storage_error_since,
                                 recovered_at=self.clock())
                    storage_error_since = None
            except OSError as error:
                if error.errno not in {errno.ENOSPC, errno.EDQUOT}:
                    raise
                if storage_error_since is None:
                    storage_error_since = self.clock()
                try:
                    os.write(2, ('[QUEUE-STORAGE-WAIT] ' + repr(error) + '\n').encode())
                except OSError:
                    pass
                time.sleep(30)
                continue
            if result is not None:
                return result
            time.sleep(POLL_SECONDS)


def run_rolling(queue, work, rows):
    return Scheduler(queue, work, rows).run()


class ExistingDownload:
    """이전 큐가 시작한 다운로드 lock을 인식하여 중복 전송을 방지함."""
    def __init__(self, scheduler, entry):
        self.scheduler = scheduler
        self.entry = entry

    def poll(self):
        s = self.scheduler
        ctx = s.context(self.entry[0])
        state = s.download_lock_state(self.entry)
        if s.q.sample_locked(state, self.entry[2]['sample_id']):
            return None
        return 0 if s.ready(self.entry) else 1
