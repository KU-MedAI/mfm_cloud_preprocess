"""서버/날짜 전용 경로 없이 handoff local run을 연속 실행함."""
from pathlib import Path
from types import SimpleNamespace
import copy
import fcntl
import hashlib
import json
import os
import signal
import subprocess
import sys
import time

from . import download_only, planner
from .queue_resources import admit, active_ram, cpu_limit, effective_memory, full_ram, live_leases, download_admit
from .ready_queue import Scheduler
from .resume import sample_environment
from .sequence_mode import completion_mode
from .state_bindings import StateBindings
from .lease_identity import lease_matches_config
from .queue_recovery import request_recovery, recovery_stage_running


def locked(state, sample):
    path = Path(state) / 'samples' / (sample + '.lock')
    if not path.exists():
        return False
    with path.open('a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            return False
        except BlockingIOError:
            return True


def read(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return {}


def execution_contract(config):
    value = copy.deepcopy(config)
    value.setdefault('cpu_allocation_mode','dynamic')
    resolution = value.pop('resource_resolution', {})
    for key in ['sample_cpus', 'local_jobs', 'download_jobs', 'elastic_max_jobs']:
        value[key] = resolution.get('requested_' + key, value.get(key))
    # 큐의 선행 전송량은 이미 실행 중인 샘플의 데이터/계산 계약을 바꾸지 않는다.
    for key in ['download_jobs', 'prefetch_samples', 'prefetch_max_bytes']:
        value.pop(key, None)
    # auto 계산은 현재 MemAvailable에 따라 달라질 수 있음. 사용자 요청이 같은지 비교함.
    return value


def capacity_increase_only(before, after):
    """활성 worker의 데이터 계약을 보존하는 CPU/RAM 상한 증가만 허용함."""
    if effective_memory(after) < effective_memory(before) or cpu_limit(after) < cpu_limit(before):
        return False
    adjusted = copy.deepcopy(before)
    for key in ['memory_limit_bytes', 'cpu_pool_limit']:
        adjusted[key] = after.get(key)
    for key in ['MFM_RAM_LIMIT_BYTES', 'MFM_CPU_BUDGET', 'MFM_WGS_MIN_ALIGN_THREADS']:
        if key in after.get('pipeline_env', {}):
            adjusted.setdefault('pipeline_env', {})[key] = after['pipeline_env'][key]
    return execution_contract(adjusted) == execution_contract(after)


class ExistingWorker:
    def __init__(self, adapter, row):
        self.adapter, self.row = adapter, row

    def poll(self):
        sample = self.row['sample_id']
        if locked(self.adapter.sample_state(sample), sample):
            return None
        record = read(self.adapter.sample_state(sample) / 'samples' / (sample + '.json'))
        return 0 if record.get('status') == 'success' and self.adapter.verified(self.row) else 1


class QueueStopped(Exception):
    pass


class Adapter:
    def __init__(self, config, rows, api):
        self.DEFER_FINALIZE = True
        self.PRIORITY_POLICY = config.get('queue_policy', {})
        self.config = copy.deepcopy(config)
        self.config['memory_limit_bytes'] = effective_memory(config)
        self.rows, self.api = rows, api
        self.root = Path(api.ROOT)
        self.state = Path(config['state_root'])
        self.bindings = StateBindings(config, rows)
        self.recovery_stage_running = recovery_stage_running
        self.BASE = self.state / 'queue'
        self.BASE.mkdir(parents=True, exist_ok=True)
        self.config_path = self.state / 'effective_config.json'
        self.selection = self.state / 'selected_samples.tsv'
        self.SHARED_CPU_POOL = Path(config.get('pipeline_env', {}).get(
            'MFM_CPU_POOL_DIR', str(self.state / 'cpu_pool' / config['run_id'])))
        self.cpus = cpu_limit(config)
        self.config.setdefault('pipeline_env', {}).update(
            MFM_CPU_POOL_DIR=str(self.SHARED_CPU_POOL), MFM_CPU_BUDGET=str(self.cpus))
        self.config['_wgs_compact_stage'] = '12' if api.ORIGINAL_LAYOUT else '10'
        if 'MFM_REALIGN_PREPARE_MAX_CPUS' not in self.config['pipeline_env']:
            for row in rows:
                if row['platform'] == 'illumina':
                    self.config.setdefault('sample_pipeline_env', {}).setdefault(row['sample_id'], {}).setdefault(
                        'MFM_REALIGN_PREPARE_MAX_CPUS', '4')
        resolution = config.get('resource_resolution', {})
        requested = resolution.get('requested_local_jobs', 'auto')
        self.MAX_ACTIVE_SAMPLES = int(config['local_jobs'])
        downloads = resolution.get('requested_download_jobs', 'auto')
        self.DOWNLOAD_TARGET = int(config['download_jobs']) if downloads != 'auto' else None
        self.PREFETCH_LIMIT = config.get('prefetch_samples')
        for label, value in [('prefetch_samples', self.PREFETCH_LIMIT), ('download_jobs', self.DOWNLOAD_TARGET)]:
            if value is not None and (isinstance(value, bool) or not isinstance(value, int) or value < 1):
                raise ValueError(label + '는 양의 정수여야 함')
        from .queue_recovery import stages
        self.RECOVERY_STAGES = stages(api.ORIGINAL_LAYOUT)
        self.download_config = download_only.configure(self.config)
        self.processes = []
        self.stop_requested = False
        self.scheduler = None
        self.previous_status = None
        self.last_event = 0
        self.completed_sizes = {}
        self.common = SimpleNamespace(
            save=api.atomic_json, completed_ids=self.completed_ids, verified=self.verified,
            files=lambda number, wave: ('local', self.selection, self.config_path), event=self.event)
        self.ExternalWorker = lambda row, state: ExistingWorker(self, row)
        self.sample_locked = lambda state, sample: locked(
            self.sample_state(sample) if Path(state) == self.state else state, sample)
        self.state_root = lambda path: self.state
        self.worker_config = self.recovery_config
        self.prepare = lambda *args: None  # 전체 선택 references는 run에서 한 번 준비함.
        self.cpu_quota = lambda: self.cpus
        self.reserved_ram = lambda rows: full_ram(rows, self.config)
        self.admission_ram = lambda rows, candidate=None: active_ram(
            self.SHARED_CPU_POOL, rows, candidate, final_stage=self.config['_wgs_compact_stage'], config=self.config)

    def sample_state(self, sample):
        return self.bindings.state(sample)

    def download_lock_state(self, row):
        return self.bindings.download_state(row['sample_id'], self.download_config['state_root'])

    def verified(self, row, state=None):
        mode = completion_mode(self.config,row)
        root = self.api.sample_root(self.config, row)
        if sample_environment(self.config,row['sample_id']).get('CLEAN_INTERMEDIATES_ON_SUCCESS')=='1':
            from .cleanup import cleanup_complete
            if not cleanup_complete(root,row['sample_id'],mode):return False
        if not self.api.processing_complete(root, row['sample_id'], mode):
            return False
        return self.policy_matches(row)

    def policy_matches(self, row):
        # 미정리 결과는 복구 대상임. 저장 정책 변경과 구분함.
        root = self.api.sample_root(self.config, row)
        from .sequence_mode import payload_matches
        if not payload_matches(self.config,row,root):return False
        expected = sample_environment(self.config, row['sample_id']).get('WINDOW_STORAGE_POLICY')
        return expected is None or read(root / 'metadata/window_policy.json').get('storage_policy') == expected

    def completed_ids(self, rows):
        return {row['sample_id'] for row in rows
                if self.verified(row) and not locked(self.sample_state(row['sample_id']), row['sample_id'])}

    def target_reached(self):
        if self.config.get('selection_target'):
            return False
        total = 0
        for row in self.rows:
            sample = row['sample_id']
            if self.verified(row) and not locked(self.sample_state(sample), sample):
                if sample not in self.completed_sizes:
                    info = read(self.state / 'samples' / (sample + '.json'))
                    amount = info.get('model_bytes') if info.get('status') == 'success' else None
                    if amount is None:
                        root = self.api.sample_root(self.config, row)
                        amount = read(root / 'metadata/cleanup_completed.json').get('model_bytes')
                        if amount is None:
                            amount = self.api.model_bytes(root)
                    self.completed_sizes[sample] = int(amount)
                total += self.completed_sizes[sample]
        return total >= self.config['model_target_bytes']

    def detached_running(self, row):
        return (not locked(self.sample_state(row['sample_id']), row['sample_id'])
                and any(lease.get('branch') in {'WGS', 'RNA', 'QC'}
                        and lease_matches_config(lease, self.config)
                        for lease in live_leases(self.SHARED_CPU_POOL).get(row['sample_id'], [])))

    def recovery_config(self, context, sample, grant):
        request = read(self.BASE / '작업재개요청.json').get('samples', {}).get(sample)
        if not request:
            return self.config_path
        if request.get('stage') not in self.RECOVERY_STAGES:
            raise ValueError('현재 배포본의 WGS 복구 단계가 아님')
        config = copy.deepcopy(self.config)
        request_recovery(config, sample, request)
        path = self.BASE / ('resume_' + sample + '.json')
        self.api.atomic_json(path, config)
        return path

    def resource_check(self, rows, candidate=None, prefetch=()):
        if candidate and self.detached_running(candidate):
            return False, '분리 실행 중인 기존 단계 완료 대기'
        target = self.scheduler.download_target() if self.scheduler else self.DOWNLOAD_TARGET or 3
        return admit(self.config, self.SHARED_CPU_POOL, rows, candidate, prefetch, target)

    def download_resource_check(self, rows, candidate=None):
        from .sequence_mode import uses_bundled_only
        local_only = candidate is not None and uses_bundled_only(self.config,candidate,self.root)
        allowed, reason = download_admit(self.config, self.SHARED_CPU_POOL,local_only=local_only)
        if not allowed or candidate is None or local_only:
            return allowed, reason
        from .download_budget import storage_admit
        pending = [entry[2] for entry, _ in self.scheduler.downloads.values()] if self.scheduler else []
        return storage_admit(self.config, self.SHARED_CPU_POOL, pending, candidate)

    def download_state(self, state, wave):
        return Path(self.download_config['state_root'])

    def download_record(self, state, wave, sample):
        directory = self.bindings.download_state(sample, self.download_config['state_root'])
        return read(directory / 'samples' / (sample + '.json'))

    def event(self, status, **details):
        now = time.time()
        value = dict(status=status, controller_pid=os.getpid(), updated=now, **details)
        self.api.atomic_json(self.BASE / 'status.json', value)
        fingerprint = json.dumps([status, details], sort_keys=True, ensure_ascii=False)
        if status != 'sample_ready_queue' or fingerprint != self.previous_status or now - self.last_event >= 30:
            with (self.BASE / 'events.jsonl').open('a') as handle:
                handle.write(json.dumps(value, ensure_ascii=False) + '\n')
            if status == 'sample_ready_queue':
                self.previous_status, self.last_event = fingerprint, now

    def maintenance(self):
        self.ensure_watchdog()
        self.bindings.mirror(self.api)
        if self.stop_requested:
            raise QueueStopped('신규 투입 중지. 실행 중 worker는 완료까지 보존함.')

    def ensure_watchdog(self):
        """구 worker를 중단하지 않고 새 정책의 RSS 감시·대기 할당을 유지함."""
        script=self.root/'scripts/resource_watchdog.py'
        if not script.is_file():
            return
        value=read(self.SHARED_CPU_POOL/'watchdog.json')
        from integrate_wgs_rna.pool_leases import process_alive
        if process_alive(value.get('pid')) and time.time()-float(value.get('updated',0))<10:
            return
        log=self.BASE/'resource_watchdog.log'
        with log.open('a') as output:
            subprocess.Popen([sys.executable,str(script),
                '--pool-dir',str(self.SHARED_CPU_POOL),'--total-cpus',str(self.cpus),
                '--owner-pid',str(os.getpid())],cwd=self.root,stdout=output,stderr=subprocess.STDOUT,
                start_new_session=True)

    def request_stop(self, signum, frame):
        self.stop_requested = True

    def start(self, kind, tag, selection, config, *, task_index=None, wave=None,
              sample_cpus=None, worker_config_path=None):
        if self.stop_requested:
            raise QueueStopped('중지 요청 이후 신규 worker·다운로드 투입을 중단함.')
        if kind == 'finish':
            target = self.config_path
        elif kind == 'download-worker':
            target = Path(self.download_config['state_root']) / 'effective_config.json'
        else:
            target = worker_config_path or self.config_path
        selection = self.selection
        if task_index is not None:
            sample = self.rows[task_index]['sample_id']
            if sample in self.bindings.items:
                recovery = read(self.BASE / '작업재개요청.json').get('samples', {}).get(sample)
                target, selection, task_index = self.bindings.launch(
                    sample, kind, recovery, self.api, self.BASE / 'launch_configs')
        command = [sys.executable, '-B', '-m', 'handoff.cli', kind, '--effective-config', str(target),
                   '--selection', str(selection), '--executor', 'local']
        if task_index is not None:
            command += ['--task-index', str(task_index)]
        env = os.environ.copy()
        env['PYTHONPATH'] = os.pathsep.join(map(str, [self.root, self.root / 'src',
                                                   self.root / 'vendor/MFM-DNARNAFM']))
        for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
            env[name] = '1'
        logs = self.BASE / 'logs'
        logs.mkdir(exist_ok=True)
        with (logs / f'{kind}_{task_index}.log').open('a') as output:
            process = subprocess.Popen(command, cwd=self.root, env=env, stdout=output, stderr=subprocess.STDOUT)
        self.processes.append(process)
        return process


def run(config, rows, api):
    selected = hashlib.sha256('\n'.join(row['sample_id'] for row in rows).encode()).hexdigest()
    if config.get('selection_fingerprint') != selected:
        raise ValueError('effective-config와 선택 목록/순서가 다름. 기존 선택 파일로 재개하세요.')
    from .acceleration import preflight
    preflight(config,rows)
    adapter = Adapter(config, rows, api)
    with (adapter.BASE / 'control.lock').open('a') as lock, adapter.bindings.control_lock():
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        # 다른 제어큐가 살아 있을 때의 설정 덮어쓰기는 위 lock에서 먼저 차단함.
        prior = read(adapter.config_path)
        binding = adapter.state / 'selection_fingerprint.txt'
        expected = adapter.config['selection_fingerprint']
        if ((prior.get('selection_fingerprint') and prior['selection_fingerprint'] != expected)
                or (binding.exists() and binding.read_text().strip() != expected)):
            raise ValueError('기존 선택 범위와 다름. 별도 state_root를 사용하세요.')
        if (prior and any(locked(adapter.sample_state(r['sample_id']), r['sample_id']) for r in rows)
                and execution_contract(prior) != execution_contract(adapter.config)
                and not capacity_increase_only(prior, adapter.config)):
            raise ValueError('실행 중 worker가 있어 설정을 변경할 수 없음. 기존 effective-config로 재개하세요.')
        for row in rows:
            expected_policy = sample_environment(adapter.config, row['sample_id']).get('WINDOW_STORAGE_POLICY')
            if expected_policy and not locked(adapter.sample_state(row['sample_id']), row['sample_id']):
                mode = completion_mode(adapter.config,row)
                root = api.sample_root(adapter.config, row)
                if api.processing_complete(root, row['sample_id'], mode) and not adapter.policy_matches(row):
                    raise ValueError('기존 완료본의 window 정책 또는 서열 저장 정책이 다름. 별도 output_root/state_root를 사용하세요: '
                                     + row['sample_id'])
        api.atomic_json(adapter.config_path, adapter.config)
        temporary = adapter.selection.with_name('.selected_samples.' + str(os.getpid()) + '.tsv')
        planner.write_tsv(temporary, rows, fields=list(rows[0]))
        os.replace(temporary, adapter.selection)
        handlers = {s: signal.signal(s, adapter.request_stop)
                    for s in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)}
        stopped = False
        try:
            api.prepare(adapter.config, rows)
            from .stage_resources import longread_memory_plan
            api.atomic_json(adapter.BASE / 'ram_admission_plan.json', {
                'schema': 'mfm_queue_memory_plan_v2', 'limit_bytes': adapter.config['memory_limit_bytes'],
                'reservation_bytes': adapter.config['memory_limit_bytes']*9//10,
                'samples': {row['sample_id']: longread_memory_plan(row,adapter.config) for row in rows}})
            download_only.initialize(adapter.download_config, rows)
            adapter.ensure_watchdog()
            adapter.scheduler = Scheduler(adapter, [(1, rows)], rows)
            adapter.event('queue_started', cpu_pool=adapter.cpus, ram_limit_bytes=adapter.config['memory_limit_bytes'],
                          ram_reservation_bytes=adapter.config['memory_limit_bytes'] * 9 // 10,
                          compute_slots=adapter.MAX_ACTIVE_SAMPLES, disk_admission_check=False)
            result = adapter.scheduler.run()
            adapter.bindings.mirror(api)
            summary = api.finish(adapter.config, rows)
            adapter.event('queue_finished', returncode=result, selected=len(rows), success=summary['success'],
                          target_reached=adapter.target_reached())
            return result
        except QueueStopped:
            stopped = True
            return 130
        finally:
            # PTY 리더 종료로 하위 worker가 HUP을 받지 않도록 실제 자식 완료까지 유지함.
            if any(process.poll() is None for process in adapter.processes):
                adapter.event('draining', stopped=stopped)
                while any(process.poll() is None for process in adapter.processes):
                    time.sleep(1)
            for number, handler in handlers.items():
                signal.signal(number, handler)


def plan(config):
    """dry-run용 계획. 경로/상태/입력을 생성하지 않음."""
    resolution = config.get('resource_resolution', {})
    requested = resolution.get('requested_local_jobs', 'auto')
    return {'mode': 'ready_first_local', 'poll_seconds': 1,
            'completion_policy': 'all_selected' if config.get('selection_target') else 'model_target_bytes',
            'cpu_pool': cpu_limit(config), 'ram_limit_bytes': effective_memory(config),
            'ram_reservation_bytes': effective_memory(config) * 9 // 10,
            'compute_slots': int(config['local_jobs']),
            'download_request': resolution.get('requested_download_jobs', 'auto'),
            'prefetch_samples': config.get('prefetch_samples', 'download_target*2'),
            'disk_admission_check': False, 'download_retries': 3,
            'queue_policy': config.get('queue_policy', {'mode': 'ready_first'}),
            'acceleration': __import__('handoff.acceleration',fromlist=['plan']).plan(config)}
