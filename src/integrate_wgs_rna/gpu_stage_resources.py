"""GPU HAP의 준비·대기·계산·출력 예약을 같은 부모 장부에서 교환한다."""
from contextlib import contextmanager
from pathlib import Path
import fcntl
import json
import os
import shutil
import time

from .memory_policy import GIB, usable_memory
from .pool_leases import read_pool_records, regular_cpu_claim
from .realign.phases import _atomic


class StageBudget:
    def __init__(self):
        value = os.environ.get('MFM_CPU_POOL_LEASE_FILE')
        self.path = Path(value) if value and os.environ.get('MFM_GPU_HAP_ID') else None
        from .stage_gpu_config import apply_stage_environment
        apply_stage_environment(os.environ)
        from .gpu_compare_resources import enabled
        self.unreserved=enabled()
        if self.unreserved:self.path=None
        self.local=bool(os.environ.get('MFM_GPU_LOCAL_LEASE_FILE'))
        if self.local:self.path=Path(os.environ['MFM_GPU_LOCAL_LEASE_FILE'])
        self.hap = os.environ.get('MFM_GPU_HAP_ID', '1')
        self.events = []
        self.last_wait_reason = {}
        self.write_paths = []
        self.disk_peak_bytes = 0
        self.cpus = max(1, int(os.environ.get('MFM_GPU_STAGE_CPUS', '1')))
        self.maximum_cpus=self.cpus
        host=int(os.environ.get('MFM_GPU_HOST_MEMORY_BYTES',str(128*GIB)))
        peak=((host+8*GIB)*11//10+GIB-1)//GIB*GIB
        self.memory = int(os.environ.get('MFM_GPU_STAGE_MEMORY_BYTES', str(peak)))
        if self.memory>usable_memory(self.path.parent.parent if self.path else None):
            raise ValueError('GPU 도구 host RAM 피크가 실제 전체 예약 예산보다 큼')
        if self.path:
            with self._lock():
                row = json.loads(self.path.read_text())
                if row.get('effective_alignment',{}).get('hap_jobs',1)>1 and row.get('hap_resources_version')!=2:
                    # 이미 실행 중인 구형 부모의 HAP 장부를 비동기 갱신하지 않는다.
                    self.path=None
                else:
                    states = self._states(row)
                    state = states[self.hap]
                    self.cpus = int(state['cpus'])
                    if 'MFM_GPU_STAGE_CPUS' not in os.environ:self.maximum_cpus=self.cpus
                    self._write(row, states)
        self.cpus=self.maximum_cpus

    @contextmanager
    def _lock(self):
        with (self.path.parent.parent / 'pool.lock').open('a+') as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            yield

    def _states(self, row):
        if row.get('hap_resources') and (self.hap in row['hap_resources'] or row.get('effective_alignment',{}).get('hap_jobs',1)>1):
            return {str(k): dict(v) for k, v in row['hap_resources'].items()}
        count = max(1, int(row.get('effective_alignment', {}).get('hap_jobs', 1)))
        keys = [str(i) for i in (1, 2)] if count == 2 else [self.hap]
        return {h: {'cpus': int(row['cpus']) // count,
                    'memory_bytes': int(row['memory_bytes']) // count,
                    'phase': 'input_prepare'} for h in keys}

    def _write(self, row, states):
        row.update(hap_resources=states, cpus=max(0, sum(s['cpus'] for s in states.values())),
                   memory_bytes=max(4 * GIB, sum(s['memory_bytes'] for s in states.values())))
        row['peak_cpus'] = max(int(row.get('peak_cpus', 0)), row['cpus'])
        _atomic(self.path, row)

    def try_phase(self, label, cpus, memory, *, write_paths=(), disk_peak_bytes=0,preserve_disk=False):
        if self.path:
            with self._lock():
                row = json.loads(self.path.read_text())
                states = self._states(row)
                peers = [v for k, v in states.items() if k != self.hap]
                others = [] if self.local else [r for p, r in read_pool_records(self.path.parent) if p != self.path]
                old = states[self.hap]
                if disk_peak_bytes:
                    root=self.path.parent.parent
                    policy_path=root/'policy.json'
                    policy=json.loads(policy_path.read_text()) if policy_path.exists() else {}
                    reserve=int(policy.get('adaptive_prepare_storage',{}).get('reserve_bytes',150_000_000_000))
                    folder=Path(write_paths[0]).parent
                    while not folder.exists():folder=folder.parent
                    device=folder.stat().st_dev
                    commitments=list(peers)
                    for r in others:
                        commitments.extend(r.get('hap_resources',{}).values())
                        commitments.extend(r.get('prepare_components',{}).values())
                    outstanding=0
                    for s in commitments:
                        paths=[Path(p) for p in s.get('write_paths',[])]
                        if not paths or not paths[0].parent.exists() or paths[0].parent.stat().st_dev!=device:continue
                        written=sum(p.stat().st_size for p in paths if p.exists())
                        outstanding+=max(0,int(s.get('disk_peak_bytes',0))-written)
                    written=sum(Path(p).stat().st_size for p in write_paths if Path(p).exists())
                    need=max(0,disk_peak_bytes-written)
                    free=shutil.disk_usage(folder).free
                    if need+outstanding+reserve>free:
                        self.last_wait_reason={'resource':'disk','free_bytes':free,'reserve_bytes':reserve,
                                               'own_remaining_bytes':need,'peer_remaining_bytes':outstanding}
                        return False
                shrinking = cpus <= old['cpus'] and memory <= old['memory_bytes']
                peer_cpu=sum(s['cpus'] for s in peers)
                peer_memory=sum(s['memory_bytes'] for s in peers)
                other_memory=sum(int(r.get('memory_bytes',0)) for r in others)
                cpu_fits=cpus+peer_cpu+sum(regular_cpu_claim(r) for r in others)<=int(row['total_cpus'])
                if not shrinking and not self.local and row.get('cpu_backfill'):
                    from .realign.idle_resources import loan_capacity
                    capacity=loan_capacity(self.path.parent.parent,int(row['total_cpus']),row,others,
                                           usable_memory(self.path.parent.parent)-other_memory-peer_memory)
                    initial=int(row['cpu_backfill'].get('cpus',0))
                    cpu_fits=capacity is not None and cpus+peer_cpu<=min(initial,capacity[0]) and memory<=capacity[1]
                if not shrinking and (not cpu_fits or
                        memory+peer_memory+other_memory>usable_memory(self.path.parent.parent)):
                    self.last_wait_reason={'resource':'cpu_or_ram','requested_cpus':cpus,'requested_memory_bytes':memory}
                    return False
                states[self.hap] = {'cpus': cpus, 'memory_bytes': memory, 'phase': label,
                                    'worker_pid': os.getpid(),
                                    'write_paths':old.get('write_paths',[]) if preserve_disk else list(map(str,write_paths)),
                                    'disk_peak_bytes':old.get('disk_peak_bytes',0) if preserve_disk else disk_peak_bytes,
                                    'wait_reason':dict(self.last_wait_reason)}
                self._write(row, states)
        event={'phase': label, 'time': time.time(), 'cpus': cpus,'memory_bytes': memory, 'hap': self.hap,'reservation_applied':not self.unreserved}
        if not self.events or any(self.events[-1][k]!=event[k] for k in ['phase','cpus','memory_bytes']):
            self.events.append(event)
        return True

    def wait(self, label='gpu_wait', *, preserve_disk=False):
        if not self.try_phase(label, 0 if label=='storage_wait' else 1, 4 * GIB,preserve_disk=preserve_disk):
            raise RuntimeError('GPU 대기 예약을 반환할 수 없음')

    def phase(self, label, cpus, memory, *, write_paths=(), disk_peak_bytes=0,preserve_disk=False):
        self.wait(label + ':waiting',preserve_disk=preserve_disk)
        deadline = time.monotonic() + float(os.environ.get('MFM_GPU_WAIT_SECONDS', '86400'))
        while not self.try_phase(label, cpus, memory,write_paths=write_paths,disk_peak_bytes=disk_peak_bytes,preserve_disk=preserve_disk):
            if self.last_wait_reason.get('resource')=='disk':self.wait('storage_wait',preserve_disk=preserve_disk)
            if time.monotonic() >= deadline:
                raise TimeoutError('GPU host 자원 대기 초과: ' + label)
            time.sleep(1)

    def available_cpus(self):
        if not self.path:return self.maximum_cpus
        with self._lock():
            row=json.loads(self.path.read_text());states=self._states(row)
            others=[] if self.local else [r for p,r in read_pool_records(self.path.parent) if p!=self.path]
            peers=[s for h,s in states.items() if h!=self.hap]
            peer_cpu=sum(s['cpus'] for s in peers)
            used=peer_cpu+sum(regular_cpu_claim(r) for r in others)
            available=int(row['total_cpus'])-used
            if not self.local and row.get('cpu_backfill'):
                from .realign.idle_resources import loan_capacity
                ram=usable_memory(self.path.parent.parent)-sum(int(r.get('memory_bytes',0)) for r in others)-sum(s['memory_bytes'] for s in peers)
                capacity=loan_capacity(self.path.parent.parent,int(row['total_cpus']),row,others,ram)
                if capacity:available=max(available,min(int(row['cpu_backfill']['cpus']),capacity[0])-peer_cpu)
            return max(0,min(self.maximum_cpus,available))

    @contextmanager
    def gpu(self, required_bytes, *, purpose, count=1, producer_cpus=0):
        from .gpu_pool import acquire
        minimum=max(min(4,self.maximum_cpus),producer_cpus+2 if producer_cpus else 1)
        if minimum>self.maximum_cpus:raise ValueError('실행 중 입력 생산자와 GPU 보조 CPU 예산이 부족함')
        if producer_cpus:
            if not self.try_phase('gpu_wait_input_active',producer_cpus,4*GIB,preserve_disk=True):
                raise RuntimeError('실행 중 FIFO 생산자의 예약을 보존할 수 없음')
        else:self.wait()
        deadline = time.monotonic() + float(os.environ.get('MFM_GPU_WAIT_SECONDS', '86400'))
        while True:
            lease = None
            try:
                lease = acquire(required_bytes, purpose=purpose, exclusive=True, timeout=0, count=count)
                target=self.available_cpus()
                if target>=minimum and self.try_phase('gpu_compute', target, self.memory,
                                  write_paths=self.write_paths,disk_peak_bytes=self.disk_peak_bytes):
                    self.cpus=target
                    break
            except TimeoutError:
                pass
            except BaseException:
                if lease is not None:
                    lease.close()
                raise
            if lease is not None:
                lease.close()
            if self.last_wait_reason.get('resource')=='disk' and not producer_cpus:self.wait('storage_wait')
            if time.monotonic() >= deadline:
                raise TimeoutError('GPU와 CPU/RAM 동시 확보 대기 초과: ' + purpose)
            time.sleep(1)
        try:
            yield lease
        finally:
            lease.close()

    def output(self):
        target=min(self.maximum_cpus,int(os.environ.get('MFM_GPU_OUTPUT_CPUS',str(self.cpus))))
        if self.path:
            with self._lock():
                row=json.loads(self.path.read_text());peers=[s for h,s in self._states(row).items() if h!=self.hap]
            # 다음 GPU 계산을 기다리는 HAP에 CPU를 남기며 기존 실행을 선점하지 않음.
            if any(s['phase'] not in {'complete','output_restore_sort','output_index','output_qc','axis_qc'} for s in peers):
                target=min(target,16)
        target=max(4,min(target,self.available_cpus()))
        from .gpu_cram import plan
        memory=max(8*GIB,plan(target,64)['estimated_peak_bytes'])
        self.phase('output_restore_sort',target,memory,preserve_disk=True)
        self.cpus=target
