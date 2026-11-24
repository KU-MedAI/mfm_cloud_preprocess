"""Short 인덱스와 FASTQ 캐시 준비를 하나의 CPU·RAM·디스크 예약 안에서 겹친다."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import os
import time
from .memory_policy import GIB
from .gpu_fastq import cached_inputs


class _InputShare:
    def __init__(self, budget, index_memory, write_path, index_disk):
        self.budget = budget
        self.index_memory = index_memory
        self.write_path, self.index_disk = write_path, index_disk
        self.cpus = self.maximum_cpus = budget.maximum_cpus - 1

    def available_cpus(self):
        return max(0, self.budget.available_cpus() - 1)

    def wait(self, *args, **kwargs):
        if not self.budget.try_phase('input_cache_wait_index_active', 1, self.index_memory + 4*GIB,
                write_paths=[self.write_path], disk_peak_bytes=self.index_disk):
            raise RuntimeError('실행 중인 BWA index의 예약을 보존할 수 없음')

    def phase(self, label, cpus, memory, **kwargs):
        if not 0 <= cpus <= self.maximum_cpus or memory > 4 * GIB:
            raise ValueError('동시 입력 준비의 CPU/RAM 배분 초과')
        deadline = time.monotonic() + float(os.environ.get('MFM_GPU_WAIT_SECONDS','86400'))
        while not self.budget.try_phase('bwa_index_and_input', cpus+1, memory+self.index_memory,
                write_paths=[self.write_path], disk_peak_bytes=kwargs.get('disk_peak_bytes',0)+self.index_disk):
            self.wait()
            if time.monotonic() >= deadline:
                raise TimeoutError('BWA index와 입력의 합산 예약 대기 초과')
            time.sleep(1)


def prepare_short(reference, work, source, key, input_root, budget, producer, expected_bytes, indexer, *, index_root=None):
    storage=Path(index_root) if index_root is not None else work
    storage.mkdir(parents=True,exist_ok=True)
    if storage.stat().st_dev != input_root.stat().st_dev:
        return None  # 공유 장부의 디스크 예약은 한 파일시스템 단위다.
    cpus = min(budget.maximum_cpus, budget.available_cpus())
    if cpus < producer + 8:
        return None
    memory = (reference.stat().st_size * 16 + 4 * GIB) * 11 // 10 + 4 * GIB
    if memory > budget.memory:
        return None
    disk_peak = max(8 * GIB, expected_bytes * 2 + 8 * GIB) + reference.stat().st_size * 16
    index_memory = memory - 4*GIB
    index_disk = reference.stat().st_size * 16
    write_path = storage / 'gpu_prepare'
    budget.phase('bwa_index_prepare', 1, memory, write_paths=[write_path], disk_peak_bytes=index_disk)
    share = _InputShare(budget, index_memory, write_path, index_disk)
    started = time.monotonic()
    times = {}

    def measure(label, fn):
        before = time.monotonic()
        try:
            return fn()
        finally:
            times[label] = time.monotonic() - before

    def prepare_inputs():
        try:
            return cached_inputs(source,key,input_root,True,share,producer,expected_bytes)
        finally:
            # 입력 생산/압축 자식이 모두 join된 뒤 index 몫만 유지함.
            share.wait()

    # index 1CPU·RAM은 유지하고 캐시 잠금 대기 중 FASTQ 몫만 반환한다.
    with ThreadPoolExecutor(2) as pool:
        index = pool.submit(measure, 'index_seconds', lambda: indexer(reference, work))
        inputs = pool.submit(measure, 'input_seconds', prepare_inputs)
        reference = index.result()
        fastqs, count, reused = inputs.result()
    return reference, fastqs, count, reused, dict(times, elapsed_seconds=time.monotonic()-started,
        cpu_budget=share.cpus+1 if not reused else 1, index_cpus=1, input_cpus=share.cpus if not reused else 0, post_input_reserved_cpus=1, reserved_memory_bytes=memory,
        expected_disk_peak_bytes=disk_peak)
