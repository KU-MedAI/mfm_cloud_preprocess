"""FASTQ 계약을 검증하며 pigz로 압축하고 동일 입력을 HAP 사이에 공유한다."""
from pathlib import Path
import fcntl
import gzip
import hashlib
import json
import os
import shutil
import signal
import subprocess
import time

from .memory_policy import GIB
from .realign.phases import _atomic

BATCH_BYTES = 4 << 20
MAX_LINE_BYTES = 64 << 20


def split_fastq(stream, paths, paired, *, cpus=1, producer_cpus=0, mirrors=(), execution=None):
    processes, files, buffers = [], [], [bytearray() for _ in paths]
    mode = os.environ.get('MFM_GPU_FASTQ_BACKEND', 'auto')
    if mode not in {'auto', 'native', 'python'}:
        raise ValueError('GPU FASTQ 검증은 auto|native|python이어야 함')
    tool = None
    if (mode == 'native' or (mode == 'auto' and paired)) and cpus - producer_cpus >= 3:
        from .record_stats_parallel import native_tool, NativeUnavailable
        try:
            tool = native_tool('gpu_fastq')
        except NativeUnavailable:
            if mode == 'native':
                raise
    elif mode == 'native':
        raise ValueError('네이티브 FASTQ 연결부의 CPU 3개가 부족함')
    # parser·feed·프레임 I/O와 pigz reader/writer를 모두 같은 grant에 포함한다.
    bridge = 3 if tool else 1
    workers = (cpus - producer_cpus - bridge - 2 * len(paths)) // len(paths)
    pigz = shutil.which('pigz') if workers >= 2 else None
    if execution is not None:
        execution.update(backend='native' if tool else 'python', cpu_budget=cpus,
                         producer_cpus=producer_cpus, bridge_cpus=bridge,
                         compression_workers_per_file=workers if pigz else 0,
                         frame_bytes=BATCH_BYTES if tool else None,
                         native_tool=str(tool) if tool else None)
    count = 0
    mate_name = None
    ok = False
    try:
        for p in paths:
            if pigz:
                out = open(p, 'wb')
                try:
                    proc = subprocess.Popen([pigz, '-1', '-p', str(workers), '-n', '-c'],
                                            stdin=subprocess.PIPE, stdout=out)
                finally:
                    out.close()
                processes.append(proc)
                files.append(proc.stdin)
            else:
                files.append(gzip.open(p, 'wb', compresslevel=1))
        if tool:
            from .gpu_fastq_native import native_frames
            count = native_frames(tool, stream, paired, files, mirrors)
        else:
            while True:
                first = stream.readline(MAX_LINE_BYTES+1)
                if not first:
                    break
                sequence, plus, quality = (stream.readline(MAX_LINE_BYTES+1) for _ in range(3))
                if any(len(v)>MAX_LINE_BYTES for v in [first,sequence,plus,quality]):
                    raise MemoryError('FASTQ 단일 행이 bounded 입력 버퍼를 초과함')
                if not first.startswith(b'@') or not plus.startswith(b'+') or not quality or len(sequence.rstrip(b'\r\n')) != len(quality.rstrip(b'\r\n')):
                    raise ValueError('GPU 입력 FASTQ 레코드 손상')
                name = first.split()[0].removesuffix(b'/1').removesuffix(b'/2')
                if paired:
                    if count % 2 == 0:
                        mate_name = name
                    elif name != mate_name:
                        raise ValueError('GPU paired FASTQ R1/R2 QNAME 불일치')
                i = count % 2 if paired else 0
                buffers[i].extend(first + sequence + plus + quality)
                if len(buffers[i]) >= BATCH_BYTES:
                    files[i].write(buffers[i])
                    if mirrors:mirrors[i].write(buffers[i])
                    buffers[i].clear()
                count += 1
            if paired and count % 2:
                raise ValueError('GPU paired FASTQ mate 누락')
            if not count:
                raise ValueError('GPU 입력 FASTQ가 비어 있음')
            for i,(f, b) in enumerate(zip(files, buffers)):
                if b:
                    f.write(b)
                    if mirrors:mirrors[i].write(b)
        ok = True
    finally:
        for f in files:
            try:
                f.close()
            except BrokenPipeError:
                ok = False
        for proc in processes:
            if not ok and proc.poll() is None:
                proc.terminate()
            if proc.wait() != 0 and ok:
                raise RuntimeError('pigz FASTQ 압축 실패')
    if not ok:
        raise RuntimeError('FASTQ 압축 pipe 실패')
    return count


def cached_inputs(source_script, key, work, paired, budget, producer_cpus, expected_bytes=0, *, streaming=None):
    directory = Path(work) / ('gpu_input_' + hashlib.sha256(key.encode()).hexdigest()[:24])
    directory.mkdir(parents=True, exist_ok=True)
    paths = [directory / f'input_R{i + 1}.fq.gz' for i in range(2 if paired else 1)]
    receipt = directory / 'complete.json'
    budget.wait('input_cache_wait')
    deadline=time.monotonic()+float(os.environ.get('MFM_GPU_WAIT_SECONDS','86400'))
    with (directory / 'input.lock').open('a+') as lock:
        while True:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if time.monotonic()>=deadline:
                    raise TimeoutError('공통 FASTQ 생성 대기 초과')
                time.sleep(1)
        if receipt.exists():
            r = json.loads(receipt.read_text())
            if r['key'] != key or r['paired'] != paired or any(
                    not p.is_file() or p.stat().st_size != v['bytes'] or p.stat().st_mtime_ns != v['mtime_ns']
                    for p, v in zip(paths, r['files'])) or len(r['files']) != len(paths):
                raise ValueError('GPU 공통 FASTQ 완료 장부와 파일 불일치')
            return paths, r['records'], True
        temporary = [p.with_suffix(p.suffix + '.part') for p in paths]
        budget.cpus=min(budget.maximum_cpus,max(producer_cpus+3,budget.available_cpus()))
        budget.phase('input_prepare', min(budget.cpus,streaming.cpus) if streaming else budget.cpus, 4 * GIB,write_paths=temporary,
                     disk_peak_bytes=max(8*GIB,int(expected_bytes)*2+8*GIB))
        _atomic(directory/'resource_plan.json',{'cpus':min(budget.cpus,streaming.cpus) if streaming else budget.cpus,'producer_cpus':producer_cpus,
                'reserved_memory_bytes':4*GIB,'batch_bytes':BATCH_BYTES,'max_line_bytes':MAX_LINE_BYTES,
                'buffer_peak_bound_bytes':MAX_LINE_BYTES*16+2*BATCH_BYTES+512*(1<<20),
                'margin_percent':10,'expected_disk_peak_bytes':max(8*GIB,int(expected_bytes)*2+8*GIB)})
        started = time.time()
        proc = subprocess.Popen(['bash', str(source_script)], stdout=subprocess.PIPE,start_new_session=True)
        try:
            stream=streaming.attach(proc) if streaming else proc.stdout
            execution = {}
            count = split_fastq(stream, temporary, paired, cpus=min(budget.cpus,streaming.cpus) if streaming else budget.cpus, producer_cpus=producer_cpus+len(streaming.mirrors) if streaming else producer_cpus, mirrors=streaming.mirrors if streaming else (), execution=execution)
            proc.stdout.close()
            if proc.wait() != 0:
                raise RuntimeError('GPU 입력 source command 실패')
            if streaming:streaming.finish_mirrors()
            for src, dst in zip(temporary, paths):
                os.replace(src, dst)
            _atomic(receipt, {'key': key, 'paired': paired, 'records': count,
                             'elapsed_seconds': time.time() - started, 'execution': execution,
                             'files': [{'bytes': p.stat().st_size, 'mtime_ns': p.stat().st_mtime_ns} for p in paths]})
        finally:
            if proc.poll() is None:
                os.killpg(proc.pid,signal.SIGTERM); proc.wait()
            proc.stdout.close()
            for p in temporary:
                p.unlink(missing_ok=True)
        return paths, count, False
