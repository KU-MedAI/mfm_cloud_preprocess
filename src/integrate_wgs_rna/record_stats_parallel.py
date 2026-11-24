"""BAM을 인덱스 없이 한 번 읽는 bounded 네이티브 집계의 빌드·자원 계획."""
from pathlib import Path
import fcntl
import hashlib
import json
import os
import shutil
import subprocess
import sys
import sysconfig
import tempfile

MIB = 1024**2
GIB = 1024**3


class NativeUnavailable(RuntimeError):
    pass


class RecordSizeLimit(RuntimeError):
    pass


def scan_plan(cpus):
    cpus = max(1, int(cpus))
    if cpus == 1:
        decode = count = 0
    elif cpus < 4:
        decode, count = 0, cpus - 1
    else:
        decode = max(1, (cpus - 2) * 2 // 3)
        count = cpus - 2 - decode
    batch, maximum = 8 * MIB, 64 * MIB
    # 생산자 1개, HTSlib I/O 1개(병렬 해독 시), 해독·집계 worker를 모두 계상한다.
    threads = 1 + int(decode > 0) + decode + count
    peak = (256 * MIB + (3 * count + 3) * (maximum + batch) + decode * 16 * MIB) * 12 // 10
    return dict(cpu_budget=cpus, native_thread_limit=threads, decode_threads=decode,
                count_workers=count, batch_bytes=batch, max_record_bytes=maximum,
                estimated_peak_rss_bytes=peak, passes=1, index_required=False,
                mode='single_pass_parallel_bam')


def native_tool(name="record_stats"):
    if name not in {"gpu_fastq","record_stats","restore_tags","trace_hash","rna_depth","rna_fastq","rna_ase_partition","rna_primary","rna_variant_support","rna_cnv_sum","rna_ase_read","pileup_oracle","pileup_counts","pileup_counts_complete","pileup_oracle_complete"}:raise ValueError("미지원 native 도구")
    compiler = shutil.which(os.environ.get('CXX','c++')) or shutil.which('g++') or shutil.which('x86_64-conda-linux-gnu-c++')
    include, library = Path(sys.prefix) / 'include', Path(sys.prefix) / 'lib'
    header, shared = include / 'htslib/bgzf.h', library / 'libhts.so'
    if not compiler or not header.exists() or not shared.exists():
        raise NativeUnavailable('C++ 컴파일러 또는 현재 Python 환경의 HTSlib 개발 파일 없음')
    source = Path(__file__).parent / 'native' / (name+'.cpp')
    signature = source.read_bytes() + str((compiler, sys.prefix, shared.resolve(), shared.stat().st_size,
                                         shared.stat().st_mtime_ns)).encode()
    if name in {'pileup_counts_complete','pileup_oracle_complete'}:
        signature += (source.parent/'pileup_complete.hpp').read_bytes()
    if name=='rna_ase_read':
        pyinclude=Path(sysconfig.get_path('include'))
        if not (pyinclude/'Python.h').is_file():raise NativeUnavailable('현재 Python 개발 header 없음')
        signature+=str((sys.version_info[:3],pyinclude)).encode()
        signature+=(source.parent/'rna_ase_fused.hpp').read_bytes()
    key = hashlib.sha256(signature).hexdigest()[:24]
    cache = Path(os.environ.get('XDG_CACHE_HOME', str(Path.home() / '.cache')))
    root = Path(os.environ.get('MFM_NATIVE_CACHE_ROOT', str(cache / 'mfm_native')))
    try:
        root.mkdir(mode=0o700, parents=True, exist_ok=True)
    except OSError as error:
        raise NativeUnavailable('네이티브 cache 생성 실패: ' + str(error)) from error
    output = root / (name+'_' + key)
    with (root / (name+'_' + key + '.lock')).open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if output.exists():
            return output
        fd, temporary = tempfile.mkstemp(prefix=name+'_build_', dir=root); os.close(fd)
        try:
            command = [compiler, '-O3', '-std=c++17', '-pthread', str(source), '-I' + str(include),
                       '-L' + str(library), '-Wl,-rpath,' + str(library), '-lhts', '-o', temporary]
            if name=='rna_fastq':command+=['-lz']
            if name=='rna_ase_read':command+=['-shared','-fPIC','-I'+str(pyinclude)]
            if name=='trace_hash':command+=['-shared','-fPIC','-lcrypto']
            if name in {'pileup_oracle','pileup_counts','pileup_counts_complete','pileup_oracle_complete'}:command+=['-shared','-fPIC']
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode:
                raise NativeUnavailable('네이티브 집계 빌드 실패: ' + result.stderr[-3000:])
            os.chmod(temporary, 0o700); os.replace(temporary, output)
        finally:
            Path(temporary).unlink(missing_ok=True)
    return output


def run_scan(path, cpus, *, tool=None, plan=None):
    plan = plan or scan_plan(cpus)
    tool = tool or native_tool()
    command = [str(tool), str(path), str(plan['decode_threads']), str(plan['count_workers']),
               str(plan['batch_bytes']), str(plan['max_record_bytes'])]
    try:
        result = subprocess.run(command, capture_output=True, text=True)
    except OSError as error:
        raise NativeUnavailable('네이티브 집계 실행 불가: ' + str(error)) from error
    if result.returncode == 75:
        raise RecordSizeLimit(result.stderr.strip())
    if result.returncode:
        raise RuntimeError('전수 병렬 집계 실패: ' + result.stderr[-4000:])
    values = json.loads(result.stdout)
    if values['execution']['native_threads'] > plan['cpu_budget']:
        raise RuntimeError('네이티브 집계 실제 스레드가 CPU grant를 초과함')
    return values
