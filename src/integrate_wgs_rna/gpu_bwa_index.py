"""같은 FASTA와 BWA 실행파일의 완료된 인덱스만 내용 지문으로 재사용함."""
from pathlib import Path
from contextlib import contextmanager
import fcntl,hashlib,json,os,shutil,subprocess,sys,tempfile,time
import pysam
from .memory_policy import GIB

FILES=['reference.fa','reference.fa.fai',*[f'reference.fa{s}' for s in ['.amb','.ann','.bwt','.pac','.sa']]]

def digest(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        while block:=stream.read(4<<20):value.update(block)
    return value.hexdigest()

def fingerprint(path):
    st=Path(path).stat()
    return st.st_dev,st.st_ino,st.st_size,st.st_mtime_ns,st.st_ctime_ns

def cache_root(work):
    # 공유 보존은 명시 경로에서만 수행하고 기본값은 기존 work 정리 대상에 둠.
    configured=os.environ.get('MFM_GPU_BWA_INDEX_ROOT','').strip()
    return Path(configured).expanduser().resolve() if configured else Path(work).resolve()/'bwa_indexes'

def identity(target,work):
    target=Path(target).resolve();stamp=fingerprint(target)
    binary=Path(os.environ.get('MFM_GPU_BWA',str(Path(__file__).resolve().parents[2]/'tools/gpu/bwa')))
    if not binary.is_file():binary=Path(shutil.which(str(binary)) or shutil.which('bwa') or binary)
    if not binary.is_file():raise RuntimeError('GPU BWA 인덱스 생성에 bwa가 필요함. 00_environment.sh install-gpu-tools 사용함')
    binary=binary.resolve();tool_stamp=fingerprint(binary)
    contract=dict(schema='mfm_gpu_bwa_index_v1',reference_sha256=digest(target),reference_bytes=stamp[2],
                  bwa_sha256=digest(binary),index_options=['index'])
    if fingerprint(target)!=stamp or fingerprint(binary)!=tool_stamp:raise ValueError('인덱스 지문 계산 중 입력 변경')
    key=hashlib.sha256(json.dumps(contract,sort_keys=True).encode()).hexdigest()
    root=cache_root(work);root.mkdir(parents=True,exist_ok=True)
    return target,binary,contract,root/key,stamp,tool_stamp

@contextmanager
def locked(entry,*,blocking=True,budget=None):
    with (entry.parent/(entry.name+'.lock')).open('a') as handle:
        try:fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            if not blocking:
                yield False
                return
            if budget is not None:budget.wait('bwa_index_cache_wait')
            fcntl.flock(handle,fcntl.LOCK_EX)
        try:yield True
        finally:fcntl.flock(handle,fcntl.LOCK_UN)

def validate(entry,contract):
    if not entry.exists():return None
    if entry.is_symlink() or not entry.is_dir():raise ValueError('BWA 캐시 경로 형식 불일치')
    marker=entry/'complete.json'
    if not marker.is_file() or marker.is_symlink():raise ValueError('완료되지 않은 BWA 캐시 재사용 거부: '+str(entry))
    row=json.loads(marker.read_text())
    if row.get('contract')!=contract or set(row.get('files',{}))!=set(FILES):raise ValueError('BWA 캐시 계약 불일치')
    for name in FILES:
        path=entry/name;item=row['files'][name]
        if path.is_symlink() or not path.is_file() or path.stat().st_size!=item['size_bytes'] or digest(path)!=item['sha256']:
            raise ValueError('BWA 캐시 무결성 실패: '+str(path))
    if row['files']['reference.fa']['sha256']!=contract['reference_sha256']:raise ValueError('BWA 캐시 FASTA 지문 불일치')
    return row

def ready_hint(target,work):
    """실행 순서 결정용 힌트임. 실제 사용 전 ensure_index/lookup의 전량 검증이 필요함."""
    entry=identity(target,work)[3]
    return (entry/'complete.json').is_file()


def receipt(metrics,entry,row,reused,started):
    result=dict(cache_directory=str(entry),reference_sha256=row['contract']['reference_sha256'],
        bwa_sha256=row['contract']['bwa_sha256'],reused=reused,
        lookup_and_prepare_seconds=time.monotonic()-started,
        index_build_resources=row.get('build_resources'),validation='full_sha256_all_7_files')
    if metrics is not None:metrics.update(result)
    print('[GPU-BWA-INDEX] '+json.dumps(result,ensure_ascii=False),flush=True)

def lookup(target,work,metrics=None):
    started=time.monotonic();target,binary,contract,entry,stamp,tool_stamp=identity(target,work)
    with locked(entry,blocking=False) as acquired:
        if not acquired:return None
        row=validate(entry,contract)
        if row is None:return None
        if fingerprint(target)!=stamp or fingerprint(binary)!=tool_stamp:raise ValueError('BWA 캐시 검증 중 입력 변경')
        receipt(metrics,entry,row,True,started)
        return entry/'reference.fa'

def ensure_index(target,work,budget=None,metrics=None):
    started=time.monotonic();target,binary,contract,entry,stamp,tool_stamp=identity(target,work)
    with locked(entry,budget=budget):
        row=validate(entry,contract)
        reused=row is not None
        if row is None:
            memory=(target.stat().st_size*16+4*GIB)*11//10
            if budget is not None:
                if memory>budget.memory:raise ValueError('BWA 인덱스 RAM 예상량이 HAP 예약 한도를 초과함')
                budget.phase('bwa_index',1,memory,write_paths=[entry.parent/(entry.name+'.reserved')],disk_peak_bytes=target.stat().st_size*16+GIB)
            with tempfile.TemporaryDirectory(prefix='.bwa_build_',dir=entry.parent) as tmp:
                stage=Path(tmp);reference=stage/'reference.fa';shutil.copyfile(target,reference)
                if digest(reference)!=contract['reference_sha256'] or fingerprint(target)!=stamp:raise ValueError('BWA FASTA 복사 중 입력 변경')
                measured=stage/'resource.json'
                runner=Path(__file__).resolve().parents[2]/'scripts/resource_runner.py'
                command=[sys.executable,'-B',str(runner),'--output',str(measured),'--address-space-bytes',str(memory),'--',str(binary),'index',str(reference)]
                env=dict(os.environ)
                # 인덱스 전용 단일 프로세스 측정이며 정렬 boundary 증설을 호출하지 않음.
                env.pop('MFM_CPU_POOL_LEASE_FILE',None)
                with (stage/'index.log').open('w') as log:
                    try:subprocess.run(command,env=env,stdout=log,stderr=log,check=True)
                    except subprocess.CalledProcessError:
                        log.flush()
                        with (stage/'index.log').open('rb') as detail:
                            detail.seek(max(0,(stage/'index.log').stat().st_size-4096))
                            print('[BWA-INDEX-FAIL] '+detail.read().decode(errors='replace'),file=sys.stderr,flush=True)
                        raise
                pysam.faidx(str(reference))
                files={}
                for name in FILES:
                    path=stage/name
                    if not path.is_file() or not path.stat().st_size:raise ValueError('BWA 생성 파일 누락: '+name)
                    files[name]=dict(size_bytes=path.stat().st_size,sha256=digest(path))
                if files['reference.fa']['sha256']!=contract['reference_sha256']:raise ValueError('BWA 생성 중 캐시 FASTA 변경')
                if fingerprint(target)!=stamp or fingerprint(binary)!=tool_stamp:raise ValueError('BWA 인덱스 생성 중 입력/도구 변경')
                row=dict(contract=contract,files=files,build_resources=json.loads(measured.read_text()),created_at=time.time())
                (stage/'complete.json').write_text(json.dumps(row,ensure_ascii=False,indent=2)+'\n')
                stage.rename(entry)
        if fingerprint(target)!=stamp or fingerprint(binary)!=tool_stamp:raise ValueError('BWA 캐시 준비 중 입력 변경')
        receipt(metrics,entry,row,reused,started)
        return entry/'reference.fa'
