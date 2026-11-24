"""기존 CPU 손상 검출 계약을 유지하면서 CRAM block CRC만 선택적으로 GPU 실행함."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,tempfile,time
from . import alignment_integrity as cpu
from .integrity_cache import atomic


def validate(path,report):
    backend=os.environ.get('MFM_IO_BACKEND','cpu')
    if backend not in {'cpu','nvcomp'}:raise ValueError('MFM_IO_BACKEND는cpu/nvcomp여야함')
    path,report=Path(path),Path(report)
    with path.open('rb') as stream:magic=stream.read(4)
    # 손상 Deflate에 대한 CPU 해독/CRC 검출 보장을 유지함. BAM은 CPU로 검사함.
    if backend=='cpu' or magic!=b'CRAM' or path.stat().st_size<int(os.environ.get('MFM_GPU_CRAM_MIN_BYTES',str(256<<20))):
        return cpu.validate(path,report)
    before=cpu.fingerprint(path)
    identity={'schema':cpu.SCHEMA,'input':before,'validator_sha256':hashlib.sha256(Path(cpu.__file__).read_bytes()).hexdigest()}
    if report.is_file() and time.time_ns()-max(before['mtime_ns'],before['ctime_ns'])>2_000_000_000:
        try:
            old=json.loads(report.read_text())
            if old.get('status')=='pass' and all(old.get(k)==v for k,v in identity.items()):return old
        except (OSError,ValueError,TypeError):pass
    payload=dict(identity,status='checking',started=time.time())
    report.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='gpu_cram_check_',dir=report.parent) as temporary:
        result_path=Path(temporary)/'result.json'
        command=[sys.executable,'-B','-m','integrate_wgs_rna.gpu_frontend','cram-crc','--path',str(path),'--report',str(result_path)]
        pool=os.environ.get('MFM_CPU_POOL_DIR')
        if pool and not os.environ.get('MFM_CPU_POOL_LEASE_FILE'):
            affinity=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else (os.cpu_count() or 1)
            total=min(affinity,int(os.environ.get('MFM_CPU_BUDGET',str(affinity))))
            if os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_CPUS_PER_TASK'):total=min(total,int(os.environ['SLURM_CPUS_PER_TASK']))
            command=[sys.executable,'-B',str(Path(__file__).resolve().parents[1]/'scripts/cpu_pool.py'),
                     '--pool-dir',pool,'--total-cpus',str(total),'--branch','IO','--stage','cram_crc',
                     '--sample',os.environ.get('SAMPLE','input')+'_GPU_IO','--audit',str(Path(pool)/'io_usage.jsonl'),'--',*command]
        process=subprocess.run(command,capture_output=True,text=True)
        result=json.loads(result_path.read_text()) if result_path.exists() else {'status':'failed','error':process.stderr[-4000:] or 'GPU 검사 자식 미완료'}
    if cpu.fingerprint(path)!=before:raise RuntimeError('무결성 검사 중 입력 파일 변경')
    payload.update(result,finished=time.time())
    atomic(report,payload)
    if process.returncode or payload['status']!='pass':
        if result.get('corruption'):raise cpu.AlignmentCorruption(result['error'],result.get('error_offset',0),result.get('error_bytes'))
        if result.get('error_type')=='ValueError':raise ValueError(result.get('error',''))
        raise RuntimeError('GPU CRAM 검사 실행 오류: '+result.get('error',''))
    return payload
