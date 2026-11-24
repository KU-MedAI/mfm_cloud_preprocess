"""동일최종BAM/CRAM의물리검증증거를다운로드와worker가공유함."""
from pathlib import Path
import fcntl,hashlib,json,os,tempfile,time
from . import alignment_integrity as core

def atomic(path,payload):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(prefix='.'+path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as handle:json.dump(payload,handle,ensure_ascii=False,indent=2);handle.write('\n');handle.flush();os.fsync(handle.fileno())
        os.replace(name,path)
    finally:Path(name).unlink(missing_ok=True)

def publish_relocated(payload,target,report):
    """검증직후의통제된rename만최종경로fingerprint로연결함."""
    current=core.fingerprint(target);prior=payload['input']
    if payload.get('status')!='pass' or payload.get('validator_sha256')!=hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest():raise ValueError('유효한물리검증증거가아님')
    if any(prior[k]!=current[k] for k in ['inode','size','mtime_ns']):raise ValueError('검증한파일과rename후파일이다름')
    if Path(prior['path'])!=Path(str(Path(target).resolve())+'.part'):raise ValueError('다운로드부분파일rename경로와다름')
    # download소유lock아래에서완료게시를늦춤. worker에서fresh-cache검사를건너뛰지않음.
    age=(time.time_ns()-max(current['mtime_ns'],current['ctime_ns']))/1e9
    if 0<=age<2:time.sleep(2.01-age)
    if core.fingerprint(target)!=current:raise ValueError('rename안정화중파일변경')
    updated=dict(payload,input=current,relocated_from=prior,relocated_at=time.time());atomic(report,updated);return updated

def validate_shared(path,report,data_root):
    path=Path(path).resolve();state=Path(data_root)/'.download_records';state.mkdir(parents=True,exist_ok=True)
    key=hashlib.sha256(str(path).encode()).hexdigest();shared=state/(key+'.integrity.json')
    with (state/(key+'.lock')).open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        from .io_integrity import validate
        try:payload=validate(path,shared)
        except BaseException:
            if shared.is_file():atomic(report,json.loads(shared.read_text()))
            raise
        mirrored=dict(payload,shared_report=str(shared));atomic(report,mirrored);return mirrored
