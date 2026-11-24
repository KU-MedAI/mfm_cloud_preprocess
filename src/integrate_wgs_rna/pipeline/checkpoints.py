"""입력·정책·출력 지문이 같은 성공 단계만 재사용함."""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import time
import uuid


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()


def file_identity(path):
    p=Path(path).resolve(strict=True)
    s=p.stat()
    if not p.is_file():raise ValueError(f'일반 파일이 아님: {p}')
    identity=[str(p),s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns]
    # 공유 파일시스템의 시간 해상도 때문에 작은 설정/인덱스의 동크기 변경을 놓치지 않음.
    if s.st_size <= 1 << 20:
        identity.append(hashlib.sha256(p.read_bytes()).hexdigest())
        after=p.stat()
        if (s.st_size,s.st_mtime_ns,s.st_ctime_ns)!=(after.st_size,after.st_mtime_ns,after.st_ctime_ns):
            raise ValueError(f'입력 지문 계산 중 파일 변경: {p}')
    return identity


def tree_identity(path, *, exclude=()):
    root=Path(path);items=[]
    if not root.is_dir():raise FileNotFoundError(root)
    for parent,dirs,files in os.walk(root,followlinks=False):
        for name in sorted(dirs+files):
            p=Path(parent)/name
            if p.is_symlink():raise ValueError(f'완료 출력에 symlink가 있음: {p}')
        for name in files:
            if name in exclude:continue
            p=Path(parent)/name;s=p.stat()
            items.append([str(p.relative_to(root)),s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns])
    return {'files':len(items),'stat_sha256':digest(sorted(items))}


def atomic_json(path,value):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_name(f'.{p.name}.{os.getpid()}.{uuid.uuid4().hex}.part')
    try:
        with tmp.open('w') as f:
            json.dump(value,f,ensure_ascii=False,sort_keys=True);f.write('\n');f.flush();os.fsync(f.fileno())
        os.replace(tmp,p)
    finally:tmp.unlink(missing_ok=True)


def run_step(marker, key, action, output_identity, *, inputs_identity=None):
    """action 결과가 pass/complete이고 작업 중 입출력이 안정적일 때만 완료 확정함."""
    marker=Path(marker);started=time.monotonic()
    before_inputs=inputs_identity() if inputs_identity else None
    identity=digest({'key':key,'inputs':before_inputs})
    try:old=json.loads(marker.read_text())
    except (OSError,ValueError):old={}
    if old.get('schema')=='mfm_pipeline_step_v1' and old.get('status')=='pass' and old.get('key')==identity:
        if old.get('output')==output_identity():
            print(f'[CHECKPOINT] reuse step={marker.stem}',flush=True)
            return old['result']
    atomic_json(marker,{'schema':'mfm_pipeline_step_v1','status':'running','key':identity,'started':time.time()})
    try:
        result=action()
        if not isinstance(result,dict):raise ValueError('단계 결과 형식 오류')
        if result.get('status') not in {'pass','complete'}:
            atomic_json(marker,{'schema':'mfm_pipeline_step_v1','status':'failed','key':identity,'result':result})
            return result
        if inputs_identity and before_inputs!=inputs_identity():raise ValueError('처리 중 입력 지문 변경')
        output=output_identity()
        atomic_json(marker,{'schema':'mfm_pipeline_step_v1','status':'pass','key':identity,
                            'output':output,'result':result,'completed':time.time()})
        return result
    except BaseException as exc:
        atomic_json(marker,{'schema':'mfm_pipeline_step_v1','status':'failed','key':identity,'error':str(exc)})
        raise
    finally:
        marker.parent.mkdir(parents=True,exist_ok=True)
        event={'step':marker.stem,'key':identity,'seconds':round(time.monotonic()-started,6),'time':time.time()}
        fd=os.open(marker.parent/'events.jsonl',os.O_WRONLY|os.O_CREAT|os.O_APPEND,0o600)
        try:os.write(fd,(json.dumps(event)+'\n').encode())
        finally:os.close(fd)
