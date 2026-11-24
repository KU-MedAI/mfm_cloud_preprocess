"""손상 입력을 보존하고 별도 복구본의 손상 구간만 원격에서 다시 확보한다."""
from __future__ import annotations
import hashlib
import json
import os
import shutil
import time
from dataclasses import asdict
from pathlib import Path
from .alignment_integrity import AlignmentCorruption, fingerprint, scan, validate
from .http_download import CHUNK_BYTES, error_offsets, identify, read_range

def repair(source,target,url,report,expected_size,expected_etag=None,max_patches=20):
    source=Path(source);target=Path(target);report=Path(report)
    before=fingerprint(source)
    remote=identify(url)
    if remote.size!=int(expected_size) or before['size']!=remote.size:
        raise ValueError('손상 파일/원격 파일/동결 catalog 크기가 다릅니다.')
    if expected_etag and remote.etag!=expected_etag:
        raise ValueError('검토한 원격 버전과 현재 버전이 다릅니다.')
    # 파일 식별/header와 EOF가 같은 원격 자원인지 추가 확인한다.
    with source.open('rb') as stream:
        for start in [0,max(0,remote.size-512)]:
            end=min(remote.size-1,start+511);stream.seek(start)
            if stream.read(end-start+1)!=read_range(remote,start,end):
                raise ValueError('손상 파일의 header/EOF와 원격 자원이 다릅니다.')
    target.parent.mkdir(parents=True,exist_ok=True)
    # 기존 손상 입력과 기존 복구 후보를 덮어쓰지 않는다.
    with source.open('rb') as src,target.open('xb') as dst:
        shutil.copyfileobj(src,dst,8<<20)
    if fingerprint(source)!=before:raise RuntimeError('복사 중 원본 파일 변경')
    payload={'schema':'mfm_alignment_range_recovery_v1','status':'repairing','started':time.time(),
             'source':before,'replacement':str(target),'remote':asdict(remote),'patches':[]}
    report.parent.mkdir(parents=True,exist_ok=True)
    def save():
        report.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n')
    def patch(start,end,reason):
        if len(payload['patches'])>=max_patches:
            raise ValueError('부분 복구 한도 초과: 전체 재확보가 필요합니다.')
        data=read_range(remote,start,end)
        with target.open('r+b') as stream:
            stream.seek(start);previous=stream.read(end-start+1)
            if data==previous:
                raise ValueError('원격과 동일한 구간에서도 손상이 발생합니다. 자동 복구 중단.')
            stream.seek(start);stream.write(data);stream.flush();os.fsync(stream.fileno())
        payload['patches'].append({'offset':start,'bytes':len(data),'reason':reason,
            'changed_bytes':sum(a!=b for a,b in zip(data,previous)),
            'before_sha256':hashlib.sha256(previous).hexdigest(),'after_sha256':hashlib.sha256(data).hexdigest()})
        save()
        print(f'[alignment-recovery] path={target} offset={start} bytes={len(data)} reason={reason}',flush=True)
    try:
        offsets=error_offsets(target)
        payload['xml_error_offsets']=offsets;save()
        # 응답 문자열이 블록 header를 손상시킨 경우에도 먼저 주변 구간을 복구한다.
        windows=sorted({offset//(1<<20)*(1<<20) for offset in offsets})
        for start in windows:patch(start,min(remote.size-1,start+(1<<20)-1),'XML 오류 응답')
        for _ in range(max_patches+1):
            try:
                stats=scan(target)
                break
            except AlignmentCorruption as error:
                start=error.offset
                length=error.length
                if not length or length>CHUNK_BYTES:
                    start=max(0,start-65536);length=1<<20
                patch(start,min(remote.size-1,start+length-1),str(error))
        else:raise ValueError('전체 무결성 검사 통과 실패')
        physical=validate(target,Path(str(report)+'.physical.json'))
        if fingerprint(source)!=before:raise RuntimeError('복구 중 원본 파일 변경')
        payload.update(status='pass',finished=time.time(),physical_validation=physical,
                       network_repair_bytes=sum(p['bytes'] for p in payload['patches']),stats=stats,
                       replacement_fingerprint=fingerprint(target))
        save();return payload
    except Exception as error:
        payload.update(status='failed',finished=time.time(),error=str(error));save();raise

