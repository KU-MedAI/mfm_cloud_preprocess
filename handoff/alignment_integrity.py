"""BAM 전체 BGZF 압축해독과 CRAM3 전체 물리 블록 CRC를 검증한다."""
from __future__ import annotations
import hashlib
import json
import os
import struct
import tempfile
import time
import zlib
from pathlib import Path

SCHEMA='mfm_alignment_physical_integrity_v1'
CRAM_EOF=bytes.fromhex('0f000000ffffffff0fe0454f4600000000010005bdd94f0001000606010001000100ee63014b')

class AlignmentCorruption(ValueError):
    def __init__(self,message,offset=0,length=None):
        super().__init__(message)
        self.offset=offset
        self.length=length

def fingerprint(path):
    p=Path(path);s=p.stat()
    return {'path':str(p.resolve()),'size':s.st_size,'mtime_ns':s.st_mtime_ns,
            'ctime_ns':s.st_ctime_ns,'inode':s.st_ino}

def _take(stream,length,wire=None):
    offset=stream.tell();data=stream.read(length)
    if len(data)!=length:
        raise AlignmentCorruption('파일 내용이 잘렸습니다.',offset,length)
    if wire is not None:wire.extend(data)
    return data

def _itf8(stream,wire):
    lead=_take(stream,1,wire)[0]
    if lead<128:return lead
    if lead<192:count,mask=1,63
    elif lead<224:count,mask=2,31
    elif lead<240:count,mask=3,15
    else:
        b=_take(stream,4,wire)
        return ((lead&15)<<28)|(b[0]<<20)|(b[1]<<12)|(b[2]<<4)|(b[3]&15)
    return ((lead&mask)<<(count*8))|int.from_bytes(_take(stream,count,wire),'big')

def _ltf8(stream,wire):
    lead=_take(stream,1,wire)[0];count=0
    while count<8 and lead&(128>>count):count+=1
    _take(stream,count,wire)

def scan_bam(path):
    size=Path(path).stat().st_size;blocks=0;decoded=0;last_empty=False
    with Path(path).open('rb',buffering=8<<20) as stream:
        while stream.tell()<size:
            start=stream.tell();header=_take(stream,12)
            if header[:4]!=b'\x1f\x8b\x08\x04':
                raise AlignmentCorruption('BGZF 헤더 오류',start)
            xlen=struct.unpack('<H',header[10:12])[0]
            if xlen>65512:raise AlignmentCorruption('BGZF extra 크기 오류',start)
            extra=_take(stream,xlen);position=0;length=None
            while position+4<=len(extra):
                tag=extra[position:position+2];n=struct.unpack('<H',extra[position+2:position+4])[0]
                if position+4+n>len(extra):raise AlignmentCorruption('BGZF extra 구조 오류',start)
                if tag==b'BC' and n==2:length=struct.unpack('<H',extra[position+4:position+6])[0]+1
                position+=4+n
            if length is None or not 12+xlen+8<=length<=65536 or start+length>size:
                raise AlignmentCorruption('BGZF 블록 크기 오류',start,length)
            block=header+extra+_take(stream,length-12-xlen)
            try:
                raw=zlib.decompress(block,wbits=31)
            except zlib.error as error:
                raise AlignmentCorruption('BGZF 압축/CRC 오류: '+str(error),start,length) from error
            if len(raw)>65536:raise AlignmentCorruption('BGZF 해독 크기 오류',start,length)
            if blocks==0 and raw[:4]!=b'BAM\x01':raise AlignmentCorruption('BAM 식별자 오류',start,length)
            blocks+=1;decoded+=len(raw);last_empty=not raw
        if not blocks or not last_empty:raise AlignmentCorruption('BAM EOF 블록 누락',max(0,size-65536))
    return {'format':'BAM/BGZF','blocks_checked':blocks,'decoded_bytes':decoded,'bytes_checked':size}

def scan_cram(path):
    size=Path(path).stat().st_size;containers=0;blocks=0
    with Path(path).open('rb',buffering=8<<20) as stream:
        definition=_take(stream,26)
        if definition[:4]!=b'CRAM' or definition[4] != 3 or definition[5] not in {0,1}:
            raise ValueError('물리 CRC 검사는 CRAM3.0/3.1만 지원합니다.')
        while stream.tell()<size:
            start=stream.tell();header=bytearray();length=struct.unpack('<i',_take(stream,4,header))[0]
            for _ in range(4):_itf8(stream,header)
            _ltf8(stream,header);_ltf8(stream,header)
            count=_itf8(stream,header);landmarks=_itf8(stream,header)
            if not 0<=landmarks<=100000:raise AlignmentCorruption('CRAM landmark 구조 오류',start)
            for _ in range(landmarks):_itf8(stream,header)
            stored=struct.unpack('<I',_take(stream,4))[0]
            if stored != zlib.crc32(header)&0xffffffff:
                raise AlignmentCorruption('CRAM container 헤더 CRC 오류',start,len(header)+4)
            end=stream.tell()+length
            if length<0 or end>size or not 0<count<=100000:
                raise AlignmentCorruption('CRAM container 크기 오류',start)
            containers+=1
            for _ in range(count):
                block_start=stream.tell();wire=bytearray()
                _take(stream,2,wire);_itf8(stream,wire)
                compressed=_itf8(stream,wire);_itf8(stream,wire)
                if compressed>64<<20 or stream.tell()+compressed+4>end:
                    raise AlignmentCorruption('CRAM block 크기 오류',block_start)
                data=_take(stream,compressed);stored=struct.unpack('<I',_take(stream,4))[0]
                calculated=zlib.crc32(data,zlib.crc32(wire))&0xffffffff
                if stored!=calculated:
                    raise AlignmentCorruption(f'CRAM block CRC 오류: 저장={stored:08x},계산={calculated:08x}',
                                              block_start,len(wire)+compressed+4)
                blocks+=1
            if stream.tell()!=end:raise AlignmentCorruption('CRAM container 경계 오류',start,length)
        stream.seek(max(0,size-len(CRAM_EOF)))
        if stream.read()!=CRAM_EOF:raise AlignmentCorruption('CRAM EOF 누락',max(0,size-len(CRAM_EOF)))
    return {'format':f'CRAM{definition[4]}.{definition[5]}','containers_checked':containers,
            'blocks_checked':blocks,'bytes_checked':size}

def scan(path):
    with Path(path).open('rb') as stream:magic=stream.read(4)
    if magic==b'CRAM':return scan_cram(path)
    if magic[:2]==b'\x1f\x8b':return scan_bam(path)
    raise AlignmentCorruption('BAM/CRAM 파일 식별 실패',0)

def validate(path,report):
    """동일 파일·동일 검사 코드 fingerprint일 때만 전체 검사 결과를 재사용한다."""
    report=Path(report);before=fingerprint(path)
    validator=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    identity={'schema':SCHEMA,'input':before,'validator_sha256':validator}
    if report.is_file() and time.time_ns()-max(before['mtime_ns'],before['ctime_ns'])>2_000_000_000:
        try:
            old=json.loads(report.read_text())
            if old.get('status')=='pass' and all(old.get(k)==v for k,v in identity.items()):return old
        except (OSError,ValueError,TypeError):pass
    payload=dict(identity,status='checking',started=time.time())
    error=None
    try:
        payload.update(scan(path),status='pass')
    except Exception as failure:
        error=failure;payload.update(status='failed',error=str(failure),
            error_offset=getattr(failure,'offset',None),error_bytes=getattr(failure,'length',None))
    if fingerprint(path)!=before:raise RuntimeError('무결성 검사 중 입력 파일 변경')
    payload['finished']=time.time();report.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile('w',dir=report.parent,delete=False,encoding='utf-8') as stream:
        json.dump(payload,stream,ensure_ascii=False,indent=2);stream.write('\n');tmp=stream.name
    os.replace(tmp,report)
    if error:raise error
    return payload

