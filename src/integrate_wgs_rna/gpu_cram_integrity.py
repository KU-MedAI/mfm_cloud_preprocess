"""CRAM3 구조·header CRC는 CPU, 압축 block의 IEEE CRC32는 GPU 배치로 검사함."""
from pathlib import Path
import struct,zlib
from .gpu_nvcomp import NVComp


def scan(path,cp):
    from handoff.alignment_integrity import _take,_itf8,_ltf8,AlignmentCorruption,CRAM_EOF
    from .gpu_io_stream import prefetch
    codec=NVComp(cp);size=Path(path).stat().st_size;summary={}
    def batches():
        containers=count=0;pending=[];pending_bytes=0
        try:
            with Path(path).open('rb',buffering=8<<20) as stream:
                definition=_take(stream,26)
                if definition[:4]!=b'CRAM' or definition[4]!=3 or definition[5] not in {0,1}:
                    raise ValueError('물리 CRC 검사는 CRAM3.0/3.1만 지원합니다.')
                while stream.tell()<size:
                    start=stream.tell();header=bytearray();length=struct.unpack('<i',_take(stream,4,header))[0]
                    for _ in range(4):_itf8(stream,header)
                    _ltf8(stream,header);_ltf8(stream,header)
                    blocks=_itf8(stream,header);landmarks=_itf8(stream,header)
                    if not 0<=landmarks<=100000:raise AlignmentCorruption('CRAM landmark 구조 오류',start)
                    for _ in range(landmarks):_itf8(stream,header)
                    stored=struct.unpack('<I',_take(stream,4))[0]
                    if stored!=zlib.crc32(header)&0xffffffff:raise AlignmentCorruption('CRAM container 헤더 CRC 오류',start,len(header)+4)
                    end=stream.tell()+length
                    if length<0 or end>size or not 0<blocks<=100000:raise AlignmentCorruption('CRAM container 크기 오류',start)
                    containers+=1
                    for _ in range(blocks):
                        block_start=stream.tell();wire=bytearray();_take(stream,2,wire);_itf8(stream,wire)
                        compressed=_itf8(stream,wire);_itf8(stream,wire)
                        if compressed>64<<20 or stream.tell()+compressed+4>end:raise AlignmentCorruption('CRAM block 크기 오류',block_start)
                        data=_take(stream,compressed);stored=struct.unpack('<I',_take(stream,4))[0]
                        value=bytes(wire)+data;pending.append((value,stored,block_start,len(value)+4));pending_bytes+=len(value)
                        count+=1
                        if pending_bytes>=32<<20 or len(pending)>=2048:
                            yield pending;pending=[];pending_bytes=0
                    if stream.tell()!=end:raise AlignmentCorruption('CRAM container 경계 오류',start,length)
                if pending:yield pending;pending=[]
                stream.seek(max(0,size-len(CRAM_EOF)))
                if stream.read()!=CRAM_EOF:raise AlignmentCorruption('CRAM EOF 누락',max(0,size-len(CRAM_EOF)))
                summary.update(format=f'CRAM{definition[4]}.{definition[5]}',containers_checked=containers,blocks_checked=count,bytes_checked=size)
        except (AlignmentCorruption,ValueError):
            if pending:yield pending
            raise
    for batch in prefetch(batches()):
        calculated=codec.crc32([r[0] for r in batch])
        for (_,expected,start,length),actual in zip(batch,calculated):
            if int(actual)!=expected:raise AlignmentCorruption(f'CRAM block CRC 오류: 저장={expected:08x},계산={int(actual):08x}',start,length)
    return summary
