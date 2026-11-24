"""기존 canonical JSON SHA256을 bounded 네이티브 배치로 계산함."""
import ctypes
from functools import lru_cache
import os
import numpy as np
from ..record_stats_parallel import native_tool, NativeUnavailable
from .trace_decode import decode_reads, read_digest, segment_digest


@lru_cache(maxsize=1)
def _library():
    try:
        lib=ctypes.CDLL(str(native_tool('trace_hash')))
    except (NativeUnavailable,OSError):
        return None
    fn=lib.trace_hash_batch
    pointer=ctypes.c_void_p
    fn.argtypes=[pointer,ctypes.c_int64,pointer,pointer,pointer,ctypes.c_int64,
                 ctypes.c_int64,ctypes.c_int,pointer]
    fn.restype=ctypes.c_int
    lib.trace_segment_hash_batch.argtypes=[pointer,pointer,ctypes.c_int64,pointer]
    lib.trace_segment_hash_batch.restype=ctypes.c_int
    return lib


def hash_backend():
    return 'python' if os.environ.get('MFM_TRACE_HASH_BACKEND')=='python' or _library() is None else 'native_openssl_sha256'


def digests(bases, offsets, metadata, cigar, *, packed):
    """metadata 열: HP, PS, start, end, type, CIGAR offset/length, query length임."""
    mode=os.environ.get('MFM_TRACE_HASH_BACKEND','auto')
    if mode not in {'auto','native','python'}:raise ValueError('trace hash backend는 auto/native/python임')
    bases=np.ascontiguousarray(bases,dtype=np.uint8)
    offsets=np.ascontiguousarray(offsets,dtype=np.int64)
    metadata=np.ascontiguousarray(metadata,dtype=np.int64)
    cigar=np.ascontiguousarray(cigar,dtype=np.int64)
    if metadata.ndim!=2 or metadata.shape[1]!=8 or offsets.shape!=(len(metadata)+1,) or bases.ndim!=1 or cigar.ndim!=1:
        raise ValueError('trace hash batch shape 오류')
    n=len(metadata);out=np.empty(n,dtype='V32')
    lib=None if mode=='python' else _library()
    if mode=='native' and lib is None:raise NativeUnavailable('trace hash native 빌드 불가')
    if lib is not None:
        rc=lib.trace_hash_batch(bases.ctypes.data,len(bases),offsets.ctypes.data,metadata.ctypes.data,
                                cigar.ctypes.data,len(cigar),n,int(packed),out.ctypes.data)
        if rc:raise ValueError('trace hash batch 범위/염기 오류: '+str(rc))
        return out
    if np.any(offsets<0) or np.any(offsets>len(bases)) or np.any(offsets[1:]<offsets[:-1]):
        raise ValueError('trace hash offsets 오류')
    if packed:sequences=decode_reads(bases,offsets,metadata[:,7])
    else:sequences=[bytes(bases[offsets[i]:offsets[i]+metadata[i,7]]).decode('ascii') for i in range(n)]
    for i,(seq,row) in enumerate(zip(sequences,metadata)):
        hp,ps,start,end,kind,c,count,length=map(int,row)
        if c<0 or count<0 or c+count>len(cigar) or length<0 or length>(int(offsets[i+1]-offsets[i])*(2 if packed else 1)):
            raise ValueError('trace hash metadata 오류')
        out[i]=np.void(read_digest(seq,hp,ps,start,end,kind,tuple(map(int,cigar[c:c+count]))))
    return out


def source_digests(canonical):
    offsets=[0];bases=bytearray();metadata=[];cigar=[]
    for seq,hp,ps,start,end,kind,ops in canonical:
        raw=seq.encode('ascii');bases.extend(raw);offsets.append(len(bases))
        metadata.append((hp,ps,start,end,kind,len(cigar),len(ops),len(raw)));cigar.extend(ops)
    return digests(np.frombuffer(bases,dtype=np.uint8),offsets,
                   np.asarray(metadata,dtype=np.int64).reshape(-1,8),cigar,packed=False)


def segment_digests(reads,fields):
    ids=np.ascontiguousarray(reads,dtype='V32');fields=np.ascontiguousarray(fields,dtype=np.int64)
    if ids.ndim!=1 or fields.shape!=(len(ids),6):raise ValueError('segment hash batch shape 오류')
    mode=os.environ.get('MFM_TRACE_HASH_BACKEND','auto')
    if mode not in {'auto','native','python'}:raise ValueError('trace hash backend 오류')
    lib=None if mode=='python' else _library();out=np.empty(len(ids),dtype='V32')
    if mode=='native' and lib is None:raise NativeUnavailable('segment hash native 빌드 불가')
    if lib is not None:
        if lib.trace_segment_hash_batch(ids.ctypes.data,fields.ctypes.data,len(ids),out.ctypes.data):raise ValueError('segment hash 실패')
    else:
        for i,(read,row) in enumerate(zip(ids,fields)):out[i]=np.void(segment_digest(read,*map(int,row)))
    return out


def add_hashes(bag,side,kind,values):
    if hasattr(bag,'add_many'):bag.add_many(side,kind,values)
    else:
        for value in values:bag.add(side,kind,bytes(value))


class SegmentBatch:
    def __init__(self,bag):self.bag=bag;self.reads=[];self.fields=[]
    def add(self,read,*fields):
        self.reads.append(read);self.fields.append(fields)
        if len(self.reads)>=8192:self.flush()
    def flush(self):
        if not self.reads:return
        values=segment_digests(np.frombuffer(b''.join(self.reads),dtype='V32'),self.fields)
        add_hashes(self.bag,'src','segment',values);self.reads.clear();self.fields.clear()
