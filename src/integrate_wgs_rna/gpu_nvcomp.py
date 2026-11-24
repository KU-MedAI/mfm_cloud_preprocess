"""nvCOMP6의 무손실 Deflate 압축·IEEE CRC32만 호출함. 미검증 압축입력 해독은 하지 않음."""
from pathlib import Path
import ctypes as C
import importlib.metadata
import os
import numpy as np

class DeflateOptions(C.Structure):
    _fields_=[('algorithm',C.c_int),('reserved',C.c_char*60)]
class CRCSpec(C.Structure):
    _fields_=[('poly',C.c_uint32),('initial',C.c_uint32),('ref_in',C.c_bool),('ref_out',C.c_bool),('xor_out',C.c_uint32),('reserved',C.c_char*16)]
class CRCConf(C.Structure):
    _fields_=[('kind',C.c_int),('bytes_per_read',C.c_int32),('blocks_per_message',C.c_int32),('reserved',C.c_char*20)]
class CRCOptions(C.Structure):
    _fields_=[('spec',CRCSpec),('conf',CRCConf),('reserved',C.c_char*64)]


def library_path():
    distribution=importlib.metadata.distribution('nvidia-libnvcomp-cu13')
    if distribution.version!='6.0.0.12':
        raise RuntimeError('검증된 nvidia-libnvcomp-cu13==6.0.0.12 필요함')
    return Path(distribution.locate_file('nvidia/libnvcomp/lib64/libnvcomp.so.6'))


class NVComp:
    def __init__(self,cp):
        self.cp=cp;self.lib=C.CDLL(str(library_path()));p=C.c_void_p;s=C.c_size_t
        signatures={
            'nvcompBatchedDeflateCompressGetMaxOutputChunkSize':[s,DeflateOptions,C.POINTER(s)],
            'nvcompBatchedDeflateCompressGetTempSize':[s,s,DeflateOptions,C.POINTER(s),s,p],
            'nvcompBatchedDeflateCompressAsync':[p,p,s,s,p,s,p,p,DeflateOptions,p,p],
            'nvcompBatchedCRC32Async':[p,p,s,p,CRCOptions,C.c_int,p,p]}
        for name,args in signatures.items():
            f=getattr(self.lib,name);f.argtypes=args;f.restype=C.c_int
        self.algorithm=int(os.environ.get('MFM_GPU_DEFLATE_ALGORITHM','2'))
        if self.algorithm not in {0,1,2,4,5}:raise ValueError('nvCOMP Deflate algorithm 오류')
        self.scratch={}
        self.pack=cp.RawKernel(r"""
extern "C" __global__ void compact_chunks(const unsigned char* input, const unsigned long long* sizes,
 const unsigned long long* offsets, unsigned char* output, const unsigned long long stride) {
 unsigned int row=blockIdx.x;
 for(unsigned long long i=threadIdx.x;i<sizes[row];i+=blockDim.x) output[offsets[row]+i]=input[row*stride+i];
}
""",'compact_chunks')
    def check(self,status):
        if status:raise RuntimeError('nvCOMP API 오류: '+str(status))
    def _status(self,status,label):
        values=self.cp.asnumpy(status)
        if np.any(values):raise RuntimeError(label+' GPU chunk 상태 오류: '+str(values[values!=0][:8]))
    def _crc(self,pointers,lengths,max_bytes):
        cp=self.cp;n=len(lengths);result=cp.zeros(n,dtype=cp.uint32);status=cp.zeros(n,dtype=cp.int32)
        big=max_bytes>65536
        conf=CRCConf(1 if big else 0,4,min(32,1<<max(0,(((max_bytes+65535)//65536)-1).bit_length())) if big else 1)
        options=CRCOptions(CRCSpec(0x04C11DB7,0xffffffff,True,True,0xffffffff),conf)
        self.check(self.lib.nvcompBatchedCRC32Async(pointers.data.ptr,lengths.data.ptr,n,result.data.ptr,options,0,status.data.ptr,cp.cuda.get_current_stream().ptr))
        cp.cuda.get_current_stream().synchronize();self._status(status,'CRC32')
        return cp.asnumpy(result)
    def crc32(self,chunks):
        if not chunks:return np.empty(0,dtype='u4')
        cp=self.cp;lengths=np.asarray([len(v) for v in chunks],dtype='u8');offsets=np.r_[np.uint64(0),np.cumsum(lengths,dtype=np.uint64)[:-1]]
        host=np.frombuffer(b''.join(chunks),dtype='u1');data=cp.asarray(host)
        points=data.data.ptr+cp.asarray(offsets);sizes=cp.asarray(lengths)
        return self._crc(points,sizes,int(lengths.max(initial=0)))
    def compress(self,chunks):
        if not chunks:return []
        if any(not 0<len(v)<=60000 for v in chunks):raise ValueError('GPU BGZF chunk 길이 오류')
        cp=self.cp;n=len(chunks);stride=(max(map(len,chunks))+7)//8*8
        host=np.zeros((n,stride),dtype='u1')
        for i,value in enumerate(chunks):host[i,:len(value)]=np.frombuffer(value,dtype='u1')
        data=cp.asarray(host);lengths=cp.asarray([len(v) for v in chunks],dtype=cp.uint64)
        points=data.data.ptr+cp.arange(n,dtype=cp.uint64)*stride
        option=DeflateOptions(self.algorithm);temp=C.c_size_t();maximum=C.c_size_t();stream=cp.cuda.get_current_stream().ptr
        self.check(self.lib.nvcompBatchedDeflateCompressGetMaxOutputChunkSize(stride,option,C.byref(maximum)))
        self.check(self.lib.nvcompBatchedDeflateCompressGetTempSize(n,stride,option,C.byref(temp),sum(map(len,chunks)),stream))
        output_stride=(maximum.value+7)//8*8
        needed=data.nbytes+n*output_stride+temp.value
        if needed>2*1024**3:raise MemoryError('GPU BGZF 배치 예상 VRAM2GiB 초과')
        output=cp.empty((n,output_stride),dtype='u1');outpoints=output.data.ptr+cp.arange(n,dtype=cp.uint64)*output_stride
        sizes=cp.zeros(n,dtype=cp.uint64);status=cp.zeros(n,dtype=cp.int32)
        scratch=cp.empty(max(1,temp.value),dtype='u1')
        self.check(self.lib.nvcompBatchedDeflateCompressAsync(points.data.ptr,lengths.data.ptr,stride,n,scratch.data.ptr,temp.value,outpoints.data.ptr,sizes.data.ptr,option,status.data.ptr,stream))
        cp.cuda.get_current_stream().synchronize();self._status(status,'Deflate')
        lengths_out=cp.asnumpy(sizes)
        if np.any(lengths_out>maximum.value) or np.any(lengths_out==0):raise RuntimeError('GPU 압축 결과 길이 오류')
        checksums=self._crc(points,lengths,max(map(len,chunks)))
        # 가변 길이를 GPU에서 모아 실제 압축 결과만 호스트로 전송함.
        gpu_sizes=cp.asarray(lengths_out);offsets=cp.cumsum(gpu_sizes,dtype=cp.uint64)-gpu_sizes
        packed=cp.empty(int(lengths_out.sum()),dtype='u1')
        self.pack((n,),(256,),(output,gpu_sizes,offsets,packed,np.uint64(output_stride)))
        raw=cp.asnumpy(packed).tobytes();result=[];position=0
        for size,crc in zip(lengths_out,checksums):
            end=position+int(size);result.append((raw[position:end],int(crc)));position=end
        return result
