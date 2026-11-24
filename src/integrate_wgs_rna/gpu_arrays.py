"""Coverage·CIGAR QC·4bit 패킹의 선택적 CuPy 계산. 저장 형식은 유지함."""
from contextlib import contextmanager
import atexit,os,resource,threading
from pathlib import Path
import numpy as np

_lease=None
_owner=None
_cp=None
_device=None
_lock=threading.RLock()
_counts={}
_address_space={}
FUSED_DEPTH_QC=True

def enabled(stage):
    raw=os.environ.get('MFM_'+stage.upper()+'_BACKEND','cpu')
    if raw not in {'cpu','cupy'}:raise ValueError(stage+' backend는 cpu|cupy여야 함')
    return raw=='cupy'

def _pci(value):
    if isinstance(value,bytes):value=value.decode()
    domain,bus,tail=value.lower().split(':')
    return f'{int(domain,16):x}:{int(bus,16):x}:{tail}'

def _close():
    global _lease
    if _lease is not None and _owner==os.getpid():
        _lease.close();_lease=None

@contextmanager
def session(stage, needed_bytes=0):
    global _lease,_owner,_cp,_device
    with _lock:
        if _owner is not None and _owner!=os.getpid():
            raise RuntimeError('CUDA 초기화 후 fork한 worker는 지원하지 않음. spawn 방식 사용 필요함')
        if _lease is None:
            from .gpu_pool import acquire,GIB
            reservation=int(os.environ.get('MFM_GPU_ARRAY_VRAM_BYTES',str(4*GIB)))
            if needed_bytes>reservation:raise MemoryError('GPU 계산 예상 VRAM이 worker 예약을 초과함')
            lease=acquire(reservation,purpose='array:'+stage)
            previous_cuda=os.environ.get('CUDA_VISIBLE_DEVICES')
            # CuPy 초기화 전에 배정된 한 장으로 좁혀 GPU0의 불필요한 context를 방지함.
            os.environ['CUDA_VISIBLE_DEVICES']=lease.device['uuid']
            soft,hard=resource.getrlimit(resource.RLIMIT_AS)
            if soft!=resource.RLIM_INFINITY and hard==resource.RLIM_INFINITY:
                resource.setrlimit(resource.RLIMIT_AS,(resource.RLIM_INFINITY,hard))
            try:
                try:import cupy as cp
                except ImportError as e:raise RuntimeError('GPU 계산은 requirements-gpu.txt 설치가 필요함') from e
                for index in range(cp.cuda.runtime.getDeviceCount()):
                    if _pci(cp.cuda.runtime.deviceGetPCIBusId(index))==_pci(lease.device['pci']):
                        _device=index;break
                else:raise RuntimeError('예약 GPU가 CUDA 공개 장치에 없음')
                with cp.cuda.Device(_device):
                    cp.get_default_memory_pool().set_limit(size=max(GIB,reservation-GIB))
                    cp.zeros(1,dtype=cp.uint8).get()
                if soft!=resource.RLIM_INFINITY:
                    vm=int(next(line.split()[1] for line in Path('/proc/self/status').read_text().splitlines()
                                if line.startswith('VmSize:')))*1024
                    target=vm+soft
                    if hard!=resource.RLIM_INFINITY:target=min(target,hard)
                    resource.setrlimit(resource.RLIMIT_AS,(target,hard))
                    _address_space.update(previous_limit_bytes=soft,cuda_initialized_vm_bytes=vm,
                                          effective_virtual_limit_bytes=target,rss_guard_unchanged=True)
            except BaseException:
                if previous_cuda is None:os.environ.pop('CUDA_VISIBLE_DEVICES',None)
                else:os.environ['CUDA_VISIBLE_DEVICES']=previous_cuda
                lease.close()
                if soft!=resource.RLIM_INFINITY:resource.setrlimit(resource.RLIMIT_AS,(soft,hard))
                raise
            _cp,_lease,_owner=cp,lease,os.getpid();atexit.register(_close)
        if needed_bytes>_lease.reserved:raise MemoryError('GPU 계산 예상 VRAM이 worker 예약을 초과함')
        _counts[stage]=_counts.get(stage,0)+1
        with _cp.cuda.Device(_device):
            yield _cp
            _cp.cuda.get_current_stream().synchronize()

def execution_metadata():
    return {'requested':{k:os.environ.get('MFM_'+k.upper()+'_BACKEND','cpu') for k in ['coverage','qc','zarr','io']},
            'gpu_uuid':_lease.device['uuid'] if _lease is not None and _owner==os.getpid() else None,
            'gpu_calls':dict(_counts),'counter_scope':'process_lifetime','worker_pid':os.getpid(),'address_space':dict(_address_space),
            'zarr_storage':'CPU Blosc/Zstd; GPU packing/reductions only'}

def depth_from_difference(diff):
    with session('coverage',int(np.asarray(diff).nbytes)*5) as cp:
        values=cp.cumsum(cp.asarray(diff[:-1],dtype=cp.int64),dtype=cp.int64)
        if bool(cp.any(values<0)) or (values.size and int(cp.max(values))>np.iinfo(np.uint32).max):
            raise OverflowError('관측 depth가 uint32 범위를 벗어남')
        return cp.asnumpy(values.astype(cp.uint32))

def depth_from_intervals(starts,ends,length,*,return_metrics=False):
    starts=np.asarray(starts,dtype=np.int64);ends=np.asarray(ends,dtype=np.int64)
    if starts.shape!=ends.shape or np.any(starts<0) or np.any(ends>length) or np.any(ends<starts):
        raise ValueError('coverage interval 범위 오류')
    with session('coverage',(int(length)+1)*40+starts.nbytes*6) as cp:
        delta=cp.zeros(int(length)+1,dtype=cp.int64)
        cp.add.at(delta,cp.asarray(starts),1);cp.add.at(delta,cp.asarray(ends),-1)
        values=cp.cumsum(delta[:-1],dtype=cp.int64)
        if values.size and int(values.max())>np.iinfo(np.uint32).max:raise OverflowError('depth uint32 초과')
        result=cp.asnumpy(values.astype(cp.uint32))
        if return_metrics:
            _counts['qc']=_counts.get('qc',0)+1
            return result,_depth_metrics(cp,values,length)
        return result

def _depth_metrics(cp,actual,actual_length):
    total=int(cp.sum(actual,dtype=cp.uint64));nonzero=int(cp.count_nonzero(actual))
    return {'depth_sum':total,'depth_mean':total/actual_length,
            'depth_zero_fraction':(actual_length-nonzero)/actual_length,'depth_max':int(actual.max()),
            'depth_nonzero_count':nonzero}

def depth_metrics(values,actual_length):
    values=np.asarray(values)
    if actual_length<=0 or actual_length>values.size:raise ValueError('관측 구간 길이 오류')
    if values.dtype.kind not in 'ui':raise ValueError('음수 관측 깊이')
    with session('qc',values[:actual_length].nbytes*3) as cp:
        actual=cp.asarray(values[:actual_length])
        if bool(cp.any(actual<0)):raise ValueError('음수 관측 깊이')
        result=_depth_metrics(cp,actual,actual_length)
        result.pop('depth_nonzero_count')
        return result

def pack_codes(raw,lookup,valid,odd_positions,invalid_policy):
    """입력은 레코드별 홀수 padding이 삽입된 바이트열임. 오류 위치도 보존함."""
    with session('zarr',np.asarray(raw).nbytes*12+4096) as cp:
        data=cp.asarray(raw);codes=cp.asarray(lookup)[data]
        invalid=~cp.asarray(valid)[codes]
        if len(odd_positions):invalid[cp.asarray(odd_positions,dtype=cp.int64)]=False
        first=-1;hist={}
        if bool(cp.any(invalid)):
            if invalid_policy=='fail':first=int(cp.flatnonzero(invalid)[0])
            else:
                values,frequencies=cp.unique(data[invalid],return_counts=True)
                for v,n in zip(cp.asnumpy(values),cp.asnumpy(frequencies)):
                    key=chr(int(v)).upper();hist[key]=hist.get(key,0)+int(n)
                codes[invalid]=int(lookup[ord('N')])
        packed=(codes[0::2]<<4)|codes[1::2]
        return cp.asnumpy(packed),first,hist

def cigar_summary(offsets,sizes,ops,read_lengths,ref_lengths,first_read):
    """GPU에서 CIGAR 유효성 및 read별 query/reference 소비 길이를 검증함."""
    size=sum(np.asarray(v).nbytes for v in [offsets,sizes,ops,read_lengths,ref_lengths])
    if np.any(np.asarray(offsets)<0) or np.any(np.asarray(offsets)+np.asarray(sizes)>len(ops)):
        raise ValueError('full validation truncated CIGAR')
    with session('qc',size*8+4096) as cp:
        starts=cp.asarray(offsets,dtype=cp.int64);counts_per=cp.asarray(sizes,dtype=cp.int64)
        ends=starts+counts_per;packed=cp.asarray(ops,dtype=cp.uint32)
        codes=packed&15;lengths=(packed>>4).astype(cp.int64)
        def first_bad(mask):
            return int(cp.flatnonzero(mask)[0]) if bool(cp.any(mask)) else None
        bad=first_bad(counts_per<=0)
        if bad is not None:raise ValueError(f'full validation empty CIGAR for read_id={first_read+bad}')
        for allowed,expected,label in [([0,1,2,3,6,7,8],None,'invalid stored CIGAR op'),
                    ([0,1,7,8],read_lengths,'CIGAR/query mismatch'),
                    ([0,2,3,7,8],ref_lengths,'CIGAR/reference mismatch')]:
            weights=(~(cp.isin(codes,cp.asarray(allowed))&(lengths>0))).astype(cp.int64) if expected is None else lengths*cp.isin(codes,cp.asarray(allowed))
            prefix=cp.concatenate([cp.zeros(1,dtype=cp.int64),cp.cumsum(weights,dtype=cp.int64)])
            observed=prefix[ends]-prefix[starts]
            bad=first_bad(observed!=cp.asarray(expected) if expected is not None else observed!=0)
            if bad is not None:raise ValueError(f'full validation {label} read_id={first_read+bad}')
    return len(sizes)
