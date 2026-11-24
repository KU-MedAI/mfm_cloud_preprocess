"""samtools의 BAM record 순서를 유지하며 BGZF 압축만 GPU에서 수행함."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,hashlib,json,os,struct,sys,time,zlib
from .gpu_nvcomp import NVComp

EOF=bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')
HEADER=bytes.fromhex('1f8b08040000000000ff060042430200')


def blocks(stream):
    count=0;last_empty=False
    while True:
        header=stream.read(12)
        if not header:break
        if len(header)!=12 or header[:4]!=b'\x1f\x8b\x08\x04':raise ValueError('BGZF 입력 헤더 손상')
        xlen=struct.unpack('<H',header[10:12])[0]
        if not 6<=xlen<=65512:raise ValueError('BGZF 입력 extra 길이 오류')
        extra=stream.read(xlen)
        if len(extra)!=xlen:raise ValueError('BGZF extra 잘림')
        offset=0;length=None
        while offset+4<=len(extra):
            kind=extra[offset:offset+2];size=struct.unpack_from('<H',extra,offset+2)[0]
            if offset+4+size>len(extra):raise ValueError('BGZF extra 구조 오류')
            if kind==b'BC' and size==2:length=struct.unpack_from('<H',extra,offset+4)[0]+1
            offset+=4+size
        if length is None or not 12+xlen+8<=length<=65536:raise ValueError('BGZF block 길이 오류')
        tail=stream.read(length-12-xlen)
        if len(tail)!=length-12-xlen:raise ValueError('BGZF block 잘림')
        raw=zlib.decompress(header+extra+tail,wbits=31)
        if len(raw)>65536:raise ValueError('BGZF 해독 크기 초과')
        if count==0 and raw[:4]!=b'BAM\x01':raise ValueError('BAM 식별자 오류')
        count+=1;last_empty=not raw
        if raw:yield raw
    if not count or not last_empty:raise ValueError('BAM EOF block 누락')


def verify_partition(rows):
    for raw,(compressed,crc) in rows:
        if zlib.decompress(compressed,wbits=-15)!=raw or zlib.crc32(raw)!=crc:
            raise ValueError('GPU 압축 결과의 CPU 독립 해독/CRC 전수 대조 불일치')


def compress_stream(source,target,codec,*,batch_chunks=512,verify_workers=1):
    from .gpu_io_stream import prefetch
    if not 1<=batch_chunks<=1024 or verify_workers not in {1,2}:raise ValueError('GPU IO 배치/CPU 검증 자원 오류')
    digest=hashlib.sha256();count=total=written=fallback=0
    def batches():
        pending=[];carry=b''
        for raw in blocks(source):
            digest.update(raw)
            data=carry+raw;position=0
            while len(data)-position>=60000:
                pending.append(data[position:position+60000]);position+=60000
                if len(pending)>=batch_chunks:yield pending;pending=[]
            carry=data[position:]
        if carry:pending.append(carry)
        if pending:yield pending
    with ThreadPoolExecutor(max_workers=verify_workers) as pool:
        waiting=[]
        def write_verified():
            nonlocal count,total,written,fallback
            rows,futures=waiting.pop(0)
            for future in futures:future.result()
            wire=bytearray()
            for raw,(compressed,crc) in rows:
                if len(compressed)>65510:
                    c=zlib.compressobj(level=0,wbits=-15);compressed=c.compress(raw)+c.flush();fallback+=1
                size=18+len(compressed)+8
                if size>65536:raise ValueError('BGZF 출력 크기 초과')
                wire+=HEADER+struct.pack('<H',size-1)+compressed+struct.pack('<II',crc,len(raw))
                count+=1;total+=len(raw)
            target.write(wire);written+=len(wire)
        for pending in prefetch(batches()):
            encoded=codec.compress(pending)
            if len(encoded)!=len(pending):raise ValueError('GPU 압축 chunk 개수 불일치')
            rows=list(zip(pending,encoded))
            futures=[pool.submit(verify_partition,rows[i::verify_workers]) for i in range(verify_workers)]
            waiting.append((rows,futures))
            if len(waiting)>=2:write_verified()
        while waiting:write_verified()
        target.write(EOF);written+=len(EOF)
    return {'schema':'mfm_gpu_bgzf_compression_v1','decoded_bytes':total,'compressed_bytes':written,'blocks_written':count+1,'raw_sha256':digest.hexdigest(),'cpu_full_decode_and_byte_compare':True,'cpu_crc_all_chunks_equal':True,'stored_fallback_blocks':fallback,'bounded_input_batches':2,'bounded_verify_batches':2,'record_order':'입력 순서 보존함; 정렬은 기존 samtools 사용함'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--metadata',type=Path,required=True)
    parser.add_argument('--verify-workers',type=int,default=1);args=parser.parse_args()
    temporary=args.output.with_name(args.output.name+'.gpu_io.part.'+str(os.getpid()))
    os.environ['MFM_IO_BACKEND']='nvcomp'
    started=time.monotonic()
    try:
        from .gpu_arrays import session,execution_metadata
        with session('io') as cp,temporary.open('wb',buffering=8<<20) as output:
            result=compress_stream(sys.stdin.buffer,output,NVComp(cp),verify_workers=args.verify_workers)
            result.update(elapsed_seconds=time.monotonic()-started,compute_execution=execution_metadata(),nvcomp_version='6.0.0.12')
        os.replace(temporary,args.output)
        args.metadata.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
        print('[GPU-BGZF] '+json.dumps(result,ensure_ascii=False),file=sys.stderr,flush=True)
    except BaseException:
        temporary.unlink(missing_ok=True);raise


if __name__=='__main__':main()
