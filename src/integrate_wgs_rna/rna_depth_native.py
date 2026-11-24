"""한 번 읽은 CIGAR 차분을 CPU/GPU로 변환하고 두 출력에 같은 buffer를 전달함."""
from pathlib import Path
import hashlib,json,os,resource,subprocess,sys,tempfile,time
import numpy as np
import zarr
from .record_stats_parallel import native_tool


def gpu_difference(source,output,length):
    # 파일 읽기는 GPU 예약 전에, 압축 저장은 이 프로세스가 종료된 뒤 부모에서 함.
    expected=(length+1)*40
    os.environ['MFM_GPU_ARRAY_VRAM_BYTES']=str(max(expected+(1<<30),int(os.environ.get('MFM_GPU_ARRAY_VRAM_BYTES',str(4<<30)))))
    # CUDA 장치 메모리 mapping도 RLIMIT_AS에 포함됨. CPU RSS 보호와 분리해 예측 VRAM만큼 주소 여유를 둠.
    reservation=int(os.environ['MFM_GPU_ARRAY_VRAM_BYTES'])
    soft,hard=resource.getrlimit(resource.RLIMIT_AS)
    if soft!=resource.RLIM_INFINITY:
        virtual_limit=soft+2*reservation
        if hard!=resource.RLIM_INFINITY:virtual_limit=min(virtual_limit,hard)
        resource.setrlimit(resource.RLIMIT_AS,(virtual_limit,hard))
        print(f'[RNA-GPU-VAS] previous={soft} gpu_reservation={reservation} adjusted={virtual_limit}; CPU RSS guard 유지함',flush=True)
    diff=np.fromfile(source,dtype='<i8')
    if diff.size!=length+1:raise ValueError('잘린 RNA 차분')
    from .gpu_arrays import depth_from_difference
    values=depth_from_difference(diff)
    if int(diff.sum(dtype=np.int64))!=0:raise ValueError('RNA 차분 불균형')
    values.astype('<u4',copy=False).tofile(output)


def contig_depth(task):
    contig,length,bam,samtools,store=task[:5];snapshot=Path(task[5]) if len(task)>5 and task[5] else None
    started=time.monotonic();tool=native_tool('rna_depth')
    from numcodecs import blosc
    blosc.set_nthreads(1)
    root=zarr.open_group(store,mode='r+');array=root[contig];chunk=int(array.chunks[0])
    temporary=tempfile.TemporaryDirectory(prefix='rna_depth_',dir=snapshot) if snapshot else tempfile.TemporaryDirectory(prefix='rna_depth_')
    try:
        folder=Path(temporary.name);depth=folder/'depth.bin'
        gpu=os.environ.get('MFM_COVERAGE_BACKEND','cpu')=='cupy'
        logfile=folder/'native.log'
        with logfile.open('wb') as log:
            def run(command,**kwargs):
                try:subprocess.run(command,stderr=log,check=True,**kwargs)
                except subprocess.CalledProcessError as error:
                    log.flush()
                    raise RuntimeError(contig+' RNA native/GPU 실패: '+logfile.read_text(errors='replace')[-8000:]) from error
            if gpu:
                diff=folder/'difference.bin'
                with diff.open('wb') as out:run([str(tool),bam,contig,str(length),'difference'],stdout=out)
                run([sys.executable,'-m','integrate_wgs_rna.rna_depth_native',str(diff),str(depth),str(length)])
                diff.unlink()
            else:
                with depth.open('wb') as out:run([str(tool),bam,contig,str(length),'depth'],stdout=out)
        if depth.stat().st_size!=length*4:raise ValueError('잘린 RNA depth')
        values=np.memmap(depth,mode='r',dtype='<u4',shape=(length,));total=nonzero=0;digest=hashlib.sha256()
        for begin in range(0,length,chunk):
            block=np.asarray(values[begin:min(length,begin+chunk)])
            array[begin:begin+len(block)]=block
            total+=int(block.sum(dtype=np.uint64));nonzero+=int(np.count_nonzero(block));digest.update(memoryview(block).cast('B'))
        del values
        if snapshot:
            path=snapshot/(hashlib.sha256(contig.encode()).hexdigest()+'.bin');os.replace(depth,path)
            stamp=path.stat();(path.with_suffix('.json')).write_text(json.dumps({'contig':contig,'path':str(path),'bytes':stamp.st_size,'mtime_ns':stamp.st_mtime_ns,'sha256':digest.hexdigest()}))
        print(f'[cov:native] {contig} mode={"gpu_prefix" if gpu else "cpu_prefix"} elapsed={time.monotonic()-started:.3f}s',flush=True)
        return total,nonzero,total,digest.hexdigest()
    finally:temporary.cleanup()


if __name__=='__main__':gpu_difference(Path(sys.argv[1]),Path(sys.argv[2]),int(sys.argv[3]))
