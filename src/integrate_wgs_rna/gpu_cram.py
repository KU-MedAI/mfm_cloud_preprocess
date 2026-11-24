"""좌표 정렬 BAM을 염색체별 복원·CRAM 인코딩 후 원래 순서로 결합함."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json, os, resource, shutil, subprocess, tempfile, time
import pysam
from .memory_policy import GIB
from .record_stats_parallel import native_tool, NativeUnavailable


def plan(cpus, references):
    # fetch/restore/encode main 세 개와 codec·I/O 두 개, 부모 두 개를 포함함.
    workers=min(int(references),max(0,(int(cpus)-2)//5))
    return {'workers':workers,'cpu_budget':int(cpus),'native_cpu_limit':workers*5+2,
            'estimated_peak_bytes':(4+3*workers)*GIB*11//10}


def encode(source, reference, destination, read_group, sample, duplicate_mode, *, cpus, work):
    started=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    source,reference,destination=map(Path,(source,reference,destination))
    samtools=shutil.which('samtools')
    if not samtools:raise RuntimeError('samtools 없음')
    destination.parent.mkdir(parents=True,exist_ok=True)
    directory=Path(tempfile.mkdtemp(prefix='cram_shards_',dir=work))
    try:
        tool=None
        if os.environ.get('MFM_GPU_TAG_RESTORE','auto')!='python':
            try:tool=native_tool('restore_tags')
            except NativeUnavailable:
                if os.environ.get('MFM_GPU_TAG_RESTORE')=='native':raise
        with pysam.AlignmentFile(str(source),'rb') as bam:
            header=bam.header.to_dict();references=list(bam.references)
        if header.get('HD',{}).get('SO')!='coordinate':raise ValueError('좌표 정렬되지 않은 BAM임')
        cfg=plan(cpus,len(references)+1)
        if tool and cfg['workers']>=2:
            if not Path(str(source)+'.bai').exists() and not Path(str(source)+'.csi').exists():
                subprocess.run([samtools,'index','-@',str(max(0,cpus-2)),str(source)],check=True)
            with pysam.AlignmentFile(str(source),'rb') as bam:
                counts={s.contig:s.mapped+s.unmapped for s in bam.get_index_statistics()}
                regions=[r for r in references if counts.get(r,0)]
                counts['*']=bam.nocoordinate
                if bam.nocoordinate:regions.append('*')
            cfg=plan(cpus,len(regions))
        if tool and cfg['workers']>=2:
            rg=dict(field.split(':',1) for field in read_group.replace('\\t','\t').split('\t') if ':' in field)
            header['RG']=[rg];header_path=directory/'header.sam'
            header_path.write_text(str(pysam.AlignmentHeader.from_dict(header)))
            tasks=list(enumerate(regions))
            def shard(task):
                index,region=task;part=directory/f'{index:05d}.cram';processes=[]
                with (directory/f'{index:05d}.log').open('wb') as log:
                    try:
                        fetch=subprocess.Popen([samtools,'view','--no-PG','-u',str(source),region],stdout=subprocess.PIPE,stderr=log);processes.append(fetch)
                        restore=subprocess.Popen([str(tool),'-','-',str(header_path),rg['ID'],duplicate_mode or 'unset','wb0','0'],stdin=fetch.stdout,stdout=subprocess.PIPE,stderr=log);processes.append(restore);fetch.stdout.close()
                        writer=subprocess.Popen([samtools,'view','--no-PG','-C','-@','1','-T',str(reference),'-o',str(part),'-'],stdin=restore.stdout,stderr=log);processes.append(writer);restore.stdout.close()
                        codes=[p.wait() for p in reversed(processes)]
                        if any(codes):raise RuntimeError('염색체 CRAM 실패: '+region+' '+str(codes))
                        subprocess.run([samtools,'quickcheck',str(part)],check=True,stderr=log)
                    finally:
                        for p in processes:
                            if p.poll() is None:p.terminate()
                        for p in processes:p.wait()
                numbers=[int(line.split('=',1)[1]) for line in (directory/f'{index:05d}.log').read_text().splitlines() if line.startswith('MFM_RESTORE_RECORDS=')]
                if numbers!=[counts[region]]:raise RuntimeError('BAM index와 실제 shard read 수 불일치: '+region+' '+str(numbers))
                return part
            with ThreadPoolExecutor(max_workers=cfg['workers']) as pool:parts=list(pool.map(shard,tasks))
            subprocess.run([samtools,'cat','--no-PG','-o',str(destination),*map(str,parts)],check=True)
            backend='native_chromosome_cram'
        else:
            decode=min(8,max(0,(cpus-3)//4));codec=max(0,cpus-decode-3)
            cfg.update(workers=1,native_cpu_limit=cpus,estimated_peak_bytes=(4*GIB+codec*(GIB//2)+decode*64*(1<<20))*11//10)
            with (directory/'serial.log').open('wb') as log:
                writer=subprocess.Popen([samtools,'view','--no-PG','-C','-@',str(codec),'-T',str(reference),'-o',str(destination),'-'],stdin=subprocess.PIPE,stderr=log)
                # fd 중복 변경 없이 기존 네이티브 복원을 독립 자식에서 실행함.
                code='from integrate_wgs_rna.gpu_wgs import restore_bam;import sys;a=sys.argv[1:6];a[3]=a[3] or None;restore_bam(*a)'
                import sys
                env=dict(os.environ,MFM_GPU_RESTORE_DECODE_THREADS=str(decode))
                try:
                    subprocess.run([sys.executable,'-c',code,str(source),read_group,sample,duplicate_mode or '', '-'],stdout=writer.stdin,env=env,stderr=log,check=True)
                except BaseException:
                    writer.terminate();writer.wait();raise
                finally:writer.stdin.close()
                if writer.wait():raise RuntimeError('직렬 CRAM 출력 실패')
            backend='serial_cram'
        subprocess.run([samtools,'quickcheck',str(destination)],check=True)
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        details={**cfg,'backend':backend,'wall_seconds':time.monotonic()-started,
                 'cpu_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
                 'order':'원래 SQ 순서와 각 contig 내부 record 순서, 마지막 unplaced를 보존함'}
        shutil.rmtree(directory)
        return details
    except BaseException as exc:
        (directory/'failure.json').write_text(json.dumps({'error':repr(exc),'source':str(source)},ensure_ascii=False))
        raise
