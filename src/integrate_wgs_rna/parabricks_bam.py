"""공유 자원 계약 안에서 Parabricks BAM/FASTQ 도구를 사용함."""
from pathlib import Path
import gzip
import hashlib
import json
import os
import shutil
import tempfile
import time

from .gpu_pool import reserved, GIB
from .parabricks_runtime import Runtime
from .memory_policy import memory_limit
from .pool_leases import cgroup_cpu_limit


def enabled():
    value=os.environ.get('MFM_BAM_BACKEND') or ('parabricks' if any(os.environ.get(k)=='parabricks' for k in ['MFM_WGS_BACKEND','MFM_RNA_BACKEND']) else 'cpu')
    if value not in {'cpu','parabricks'}:raise ValueError('MFM_BAM_BACKEND는 cpu/parabricks여야 함')
    return value=='parabricks'


def resource_plan():
    affinity=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else os.cpu_count() or 1
    caps=[affinity,cgroup_cpu_limit() or affinity]
    for key in ['MFM_CPU_BUDGET','MFM_CPU_POOL_GRANTED','MFM_SAMPLE_CPU_BUDGET','SLURM_CPUS_PER_TASK']:
        if os.environ.get(key):caps.append(int(os.environ[key]))
    cpu=min(caps);ram=memory_limit()*9//10
    value=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if value:
        row=json.loads(Path(value).read_text());cpu=min(cpu,int(row['cpus']));ram=min(ram,int(row['memory_bytes']))
    if cpu<1 or ram<2*GIB:raise ValueError('Parabricks BAM 도구의 최소 CPU/RAM 예약 부족')
    return {'cpus':cpu,'host_memory_bytes':ram,'mem_limit_gb':max(1,(ram-GIB)//10**9)}


def run_tool(args, paths, log, plan, *, gpu):
    runtime=Runtime();runtime.validate()
    if gpu:
        with reserved(int(os.environ.get('MFM_GPU_BAM_VRAM_BYTES',str(24*GIB))),
                      purpose='parabricks_'+args[0],exclusive=True) as device:
            result=runtime.run(args,device,paths,host_limit=plan['host_memory_bytes'],
                               cpu_limit=plan['cpus'],log=log)
    else:
        result=runtime.run(args,None,paths,host_limit=plan['host_memory_bytes'],
                           cpu_limit=plan['cpus'],log=log)
    return {**result,'resource_plan':plan,'gpu_kernels_requested':gpu,'tool':args[0]}


def unpack_fastq(source,destination,mate=None):
    from .rna_fastq_native import process
    return process(source,destination,mate)['records']


def rna_fastq(source, reference, parts, log):
    plan=resource_plan()
    with tempfile.TemporaryDirectory(prefix='parabricks_fastq_',dir=Path(parts['r1']).parent) as tmp:
        root=Path(tmp);prefix=root/'reads'
        args=['bam2fq','--in-bam',str(source),'--out-prefix',str(prefix),
              '--out-suffixO','_orphan1.fq.gz','--out-suffixO2','_orphan2.fq.gz',
              '--out-suffixS','_single.fq.gz','--num-threads',str(plan['cpus']),'--tmp-dir',str(root)]
        if reference:args+=['--ref',str(reference)]
        info=run_tool(args,[source,root,*([reference] if reference else [])],log,plan,gpu=False)
        from .rna_fastq_native import batch
        first=root/'orphan1.fq';second=root/'orphan2.fq'
        tasks=[(root/'reads_1.fastq.gz',parts['r1'],1),(root/'reads_2.fastq.gz',parts['r2'],2),
               (root/'reads_single.fq.gz',parts['un'],None),(root/'reads_orphan1.fq.gz',first,1),(root/'reads_orphan2.fq.gz',second,2)]
        for task in tasks[:2]:
            if not task[0].is_file():raise ValueError('Parabricks paired FASTQ 누락: '+str(task[0]))
        stats=batch(tasks,plan['cpus'])
        with Path(parts['sg']).open('wb') as out:
            for path in [first,second]:
                if path.exists():
                    with path.open('rb') as stream:shutil.copyfileobj(stream,out,16<<20)
        info['fastq_stats']=dict(zip(['r1','r2','un'],stats[:3]))
        info['fastq_stats']['sg']={'records':stats[3]['records']+stats[4]['records'],'max_read_length':max(stats[3]['max_read_length'],stats[4]['max_read_length'])}
    return info


def alignment_reference(alignment, work):
    """정렬에 사용한 HAP·흡수 contig 전체 reference를 확인함."""
    index=json.loads(Path(alignment['index_json']).read_text())
    contract=Path(index['index_dir'])/'gpu_index_contract.json'
    if contract.is_file():
        return Path(json.loads(contract.read_text())['reference'])
    sources=[Path(index['target_fasta'])]
    if index.get('addnl_fasta'):sources.append(Path(index['addnl_fasta']))
    reference=Path(work)/'parabricks_markdup_reference.fa'
    temporary=reference.with_suffix('.fa.part')
    with temporary.open('wb') as out:
        for path in sources:
            with path.open('rb') as source:shutil.copyfileobj(source,out,16<<20)
    os.replace(temporary,reference)
    import pysam
    pysam.faidx(str(reference))
    return reference


def normalize_markdup_header(output, work, log):
    """Parabricks GPU writer의 중복 HD만 제거하고 BAM read 본문을 보존함."""
    import pysam
    import subprocess
    with pysam.AlignmentFile(str(output)) as handle:lines=handle.text.splitlines()
    hd=[line for line in lines if line.startswith('@HD\t')]
    coordinate=[line for line in hd if 'SO:coordinate' in line.split('\t')]
    if not coordinate:raise ValueError('Parabricks markdup 출력이 coordinate 정렬이 아님')
    if len(hd)==1:return {'duplicate_hd_lines_removed':0}
    header=Path(work)/'coordinate.header.sam';part=Path(work)/'reheader.bam'
    header.write_text(coordinate[-1]+'\n'+'\n'.join(line for line in lines if not line.startswith('@HD\t'))+'\n')
    with part.open('wb') as out,Path(log).open('ab') as err:
        subprocess.run([os.environ.get('SAMTOOLS','samtools'),'reheader','-P',str(header),str(output)],stdout=out,stderr=err,check=True)
    pysam.quickcheck(str(part))
    with pysam.AlignmentFile(str(part)) as handle:
        if handle.header.to_dict().get('HD',{}).get('SO')!='coordinate':raise ValueError('Parabricks HD 정상화 실패')
    os.replace(part,output)
    return {'duplicate_hd_lines_removed':len(hd)-1,'read_body':'samtools reheader로 보존함'}


def check_markdup_validation():
    validation=os.environ.get('MFM_PARABRICKS_MARKDUP_VALIDATION')
    if validation:
        deadline=time.monotonic()+float(os.environ.get('MFM_GPU_WAIT_SECONDS','86400'))
        while True:
            checked=json.loads(Path(validation).read_text()) if Path(validation).is_file() else {}
            if checked.get('status') not in {None,'pending','running'}:break
            if time.monotonic()>=deadline:raise TimeoutError('Parabricks GPU 중복표시 동등성 검증 대기 시간 초과')
            time.sleep(1)
        if checked.get('status')!='pass' or checked.get('backend_sha256')!=hashlib.sha256(Path(__file__).read_bytes()).hexdigest():
            raise RuntimeError('Parabricks GPU 중복표시 동등성 검증 미통과 또는 코드 변경: '+validation)


def mark_duplicates(source, reference, output, metrics, work, log):
    check_markdup_validation()
    plan=resource_plan()
    if plan['cpus']<2:raise ValueError('Parabricks sort/markdup에는 최소2CPU가 필요함')
    worker=max(1,plan['cpus']//2);zip_threads=plan['cpus']-worker
    import pysam
    with pysam.FastaFile(str(reference)) as fasta,pysam.AlignmentFile(str(source)) as bam:
        lengths=dict(zip(fasta.references,fasta.lengths))
        if any(lengths.get(name)!=length for name,length in zip(bam.references,bam.lengths)):
            raise ValueError('Parabricks MarkDuplicates reference가 입력 BAM contig와 다름')
    with tempfile.TemporaryDirectory(prefix='parabricks_markdup_',dir=work) as tmp:
        root=Path(tmp);queryname=root/'queryname.bam'
        common=['--ref',str(reference),'--num-zip-threads',str(zip_threads),
                '--mem-limit',str(plan['mem_limit_gb']),'--tmp-dir',str(root),'--gpusort','--gpuwrite']
        sort=run_tool(['bamsort','--in-bam',str(source),'--out-bam',str(queryname),
                       '--sort-order','queryname','--num-sort-threads',str(worker),*common],
                      [source,reference,root],log,plan,gpu=True)
        mark=run_tool(['markdup','--in-bam',str(queryname),'--out-bam',str(output),
                       '--out-duplicate-metrics',str(metrics),'--optical-duplicate-pixel-distance','100',
                       '--num-worker-threads',str(worker),*common],
                      [queryname,reference,output,metrics,root],log,plan,gpu=True)
        import pysam
        pysam.quickcheck(str(output))
        header_repair=normalize_markdup_header(output,root,log)
        with pysam.AlignmentFile(str(output)) as handle:
            if handle.header.to_dict().get('HD',{}).get('SO')!='coordinate':
                raise ValueError('Parabricks markdup 출력이 coordinate 정렬이 아님')
        if not Path(metrics).is_file() or not Path(metrics).stat().st_size:
            raise ValueError('Parabricks duplicate metrics 누락')
    return {'backend':'parabricks','sort':sort,'markdup':mark,'remove_duplicates':False,
            'baseline_sort_order':'coordinate','optical_duplicate_pixel_distance':100,'header_repair':header_repair,
            'duplicate_metrics_sort_note':'Parabricks4.8은 metrics 요청 시 sorting phaseII를 CPU로 수행함'}
