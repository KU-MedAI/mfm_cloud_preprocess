"""기존 public VCF 변환을 염색체별로 병렬 실행하고 전역 정렬한다."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json
import fcntl
import os
import tempfile
import time
from .input_parallel import cpu_budget
from .memory_policy import GIB, memory_limit
from .cohort_vcf import run_command


def estimated_memory(workers):
    return ((2*GIB + workers*2*GIB)*6 + 4)//5


def entry_memory(cpus, chromosomes=22):
    """입장 전에는 가능한 worker 상한을 예약하고 실제 grant로 추가 감축함."""
    return estimated_memory(max(1, min(chromosomes, max(1, (cpus-1)//2))))


def reserve_validation_memory():
    value=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if not value:return
    row=json.loads(Path(value).read_text())
    reserve_stage_memory({'phase':'tag_and_validation', 'cpu_budget':min(2,int(row['cpus'])),
                          'estimated_peak_bytes':4*GIB})


def reserve_stage_memory(plan, lease_path=None):
    """단계 worker·도구 버퍼로 계산한 피크로만 감축하며 CPU와 전역 보호를 유지한다."""
    value=lease_path or os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if not value:return
    path=Path(value);root=path.parent.parent
    with (root/'pool.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        row=json.loads(path.read_text());required=int(plan['estimated_peak_bytes'])
        if row.get('branch')!='WGS' or row.get('stage')!='03':return
        active_cpu=int(plan.get('active_cpus',plan['cpu_budget']))
        if required>int(row['memory_bytes']) or active_cpu>int(row['cpus']):
            raise ValueError('VCF 계획이 현재 CPU/RAM 예약을 초과함')
        row.setdefault('public_vcf_initial_memory_bytes',row['memory_bytes'])
        row.setdefault('public_vcf_initial_cpus',row['cpus'])
        row.update(memory_bytes=required,cpus=active_cpu,public_vcf_memory_plan=plan)
        fd,temp=tempfile.mkstemp(prefix='.'+path.name,dir=path.parent)
        try:
            with os.fdopen(fd,'w') as out:json.dump(row,out,ensure_ascii=False);out.write('\n')
            os.replace(temp,path)
        finally:Path(temp).unlink(missing_ok=True)


def resource_plan(chromosomes):
    cpus=cpu_budget();ram=memory_limit()*9//10
    lease=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if lease:
        row=json.loads(Path(lease).read_text());cpus=min(cpus,int(row['cpus']));ram=min(ram,int(row['memory_bytes']))
    requested=os.environ.get('MFM_PUBLIC_VCF_JOBS','auto')
    if requested!='auto' and (not requested.isdigit() or int(requested)<1):
        raise ValueError('MFM_PUBLIC_VCF_JOBS는 auto 또는 양수여야 함')
    # worker당 native pipe 최대2CPU, 정렬/FASTA/BCF버퍼2GiB, 부모2GiB+20%여유.
    memory_workers=max(0,(ram*5//6-2*GIB)//(2*GIB))
    if not memory_workers:raise ValueError('public VCF 변환의 최소 RAM 예산 부족')
    workers=max(1,min(chromosomes,max(1,(cpus-1)//2),memory_workers))
    if requested!='auto':workers=min(workers,int(requested))
    return {'cpu_budget':cpus,'active_cpus':min(cpus,workers*2+1),'ram_budget_bytes':ram,'workers':workers,
            'estimated_peak_bytes':estimated_memory(workers),'chromosomes':chromosomes}


def build_parallel(config,sources,reference,chroms,policy,output,serial):
    plan=resource_plan(len(chroms));started=time.monotonic()
    reserve_stage_memory(plan)
    print('[PUBLIC-VCF-PARALLEL] '+json.dumps({'phase':'plan',**plan}),flush=True)
    def merge_resources():
        reserve_stage_memory({'phase':'global_merge','cpu_budget':min(2,plan['active_cpus']),
                              'estimated_peak_bytes':min(4*GIB,plan['estimated_peak_bytes'])})
    if plan['workers']==1:
        result=serial(config,sources,reference,chroms,policy,output)
        merge_resources()
        return result
    timings=[]
    try:
        with tempfile.TemporaryDirectory(prefix='.public_vcf_parallel_',dir=output.parent) as temporary:
            root=Path(temporary)
            def chromosome(item):
                index,chrom=item;folder=root/str(index);folder.mkdir();path=folder/'sequence.vcf.gz';start=time.monotonic()
                serial(config,sources,reference,[chrom],policy,path)
                return path,(start,time.monotonic())
            with ThreadPoolExecutor(max_workers=plan['workers']) as executor:
                parts=list(executor.map(chromosome,enumerate(chroms)))
            paths=[path for path,_ in parts];timings=[period for _,period in parts]
            # 변환 worker가 모두 join된 뒤 직렬 concat/sort/index의 몫만 유지한다.
            merge_resources()
            combined=root/'combined.bcf';ordered=root/'ordered.vcf.gz'
            run_command([config.bcftools_bin,'concat','-a','-Ob','-o',str(combined),*map(str,paths)])
            run_command([config.bcftools_bin,'sort','-Oz','-o',str(ordered),str(combined)])
            # 실패한 index와 완료본을 혼동하지 않도록 임시 파일에서 index까지 완성한다.
            run_command([config.bcftools_bin,'index','-t',str(ordered)])
            os.replace(ordered,output);os.replace(Path(str(ordered)+'.tbi'),Path(str(output)+'.tbi'))
    finally:
        events=sorted([(a,1) for a,b in timings]+[(b,-1) for a,b in timings]);active=peak=0
        for _,change in events:active+=change;peak=max(peak,active)
        print('[PUBLIC-VCF-PARALLEL] '+json.dumps({'phase':'complete',**plan,
            'wall_seconds':time.monotonic()-started,'completed_chromosomes':len(timings),
            'peak_active_chromosomes':peak}),flush=True)
