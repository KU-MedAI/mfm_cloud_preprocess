"""RNA용 WGS 지지수·CNV native 도구와 CPU/RAM 예산을 공유함."""
from pathlib import Path
import json,os,subprocess,tempfile
from .record_stats_parallel import NativeUnavailable,native_tool
from .gpu_compare_resources import cpu_limit
from .memory_policy import GIB,memory_limit


def select(variable,tool):
    mode=os.environ.get(variable,'auto')
    if mode not in {'auto','native','legacy'}:raise ValueError(variable+'는 auto/native/legacy임')
    if mode=='legacy':
        print('[RNA-WGS-BACKEND] '+variable+'=legacy',flush=True);return None
    try:
        path=str(native_tool(tool))
        print('[RNA-WGS-BACKEND] '+json.dumps({'selector':variable,'backend':'native','tool':path}),flush=True)
        return path
    except NativeUnavailable as error:
        if mode=='native':raise
        print('[RNA-WGS] native 빌드 불가, legacy 사용: '+str(error),flush=True)
        return None


def cpu_upper(initial):
    upper=cpu_limit(int(os.environ.get('MFM_CPU_BUDGET',str(initial))))
    ceiling=os.environ.get('MFM_CPU_POOL_STAGE_MAX')
    dynamic=os.environ.get('MFM_RNA_DAG_DYNAMIC_CPU_MAX')
    lease=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    role=os.environ.get('MFM_RNA_DAG_ROLE')
    elastic=(os.environ.get('MFM_CPU_UTILIZATION','maximum')=='maximum' and os.environ.get('MFM_CPU_POOL_PROFILE','throughput')=='throughput')
    if elastic and dynamic and lease and role in {'variants','cnv','ase'}:
        row=json.loads(Path(lease).read_text())
        legacy=os.environ.get('MFM_CPU_POOL_LAYOUT')=='legacy'
        expected={'variants':'18' if legacy else '20','cnv':'20' if legacy else '22','ase':'19' if legacy else '21'}[role]
        if row.get('branch')=='RNA' and str(row.get('stage'))==expected and not row.get('cpu_backfill'):
            ceiling=min(int(dynamic),int(row['total_cpus']))
    if ceiling is not None:upper=min(upper,max(1,int(ceiling)))
    return upper


def plan(tasks,*,cpu_request=None,requested='auto',decode_threads=None):
    initial=cpu_limit(int(os.environ.get('MFM_SAMPLE_CPU_BUDGET',os.environ.get('MFM_CPU_BUDGET',str(cpu_request or 68)))))
    upper=cpu_upper(initial)
    if cpu_request is not None:
        if os.environ.get('MFM_DYNAMIC_CPU_POOL_ACTIVE')!='1' or os.environ.get('MFM_CPU_POOL_STAGE','') in os.environ.get('MFM_DYNAMIC_RNA_FIXED_STAGES','').split(','):
            upper=min(upper,max(1,int(cpu_request)))
    # HTSlib 병렬 해독에는 decoder 외에 I/O thread 한 개가 추가됨.
    decode=(1 if upper>=6 else 0) if decode_threads is None else int(decode_threads)
    if decode<0:raise ValueError('해독 thread는 음수가 될 수 없음')
    cost=1+decode+int(decode>0)
    workers=min(tasks,max(1,(upper-1)//cost))
    if requested!='auto':workers=min(workers,max(1,int(requested)))
    per_worker=4*GIB;parent=4*GIB;memory=min(1_650_000_000_000,memory_limit())*9//10
    lease=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if lease:
        from .memory_policy import usable_memory
        from .pool_leases import read_pool_records
        path=Path(lease);root=path.parent.parent
        others=sum(int(row.get('memory_bytes',0)) for other,row in read_pool_records(root/'leases') if other!=path)
        memory=min(memory,max(0,usable_memory(root)-others))
    workers=min(workers,(memory-parent)//per_worker)
    if workers<1:raise ValueError('RNA WGS 검증 한 작업의 예상 RAM이 예산을 초과함')
    info=dict(workers_upper=workers,workers_initial=min(workers,max(1,(initial-1)//cost)),
              decode_threads=decode,cpu_per_worker=cost,cpu_upper=upper,parent_cpus=int(upper>1),
              worker_address_space_bytes=per_worker,parent_bytes=parent,
              shared_ram_budget_bytes=memory,peak_ram_upper_bytes=parent+workers*per_worker)
    print('[RNA-WGS-PREFLIGHT] '+json.dumps(info),flush=True)
    return info


def count_support(tool,bam,reference,contig,groups,mapq,decode=0):
    rows=[(region,v) for region,variants in groups for v in variants]
    if not rows:return
    with tempfile.TemporaryDirectory(prefix='rna_variants_') as raw:
        folder=Path(raw);source=folder/'variants.tsv';target=folder/'counts.tsv'
        with source.open('w') as out:
            for (start,end),v in rows:out.write(f"{start}\t{end}\t{v['p0']}\t{v['ref']}\t{v['alt']}\t{v['kind']}\n")
        with target.open('w') as out:
            result=subprocess.run([tool,str(bam),str(reference or '-'),contig,str(source),str(mapq),str(decode)],stdout=out,stderr=subprocess.PIPE,text=True)
        if result.returncode:raise RuntimeError(contig+' native variant 지원도 실패: '+result.stderr[-4000:])
        fields=['wgs_depth','wgs_alt','hp1_alt','hp1_depth','hp2_alt','hp2_depth'];n=0
        with target.open() as stream:
            for line in stream:
                values=list(map(int,line.rstrip().split('\t')))
                if n>=len(rows) or len(values)!=6 or min(values)<0:raise ValueError('variant 지원도 출력 손상')
                if values[1]>values[0] or values[2]>values[3] or values[4]>values[5]:raise ValueError('variant 지원도 범위 오류')
                rows[n][1].update(zip(fields,values));n+=1
        if n!=len(rows):raise ValueError('variant 지원도 행 수 불일치')
        print('[variant:native] '+result.stderr.strip(),flush=True)
