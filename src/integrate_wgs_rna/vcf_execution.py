"""다인 패널 concat 없이 개인 변이를 먼저 추출하고 검증된 출력만 재사용함."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json
import os
import tempfile
import threading
from .memory_policy import memory_limit
from .pipeline.checkpoints import atomic_json,digest,file_identity


SAMPLE_CONTEXT=threading.local()


def sample_parallelism(requested, count):
    plan=resource_plan(1)
    minimum=int((2**30+768*2**20)*1.1)
    return max(1,min(requested,count,plan["cpu_budget"],plan["ram_budget_bytes"]//minimum))


def resource_plan(samples=1):
    caps=[len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else os.cpu_count() or 1]
    lease={}
    if os.environ.get('MFM_CPU_POOL_LEASE_FILE'):
        lease=json.loads(Path(os.environ['MFM_CPU_POOL_LEASE_FILE']).read_text())
    for key in ['MFM_CPU_BUDGET','SLURM_CPUS_PER_TASK']:
        if os.environ.get(key):caps.append(int(os.environ[key]))
    if lease.get('cpus'):caps.append(int(lease['cpus']))
    else:
        for key in ['MFM_CPU_POOL_GRANTED','MFM_SAMPLE_CPU_BUDGET']:
            if os.environ.get(key):caps.append(int(os.environ[key]))
    cpu=max(1,min(caps)//samples)
    ram=memory_limit()*9//10
    if lease.get('memory_bytes'):ram=min(ram,int(lease['memory_bytes']))
    for name,factor in [('SLURM_MEM_PER_NODE',1),('SLURM_MEM_PER_CPU',min(caps))]:
        if int(os.environ.get(name,'0'))>0:ram=min(ram,int(os.environ[name])*2**20*factor*9//10)
    ram//=samples
    # 한 염색체는 최대3개 단일코어 도구를 pipe로 실행함.
    memory_workers=max(1,(ram*10//11-2**30)//(768*2**20))
    workers=min(max(1,cpu//3),memory_workers)
    selection=min(max(1,cpu-1),memory_workers)
    if ram<(2**30+768*2**20)*11//10:raise ValueError('VCF 추출 RAM 예약 부족')
    return {'cpu_budget':cpu,'ram_budget_bytes':ram,'chromosome_workers':workers,
            'selection_workers':selection,
            'estimated_peak_bytes':int((2**30+max(workers,selection)*768*2**20)*1.1)}


def compression_threads(config, cpu_budget=None):
    """염색체 worker 종료 뒤 압축은 main1CPU를 제외한 예산을 활용함."""
    if config.concat_threads < 0:raise ValueError('concat_threads는0이상이어야 함')
    budget=resource_plan()['cpu_budget'] if cpu_budget is None else cpu_budget
    return max(0,min(config.concat_threads or max(0,budget-1),budget-1))


def panels(config):
    from .sample_vcf import resolve_common_vcf_path
    from .cohort_vcf import expand_chroms
    return [resolve_common_vcf_path(config.common_vcf_dir.resolve(),config.cohort,c,config.common_vcf_pattern)
            for c in expand_chroms(config.chroms)]


def source_signature(config,concat=None):
    paths=[concat] if concat is not None else panels(config)
    expanded=[]
    for p in paths:
        expanded.append(p)
        index=next((Path(str(p)+s) for s in ['.tbi','.csi'] if Path(str(p)+s).is_file()),None)
        if index is None:raise FileNotFoundError(f'VCF index 없음: {p}')
        expanded.append(index)
    return {'sources':[file_identity(p) for p in expanded],'policy':config.sequence_variant_policy,
            'sample_mode':config.extraction_mode}


def content(path):
    before=file_identity(path);h=hashlib.sha256()
    with Path(path).open('rb') as f:
        while block:=f.read(8<<20):h.update(block)
    if before!=file_identity(path):raise ValueError('VCF 검사 중 파일 변경')
    return {'stat':before,'sha256':h.hexdigest()}


def output_identity(seq,tag):
    output={}
    for p in [seq,tag]:
        output[p.name]=content(p)
        indexes=[Path(str(p)+s) for s in ['.tbi','.csi'] if Path(str(p)+s).is_file()]
        if not indexes:raise FileNotFoundError(p)
        for i in indexes:output[i.name]=content(i)
    return output


def validation_revision():
    root=Path(__file__).parent
    return digest([(root/'sample_vcf.py').read_text(),Path(__file__).read_text()])


def validated_metrics(config,sample,seq,tag,source,*,fresh=False):
    from .sample_vcf import validate_seq_vcf,validate_tag_vcf
    marker=seq.parent/f'.{sample}.vcf_validation.json'
    try:old=json.loads(marker.read_text())
    except (OSError,ValueError):old={}
    revision=validation_revision()
    observed=None
    if not fresh and old.get('schema')=='mfm_sample_vcf_validation_v1' and old.get('status')=='pass':
        if old['source']!=source:raise RuntimeError('VCF 입력/정책 변경: --overwrite로 재생성 필요')
        observed=output_identity(seq,tag)
        if {k:v['sha256'] for k,v in observed.items()}!={k:v['sha256'] for k,v in old['outputs'].items()}:
            raise RuntimeError('검증 이후 VCF/index 내용 변경: 완료 결과 재사용 차단')
        if old['revision']==revision and old.get('fast_validate')==config.fast_validate:
            print(f'[VCF-CACHE] reuse validated sample={sample}',flush=True)
            return old['seq_metrics'],old['tag_metrics']
    before=output_identity(seq,tag)
    a,b=validate_seq_vcf(config,seq,config.sequence_variant_policy),validate_tag_vcf(config,tag)
    if before!=output_identity(seq,tag):raise RuntimeError('VCF 검증 중 출력 변경')
    atomic_json(marker,{'schema':'mfm_sample_vcf_validation_v1','status':'pass','source':source,
                        'source_verified':fresh or old.get('source_verified',False),'revision':revision,'fast_validate':config.fast_validate,
                        'outputs':before,'seq_metrics':a,'tag_metrics':b})
    return a,b


def worker_affinity(sample,cpu_budget,slot):
    available=sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else []
    if available and os.environ.get('MFM_WGS_PAYLOAD')=='sequence_only':
        # 여러 샘플의 worker가 같은 첫 CPU에 모두 몰리지 않도록 분산함.
        shift=int.from_bytes(hashlib.sha256(sample.encode()).digest()[:8],'little')%len(available)
        available=available[shift:]+available[:shift]
    return available[slot*cpu_budget:(slot+1)*cpu_budget]


def extract_chromosomes(config,sample,seq,tag,vtypes):
    from .sample_vcf import extract_sequence_vcf,extract_tag_vcf,run_command
    sources=panels(config)
    from .cohort_vcf import load_sample_ids
    parallel_samples=sample_parallelism(config.jobs,len(load_sample_ids(config.sample_list,config.samples)))
    plan=resource_plan(max(1,parallel_samples))
    workers=min(len(sources),plan['chromosome_workers'])
    selectors=min(len(sources),plan['selection_workers'])
    peak=int((2**30+max(workers,selectors)*768*2**20)*1.1)
    # 다인 호출은 다른 sample의 선택 단계가 남을 수 있어 공통 예약을 줄이지 않음.
    reduce_memory=os.environ.get('MFM_WGS_PAYLOAD')=='sequence_only' and parallel_samples==1
    if reduce_memory:
        from .public_vcf_parallel import reserve_stage_memory
        # 최종 concat·변이 검증까지 포함하는4GiB 하한을 유지함.
        reserve_stage_memory({**plan,'phase':'sample_extract','estimated_peak_bytes':max(4<<30,peak*parallel_samples)})
    print('[VCF-PLAN] '+json.dumps({**plan,'chromosome_workers':workers,'selection_workers':selectors,'estimated_peak_bytes':peak,'sample':sample}),flush=True)
    # 개인 nonref만 BCF로 한 번 추출하여 원 패널의 두 번째 순회를 피함.
    with tempfile.TemporaryDirectory(prefix=f'.{sample}.vcf.',dir=seq.parent) as tmp:
        temp=Path(tmp)
        def select(item):
            i,p=item;selected=(temp/f'{i:02d}').with_suffix('.bcf')
            run_command([config.bcftools_bin,'view','-s',sample,'-c','1','-Ob','-o',str(selected),str(p)])
            return i
        def convert(i):
            base=temp/f'{i:02d}';selected=base.with_suffix('.bcf')
            a,b=base.with_suffix('.seq.vcf.gz'),base.with_suffix('.tag.vcf.gz')
            extract_sequence_vcf(config,sample,selected,a,vtypes)
            extract_tag_vcf(config,sample,selected,b)
            selected.unlink()
            return a,b
        slot=getattr(SAMPLE_CONTEXT,'slot',0)
        available=worker_affinity(sample,plan['cpu_budget'],slot)
        def pinned_map(operation,items,count,width):
            grants=[available[i*width:(i+1)*width] or available for i in range(count)]
            lock=threading.Lock()
            def pin_worker():
                if grants and hasattr(os,'sched_setaffinity'):
                    with lock:grant=grants.pop(0)
                    os.sched_setaffinity(0,grant)
            with ThreadPoolExecutor(max_workers=count,initializer=pin_worker) as pool:
                return list(pool.map(operation,items))
        # 긴1CPU 패널 선택에3CPU 슬롯을 묶어 두지 않음. 변환 순서는 그대로 유지함.
        selected=pinned_map(select,enumerate(sources),selectors,1)
        if reduce_memory:
            reserve_stage_memory({**plan,'phase':'sample_convert',
                'estimated_peak_bytes':max(4<<30,int((2**30+workers*768*2**20)*1.1))})
        parts=pinned_map(convert,selected,workers,3)
        if reduce_memory:
            reserve_stage_memory({**plan,'phase':'sample_finalize','estimated_peak_bytes':4<<30})
        threads=compression_threads(config,plan['cpu_budget'])
        print('[VCF-FINALIZE] '+json.dumps(dict(sample=sample,cpu_budget=plan['cpu_budget'],
              compression_threads=threads,index_threads=threads)),flush=True)
        for i,output in enumerate([seq,tag]):
            part=output.with_suffix(output.suffix+'.part')
            run_command([config.bcftools_bin,'concat','-Oz','--threads',str(threads),'-o',str(part),*[str(x[i]) for x in parts]])
            part.replace(output)
            run_command([config.bcftools_bin,'index','--threads',str(threads),'-t',str(output)])
