"""전수 flag/HP/PS 카운트를 네이티브 samtools로 계산함. 인덱스에 의존하지 않음."""
import json,os,shutil,subprocess,time
from pathlib import Path
from .input_parallel import cpu_budget
from .cpu_executor import WorkerBudget
from concurrent.futures import ThreadPoolExecutor
from .memory_policy import GIB


def fingerprint(path):
    p=Path(path).resolve();s=p.stat()
    return str(p),s.st_size,s.st_mtime_ns,s.st_ctime_ns,s.st_ino


def _run(command):
    result=subprocess.run(command,capture_output=True,text=True)
    if result.returncode:raise RuntimeError('전수 레코드 집계 실패: '+str(command)+': '+result.stderr[-4000:])
    return result.stdout


def _legacy_record_stats(path,reference=None,*,samtools=None,workers=None):
    path=Path(path);before=fingerprint(path)
    samtools=samtools or os.environ.get('SAMTOOLS_BIN') or shutil.which('samtools')
    if not samtools:raise FileNotFoundError('samtools 없음')
    total=max(1,min(cpu_budget(),workers or cpu_budget()))
    lease=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if lease:
        total=min(total,max(1,int(json.loads(Path(lease).read_text())['cpus'])))
    concurrent=min(3,total);cost=max(1,total//concurrent);additional=cost-1
    ref=[]
    if path.suffix.lower()=='.cram':
        if reference is None:raise ValueError('CRAM decode reference가 없음')
        ref=['--input-fmt-option','reference='+str(reference)]
    commands=[[samtools,'flagstat','-@',str(additional),'-O','json',*ref,str(path)],
              [samtools,'view','-@',str(additional),'-c','-d','HP',*ref,str(path)],
              [samtools,'view','-@',str(additional),'-c','-d','PS',*ref,str(path)]]
    started=time.monotonic()
    print('[RECORD-STATS-PLAN] '+json.dumps({'input':str(path),'cpu_budget':total,'concurrent_scans':concurrent,'cpus_per_scan':cost,'reserved_peak_bytes':4*GIB+concurrent*8*GIB,'passes':3,'index_required':False}),flush=True)
    # Stage01의 나머지 sort/index도 같은 grant를 쓰므로 집계 내부에서 caller CPU를 반환하지 않음.
    lease=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if lease:
        prior=json.loads(Path(lease).read_text())
        budget=WorkerBudget(total,total,False)
        budget.base_memory=max(int(prior.get('memory_bytes',0)),4*GIB+concurrent*8*GIB)
        admitted=budget.capacity(total,grow=False)
        while admitted<total:
            time.sleep(.2);admitted=budget.capacity(total)
    with ThreadPoolExecutor(max_workers=concurrent) as pool:
        values=list(pool.map(_run,commands))
    flagstat=json.loads(values[0]);groups=[flagstat[k] for k in ['QC-passed reads','QC-failed reads']]
    fields={'records':'total','mapped_records':'mapped','primary_records':'primary',
            'secondary_records':'secondary','supplementary_records':'supplementary',
            'duplicate_records':'duplicates','primary_duplicate_records':'primary duplicates'}
    result={key:sum(int(g[field]) for g in groups) for key,field in fields.items()}
    result.update(hp_tagged_records=int(values[1]),ps_tagged_records=int(values[2]))
    if fingerprint(path)!=before:raise RuntimeError('전수 집계 중 alignment가 변경됨')
    if any(v<0 or v>result['records'] for v in result.values()):raise ValueError('전수 집계 카운트 범위 오류')
    if sum(result[k] for k in ['primary_records','secondary_records','supplementary_records'])!=result['records']:raise ValueError('전수 집계 flag 합계 오류')
    print('[RECORD-STATS-COMPLETE] '+json.dumps({'seconds':time.monotonic()-started,'counts':result}),flush=True)
    return result


def native_record_stats(path,reference=None,*,samtools=None,workers=None):
    """BAM은 index-free 한 번 집계한다. CRAM과 개발 도구 없는 환경은 기존 경로를 쓴다."""
    from .record_stats_parallel import scan_plan, native_tool, run_scan, NativeUnavailable, RecordSizeLimit
    mode=os.environ.get('MFM_RECORD_STATS_SCAN_MODE','auto')
    if mode not in {'auto','legacy','single_pass'}:
        raise ValueError('MFM_RECORD_STATS_SCAN_MODE는 auto/legacy/single_pass여야 함')
    if mode=='legacy' or Path(path).suffix.lower()!='.bam':
        if mode=='single_pass':raise ValueError('single_pass는 BAM 입력에만 적용 가능함')
        return _legacy_record_stats(path,reference,samtools=samtools,workers=workers)
    try:
        tool=native_tool()
    except NativeUnavailable as error:
        if mode=='single_pass':raise
        print('[RECORD-STATS-FALLBACK] '+str(error),flush=True)
        return _legacy_record_stats(path,reference,samtools=samtools,workers=workers)
    before=fingerprint(path)
    total=max(1,min(cpu_budget(),workers or cpu_budget()))
    lease=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if lease:
        total=min(total,max(1,int(json.loads(Path(lease).read_text())['cpus'])))
    plan=scan_plan(total)
    from .memory_policy import memory_limit
    if plan['estimated_peak_rss_bytes']>memory_limit()*9//10:
        raise ValueError('전수 집계 예상 RAM이 운영 예산보다 큼')
    lease=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if lease:
        prior=json.loads(Path(lease).read_text());budget=WorkerBudget(total,total,False)
        budget.base_memory=max(int(prior.get('memory_bytes',0)),plan['estimated_peak_rss_bytes'])
        admitted=budget.capacity(total,grow=False)
        while admitted<total:
            time.sleep(.2);admitted=budget.capacity(total)
    print('[RECORD-STATS-PLAN] '+json.dumps({'input':str(path),**plan}),flush=True)
    try:
        values=run_scan(path,total,tool=tool,plan=plan)
    except (RecordSizeLimit,NativeUnavailable) as error:
        if mode=='single_pass':raise
        print('[RECORD-STATS-FALLBACK] '+str(error),flush=True)
        return _legacy_record_stats(path,reference,samtools=samtools,workers=workers)
    result=values['counts']
    if fingerprint(path)!=before:raise RuntimeError('전수 집계 중 alignment가 변경됨')
    if any(v<0 or v>result['records'] for v in result.values()):raise ValueError('전수 집계 카운트 범위 오류')
    if sum(result[k] for k in ['primary_records','secondary_records','supplementary_records'])!=result['records']:
        raise ValueError('전수 집계 flag 합계 오류')
    print('[RECORD-STATS-COMPLETE] '+json.dumps({'counts':result,**values['execution']}),flush=True)
    return result
