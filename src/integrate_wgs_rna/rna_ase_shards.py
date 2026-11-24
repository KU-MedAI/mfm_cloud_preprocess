"""QNAME를 보존한 binary shard에서 기존 fragment/variant 배정 함수를 병렬 호출함."""
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import get_context
from pathlib import Path
from functools import lru_cache
import importlib.util,importlib.machinery
import collections,json,os,resource,shutil,struct,subprocess,tempfile,time,hashlib,inspect
from .cpu_executor import WorkerBudget,cpu_results
from .memory_policy import GIB,memory_limit
from .record_stats_parallel import native_tool,NativeUnavailable
_CONTEXT=None
_FUSED_CONTEXT=None
_FIXED=struct.Struct('<IHBBqiI')


def selected(paths,workers):
    mode=os.environ.get('RNA_ALLELIC_INPUT_BACKEND','auto')
    if mode not in {'auto','partitioned','legacy'}:raise ValueError('RNA_ALLELIC_INPUT_BACKEND 오류')
    if mode=='legacy' or (mode=='auto' and (sum(Path(p).stat().st_size for p in paths)<64<<20)):return False
    try:native_tool('rna_ase_partition')
    except NativeUnavailable:
        if mode=='partitioned':raise
        return False
    return True


def _exact(handle,n):
    data=handle.read(n)
    if len(data)!=n:raise ValueError('잘린 ASE shard')
    return data


@lru_cache(maxsize=1)
def _native_reader():
    path=native_tool('rna_ase_read');name='rna_ase_read'
    loader=importlib.machinery.ExtensionFileLoader(name,str(path))
    spec=importlib.util.spec_from_file_location(name,path,loader=loader)
    module=importlib.util.module_from_spec(spec);loader.exec_module(module)
    return module.read_shard


def shard_reader():
    mode=os.environ.get('RNA_ASE_SHARD_READER','auto')
    if mode not in {'auto','python','native'}:raise ValueError('RNA_ASE_SHARD_READER는 auto/python/native임')
    if mode=='python':return _read_shard_python
    try:return _native_reader()
    except NativeUnavailable:
        if mode=='native':raise
        return _read_shard_python


def read_shard(path,mapq,duplicates):
    return shard_reader()(str(path),mapq,duplicates)


def _read_shard_python(path,mapq,duplicates):
    uniq={};mmap={};dupnames=set();dupcount=0
    with Path(path).open('rb',buffering=1<<20) as stream:
        if _exact(stream,8)!=b'MFMASE1\0':raise ValueError('ASE shard magic 오류')
        n=struct.unpack('<I',_exact(stream,4))[0]
        if n>1_000_000:raise ValueError('ASE reference 개수 오류')
        references=[]
        for _ in range(n):
            size=struct.unpack('<I',_exact(stream,4))[0]
            if size>1<<20:raise ValueError('ASE reference 이름 한도 초과')
            references.append(_exact(stream,size).decode())
        while raw:=stream.read(_FIXED.size):
            if len(raw)!=_FIXED.size:raise ValueError('잘린 ASE record')
            size,flag,mq,has_as,score,tid,blocks_n=_FIXED.unpack(raw)
            if size>1<<20 or blocks_n>1<<20 or tid<0 or tid>=len(references):raise ValueError('ASE shard record 범위 오류')
            name=_exact(stream,size).decode();raw_blocks=_exact(stream,16*blocks_n)
            if flag&1024:
                dupcount+=1;dupnames.add(name)
                if not duplicates:continue
            if not has_as or (mq==0 and mq<mapq):continue
            table=uniq if mq>=mapq else mmap
            mate=0 if flag&64 else 1;slot=table.setdefault(name,[None,None]);previous=slot[mate]
            if previous is None or score>previous[0]:slot[mate]=(score,references[tid],list(struct.iter_unpack('<qq',raw_blocks)))
    return uniq,mmap,dupcount,dupnames


_FUSED_POLICY_SHA256={'e2d1b4fbfebcb5238fb18001deb116ce29ebd1d5e258b777e8fb686d8d6e5141', '611ec8bb352f51de1a470e6bff1e81ba46ee095a2bff937aec7274c1ae58d6b8'}  # 전량 검증한 Python 정책 지문만 native 실행을 허용함.


def policy_fingerprint(assign,drop,length,evidence):
    genes=assign.__globals__['genes_of']
    indexed=genes.__globals__['genes_for_indexed']
    split=assign.__globals__['split_bucket']
    term=length.__globals__['star_length_term']
    pieces=[inspect.getsource(x) for x in [assign,drop,length,type(evidence),genes,indexed,split,term,_read_shard_python]]
    pieces.append(repr((indexed.__defaults__,assign.__globals__['MATES'])))
    return hashlib.sha256('\n'.join(pieces).encode()).hexdigest()


def fused_context(assign,drop,length,ex1,ex2,both,evidence,mapq,duplicates,scale,policy):
    mode=os.environ.get('RNA_ASE_ASSIGN_BACKEND','auto')
    if mode not in {'auto','tables','native'}:raise ValueError('RNA_ASE_ASSIGN_BACKEND는 auto/tables/native임')
    if mode=='tables':return None
    if policy!='validated_hap_v1' or policy_fingerprint(assign,drop,length,evidence) not in _FUSED_POLICY_SHA256:
        if mode=='native':raise ValueError('native ASE의 검증된 Python 정책 지문이 아님')
        return None
    try:
        from rna_allelic_common import split_bucket
        return _native_reader().__self__.prepare_assign(ex1,ex2,evidence.idx,both,mapq,duplicates,scale,split_bucket)
    except (NativeUnavailable,OverflowError):
        if mode=='native':raise
        return None


def _part_task(paths):
    assign,drop,length,ex1,ex2,both,evidence,mapq,duplicates,scale=_CONTEXT
    started=time.monotonic();cpu=time.process_time()
    if _FUSED_CONTEXT is not None:
        try:
            value=_native_reader().__self__.assign_shards(str(paths[0]),str(paths[1]),_FUSED_CONTEXT)
        except OverflowError:
            if os.environ.get('RNA_ASE_ASSIGN_BACKEND')=='native':raise
        else:
            *values,tm=value
            tm.update(cpu_seconds=time.process_time()-cpu,pid=os.getpid(),started=started,finished=time.monotonic())
            return (*values,tm)
    A,Am,d1,n1=read_shard(paths[0],mapq,duplicates);B,Bm,d2,n2=read_shard(paths[1],mapq,duplicates)
    read_cpu=time.process_time()-cpu;filter_started=time.monotonic();filter_cpu=time.process_time()
    read_seconds=filter_started-started;union=n1|n2
    removed=drop(union,A,Am,B,Bm) if not duplicates else 0
    corrected=sum(length(t,scale) for t in (A,Am,B,Bm))
    filter_seconds=time.monotonic()-filter_started;filter_cpu=time.process_time()-filter_cpu
    assign_started=time.monotonic();assign_cpu=time.process_time()
    counts,stats=assign(A,B,ex1,ex2,both,evidence);mm,mmstats=assign(Am,Bm,ex1,ex2,both,evidence)
    audit=dict(dup1=d1,dup2=d2,dup_union=len(union),dup_inconsistent=len(n1^n2),removed=removed,corrected=corrected,total=len(set(A)|set(B)))
    return counts,stats,mm,mmstats,audit,{'read_seconds':read_seconds,'shard_read_cpu_seconds':read_cpu,'filter_seconds':filter_seconds,'filter_cpu_seconds':filter_cpu,'assign_seconds':time.monotonic()-assign_started,'assign_cpu_seconds':time.process_time()-assign_cpu,'cpu_seconds':time.process_time()-cpu,'pid':os.getpid(),'started':started,'finished':time.monotonic()}


def _merge_counts(dest,value):
    for key in ['a','b','both_a','both_b']:dest[key].update(value[key])
    for key in ['split_a','split_b']:
        for half in [0,1]:dest[key][half].update(value[key][half])


def spool_sizes(tasks):
    files=[path for pair in tasks for path in pair[0]]
    stats=[path.stat() for path in files]
    logical=sum(s.st_size for s in stats)
    allocated=sum(s.st_blocks*512 for s in stats)
    return {'logical_bytes':logical,'allocated_bytes':allocated,'accounted_bytes':max(logical,allocated),'files':len(files)}


def assign_from_bams(paths,assign,drop,length,ex1,ex2,both,evidence,mapq,duplicates,scale,workers,initial_workers=None,*,native_policy=None):
    global _CONTEXT,_FUSED_CONTEXT
    result=None
    try:
        result=_assign_temporary(paths,assign,drop,length,ex1,ex2,both,evidence,mapq,duplicates,scale,workers,initial_workers,native_policy=native_policy)
        return result
    finally:
        # 모든 worker join과 TemporaryDirectory 정리가 끝난 뒤에만 tmpfs 예약을 반환함.
        _CONTEXT=None;_FUSED_CONTEXT=None
        rss=int(Path('/proc/self/statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE')
        release=WorkerBudget(1,1,True);release.base_memory=max(4*GIB,int(rss*1.1));release.capacity(0,grow=False)
        if result is not None:result[-1].update(spool_released_after_cleanup=True,post_cleanup_parent_bytes=release.base_memory)


def _assign_temporary(paths,assign,drop,length,ex1,ex2,both,evidence,mapq,duplicates,scale,workers,initial_workers=None,*,native_policy=None):
    global _CONTEXT,_FUSED_CONTEXT
    reader=shard_reader();reader_backend='native' if reader is not _read_shard_python else 'python'
    started=time.monotonic();raw=sum(Path(p).stat().st_size for p in paths)
    shards=min(128,max(2,workers*4));spool_estimate=raw*6+4*GIB
    scratch=Path(os.environ.get('TMPDIR',tempfile.gettempdir()))
    if shutil.disk_usage(scratch).free<spool_estimate:raise ValueError('ASE shard 임시공간 예상량 부족')
    print('[ASE-PREFLIGHT] '+json.dumps({'input_bytes':raw,'scratch_estimate_bytes':spool_estimate,'shards':shards,'cpu_upper':workers,'ram_limit_bytes':memory_limit()}),flush=True)
    with tempfile.TemporaryDirectory(prefix='rna_ase_',dir=scratch) as tmp:
        folder=Path(tmp);tool=native_tool('rna_ase_partition');prefixes=[folder/'H1',folder/'H2']
        input_cpus=workers;budget=WorkerBudget(min(workers,initial_workers or workers),workers,True,memory_per_worker=64<<20)
        budget.base_memory=4*GIB+spool_estimate
        grant=0
        while grant<1:grant=budget.capacity(input_cpus);time.sleep(.02 if grant else .2)
        input_jobs=2 if grant>=3 else 1
        decode_threads=max(0,(grant-1)//input_jobs-2)
        before=resource.getrusage(resource.RUSAGE_CHILDREN)
        print('[ASE-DECODE] '+json.dumps({'granted_cpu_upper':grant,'bam_jobs':input_jobs,'bgzf_threads_per_bam':decode_threads,'predicted_ram_bytes':budget.base_memory+grant*(64<<20)}),flush=True)
        def partition(task):
            path,prefix=task;result=subprocess.run([str(tool),str(path),str(prefix),str(shards),str(int(duplicates)),str(decode_threads),str(mapq)],capture_output=True,text=True)
            if result.returncode:raise ValueError('ASE 입력 분할 실패: '+result.stderr[-2000:])
            return json.loads(result.stdout)
        with ThreadPoolExecutor(max_workers=input_jobs) as pool:input_stats=list(pool.map(partition,zip(paths,prefixes)))
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        input_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime
        input_seconds=time.monotonic()-started
        tasks=[((Path(str(prefixes[0])+f'.{i}.bin'),Path(str(prefixes[1])+f'.{i}.bin')),) for i in range(shards)]
        spool=spool_sizes(tasks)
        if spool['accounted_bytes']>spool_estimate:raise MemoryError('ASE 실제 spool이 사전 RAM 예상을 초과함. 배정 전 중단·정리함')
        largest=max(sum(p.stat().st_size for p in pair[0]) for pair in tasks)
        vm,rss=[int(x)*os.sysconf('SC_PAGE_SIZE') for x in Path('/proc/self/statm').read_text().split()[:2]]
        per_worker=max(2*GIB,vm+largest*32+GIB);base=max(4*GIB,int(rss*1.1))+spool['accounted_bytes']
        available=memory_limit()*9//10
        lease=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
        if lease:
            from .memory_policy import usable_memory
            from .pool_leases import read_pool_records
            path=Path(lease);pool_root=path.parent.parent
            others=sum(int(row.get('memory_bytes',0)) for other,row in read_pool_records(pool_root/'leases') if other!=path)
            available=min(available,max(0,usable_memory(pool_root)-others))
        cap=min(max(1,workers-1),shards,max(0,(available-base)//per_worker))
        if cap<1:raise ValueError('ASE shard 한 작업의 CPU/RAM 예산 부족')
        _CONTEXT=(assign,drop,length,ex1,ex2,both,evidence,mapq,duplicates,scale)
        _FUSED_CONTEXT=fused_context(assign,drop,length,ex1,ex2,both,evidence,mapq,duplicates,scale,native_policy)
        counts,stats=assign({}, {},ex1,ex2,both,evidence);mm,mmstats=assign({}, {},ex1,ex2,both,evidence);audit=collections.Counter();timing=collections.Counter();traces=[]
        native_used=_FUSED_CONTEXT is not None
        began=time.monotonic();old_base=os.environ.get('MFM_MEMORY_STAGE_BASE_BYTES');os.environ['MFM_MEMORY_STAGE_BASE_BYTES']=str(base)
        results=cpu_results(_part_task,tasks,min(grant,cap),elastic=True,mp_context=get_context('fork'),max_workers=cap,memory_per_worker_bytes=per_worker,retain_worker_memory=True,parent_cpus=int(workers>1))
        try:
            for _index,value,error in results:
                if error:raise error
                cnt,st,mc,ms,au,tm=value;_merge_counts(counts,cnt);_merge_counts(mm,mc)
                for key,v in st.items():stats[key]+=v
                for key,v in ms.items():mmstats[key]+=v
                audit.update(au);timing.update({k:tm[k] for k in ['cpu_seconds','read_seconds','shard_read_cpu_seconds','filter_cpu_seconds','assign_cpu_seconds']});traces.append(tm)
        finally:
            results.close()
            _CONTEXT=None;_FUSED_CONTEXT=None
            if old_base is None:os.environ.pop('MFM_MEMORY_STAGE_BASE_BYTES',None)
            else:os.environ['MFM_MEMORY_STAGE_BASE_BYTES']=old_base
        active=peak=0
        for _,change in sorted([(v['started'],1) for v in traces]+[(v['finished'],-1) for v in traces]):
            active+=change;peak=max(peak,active)
        execution={'assignment_backend':'native' if native_used else 'tables','shard_reader':reader_backend,'shard_read_cpu_seconds':timing['shard_read_cpu_seconds'],'filter_cpu_seconds':timing['filter_cpu_seconds'],'assign_cpu_seconds':timing['assign_cpu_seconds'],'mode':'qname_shards','workers':cap,'shards':shards,'seconds':time.monotonic()-began,'parallel_observed':peak>1,'max_active_workers':peak,'worker_cpu_seconds':timing['cpu_seconds'],'shard_read_seconds_sum':timing['read_seconds'],'worker_address_space_bytes':per_worker,'parent_reserved_bytes':base,'shared_ram_budget_bytes':available,'spool_actual':spool,'serial_equal':None}
        inputs={'mode':'native_qname_partition','workers':input_jobs,'granted_cpu_upper':grant,'decode_threads_per_bam':decode_threads,'cpu_seconds':input_cpu_seconds,'seconds':input_seconds,'input_bytes':raw,'spool_estimate_bytes':spool_estimate,'spool_actual':spool,'records':input_stats}
        print('[ASE-PARTITION] '+json.dumps({'input':inputs,'assignment':execution},ensure_ascii=False),flush=True)
        return counts,stats,mm,mmstats,dict(audit),inputs,execution


def verify_legacy(paths,assign,drop,length,ex1,ex2,both,evidence,mapq,duplicates,scale,expected):
    from .allelic_input import read_pair_as
    started=time.monotonic()
    pair,inputs=read_pair_as(*paths,mapq,duplicates,workers=1)
    (A,Am,d1,n1),(B,Bm,d2,n2)=pair;union=n1|n2
    began=time.monotonic();removed=drop(union,A,Am,B,Bm) if not duplicates else 0
    corrected=sum(length(t,scale) for t in (A,Am,B,Bm))
    cnt,st=assign(A,B,ex1,ex2,both,evidence);mm,ms=assign(Am,Bm,ex1,ex2,both,evidence)
    audit=dict(dup1=d1,dup2=d2,dup_union=len(union),dup_inconsistent=len(n1^n2),removed=removed,corrected=corrected,total=len(set(A)|set(B)))
    if (cnt,st,mm,ms,audit)!=expected:raise ValueError('ASE native/legacy 전체 count·통계·duplicate audit 불일치')
    return dict(input_seconds=inputs['seconds'],assignment_seconds=time.monotonic()-began,total_seconds=time.monotonic()-started,serial_equal=True,scope='all_counts_statistics_and_duplicate_audit')
