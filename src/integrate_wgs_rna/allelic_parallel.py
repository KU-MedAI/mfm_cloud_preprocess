"""읽기 전용 fragment table을 fork로 공유하고 작은 count만 반환함."""
import json,os,sys,threading,time,resource
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
from multiprocessing import get_all_start_methods,get_context
from .cpu_executor import WorkerBudget
from .memory_policy import GIB,memory_limit,limit_worker_address_space
_CONTEXT=None


def _assign_keys(keys):
    function,A,B,ex1,ex2,both,evidence=_CONTEXT
    started=time.monotonic();cpu=time.process_time()
    counts,stats=function(A,B,ex1,ex2,both,evidence,keys=keys)
    return counts,stats,{'pid':os.getpid(),'started':started,'finished':time.monotonic(),'cpu_seconds':time.process_time()-cpu}


def execution_summary(intervals):
    events=[];cpu=0.0
    for row in intervals:
        if row['finished']>row['started']:
            events.extend([(row['started'],1),(row['finished'],-1)])
        cpu+=row['cpu_seconds']
    active=peak=0
    for _at,change in sorted(events):active+=change;peak=max(peak,active)
    return {'peak_active_workers':peak,'parallel_observed':peak>1,'worker_cpu_seconds':cpu}


def _memory_plan():
    values=[int(x) for x in Path('/proc/self/statm').read_text().split()]
    page=os.sysconf('SC_PAGE_SIZE');vms,rss=values[0]*page,values[1]*page
    ceiling=max(vms+GIB,int(rss*1.3)+256*1024**2)
    return rss,ceiling


def assign_partitioned(function,A,B,ex1,ex2,both,evidence,*,workers=1,elastic=False,benchmark=False,batch_size=50000):
    global _CONTEXT
    if batch_size<1:raise ValueError('ASE batch_size는 양수여야 함')
    started=time.monotonic();parent_cpu_started=time.process_time();keys=list(set(A)|set(B));tasks=[(keys[i:i+batch_size],) for i in range(0,len(keys),batch_size)]
    enabled=(sys.platform.startswith('linux') and workers>1 and len(tasks)>1 and 'fork' in get_all_start_methods() and threading.active_count()==1)
    if not enabled:
        result=function(A,B,ex1,ex2,both,evidence)
        return result,{'mode':'serial','workers':1,'seconds':time.monotonic()-started,
                       'serial_equal':True,'speedup':1.0,'parallel_observed':False,'reason':'작은 입력/단일 CPU 또는 안전한 fork 조건 없음'}
    rss,ceiling=_memory_plan();base=max(4*GIB,int(rss*1.1))
    if ceiling+base>memory_limit()*9//10:raise ValueError('ASE fragment 한 작업의 예상 RAM이 예산보다 큼')
    cap=min(int(workers),len(tasks),max(1,int((memory_limit()*9//10-base)//ceiling)))
    report={'mode':'fork_fragments','workers':cap,'tasks':len(tasks),'fragments':len(keys),
            'parent_rss_bytes':rss,'worker_address_space_bytes':ceiling,'parent_reserved_bytes':base,
            'memory_limit_bytes':memory_limit(),'scratch_bytes':0,
            'memory_note':'fork 공유 페이지도 RSS 예약에서는 중복 계상하는 보수적 추정임'}
    print('[ASE-ASSIGN-PLAN] '+json.dumps(report,ensure_ascii=False),flush=True)
    merged,stats=function({}, {},ex1,ex2,both,evidence)
    old_base=os.environ.get('MFM_MEMORY_STAGE_BASE_BYTES');os.environ['MFM_MEMORY_STAGE_BASE_BYTES']=str(base)
    _CONTEXT=(function,A,B,ex1,ex2,both,evidence)
    reservation=WorkerBudget(cap,cap,elastic,memory_per_worker=ceiling)
    admitted=reservation.capacity(len(tasks))
    while admitted<1:
        time.sleep(.2);admitted=reservation.capacity(len(tasks))
    cap=admitted;report['workers']=cap
    # fork로 생성한 유휴 worker의 공유 RSS도 pool 종료까지 예약함.
    budget=WorkerBudget(cap,cap,False,memory_per_worker=0)
    budget.base_memory=base+cap*ceiling
    budget.capacity(len(tasks),grow=False)
    intervals=[]
    child_before=resource.getrusage(resource.RUSAGE_CHILDREN)
    def merge(value):
        counts,part_stats,execution=value
        intervals.append(execution)
        for key in ['a','b','both_a','both_b']:merged[key].update(counts[key])
        for key in ['split_a','split_b']:
            for half in [0,1]:merged[key][half].update(counts[key][half])
        for key,value in part_stats.items():stats[key]+=value
    try:
        with ProcessPoolExecutor(max_workers=cap,mp_context=get_context('fork'),
              initializer=limit_worker_address_space,initargs=(ceiling,)) as executor:
            pending=iter(tasks);futures=set();submitted=completed=0
            while completed<len(tasks):
                while submitted<len(tasks) and len(futures)<cap:
                    futures.add(executor.submit(_assign_keys,*next(pending)));submitted+=1
                done,_=wait(futures,return_when=FIRST_COMPLETED)
                for future in done:
                    futures.remove(future);merge(future.result());completed+=1
                budget.capacity(len(tasks)-completed,grow=False)

    finally:
        budget.base_memory=base
        budget.capacity(0,grow=False)
        _CONTEXT=None
        if old_base is None:os.environ.pop('MFM_MEMORY_STAGE_BASE_BYTES',None)
        else:os.environ['MFM_MEMORY_STAGE_BASE_BYTES']=old_base
    result=(merged,stats);report['seconds']=time.monotonic()-started
    child_after=resource.getrusage(resource.RUSAGE_CHILDREN)
    report.update(execution_summary(intervals))
    report['children_cpu_seconds']=max(0.0,child_after.ru_utime+child_after.ru_stime-child_before.ru_utime-child_before.ru_stime)
    report['parent_cpu_seconds']=time.process_time()-parent_cpu_started
    report['actual_cpu_seconds']=report['children_cpu_seconds']+report['parent_cpu_seconds']
    report['average_used_cores']=report['actual_cpu_seconds']/max(.000001,report['seconds'])
    if benchmark:
        check=time.monotonic();serial=function(A,B,ex1,ex2,both,evidence)
        serial_seconds=time.monotonic()-check
        if serial!=result:raise ValueError('ASE 직렬/병렬 fragment count/stat 전체 불일치')
        report.update(serial_equal=True,serial_seconds=serial_seconds)
        report['speedup']=report['serial_seconds']/max(.000001,report['seconds'])
    print('[ASE-ASSIGN-COMPLETE] '+json.dumps(report,ensure_ascii=False),flush=True)
    return result,report
