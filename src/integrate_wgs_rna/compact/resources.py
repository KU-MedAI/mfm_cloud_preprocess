"""compact의 해독·압축 worker를 실제 CPU 한도 안에서 배분함."""
import os,json
from pathlib import Path
from ..pool_leases import cgroup_cpu_limit,pool_cpu_limit
from ..memory_policy import GIB


def plan(workers, tasks, *, hap_split=False, available_cpus=None):
    caps=[len(os.sched_getaffinity(0)),cgroup_cpu_limit() or len(os.sched_getaffinity(0))]
    keys=['MFM_CPU_BUDGET','SLURM_CPUS_PER_TASK']
    if available_cpus is None:keys.append('MFM_CPU_POOL_GRANTED')
    else:caps.append(int(available_cpus))
    for key in keys:
        if os.environ.get(key):caps.append(int(os.environ[key]))
    lease=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if lease:
        path=Path(lease);row=json.loads(path.read_text())
        caps.append(pool_cpu_limit(path.parent.parent,int(row.get('total_cpus',min(caps)))))
    total=max(1,min(caps))
    base_blosc=max(1,int(os.environ.get('MFM_BLOSC_THREADS','1')))
    compression=os.environ.get('MFM_BLOSC_THREADS_REQUESTED',str(base_blosc))
    if compression!='auto' and (not compression.isdigit() or int(compression)<1):
        raise ValueError('MFM_BLOSC_THREADS_REQUESTED는 auto 또는 양수여야 함')
    gpu=any(os.environ.get('MFM_'+stage+'_BACKEND')=='cupy' for stage in ['ZARR','QC','COVERAGE'])
    native_auto=compression=='auto' and not gpu
    blosc=1 if native_auto else base_blosc if compression=='auto' else int(compression)
    requested=os.environ.get('MFM_COMPACT_HAP_JOBS','auto')
    if requested not in {'auto','1','2'}:raise ValueError('HAP 병렬 수는 auto|1|2이어야 함')
    haps=2 if hap_split and (requested=='2' or total>=4 and requested!='1') else 1
    count=max(1,min(workers,tasks,max(1,total//(2*haps if haps>1 else blosc))))
    while True:
        share=max(1,total//(count*haps))
        if native_auto:
            default=share if haps>1 else max(1,(share+1)//2)
        else:
            default=share if haps>1 else 1+max(0,total//count-blosc)
        threads=int(os.environ.get('MFM_COMPACT_ALIGNMENT_THREADS',str(default)))
        if threads<1:raise ValueError('MFM_COMPACT_ALIGNMENT_THREADS는 양수여야 함')
        if native_auto:blosc=share if haps>1 else max(1,share-threads+1)
        # 각 HAP의 read 생성과 append/압축은 순서대로 수행된다.
        hap_cost=max(blosc,threads) if haps>1 else blosc+threads-1
        cost=haps*hap_cost
        if cost>total:
            if hap_split and requested=='auto':
                return plan(workers,tasks,hap_split=False,available_cpus=available_cpus)
            raise ValueError('compact worker 한 개의 CPU가 실제 grant를 초과함')
        allowed=min(count,max(1,total//cost))
        if not native_auto or allowed==count:
            count=allowed
            break
        count=allowed
    return {'gpu':gpu,'workers':count,'cpu_per_worker':cost,'alignment_threads':threads,
            'cpu_limit':total,'blosc_threads':blosc,'hap_jobs':haps,'hap_cpu_cost':hap_cost,
            'extra_decode_buffers_bytes':haps*(threads-1)*32*(1<<20),
            'extra_codec_buffers_bytes':haps*(blosc-1)*32*(1<<20),
            'compression_thread_request':compression,
            'recycle_tasks':1 if gpu else None,'scope':'정상 전처리/benchmark 실제 CPU quota를 함께 우선함'}


def fit_memory(workers, tasks, forecast, *, hap_split=False, available_cpus=None, available_memory=None):
    """입장 전 RAM worker 상한에 맞춰 해독 스레드와 예약도 함께 재산정함."""
    per_hap_ceiling = forecast['worker_ceiling_bytes']
    per_hap_reservation = forecast['worker_reservation_bytes']
    per_hap_peak = forecast['estimated_peak_per_worker_bytes']
    usable = forecast['memory_limit_bytes'] * 9 // 10
    if available_memory is not None:usable=min(usable,int(available_memory))
    usable-=forecast.get('stage_base_memory_bytes',4*GIB)
    candidate = max(1, workers)
    while True:
        execution = plan(candidate, tasks, hap_split=hap_split,available_cpus=available_cpus)
        haps = execution['hap_jobs']
        parent = 2 * GIB if haps > 1 else 0
        decode = execution['extra_decode_buffers_bytes']+execution.get('extra_codec_buffers_bytes',0)
        reservation = per_hap_reservation * haps + parent + decode
        peak = per_hap_peak * haps + parent + decode
        ceiling = max(per_hap_ceiling * haps + parent + decode, reservation)
        cap = max(0, usable // reservation)
        if cap < 1 or peak > reservation:
            if hap_split and os.environ.get('MFM_COMPACT_HAP_JOBS','auto')=='auto':
                return fit_memory(workers,tasks,forecast,hap_split=False,
                                  available_cpus=available_cpus,available_memory=available_memory)
            raise MemoryError('해독 버퍼를 포함한 compact worker RAM 예약이 부족함')
        if execution['workers'] <= cap:
            execution['memory_worker_cap'] = cap
            forecast.update(worker_reservation_bytes=reservation, worker_ceiling_bytes=ceiling,
                            worker_address_space_bytes=ceiling,
                            per_hap_address_space_bytes=per_hap_ceiling+decode//haps,
                            estimated_peak_per_worker_bytes=peak,
                            single_stage_memory_worker_cap=cap,
                            execution_resources=execution)
            return execution
        candidate = min(execution['workers'] - 1, cap)


def reserve_shared(workers, tasks, forecast, *, hap_split=False):
    """해싱이 끝난 뒤 같은 잠금에서 현재 공유 CPU/RAM으로 계획하고 예약함."""
    from ..cpu_executor import elastic_cpu_enabled
    from ..memory_policy import memory_active, usable_memory
    lease=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if not lease or not elastic_cpu_enabled('MFM_COMPACT_JOBS_REQUESTED',thread_cost_accounted=True) or not memory_active(Path(lease).parent.parent):
        return fit_memory(workers,tasks,forecast,hap_split=hap_split)
    import fcntl,time,uuid
    from ..pool_leases import prune_leases,read_pool_records,regular_cpu_claim,fitting_waiters,live_waiters
    path=Path(lease);root=path.parent.parent
    # auto의 최초 작은 grant를 전체 계산 worker 상한으로 남기지 않는다.
    workers=max(1,tasks)
    # 전체 할당으로도 불가능한 설정은 자원대기로 숨기지 않는다.
    initial=json.loads(path.read_text())
    fit_memory(workers,tasks,dict(forecast),hap_split=hap_split,
               available_cpus=pool_cpu_limit(root,int(initial['total_cpus'])))
    queued=time.time()
    base=max(4*GIB,int(os.environ.get('MFM_MEMORY_STAGE_BASE_BYTES',4*GIB)))
    while True:
        execution=None
        with (root/'pool.lock').open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            prune_leases(root)
            row=json.loads(path.read_text())
            if row.get('cpu_backfill'):
                return fit_memory(workers,tasks,forecast,hap_split=hap_split)
            total=pool_cpu_limit(root,int(row['total_cpus']))
            others=[r for p,r in read_pool_records(root/'leases') if p!=path]
            cpus=max(0,total-sum(regular_cpu_claim(r) for r in others))
            ram=max(0,usable_memory(root)-sum(int(r.get('memory_bytes',0)) for r in others))
            rss=int(Path('/proc/self/statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE')
            base=max(base,(rss*11+9)//10)
            policy=json.loads((root/'policy.json').read_text()) if (root/'policy.json').exists() else {}
            waiting=fitting_waiters(sorted(live_waiters(root/'waiters'),key=lambda p:str(p[0])),
                cpus,max(0,ram-base),stage_limits=policy.get('stage_concurrency'),active=[*others,row],
                prepare_read_only=policy.get('prepare_phased') is True)
            budget_cpu=max(0,cpus-sum(int(r.get('min_cpus',1)) for _,r in waiting))
            budget_ram=max(0,ram-sum(int(r.get('memory_bytes',0)) for _,r in waiting))
            candidate=dict(forecast,stage_base_memory_bytes=base)
            if budget_cpu>0 and budget_ram>=base:
                try:
                    execution=fit_memory(workers,tasks,candidate,hap_split=hap_split,
                                         available_cpus=budget_cpu,available_memory=budget_ram)
                except (MemoryError,ValueError):
                    pass
            if execution is None:
                row.update(cpus=0,memory_bytes=base)
            else:
                row.update(cpus=execution['workers']*execution['cpu_per_worker'],
                    memory_bytes=base+execution['workers']*candidate['worker_reservation_bytes'])
                row['peak_cpus']=max(int(row.get('peak_cpus',0)),row['cpus'])
            details={'available_cpus':cpus,'available_memory_bytes':ram,
                     'cpu_after_waiters':budget_cpu,'memory_after_waiters_bytes':budget_ram,
                     'priority_waiters':len(waiting),'reserved_cpus':row['cpus'],
                     'reserved_memory_bytes':row['memory_bytes'],'waiting':execution is None,
                     'queue_wait_seconds':time.time()-queued}
            row['compact_entry_allocation']=details
            temporary=path.with_name(path.name+'.'+uuid.uuid4().hex+'.tmp')
            temporary.write_text(json.dumps(row,ensure_ascii=False)+'\n');os.replace(temporary,path)
            if execution is not None:
                forecast.update(candidate)
                forecast['shared_pool_allocation']=details
                print('[COMPACT-RESOURCES] '+json.dumps(details,ensure_ascii=False),flush=True)
                return execution
        time.sleep(.2)
