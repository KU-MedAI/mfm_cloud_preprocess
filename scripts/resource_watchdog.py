#!/usr/bin/env python3
"""공유 풀 하나의 최신 CPU/RAM 정책·할당·RSS 감시를 계속 유지한다."""
from pathlib import Path
import argparse
import fcntl
import json
import os
import socket
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src')]
from scripts.cpu_pool import allocate, atomic_json
from integrate_wgs_rna.memory_policy import protect_pool_memory, memory_limit, usable_memory, available_memory
from integrate_wgs_rna.pool_leases import process_alive, pool_cpu_limit, regular_cpu_claim, logical_cpu_claim, read_pool_records
from integrate_wgs_rna.pool_backfill import CpuObserver, lend_idle_cpus
from integrate_wgs_rna.cpu_loan_policy import available_regular_cpus
from integrate_wgs_rna.resource_fairness import priority_request
from integrate_wgs_rna.storage_admission import refresh_storage_capacity

from integrate_wgs_rna.pool_cpu_observer import PreprocessingCpuObserver

OBSERVER=CpuObserver()
PREPROCESSING_OBSERVER=PreprocessingCpuObserver()
LAST_METRICS=0


def process_start(pid):
    raw=(Path('/proc')/str(pid)/'stat').read_text()
    return int(raw[raw.rfind(')')+2:].split()[19])


def active_owner(record):
    if record.get('host',socket.gethostname())!=socket.gethostname():
        return True
    pid=record.get('pid')
    if not process_alive(pid):return False
    try:
        return record.get('started') is None or int(record['started'])==process_start(pid)
    except FileNotFoundError:
        return False
    except PermissionError:
        return True


def tick(root, requested, owner, *, exclusive=False):
    global LAST_METRICS
    started=time.monotonic()
    metrics={}
    # /proc 전체 순회는 다른 할당기를 막지 않도록 풀 잠금 밖에서 수행함.
    configured=json.loads((root/'policy.json').read_text()) if (root/'policy.json').exists() else {}
    preprocessing=PREPROCESSING_OBSERVER.observe(pool_cpu_limit(root,requested),configured.get('cpu_backfill',{}))
    metrics['preprocessing_observe_seconds']=time.monotonic()-started
    lock_started=time.monotonic()
    with (root/'pool.lock').open('a+') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX)
        locked=time.monotonic()
        metrics['lock_wait_seconds']=locked-lock_started
        prior=json.loads((root/'watchdog.json').read_text()) if (root/'watchdog.json').exists() else {}
        current=prior.get('watchdog_owner') or {'pid':prior.get('pid')}
        if exclusive and current.get('pid')!=os.getpid() and active_owner(current):
            return None
        total=pool_cpu_limit(root,requested)
        policy_path=root/'policy.json'
        policy=json.loads(policy_path.read_text()) if policy_path.exists() else {}
        storage=refresh_storage_capacity(policy,total)
        if storage and storage['policy_changed']:
            atomic_json(policy_path,policy)
        phase_started=time.monotonic()
        allocate(root,total)
        metrics['allocate_seconds']=time.monotonic()-phase_started
        phase_started=time.monotonic()
        guard=protect_pool_memory(root,force=True)
        metrics['rss_guard_seconds']=time.monotonic()-phase_started
        phase_started=time.monotonic()
        # 잠금 대기로 오래된 표본만 실제 재수집함. timestamp만 갱신하지 않음.
        if time.monotonic()-preprocessing.get('sampled_at',0)>2:
            preprocessing=PREPROCESSING_OBSERVER.observe(total,policy.get('cpu_backfill',{}))
        metrics['preprocessing_refresh_seconds']=time.monotonic()-phase_started
        observed=OBSERVER.observe(total)
        phase_started=time.monotonic()
        backfill=lend_idle_cpus(root,total,observed,OBSERVER,guard,preprocessing=preprocessing)
        metrics['lend_seconds']=time.monotonic()-phase_started
        rows=[row for _,row in read_pool_records(root/'leases')]
        available_cpu=available_regular_cpus(total,policy.get('cpu_backfill',{}),rows)
        available_ram=available_memory(root,rows)
        waits=[]
        for path,row in sorted(read_pool_records(root/'waiters'),key=lambda item:str(item[0])):
            reasons=[]
            if int(row.get('min_cpus',1))>available_cpu:reasons.append('기본 CPU 예약 부족')
            if int(row.get('memory_bytes',0))>available_ram:reasons.append('RAM 예약 부족')
            label=str(row.get('branch','')).upper()+':'+str(row.get('stage',''))
            limit=policy.get('stage_concurrency',{}).get(label)
            count=sum(1 for active in rows if str(active.get('branch','')).upper()+':'+str(active.get('stage',''))==label
                      and not active.get('stage_storage_wait') and not active.get('storage_paused'))
            if limit is not None and count>=int(limit) and not (label=='WGS:01' and policy.get('prepare_phased')):
                reasons.append('현재 디스크 여유에 따른 단계 동시성 대기')
            waits.append({'sample':row.get('sample'),'branch':row.get('branch'),'stage':row.get('stage'),
                          'minimum_cpus':row.get('min_cpus'),'memory_bytes':row.get('memory_bytes'),
                          'reasons':reasons})
        internal_waits=[{'sample':r.get('sample'),'branch':r.get('branch'),'stage':r.get('stage'),
                         'phase':r.get('realign_phase'),**r.get('realign_wait',{})}
                        for r in rows if r.get('realign_wait') or str(r.get('realign_phase','')).endswith(':waiting')]
        state={'schema':'mfm_resource_watchdog_v1','pid':os.getpid(),'owner_pid':owner,
               'updated':time.time(),'cpu_pool':total,'memory_limit_bytes':memory_limit(root),
               'reservation_bytes':usable_memory(root),'reserved_cpus':sum(logical_cpu_claim(r) for r in rows),
               'regular_reserved_cpus':sum(regular_cpu_claim(r) for r in rows),
               'backfill_reserved_cpus':sum(logical_cpu_claim(r) for r in rows if r.get('cpu_backfill')),
               'observed_cgroup_used_cpus':observed.get('used_cpus'),
               'observed_preprocessing_used_cpus':preprocessing.get('used_cpus'),
               'cpu_backfill':backfill,
               'prepare_storage':storage,
               'prepare_read_only_active':sum(1 for r in rows if r.get('prepare_phase')=='read_only_qc'),
               'prepare_write_waiting':[r.get('sample') for r in rows if r.get('stage_storage_wait')],
               'waiting_resources':waits,
               'internal_waiting_resources':internal_waits,
               'protected_resource_request':priority_request(root,total,rows),
               'reserved_memory_bytes':sum(int(r.get('memory_bytes',0)) for r in rows),
               'rss_bytes':guard['rss_bytes'],'rss_sampled_at':guard.get('time'),'guard_action':guard['action'],
               'active_leases':len(rows),'waiting':len(list((root/'waiters').glob('*.json')))}
        if exclusive:
            state['watchdog_owner']={'host':socket.gethostname(),'pid':os.getpid(),
                                     'started':process_start(os.getpid())}
        elif prior.get('watchdog_owner'):
            state['watchdog_owner']=prior['watchdog_owner']
        metrics.update(lock_hold_seconds=time.monotonic()-locked,
                       tick_seconds=time.monotonic()-started,
                       preprocessing_sample_age_seconds=(time.monotonic()-preprocessing['sampled_at']
                           if 'sampled_at' in preprocessing else None),
                       scope='최종 watchdog JSON 쓰기 직전까지의 실측')
        state['tick_metrics']=metrics
        atomic_json(root/'watchdog.json',state)
    if exclusive and time.monotonic()-LAST_METRICS>=10:
        LAST_METRICS=time.monotonic()
        record={key:state.get(key) for key in ['updated','cpu_pool','reserved_cpus','regular_reserved_cpus',
            'backfill_reserved_cpus','observed_cgroup_used_cpus','observed_preprocessing_used_cpus',
            'reserved_memory_bytes','rss_bytes','rss_sampled_at','guard_action','active_leases','waiting',
            'tick_metrics','protected_resource_request']}
        try:
            with (root/'resource_metrics.jsonl').open('a') as output:
                output.write(json.dumps(record,ensure_ascii=False)+'\n')
        except OSError as exc:
            print('[RESOURCE-METRICS] 기록 실패: '+str(exc),file=sys.stderr,flush=True)
    return state


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pool-dir',type=Path,required=True)
    p.add_argument('--total-cpus',type=int,required=True)
    p.add_argument('--owner-pid',type=int,default=0)
    p.add_argument('--once',action='store_true')
    a=p.parse_args();root=a.pool_dir
    for name in ['leases','waiters']: (root/name).mkdir(parents=True,exist_ok=True)
    started=False
    while True:
        state=tick(root,a.total_cpus,a.owner_pid,exclusive=True)
        if state is None:return 0
        if not started:
            print('[RESOURCE-WATCHDOG] started pid='+str(os.getpid()),flush=True)
            started=True
        if a.once or (not process_alive(a.owner_pid) and not state['active_leases'] and not state['waiting']):
            return 0
        time.sleep(.5)


if __name__=='__main__':raise SystemExit(main())
