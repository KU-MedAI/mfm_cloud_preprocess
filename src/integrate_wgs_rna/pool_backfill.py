"""실제 유휴 CPU를 정책 우선순위로 대여한다. RAM 예약은 유지한다."""
from collections import deque
from pathlib import Path
import json
import math
import os
import tempfile
import time

from .memory_policy import usable_memory, available_memory
from .cpu_loan_policy import loan_limits, available_loan_cpus
from .pool_leases import fitting_waiters, process_alive, read_pool_records, logical_cpu_claim
from .cgroup_memory import pressure_usage


def cgroup_sample():
    """CPU quota를 실제로 강제하는 가장 작은 cgroup의 누계를 읽는다."""
    relative = next(line.split('::', 1)[1].strip().lstrip('/')
                    for line in Path('/proc/self/cgroup').read_text().splitlines() if '::' in line)
    mount = Path('/sys/fs/cgroup');leaf = mount/relative
    cpu = [];memory = []
    for parent in [leaf, *leaf.parents]:
        if parent != mount and mount not in parent.parents:
            break
        path = parent/'cpu.max'
        raw = path.read_text().split() if path.exists() else ['max']
        if raw[0] != 'max':
            cpu.append((int(raw[0])/int(raw[1]), parent))
        path = parent/'memory.max'
        raw = path.read_text().strip() if path.exists() else 'max'
        if raw != 'max':
            memory.append((int(raw), parent))
    if not cpu:
        return None
    quota, parent = min(cpu, key=lambda value:value[0])
    stat = dict(line.split() for line in (parent/'cpu.stat').read_text().splitlines())
    maximum, folder = min(memory, key=lambda value:value[0]) if memory else (None, None)
    current = int((folder/'memory.current').read_text()) if folder else None
    return {'quota':quota, 'usage_usec':int(stat['usage_usec']),
            'memory_current':current,
            'memory_working_bytes':pressure_usage(current, maximum, folder),
            'memory_max':maximum, 'cpu_cgroup':str(parent)}


class CpuObserver:
    def __init__(self):
        self.history = deque()
        self.cooldown_until = 0

    def observe(self, total, *, window=10, now=None, reading=None):
        explicit_time = now is not None
        now = time.monotonic() if now is None else now
        try:
            sample = cgroup_sample() if reading is None else reading
        except (OSError, ValueError, StopIteration):
            sample = None
        if not explicit_time:
            now = time.monotonic()
        if sample is None or abs(sample['quota']-total) > .001:
            self.history.clear()
            return {'ready':False,'reason':'실제 CPU quota와 풀 상한 일치 확인 필요'}
        if self.history and sample['usage_usec'] < self.history[-1][1]:
            self.history.clear()
        self.history.append((now, sample['usage_usec']))
        while len(self.history)>2 and self.history[1][0] <= now-window:
            self.history.popleft()
        elapsed = now-self.history[0][0]
        state = dict(sample, ready=False, observation_seconds=elapsed, sampled_at=now)
        if elapsed < window or len(self.history)<3:
            state['reason']='CPU 유휴 관측 중'
            return state
        # 2초 평균의 최근 최대값을 사용하여 짧은 직렬 구간을 여유로 오인하지 않는다.
        points = list(self.history);rates=[]
        for index,(end,usage) in enumerate(points):
            prior=[p for p in points[:index] if end-p[0]>=2]
            if prior:
                start,begin = prior[-1]
                rates.append((usage-begin)/1e6/(end-start))
        used=(sample['usage_usec']-points[0][1])/1e6/elapsed
        peak=max([used,*rates])
        state.update(ready=True, used_cpus=round(used,3), peak_used_cpus=round(peak,3),
                     idle_cpus=max(0,math.floor(total-peak-1)),reason='실제 유휴 CPU 관측 완료')
        return state


def _atomic(path, row):
    fd,name=tempfile.mkstemp(prefix='.'+path.name,dir=path.parent)
    try:
        with os.fdopen(fd,'w') as handle:
            json.dump(row,handle,ensure_ascii=False);handle.write('\n')
        os.replace(name,path)
    finally:
        Path(name).unlink(missing_ok=True)


def eligible_waiter(row):
    """고정 단계와 피크 재예약을 선언한 개인정렬의 최소 준비 진입만 대여한다."""
    if row.get('fixed_sample_token'):
        return False
    if not process_alive(row.get('pid')):
        return False
    try:
        raw=(Path('/proc')/str(row['pid'])/'environ').read_bytes().split(b'\0')
        env=dict(value.decode().split('=',1) for value in raw if b'=' in value)
    except (OSError,UnicodeError):
        return False
    if row.get('cpu_utilization',env.get('MFM_CPU_UTILIZATION','maximum')) != 'maximum':
        return False
    branch=row['branch'].upper();stage=row['stage']
    legacy=row.get('layout',env.get('MFM_CPU_POOL_LAYOUT'))=='legacy'
    if branch=='WGS' and stage==('08' if legacy else '10'):
        try:
            plan=json.loads(env.get('MFM_LONGREAD_MEMORY_PLAN') or '{}')
        except ValueError:
            return False
        return (isinstance(plan,dict) and isinstance(plan.get('stage_estimates_bytes'),dict)
                and env.get('MFM_WGS_BACKEND','cpu')=='cpu'
                and env.get('MFM_REALIGN_PHASE_RESOURCES','auto')!='0'
                and plan.get('adaptive_realign') is True
                and plan.get('stage_estimates_bytes',{}).get('personal_realign_entry')==4*1024**3
                and int(row.get('memory_bytes',0))==4*1024**3)
    return (branch=='WGS' and stage in {'01','02','03','06','23','24'}) or (branch=='RNA' and stage in
            ({'11','12','13','14','15'} if legacy else {'13','14','15','16','17'}))


def in_cpu_group(pid,folder):
    try:
        relative=next(line.split('::',1)[1].strip().lstrip('/')
                      for line in (Path('/proc')/str(pid)/'cgroup').read_text().splitlines() if '::' in line)
        leaf=Path('/sys/fs/cgroup')/relative;parent=Path(folder)
        return parent==leaf or parent in leaf.parents
    except (OSError,StopIteration):
        return False



def lower_process_tree_priority(manager):
    """계산을 멈추지 않고 manager의 기존 자식·스레드와 미래 자식을 nice19로 낮춘다."""
    try:
        os.setpriority(os.PRIO_PROCESS,manager,19)
        for _ in range(3):
            parents={}
            for folder in Path('/proc').iterdir():
                if not folder.name.isdigit():continue
                try:
                    raw=(folder/'stat').read_text();fields=raw[raw.rfind(')')+2:].split()
                    parents[int(folder.name)]=int(fields[1])
                except (OSError,ValueError):pass
            tree={manager}
            while True:
                expanded=tree|{pid for pid,parent in parents.items() if parent in tree}
                if expanded==tree:break
                tree=expanded
            for pid in tree:
                for task in (Path('/proc')/str(pid)/'task').glob('*'):
                    try:os.setpriority(os.PRIO_PROCESS,int(task.name),19)
                    except ProcessLookupError:pass
        return process_alive(manager) and all(
            os.getpriority(os.PRIO_PROCESS,int(task.name))==19
            for pid in tree for task in (Path('/proc')/str(pid)/'task').glob('*'))
    except (OSError,ValueError):
        return False


def relend_fixed_prepare(root,leases,budget,loan_capacity,observed):
    """고정 자원 BAM 준비의 CPU/RAM 수는 유지하고 정상 예약만 반환한다. pool.lock 필요."""
    converted=[]
    for row in leases:
        cpus=int(row.get('cpus',0));pid=int(row.get('manager_pid') or 0)
        if (row.get('branch')!='WGS' or row.get('stage')!='01'
                or row.get('cpu_backfill') or int(row.get('memory_bytes',0))<1
                or row.get('stage_storage_wait') or not 0<cpus<=min(budget,loan_capacity)):
            continue
        ticket=row.get('ticket');path=root/'leases'/((ticket or '')+'.json')
        if not ticket or not path.is_file() or not process_alive(pid):continue
        if not eligible_waiter(dict(row,pid=pid)):continue
        if not in_cpu_group(pid,observed['cpu_cgroup']) or not lower_process_tree_priority(pid):continue
        row['cpu_backfill']={'schema':'mfm_idle_cpu_loan_v1','cpus':cpus,'nice':19,
            'observed_used_cpus':observed['used_cpus'],'observed_peak_cpus':observed['peak_used_cpus'],
            'cpu_cgroup':observed['cpu_cgroup'],'converted_fixed_prepare':True}
        _atomic(path,row)
        with (root/'cpu_backfill_events.jsonl').open('a') as handle:
            handle.write(json.dumps(row,ensure_ascii=False)+'\n')
        converted.append({'sample':row['sample'],'cpus':cpus,'memory_bytes':row['memory_bytes']})
        budget-=cpus;loan_capacity-=cpus
    return converted


def lend_idle_cpus(root,total,observed,observer,guard,*,preprocessing=None):
    """pool.lock 안에서 호출한다. 기본 예약과 별도로 정책 상한 안에서 대여한다."""
    try:
        policy=json.loads((root/'policy.json').read_text())
    except FileNotFoundError:
        policy={}
    setting=policy.get('cpu_backfill',{})
    requested_nice=max(0,min(19,int(setting.get('nice',19))))
    report={'enabled':setting.get('enabled') is True,'launched':[],**observed}
    if not report['enabled'] or not observed.get('ready'):
        return report
    if guard.get('action') != 'continue':
        report['reason']='RSS 보호 중';return report
    now=time.monotonic()
    events=root/'cpu_backfill_events.jsonl'
    if now<observer.cooldown_until or (events.exists() and time.time()-events.stat().st_mtime<10):
        report['reason']='새 대여 작업의 CPU 사용 안정화 관측 중';return report
    leases=[row for _,row in read_pool_records(root/'leases')]
    borrowed=sum(logical_cpu_claim(r) for r in leases if r.get('cpu_backfill'))
    multiplier=max(1,min(2,float(setting.get('grant_multiplier',1))))
    claim_max,loan_max=loan_limits(total,setting)
    idle=int(observed['idle_cpus'])
    if setting.get('use_average_idle') is True:
        # hard quota는 전체 CPU 한도를 유지하며 짧은 peak가 전체 대여를 막지 않게 한다.
        idle=max(0,math.floor(total-float(observed['used_cpus'])-1))
    if (setting.get('relend_fixed_prepare') is True and idle>0
            and next((root/'waiters').glob('*.json'),None) is not None):
        converted=relend_fixed_prepare(root,leases,int(idle*multiplier),
                                   max(0,loan_max-borrowed),observed)
        report['converted_fixed_prepare']=converted
        borrowed+=sum(row['cpus'] for row in converted)
        if converted:
            observer.cooldown_until=now+10
            report['reason']='고정 BAM 준비 계산·RAM을 유지하고 정상 CPU 예약을 반환함'
            return report
    free_cpu=min(int(idle*multiplier),available_loan_cpus(total,setting,leases))
    # 다른 프로젝트 설정을 바꾸지 않고, 명시한 전처리 트리의 예약 미사용분만 재사용함.
    from .pool_cpu_observer import reuse_capacity
    reuse = reuse_capacity(total,setting,preprocessing)
    if setting.get('reservation_reuse',{}).get('enabled') is True:
        report['reservation_reuse']=preprocessing or {'ready':False}
        free_cpu=max(free_cpu,min(reuse,available_loan_cpus(total,setting,leases)))
    report.update(logical_cpu_ceiling=claim_max,loan_cpu_ceiling=loan_max,
                  allocation_basis='reservation_reuse' if reuse>int(idle*multiplier) else 'cgroup_idle')
    free_ram=available_memory(root,leases,observed=observed,rss_bytes=int(guard['rss_bytes']))
    report.update(available_idle_cpus=free_cpu,available_memory_bytes=free_ram)
    preferred={sample:index for index,sample in enumerate(setting.get('preferred_samples',[]))}
    rows=sorted(read_pool_records(root/'waiters'),key=lambda item:(preferred.get(item[1].get('sample'),len(preferred)),str(item[0])))
    rows=[(p,r) for p,r in rows if time.time()-r.get('queued_at',time.time())>=.35
          and eligible_waiter(r) and int(r.get('memory_bytes',0))>0
          and (r.get('branch')!='WGS' or r.get('stage') not in {'02','03'} or setting.get('validation_phases') is True)
          and (r.get('branch')!='WGS' or r.get('stage')!='06' or setting.get('haplotag_phases') is True)]
    rows=[(p,r) for p,r in rows if r.get('branch')!='WGS' or r.get('stage') not in {'08','10'}
          or setting.get('realign_phases') is True]
    candidates=fitting_waiters(rows,free_cpu,free_ram,
                               stage_limits=policy.get('stage_concurrency'),active=leases,
                               prepare_read_only=policy.get('prepare_phased') is True)
    for index,(path,row) in enumerate(candidates):
        grant=min(int(row['max_cpus']),free_cpu,max(int(row['min_cpus']),free_cpu//(len(candidates)-index)))
        if grant<int(row['min_cpus']):
            continue
        if row['branch']=='WGS' and row['stage'] in {'08','10'}:
            grant=int(row['min_cpus'])
        pid=int(row['pid'])
        if not in_cpu_group(pid,observed['cpu_cgroup']):
            continue
        try:
            # 대기 manager가 자식을 만들기 전에 하향한다. 자식·스레드가 상속한다.
            # 이미 낮아진 우선순위는 권한 없이 복구할 수 없으므로 현재 값을 유지한다.
            nice=max(requested_nice,os.getpriority(os.PRIO_PROCESS,pid))
            os.setpriority(os.PRIO_PROCESS,pid,nice)
            if os.getpriority(os.PRIO_PROCESS,pid)!=nice:
                continue
        except (OSError,ProcessLookupError):
            continue
        loan={'schema':'mfm_idle_cpu_loan_v1','cpus':grant,'nice':nice,
              'observed_used_cpus':observed['used_cpus'],
              'observed_peak_cpus':observed['peak_used_cpus'],
              'cpu_cgroup':observed['cpu_cgroup'],
              'allocation_basis':report['allocation_basis'],
              'preprocessing_used_cpus':(preprocessing or {}).get('used_cpus')}
        if row['branch']=='WGS' and row['stage'] in {'08','10'}:
            loan['kind']='realign_entry'
        if row['branch']=='WGS' and row['stage']=='02':
            loan['phase_reacquire_cpus']=grant
        lease=dict(ticket=path.stem,manager_pid=pid,child_pid=None,cpus=grant,total_cpus=total,
                   branch=row['branch'],stage=row['stage'],sample=row['sample'],
                   acquired_at=time.time(),memory_bytes=int(row['memory_bytes']),cpu_backfill=loan)
        _atomic(root/'leases'/path.name,lease);path.unlink(missing_ok=True)
        with (root/'cpu_backfill_events.jsonl').open('a') as handle:
            handle.write(json.dumps(lease,ensure_ascii=False)+'\n')
        report['launched'].append({'sample':row['sample'],'branch':row['branch'],'stage':row['stage'],
                                  'cpus':grant,'nice':nice,'memory_bytes':row['memory_bytes']})
        free_cpu-=grant;free_ram-=int(row['memory_bytes'])
    if report['launched']:
        observer.cooldown_until=now+10
        report['reason']=('전처리 예약 미사용분을 한도 내 재사용함'
                          if report['allocation_basis']=='reservation_reuse' else '유휴 CPU를 정책 우선순위로 대여함')
    elif not candidates:
        report['reason']='유휴 CPU·RAM 안에 들어가는 고정 자원 단계 대기 없음'
    return report
