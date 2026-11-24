"""동시HAP의실제하위단계요청을합산해정렬중인다른HAP몫을보존함."""
from pathlib import Path
import fcntl,json,threading,time,math,os
from .phases import _atomic
from ..memory_policy import GIB
from ..pool_leases import prune_leases, read_pool_records, regular_cpu_claim
from ..cpu_loan_policy import available_regular_cpus
from ..memory_policy import available_memory
from ..resource_fairness import fair_cpu_capacity
from .idle_resources import loan_capacity, commit_loan

class HaplotypeResources:
    def __init__(self,phase,config,haps,local_root=None):
        self.phase=phase;self.lock=threading.RLock();self.events={h:[] for h in haps}
        self.path=phase.path
        self.local=not self.path and local_root is not None
        if self.local:
            self.path=Path(local_root)/'leases/hap.json';self.path.parent.mkdir(parents=True,exist_ok=True)
            _atomic(self.path,{'cpus':len(haps),'memory_bytes':len(haps)*4*GIB,'total_cpus':phase.total})
        if self.path:
            row=json.loads(self.path.read_text());cpu=int(row['cpus'])//len(haps);ram=math.ceil(int(row['memory_bytes'])/len(haps))
        else:
            cpu=phase.total//len(haps);ram=phase.ceiling()//len(haps)
        if os.environ.get('MFM_WGS_BACKEND')=='parabricks':cpu,ram=1,4*GIB
        self.states={h:{'cpus':cpu,'memory_bytes':ram,'phase':'align_restore_sort_qc'} for h in haps}
        if self.path:
            with (self.path.parent.parent/'pool.lock').open('a+') as handle:
                fcntl.flock(handle,fcntl.LOCK_EX)
                row=json.loads(self.path.read_text())
                row['hap_resources']={str(h):s for h,s in self.states.items()}
                row['hap_resources_version']=2
                row.update(cpus=sum(s['cpus'] for s in self.states.values()),memory_bytes=max(4*GIB,sum(s['memory_bytes'] for s in self.states.values())))
                _atomic(self.path,row)
    def refresh(self,row):
        # GPU 자식의 대기/계산 예약 변경을 부모의 오래된 상태로 덮지 않는다.
        for h,state in row.get('hap_resources',{}).items():
            if int(h) in self.states:self.states[int(h)]=dict(state)
    def controller(self,hap):return HapController(self,hap)
    def reserve(self,hap,label,maximum,minimum,memory_per_unit,base_memory):
        started=time.time()
        while True:
            with self.lock:
                path=self.path
                handle=(path.parent.parent/'pool.lock').open('a+') if path else None
                try:
                    if handle:
                        fcntl.flock(handle,fcntl.LOCK_EX)
                        if not self.local:prune_leases(path.parent.parent)
                        row=json.loads(path.read_text());self.refresh(row)
                        others=[] if self.local else [r for p,r in read_pool_records(path.parent) if p!=path]
                    else:row={};others=[]
                    siblings=[s for h,s in self.states.items() if h!=hap]
                    if path and not self.local:
                        root=path.parent.parent
                        policy=json.loads((root/'policy.json').read_text()) if (root/'policy.json').exists() else {}
                        cpu_room=available_regular_cpus(self.phase.total,policy.get('cpu_backfill',{}),others)
                        cpu_room=fair_cpu_capacity(root,self.phase.total,[*others,row],cpu_room,request=row)
                        memory_room=available_memory(root,others,ceiling_bytes=self.phase.ceiling())
                    else:
                        cpu_room=self.phase.total;memory_room=self.phase.ceiling()
                    cpu=max(0,cpu_room-sum(s['cpus'] for s in siblings))
                    ram_other=sum(s['memory_bytes'] for s in siblings)
                    ram=max(0,memory_room-ram_other-base_memory)
                    units=min(maximum,cpu,ram//memory_per_unit if memory_per_unit else maximum)
                    enough=units>=minimum and memory_room-ram_other>=base_memory
                    loan=None
                    if path and not self.local and units<maximum:
                        other_haps=[s for h,s in self.states.items() if h!=hap]
                        capacity=loan_capacity(path.parent.parent,self.phase.total,row,others,
                            memory_room)
                        if capacity is not None:
                            loan_cpu,loan_ram,candidate=capacity
                            available_cpu=max(0,loan_cpu-sum(s['cpus'] for s in other_haps))
                            available_ram=loan_ram-sum(s['memory_bytes'] for s in other_haps)
                            loan_units=min(maximum,available_cpu,max(0,available_ram-base_memory)//memory_per_unit if memory_per_unit else maximum)
                            if loan_units>=minimum and loan_units>units and available_ram>=base_memory:
                                units,enough,loan=loan_units,True,candidate
                    prior_wait=self.states[hap].get('wait',{})
                    self.states[hap]={'cpus':units if enough else 0,'memory_bytes':base_memory+units*memory_per_unit if enough else 4*GIB,'phase':label if enough else label+':waiting'}
                    if not enough:
                        self.states[hap]['wait']={'queued_at':prior_wait.get('queued_at',started),
                            'minimum_cpus':sum(s['cpus'] for s in siblings)+minimum,
                            'memory_bytes':ram_other+base_memory+minimum*memory_per_unit}

                    if path:
                        row.update(cpus=max(1,sum(s['cpus'] for s in self.states.values())),memory_bytes=max(4*GIB,sum(s['memory_bytes'] for s in self.states.values())),hap_resources={str(h):s for h,s in self.states.items()})
                        pending=[s['wait'] for s in self.states.values() if s.get('wait')]
                        if pending:row['resource_wait']=min(pending,key=lambda r:r['queued_at'])
                        else:row.pop('resource_wait',None)
                        if enough and not self.local:commit_loan(path.parent.parent,row,loan)
                        _atomic(path,row)
                    if enough:
                        event={'phase':label,'hap':hap,'cpus':units,'memory_bytes':self.states[hap]['memory_bytes'],'queue_wait_seconds':time.time()-started};self.events[hap].append(event);return units
                finally:
                    if handle:handle.close()
            time.sleep(.5)
    def release(self,hap):
        with self.lock:
            if self.path:
                path=self.path
                with (path.parent.parent/'pool.lock').open('a+') as handle:
                    fcntl.flock(handle,fcntl.LOCK_EX);row=json.loads(path.read_text());self.refresh(row)
                    self.states[hap]={'cpus':0,'memory_bytes':0,'phase':'complete'}
                    row.update(cpus=max(1,sum(s['cpus'] for s in self.states.values())),memory_bytes=max(4*GIB,sum(s['memory_bytes'] for s in self.states.values())),hap_resources={str(h):s for h,s in self.states.items()});_atomic(path,row)
            else:self.states[hap]={'cpus':0,'memory_bytes':0,'phase':'complete'}

class HapController:
    def __init__(self,owner,hap):self.owner=owner;self.hap=hap;self.total=owner.phase.total
    def reserve(self,label,*,max_units,min_units=1,memory_per_unit=0,base_memory=4*GIB):return self.owner.reserve(self.hap,label,max_units,min_units,memory_per_unit,base_memory)
    @property
    def resource_limits(self):
        if self.owner.path:self.owner.refresh(json.loads(self.owner.path.read_text()))
        return dict(self.owner.states[self.hap])
