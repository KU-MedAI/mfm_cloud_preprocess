"""trace 임시 multiset·재생파일의 공간과 RAM scratch 사용을 사전 추정함."""
from pathlib import Path
import shutil


def is_memory_filesystem(path):
    path=Path(path).resolve();best=('',False)
    try:lines=Path('/proc/self/mountinfo').read_text().splitlines()
    except OSError:return False
    for line in lines:
        left,right=line.split(' - ',1);mount=left.split()[4].replace('\\040',' ');kind=right.split()[0]
        if path.is_relative_to(mount) and len(mount)>len(best[0]):best=(mount,kind in {'tmpfs','ramfs'})
    return best[1]


def scratch_estimate(reads,segments,cigar_ops):
    replay=reads*256+segments*64+cigar_ops*8
    raw_digest=(reads+segments)*2*32
    return {'replay_estimate_bytes':replay,'partitioned_estimate_bytes':(replay+3*raw_digest)*11//10,
            'sqlite_estimate_bytes':(replay+(reads+segments)*112)*11//10}


def select_backend(scratch,estimates,requested,*,workers=None,elastic=False):
    """전체 누적량 대신 실제 동시 실행 상한의 scratch 피크를 계획함."""
    if requested not in {'auto','partitioned','sqlite'}:raise ValueError('multiset backend 오류')
    if workers is not None and (isinstance(workers,bool) or int(workers)<1):
        raise ValueError('trace workers는 양의 정수여야 함')
    free=shutil.disk_usage(scratch).free
    budget=free*9//10
    initial=min(len(estimates),max(1,int(workers or len(estimates))))
    maximum=len(estimates) if elastic else initial
    prefixes={}
    capacities={}
    for mode in ['partitioned','sqlite']:
        sizes=sorted((int(e[mode+'_estimate_bytes']) for e in estimates),reverse=True)
        values=[0]
        for size in sizes:values.append(values[-1]+size)
        prefixes[mode]=values
        capacities[mode]=max((n for n in range(1,maximum+1) if values[n]<=budget),default=0)
    selected=requested
    if selected=='auto':
        if capacities['partitioned']>=initial:selected='partitioned'
        elif capacities['sqlite']>=initial:selected='sqlite'
        else:selected=max(['partitioned','sqlite'],key=lambda mode:capacities[mode])
    capacity=capacities[selected]
    if capacity<1:
        minimum=prefixes[selected][1] if estimates else 0
        raise ValueError(f'trace scratch 단일 작업 예상량도 여유 공간보다 큼: backend={selected} needed={minimum} free={free} usable={budget}')
    return {'backend':selected,'requested_backend':requested,'free_bytes':free,
            'estimated_scratch_bytes':prefixes[selected][capacity],
            'all_tasks_scratch_bytes':prefixes[selected][-1],
            'requested_workers':initial,'max_concurrent_workers':capacity,
            'elastic_requested':bool(elastic),'memory_filesystem':is_memory_filesystem(scratch),
            'scope':'동시 실행 가능한 가장 큰 HAP 작업들의 피크 추정임. 종료된 작업 임시는 재사용함. 다른 샘플의 디스크 예약을 강제하지 않음.'}
