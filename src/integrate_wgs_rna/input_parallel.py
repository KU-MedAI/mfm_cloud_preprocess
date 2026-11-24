"""입력 검증의 외부 명령·해시 작업을 CPU/RAM 예약 안에서 병렬 실행함."""
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import os,time
from .cpu_executor import WorkerBudget
from .memory_policy import GIB,memory_limit


def cpu_budget():
    requested=int(os.environ.get('MFM_CPU_BUDGET',os.environ.get('MFM_SAMPLE_CPU_BUDGET','1')))
    granted=int(os.environ.get('MFM_CPU_POOL_GRANTED',requested))
    return max(1,min(requested,granted,len(os.sched_getaffinity(0))))


def parallel_calls(function,tasks,workers,*,cpu_cost=1,memory_per_task=2*GIB):
    if not tasks:return []
    count=min(len(tasks),max(1,workers))
    reserve=WorkerBudget(count,count,True,memory_per_task,cpus_per_worker=cpu_cost)
    capacity=reserve.capacity(len(tasks))
    while capacity<1:
        time.sleep(.2);capacity=reserve.capacity(len(tasks))
    count=capacity
    # 실행이 끝난 thread의 외부 명령은 종료되지만 결과 처리·유휴 thread 여유도 끝까지 예약함.
    budget=WorkerBudget(count,count,False,0,cpus_per_worker=cpu_cost)
    budget.base_memory=reserve.base_memory+count*memory_per_task
    if budget.base_memory>memory_limit()*9//10:raise ValueError('입력 검증 RAM 계획이 운영 예산을 초과함')
    budget.capacity(len(tasks),grow=False)
    results=[None]*len(tasks);pending=iter(enumerate(tasks));submitted=completed=0;futures={}
    try:
        with ThreadPoolExecutor(max_workers=count) as executor:
            while completed<len(tasks):
                while submitted<len(tasks) and len(futures)<count:
                    i,args=next(pending);futures[executor.submit(function,*args)]=i;submitted+=1
                done,_=wait(futures,return_when=FIRST_COMPLETED)
                for future in done:
                    i=futures.pop(future);results[i]=future.result();completed+=1
                budget.capacity(len(tasks)-completed,grow=False)
        return results
    finally:
        budget.base_memory=reserve.base_memory;budget.capacity(0,grow=False)
