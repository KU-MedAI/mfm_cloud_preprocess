"""원래 영역·행은 자르지 않고 인접 항목을 묶으며 원래 행 번호를 보존함."""
import os,resource,time
from .gpu_compare_resources import cpu_limit


def task_layout():
    mode=os.environ.get('RNA_WGS_TASK_LAYOUT','auto')
    if mode not in {'auto','contig','chunks'}:raise ValueError('RNA_WGS_TASK_LAYOUT는 auto/contig/chunks임')
    if mode=='auto':mode='chunks' if cpu_limit(int(os.environ.get('MFM_CPU_BUDGET','68')))>=8 else 'contig'
    return mode


def partition_rows(rows,bounds):
    if not rows:return []
    if task_layout()=='contig':return [list(range(len(rows)))]
    span=int(os.environ.get('RNA_WGS_CHUNK_BP','16000000'))
    if span<1:raise ValueError('RNA_WGS_CHUNK_BP는 양수임')
    ordered=sorted(range(len(rows)),key=lambda i:(bounds(rows[i])[0],bounds(rows[i])[1],i))
    result=[];chunk=[];left=right=0
    for i in ordered:
        start,end=bounds(rows[i])
        if start<0 or end<start:raise ValueError('분할 입력 좌표 오류')
        if chunk and max(right,end)-left>span:
            result.append(chunk);chunk=[]
        if not chunk:left=start
        right=max(right,end) if chunk else end
        chunk.append(i)
    if chunk:result.append(chunk)
    return result


def cpu_seconds():
    child=resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time()+child.ru_utime+child.ru_stime


def timing(started,cpu,**fields):
    return dict(**fields,pid=os.getpid(),wall_seconds=time.monotonic()-started,cpu_seconds=cpu_seconds()-cpu)
