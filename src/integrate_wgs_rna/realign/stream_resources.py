"""minimap2 종료가 확인된 뒤 CPU만 반환함. 실행 중 sort RAM은 유지함."""
from pathlib import Path
import argparse,fcntl,json,os,time
from .phases import _atomic


def mapper_finished(hap, sort_threads, lease=None):
    value=lease or os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if not value:return {'applied':False,'reason':'독립 실행으로 공유 예약 없음'}
    path=Path(value)
    with (path.parent.parent/'pool.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        row=json.loads(path.read_text())
        if row.get('realign_phase')!='align_restore_sort_qc':
            return {'applied':False,'reason':'정렬 단계가 아님'}
        before=int(row['cpus']);memory=int(row.get('memory_bytes',0))
        tail=max(1,int(sort_threads)+2)  # sort main·압축 보조 thread와 부모 비용임.
        states=row.get('hap_resources')
        if states:
            state=states[str(hap)]
            state.update(cpus=min(int(state['cpus']),tail),phase='sort_output_tail')
            row['cpus']=max(1,sum(int(s['cpus']) for s in states.values()))
        else:
            row['cpus']=min(before,tail)
        if row.get('cpu_backfill'):row['cpu_backfill']['cpus']=row['cpus']
        event=dict(applied=True,hap=int(hap),time=time.time(),cpus_before=before,
                   cpus_after=row['cpus'],memory_bytes=memory)
        row.setdefault('stream_tail_events',[]).append(event)
        _atomic(path,row)
        return event


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--hap',type=int,choices=[1,2],required=True)
    p.add_argument('--sort-threads',type=int,required=True);a=p.parse_args()
    if a.sort_threads<0:p.error('sort thread는 0 이상이어야 함')
    print('[REALIGN-STREAM] '+json.dumps(mapper_finished(a.hap,a.sort_threads),ensure_ascii=False),flush=True)

if __name__=='__main__':main()
