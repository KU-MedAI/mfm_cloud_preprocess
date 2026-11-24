"""기관배정과배치인원을보존하며선택기관의long/pair를첫배치에모음."""
from collections import Counter
import copy

def priority(row):
    long=row.get('read_type')=='long_read'
    return (not long, row.get('data_group')!='wgs_rna_pair' if long else False,
            int(row.get('known_input_total_bytes') or 0),row['sample_id'])

def rebatch(rows,lab='MedAI'):
    result=copy.deepcopy(rows)
    chosen=[r for r in result if r['assigned_lab']==lab]
    if not chosen:raise ValueError('기관 샘플이 없음')
    if len({r['sample_id'] for r in result})!=len(result):raise ValueError('중복 sample ID')
    batches=sorted({r['batch_id'] for r in chosen});first=batches[0]
    counts=Counter(r['batch_id'] for r in chosen)
    long=sorted([r for r in chosen if r['read_type']=='long_read'],key=priority)
    if len(long)>counts[first]:raise ValueError('첫 배치 인원보다 long이 많음; 배치 설계 재검토 필요')
    buckets={b:[] for b in batches};buckets[first]=long
    displaced=[]
    for batch in batches:
        short=sorted([r for r in chosen if r['batch_id']==batch and r['read_type']!='long_read'],key=priority)
        keep=counts[batch]-len(buckets[batch]);buckets[batch]+=short[:keep];displaced+=short[keep:]
    for batch in batches:
        needed=counts[batch]-len(buckets[batch]);buckets[batch]+=displaced[:needed];displaced=displaced[needed:]
    if displaced:raise ValueError('남는 샘플 있음')
    ordered=[]
    for batch in batches:
        bucket=buckets[batch]
        if len(bucket)!=counts[batch]:raise ValueError('배치 인원 불일치')
        if len(bucket)>min(int(r.get('batch_max_samples') or 220) for r in bucket):raise ValueError('배치 최대 인원 초과')
        if sum(int(r.get('estimated_model_bytes') or 0) for r in bucket)>min(int(r.get('batch_target_bytes') or 10**13) for r in bucket):raise ValueError('배치 예상 모델 용량 상한 초과')
        for row in sorted(bucket,key=priority):
            row['batch_id']=batch;row['batching_strategy']='long_pair_first_v1';row['batch_assignment_seed']='mfm-long-pair-first-261027';ordered.append(row)
    untouched=[r for r in result if r['assigned_lab']!=lab]
    assert {r['sample_id'] for r in result}=={r['sample_id'] for r in ordered+untouched}
    return ordered+untouched
