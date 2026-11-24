"""고정 목록을 층화 분담하고 모델 입력 용량별 실행 배치를 만든다."""
from __future__ import annotations
import csv, hashlib, json, math
from collections import Counter, defaultdict
from pathlib import Path

BATCH_ASSIGNMENT_SEED = 'mfm-stratified-batches-v2'
ASSIGNMENT_SEED = 'mfm-tgil-medai-20261001'
# 연구실 표시명 변경이 기존 sample 소유권을 바꾸지 않도록 최초 분담 순서를
# 만든 byte salt만 불투명 값으로 보존한다. 출력 metadata에는 새 seed 이름만 쓴다.
_ASSIGNMENT_RANK_SALT = bytes.fromhex('6d666d2d6b696d2d3230323631303031')
_BATCH_LAB_RANK_SALTS = {
    'TGIL': bytes.fromhex('eab980ec97b0eab5acec8ba4'),
    'MedAI': bytes.fromhex('ec9ab0eba6acec97b0eab5acec8ba4'),
}

def read_tsv(path):
    with Path(path).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(path, rows, fields=None):
    path=Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None: fields=list(rows[0]) if rows else []
    with path.open('w', encoding='utf-8', newline='') as f:
        w=csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

def assign(rows, fraction=.3, seed=ASSIGNMENT_SEED, reserved_medai=None):
    if not 0 < fraction < 1: raise ValueError('분담 비율은 0과 1 사이여야 합니다.')
    ids=[r['sample_id'] for r in rows]
    if len(set(ids)) != len(ids): raise ValueError('sample_id 중복')
    groups=defaultdict(list)
    for r in rows:
        if r['platform']=='clr': raise ValueError('CLR 샘플은 전달본에서 제외해야 합니다.')
        groups[(r['project_id'],r['platform'],r['data_group'])].append(dict(r))
    reserved=set(reserved_medai or [])
    if not reserved.issubset(set(ids)):raise ValueError('분담 예약 샘플이 목록에 없습니다.')
    minima={g:sum(r['sample_id'] in reserved for r in rs) for g,rs in groups.items()}
    quotas={g:max(math.floor(len(rs)*fraction),minima[g]) for g,rs in groups.items()}
    target=math.floor(len(rows)*fraction+.5)
    order=sorted(groups,key=lambda g:(-(len(groups[g])*fraction-quotas[g]),g))
    missing=target-sum(quotas.values())
    if missing>=0:
        for g in order[:missing]:quotas[g]+=1
    else:
        for _ in range(-missing):
            candidates=[g for g in groups if quotas[g]>minima[g]]
            if not candidates:raise ValueError('예약 샘플 수가 MedAI 할당을 초과합니다.')
            g=max(candidates,key=lambda g:(quotas[g]-len(groups[g])*fraction,g));quotas[g]-=1
    out=[]
    rank_salt=_ASSIGNMENT_RANK_SALT if seed==ASSIGNMENT_SEED else seed.encode()
    for g in sorted(groups):
        ranked=sorted(groups[g],key=lambda r:(r['sample_id'] not in reserved,
            hashlib.sha256(rank_salt+b':'+r['sample_id'].encode()).hexdigest()))
        for i,r in enumerate(ranked):
            r.update(assigned_lab='MedAI' if i<quotas[g] else 'TGIL',assignment_seed=seed)
            out.append(r)
    return sorted(out,key=lambda r:r['sample_id'])

def batches(rows, max_samples=220, target_bytes=10_000_000_000_000,
            seed=BATCH_ASSIGNMENT_SEED):
    """연구실 안에서 project/platform/WGS-RNA 층과 용량을 균형화한다.

    예전 구현은 project/sample 순서대로 앞 batch부터 채워 long-read와 RNA pair가
    일부 batch에 몰렸다. 현재 구현은 연구실별 batch 수를 먼저 고정한 뒤 각 층의
    큰 sample부터 해당 층 수·총 byte·총 sample 수가 작은 batch에 배치한다.
    """
    if max_samples < 1 or target_bytes < 1: raise ValueError('샘플 수와 용량 목표는 양수여야 합니다.')
    planned=[]
    for lab in sorted({r['assigned_lab'] for r in rows}):
        selected=[dict(r) for r in rows if r['assigned_lab']==lab]
        estimates=[int(r['estimated_model_bytes']) for r in selected]
        if any(x<=0 for x in estimates):
            bad=next(r['sample_id'] for r,x in zip(selected,estimates) if x<=0)
            raise ValueError('샘플 용량 추정이 누락되었습니다: '+bad)
        batch_count=max(math.ceil(len(selected)/max_samples),math.ceil(sum(estimates)/target_bytes))
        bins=[[] for _ in range(batch_count)]
        counts=[0]*batch_count; sizes=[0]*batch_count
        strata_counts=[Counter() for _ in range(batch_count)]
        groups=defaultdict(list)
        for row in selected:
            groups[(row['project_id'],row['platform'],row['data_group'])].append(row)
        for stratum in sorted(groups,key=lambda g:(-len(groups[g]),g)):
            lab_salt=_BATCH_LAB_RANK_SALTS.get(lab,lab.encode())
            token=seed.encode()+b':'+lab_salt+b':'+':'.join(stratum).encode()
            offset=int(hashlib.sha256(token).hexdigest()[:16],16)%batch_count
            ranked=sorted(groups[stratum],key=lambda r:(
                -int(r['estimated_model_bytes']),
                hashlib.sha256((seed+':'+r['sample_id']).encode()).hexdigest()))
            for row in ranked:
                estimate=int(row['estimated_model_bytes'])
                candidates=[i for i in range(batch_count)
                    if counts[i]<max_samples and sizes[i]+estimate<=target_bytes]
                if not candidates:
                    raise ValueError(f'층화 batch 상한 안에 배치할 수 없음: lab={lab} sample={row["sample_id"]}')
                index=min(candidates,key=lambda i:(
                    strata_counts[i][stratum], sizes[i]/target_bytes,
                    counts[i]/max_samples, (i-offset)%batch_count))
                bins[index].append(row); counts[index]+=1; sizes[index]+=estimate
                strata_counts[index][stratum]+=1
        if any(not bucket for bucket in bins):raise ValueError(f'빈 batch 생성: lab={lab}')
        for stratum in groups:
            spread=[c[stratum] for c in strata_counts]
            if max(spread)-min(spread)>1:
                raise ValueError(f'층화 편차가 1을 초과: lab={lab} stratum={stratum} counts={spread}')
        for index,bucket in enumerate(bins,1):
            for row in sorted(bucket,key=lambda r:r['sample_id']):
                row.update(batch_id=f'{index:03d}',batch_max_samples=max_samples,
                           batch_target_bytes=target_bytes,
                           batching_strategy='deterministic_stratified_best_fit_v2',
                           batch_assignment_seed=seed,
                           batch_stratum='|'.join((row['project_id'],row['platform'],row['data_group'])))
                planned.append(row)
    return planned

def batch_composition(rows):
    """동결 목록과 함께 배포할 batch별 구성 감사표를 만든다."""
    out=[]
    for (lab,batch),selected in sorted(
            ((key,[r for r in rows if (r['assigned_lab'],r['batch_id'])==key])
             for key in {(r['assigned_lab'],r['batch_id']) for r in rows})):
        projects=Counter(r['project_id'] for r in selected)
        platforms=Counter(r['platform'] for r in selected)
        groups=Counter(r['data_group'] for r in selected)
        types=Counter(r['read_type'] for r in selected)
        out.append({'assigned_lab':lab,'batch_id':batch,'sample_count':len(selected),
          'estimated_model_bytes':sum(int(r['estimated_model_bytes']) for r in selected),
          'estimated_model_tb':f'{sum(int(r["estimated_model_bytes"]) for r in selected)/1e12:.3f}',
          'short_read_count':types['short_read'],'long_read_count':types['long_read'],
          'wgs_count':groups['wgs'],'wgs_rna_pair_count':groups['wgs_rna_pair'],
          'project_counts':';'.join(f'{k}:{v}' for k,v in sorted(projects.items())),
          'platform_counts':';'.join(f'{k}:{v}' for k,v in sorted(platforms.items()))})
    return out

def build(root, out=None, calibration=None):
    root=Path(root); out=Path(out or root/'manifests/분담')
    rows=read_tsv(root/'manifests'/('handoff_preprocess_manifest.tsv' if (root/'09_personal_fasta.sh').is_file() else 'preprocess_manifest.tsv'))
    cap={r['sample_id']:r for r in read_tsv(root/'manifests/capacity_catalog_260929_complete/sample_capacity_manifest.tsv')}
    calibration=calibration or {}
    for r in rows:
        c=cap[r['sample_id']]
        measured=calibration.get(r['platform'],{}).get(r['data_group'])
        r['estimated_model_bytes']=int(measured or c['estimated_compact_total_bytes_high'])
        r['capacity_basis']='동일 플랫폼/분기 파일럿 실측의 상한' if measured else '기존 용량 목록 상한 추정'
        r['capacity_confidence']='파일럿 한정' if measured else '잠정: 같은 형식 실측 부족'
    pilots={'HG00096','HG03492','NA20503','NA12878','HG00281'}
    planned=batches(assign(rows,reserved_medai=pilots & {r['sample_id'] for r in rows}))
    write_tsv(out/'전체_분담배치.tsv',planned)
    for lab,slug in [('MedAI','medai'),('TGIL','tgil')]:
        selected=[r for r in planned if r['assigned_lab']==lab]
        write_tsv(out/f'{slug}_samples.tsv',selected)
        for batch in sorted({r['batch_id'] for r in selected}):
            write_tsv(out/slug/f'batch_{batch}.tsv',[r for r in selected if r['batch_id']==batch])
    composition=batch_composition(planned)
    write_tsv(out/'batch_composition.tsv',composition)
    report={'sample_count':len(rows),'lab_counts':dict(Counter(r['assigned_lab'] for r in planned)),
        'max_samples':220,'target_model_bytes':10_000_000_000_000,
        'batch_counts':{lab:len({r['batch_id'] for r in planned if r['assigned_lab']==lab}) for lab in ['MedAI','TGIL']},
        'batching_strategy':'deterministic_stratified_best_fit_v2',
        'batch_assignment_seed':BATCH_ASSIGNMENT_SEED,
        'batch_strata':['project_id','platform','data_group'],
        'batch_composition':'batch_composition.tsv',
        'capacity_note':'10TB는 모델 입력의 목표다. 다운로드·중간 산출물·동시 진행 작업은 별도 공간이 필요하다.'}
    (out/'분담요약.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    return report
