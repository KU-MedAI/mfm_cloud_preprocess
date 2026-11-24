"""완료 모델 입력의 DNA N/P·WGS/RNA 축을 읽기 전용으로 병렬 검증함."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import csv,json,multiprocessing,os,resource,time
import numpy as np
from numcodecs import blosc
from .memory_policy import limit_worker_address_space
from .pipeline.checkpoints import atomic_json
from .vcf_execution import resource_plan
from .window_reader import WindowReader

_READER=None


def _init(root,ceiling):
    global _READER
    limit_worker_address_space(ceiling);blosc.set_nthreads(1)
    _READER=WindowReader(root,max_cached_contigs=1)


def _axis(task):
    contig,hap,rows,maximum=task;start=time.process_time()
    lookup={(int(r['raw_start']),int(r['actual_length'])):r for r in rows}
    try:
        if len(lookup)!=len(rows):raise ValueError('후보 원장의 중복 개인 좌표')
        n=_READER.window_count(contig,hap)
        if n!=len(rows):raise ValueError('원장/저장 window 수 불일치')
        indices=range(n) if not maximum else np.linspace(0,n-1,min(n,maximum),dtype=int).tolist()
        totals={'windows':0,'actual_bases':0,'n_bases':0,'padding_bases':0,'n_windows':0,'padding_windows':0,'rna_windows':0}
        for i in indices:
            data=_READER.load(contig,hap,i);length=data['actual_length']
            row=lookup[(data['raw_start'],length)]
            actual_n=int(np.count_nonzero(data['dna_codes']==4));padding=int(np.count_nonzero(data['dna_codes']==5))
            if actual_n!=int(row['n_count']) or abs(actual_n/length-float(row['n_fraction']))>1e-12:
                raise ValueError(f'N 원값 불일치: 저장 순번{i}')
            if int(data['dna_target_mask'].sum())!=length-actual_n or int(data['observation_mask'].sum())!=length:
                raise ValueError(f'N/P 정답 mask 불일치: 저장 순번{i}')
            if padding!=len(data['dna_codes'])-length:raise ValueError('padding 길이 불일치')
            if data['rna_depth'] is not None:
                if np.any(data['rna_depth'][length:]):raise ValueError('RNA padding 관측값 오류')
                totals['rna_windows']+=1
            totals['windows']+=1;totals['actual_bases']+=length;totals['n_bases']+=actual_n;totals['padding_bases']+=padding
            totals['n_windows']+=bool(actual_n);totals['padding_windows']+=bool(padding)
        return {'contig':contig,'hap':hap,'status':'pass',**totals,'cpu_seconds':time.process_time()-start,
                'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'cache':_READER.cache_info()}
    except Exception as exc:return {'contig':contig,'hap':hap,'status':'fail','error':f'{type(exc).__name__}: {exc}'}


def audit_sample(root,output,*,workers=None,max_per_axis=0):
    root=Path(root).resolve();output=Path(output).resolve()
    if output.is_relative_to(root):raise ValueError('완료 패키지를 보존하도록 보고서는 sample 바깥에 저장해야 함')
    if max_per_axis<0 or (workers is not None and workers<1):raise ValueError('worker/표본 한도 오류')
    parent_reader=WindowReader(root);parent_reader.close()
    axes={};candidates=0
    with (root/'metadata/wgs_windows.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            candidates+=1
            if row['kept']=='1':axes.setdefault((row['contig'],row['haplotype']),[]).append(row)
    caps=resource_plan(1);per_worker=512*1024**2;parent=1024**3
    selected=min(workers or caps['cpu_budget'],caps['cpu_budget'],len(axes),max(0,(caps['ram_budget_bytes']*10//11-parent)//per_worker))
    if selected<1:raise ValueError('모델 입력 검증 RAM/작업 예산 부족')
    plan={'cpu_budget':caps['cpu_budget'],'workers_used':selected,'tasks':len(axes),
          'ram_budget_bytes':caps['ram_budget_bytes'],'estimated_peak_rss_bytes':(parent+selected*per_worker)*11//10,
          'per_worker_metadata_bytes':64<<20,'per_worker_decoded_cache_bytes':64<<20,'worker_address_space_bytes':min(caps['ram_budget_bytes'],4*1024**3),
          'candidate_windows':candidates,'stored_windows':sum(map(len,axes.values()))}
    print('[INPUT-AUDIT] '+json.dumps(plan),flush=True);started=time.monotonic()
    tasks=[(c,h,rows,max_per_axis) for (c,h),rows in sorted(axes.items())]
    with ProcessPoolExecutor(max_workers=selected,mp_context=multiprocessing.get_context('spawn'),
                             initializer=_init,initargs=(str(root),plan['worker_address_space_bytes'])) as pool:
        results=list(pool.map(_axis,tasks))
    keys=('windows','actual_bases','n_bases','padding_bases','n_windows','padding_windows','rna_windows')
    failed=[r for r in results if r['status']!='pass']
    report={'schema':'mfm_model_input_audit_v1','status':'fail' if failed else 'pass','sample':root.name,
            'root':str(root),'scope':'모든 저장 window' if not max_per_axis else '축별 표본',
            'resources':plan,'wall_seconds':time.monotonic()-started,
            'totals':{k:sum(r.get(k,0) for r in results) for k in keys},'axes':results,
            'errors':failed,'read_only':True}
    atomic_json(output,report);return report
