"""분담 manifest 기반 단일 진입점: 환경 준비 이후 다운로드·실행·요약."""
from __future__ import annotations
import argparse,concurrent.futures,fcntl,hashlib,json,os,shutil,subprocess,sys,time,tempfile
from collections import Counter
from pathlib import Path
from . import planner,slurm,bundled_assets,resources
from .storage_paths import storage_defaults
from . import acceleration
from .resources import cpu_budget, memory_budget
from .cleanup import cleanup_complete,model_bytes
from .resume import sample_environment,sample_downloads,rna_relative,alignment_relative,require_blacklist
from .rna_integrity import prepare_rna
from .sequence_mode import completion_mode,rna_required
from .downloads import sample_resources,fetch,fetch_or_index,remove_owned,safe_path

ROOT=Path(__file__).resolve().parents[1]
ORIGINAL_LAYOUT=(ROOT/"09_personal_fasta.sh").is_file()
CLOUD_TARGETS=('azure_poc','azure_main','aws_main','server')

def processing_complete(root,sample,mode,require_cleanup=False):
    from .archive_receipt import read as read_archive
    archived=read_archive(root)
    if archived is not None:return archived['sample']==sample and archived['mode']==mode
    if require_cleanup or not ORIGINAL_LAYOUT:return cleanup_complete(root,sample,mode)
    root=Path(root)
    if (root/'metadata/cleanup_completed.json').exists() or (root/'metadata/cleanup_in_progress.json').exists():
        return cleanup_complete(root,sample,mode)
    try:
        report=json.loads((Path(root)/"qc/final_package_validation.json").read_text())
        return (report.get("schema")=="mfm_final_package_validation_v1"
                and report.get("sample")==sample and report.get("mode")==mode
                and report.get("status")=="pass" and report.get("failed_checks")==[]
                and (Path(root)/"individual_fasta/permanent_assets.json").is_file()
                and (report.get("wgs_payload") in {"sequence_only","pileup_feature"} or (Path(root)/"realigned").is_dir()))
    except (OSError,ValueError,TypeError):return False

def atomic_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile('w',dir=path.parent,delete=False,encoding='utf-8') as f:
        json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n');temp=Path(f.name)
    os.replace(temp,path)

def load_config(path=None, resource_overrides=None):
    path=Path(path or ROOT/'handoff_config.json')
    if not path.exists():path=ROOT/'handoff_config.example.json'
    c=storage_defaults(json.loads(path.read_text()), ROOT);base=path.resolve().parent
    c.setdefault('run_id','tgil_handoff_'+time.strftime('%Y%m%d'))
    for name in ['data_root','output_root','state_root']:
        c[name]=str((base/c[name]).resolve())
    c['package_root']=str(ROOT)
    c['bundled_assets_root']=str((base/c.get('bundled_assets_root',str(ROOT/'assets'))).resolve())
    c.setdefault('local_jobs','auto');c.setdefault('sample_cpus','auto')
    c.setdefault('model_target_bytes',10_000_000_000_000)
    c.setdefault('queue_policy', {'mode':'long_read_first','pair_first':True})
    c.setdefault('remove_owned_raw_inputs',not ORIGINAL_LAYOUT);c.setdefault('cleanup_shared_fasta',not ORIGINAL_LAYOUT)
    if ORIGINAL_LAYOUT:
        policy=c.get('pipeline_env',{})
        model_only=(str(policy.get('MFM_MODEL_ONLY_CLEANUP','0'))=='1'
                    and str(policy.get('CLEAN_INTERMEDIATES_ON_SUCCESS','0'))=='1')
        c['remove_owned_raw_inputs']=bool(c['remove_owned_raw_inputs']) and model_only
        c['cleanup_shared_fasta']=False
    c=acceleration.apply(resources.resolve(c,resource_overrides))
    # Slurm array의 CPU 요청과 실제 sample 내부 예산은 항상 같다.
    if os.environ.get('SLURM_JOB_ID'):
        allocated = int(os.environ.get('SLURM_CPUS_PER_TASK', c['slurm']['cpus']))
        if c.get('cpu_allocation_mode') == 'fixed_sample':
            if c['sample_cpus'] > allocated:
                raise ValueError('고정 sample_cpus가 실제 Slurm 할당보다 큽니다.')
        else:
            c['sample_cpus'] = allocated
    roots=[Path(c[x]) for x in ['data_root','output_root','state_root']]
    if len(set(roots))!=3 or any(p==Path('/') for p in roots):
        raise ValueError('전달본의 input/output/state는 각각 별도 경로여야 합니다.')
    protected=Path(c['bundled_assets_root'])
    if any(p==protected or protected in p.parents for p in roots):
        raise ValueError('동봉 자산 경로에 data/output/state를 지정할 수 없음')
    for n in ['sample_cpus','local_jobs','model_target_bytes']:
        if not isinstance(c[n],int) or c[n]<1:raise ValueError(n+'은 양의 정수여야 합니다.')
    if int(c['slurm']['concurrency'])<1:raise ValueError('Slurm 동시 작업 수는 양수여야 합니다.')
    return c

def selection(args):
    effective=json.loads(Path(args.effective_config).read_text()) if getattr(args,'effective_config',None) else {}
    target=getattr(args,'target',None) or effective.get('selection_target')
    if target:
        if target not in CLOUD_TARGETS:raise ValueError('알 수 없는 클라우드 대상: '+str(target))
        if args.lab or args.batch or (args.selection and not args.effective_config):
            raise ValueError('--target은 --lab/--batch/--selection과 함께 지정하지 않습니다.')
        if args.effective_config and effective.get('selection_target')!=target:
            raise ValueError('저장된 실행과 --target이 다릅니다.')
        if args.effective_config and not args.selection:
            raise ValueError('저장된 클라우드 실행에는 --selection 스냅샷이 필요합니다.')
        allowed=planner.read_tsv(ROOT/'manifests/클라우드_분담'/f'{target}.tsv')
        allowed_ids={r['sample_id'] for r in allowed}
        if args.selection:
            rows=planner.read_tsv(args.selection)
            expected=set(effective.get('selection_scope_samples') or allowed_ids)
            fingerprint=hashlib.sha256('\n'.join(r['sample_id'] for r in rows).encode()).hexdigest()
            if {r['sample_id'] for r in rows}!=expected or not expected<=allowed_ids or fingerprint!=effective.get('selection_fingerprint'):
                raise ValueError('클라우드 대상의 저장된 선택 목록/순서가 다릅니다.')
        else:
            rows=allowed
            if args.sample and getattr(args,'command',None)!='resume':
                if len(set(args.sample))!=len(args.sample) or set(args.sample)-allowed_ids:
                    raise ValueError('해당 클라우드 대상에 없는 샘플 또는 중복 샘플입니다.')
                rows=[r for r in rows if r['sample_id'] in args.sample]
    elif args.selection:
        rows=planner.read_tsv(args.selection)
    else:
        if not args.lab:raise ValueError('--lab tgil|medai를 지정하세요. 전체3212 실행은 기본 동작이 아닙니다.')
        if not args.batch and not args.sample:raise ValueError('--batch 001 또는 --sample ID를 지정하세요.')
        path=ROOT/'manifests/분담'/f'{args.lab}_samples.tsv'
        rows=planner.read_tsv(path)
        if args.batch:rows=[r for r in rows if r['batch_id']==args.batch.zfill(3)]
        if args.sample:rows=[r for r in rows if r['sample_id'] in args.sample]
        if args.sample and {r['sample_id'] for r in rows}!=set(args.sample):raise ValueError('선택한 lab/batch에 속하지 않은 샘플이 있습니다.')
    if not rows or (not target and len(rows)>220):raise ValueError('기존 배치 선택은 1~220개, 클라우드 대상 선택은 비어 있지 않아야 합니다.')
    if len({r['sample_id'] for r in rows})!=len(rows):raise ValueError('샘플 중복')
    if not target and len({(r['assigned_lab'],r['batch_id']) for r in rows})!=1:raise ValueError('한 실행은 한 lab/batch만 허용합니다.')
    from .assignment_snapshot import master_path
    master={r['sample_id']:r for r in planner.read_tsv(master_path(ROOT,args))}
    for r in rows:
        if r['sample_id'] not in master or any(r.get(k)!=master[r['sample_id']].get(k) for k in ['assigned_lab','batch_id','project_id','platform']):
            raise ValueError('배정 목록과 선택 목록이 다릅니다: '+r['sample_id'])
    # 외부 selection의 URL/원본 경로는 사용하지 않고 동결 master 행을 사용한다.
    return [master[r['sample_id']] for r in rows]

def batch_config(c,rows):
    c=dict(c)
    fingerprint=hashlib.sha256('\n'.join(r['sample_id'] for r in rows).encode()).hexdigest()
    if c.get('selection_target'):
        key='cloud_'+c['selection_target']
        if c.get('selection_scope_samples'):
            key+='_'+rows[0]['sample_id'] if len(rows)==1 else '_subset_'+fingerprint[:12]
    else:
        key=('tgil' if rows[0]['assigned_lab']=='TGIL' else 'medai')+'_'+rows[0]['batch_id']
    if c.get('output_layout') != 'run_root':
        c['output_root']=str(Path(c['output_root'])/key)
    c['state_root']=str(Path(c['state_root'])/key)
    c['batch_key']=key;c['selection_fingerprint']=fingerprint
    return c

def prepare(c,rows):
    from .sequence_mode import validate_selection
    validate_selection(c,rows)
    require_blacklist(c,ROOT,rows)
    from .references import prepare as prepare_refs
    data=Path(c['data_root']);state=Path(c['state_root']);out=Path(c['output_root'])
    for p in [data,state,out]:p.mkdir(parents=True,exist_ok=True)
    marker=data/'.mfm_handoff_store.json'
    if marker.exists():
        if json.loads(marker.read_text()).get('schema')!='mfm_handoff_store_v1':raise ValueError('입력 저장소 소유권 마커가 다릅니다.')
    else:marker.write_text(json.dumps({'schema':'mfm_handoff_store_v1','package':str(ROOT),'created':time.time()})+'\n')
    binding=state/'selection_fingerprint.txt'
    if binding.exists() and binding.read_text().strip()!=c['selection_fingerprint']:raise ValueError('기존 배치와 선택 샘플이 다릅니다. 별도 state_root를 사용하세요.')
    binding.write_text(c['selection_fingerprint']+'\n')
    lease_dir=data/'references';lease_dir.mkdir(exist_ok=True)
    with (lease_dir/'active_samples.lock').open('a') as lease:
        fcntl.flock(lease,fcntl.LOCK_SH)
        from .references import needs_hgsvc_reference
        refs=prepare_refs(c,ROOT,any(rna_required(c,r) for r in rows),any(r['platform']=='illumina' for r in rows),
            needs_hgsvc=any(needs_hgsvc_reference(r) for r in rows))
        atomic_json(state/'references.json',refs)
    return refs

def bind_manifest(c,row,resources,refs):
    from .references import exact_cram_reference,variant_reference
    state=Path(c['state_root']);data=Path(c['data_root']);sample=row['sample_id'];r=dict(row)
    rf=[dict(x) for x in planner.read_tsv(ROOT/'manifests'/('handoff_read_files.tsv' if ORIGINAL_LAYOUT else 'read_files.tsv')) if x['sample_id']==sample]
    for x in rf:
        x['path']=str(safe_path(data,x['path']))
        x['source_reference_fasta']=refs['ncbi'] if x['alignment_state']=='aligned' else ''
    from .sequence_mode import needs_wgs,runtime_environment
    use_wgs=needs_wgs(c,row)
    if not use_wgs:
        r['read_files_manifest']='';relative=alignment_relative(c,row)
        r['alignment_path']=str(safe_path(data,relative)) if relative else ''
        r['source_decode_reference_fasta']=refs['autosomes']
    elif rf:
        file_manifest=state/'sample_manifests'/sample/'read_files.tsv';planner.write_tsv(file_manifest,rf)
        r['read_files_manifest']=str(file_manifest);r['alignment_path']='';r['source_decode_reference_fasta']=refs['ncbi']
    else:
        r['read_files_manifest']='';r['alignment_path']=str(safe_path(data,alignment_relative(c,row)))
        # CRAM M5를 검증하고 후보별 일치하는 전체 contig만 사용한다. HG00281 복구 정책과 같다.
        if r['alignment_path'].endswith('.cram'):
            candidates=[refs.get('illumina',str(data/'references/GRCh38_full_analysis_set_plus_decoy_hla.fa')),refs['broad'],refs['ncbi']]
            candidates+=c.get('extra_decode_references',[])
            r['source_decode_reference_fasta']=str(exact_cram_reference(r['alignment_path'],candidates,
                data/'references/cram_m5.fa',c.get('samtools','samtools')))
        else:r['source_decode_reference_fasta']=refs['ncbi']
    sequence_reference=variant_reference(row,refs)
    reference_policy=sample_environment(c,sample).get('MFM_REFERENCE_IUPAC_POLICY','preserve')
    if reference_policy not in {'preserve','to_n'}:raise ValueError('MFM_REFERENCE_IUPAC_POLICY는 preserve/to_n임')
    if reference_policy=='to_n':
        from .reference_iupac import normalized_reference
        sequence_reference=normalized_reference(sequence_reference,data/'references'/(Path(sequence_reference).stem+'.iupac_to_n.fa'))
    r['analysis_reference_fasta']=r['source_decode_reference_fasta'] if r['platform']=='illumina' else sequence_reference
    r['wgs_payload']=runtime_environment(c,row).get('MFM_WGS_PAYLOAD','full')
    r['vcf_reference_fasta']=sequence_reference;r['rna_reference_fasta']=refs['broad']
    relative_rna=rna_relative(c,row)
    r['rna_path']=str(safe_path(data,relative_rna)) if relative_rna else ''
    if r['public_vcf_path']:
        r['public_vcf_path']=';'.join(str(safe_path(data,p)) for p in r['public_vcf_path'].split(';'))
    r['common_vcf_dir']=str(bundled_assets.panel_directory(c,ROOT));r['persample_vcf_dir']=str(Path(c['output_root'])/'_work'/c['run_id']/sample/'persample_vcf')
    r['final_root']=str(Path(c['output_root'])/r['data_group']/r['project_id'])
    r['work_root']=str(Path(c['output_root'])/'_work'/c['run_id']/r['data_group']/r['project_id']/sample)
    r['raw_root']=str(Path(r['alignment_path']).parent) if r['alignment_path'] else str(data/'wgs'/r['project_id']/sample)
    p=state/'sample_manifests'/sample/'preprocess_manifest.tsv';planner.write_tsv(p,[r]);return p,r

def sample_root(c,row):return Path(c['output_root'])/row['data_group']/row['project_id']/row['sample_id']

def claim(c,sample):
    if c.get('selection_target'):return True
    state=Path(c['state_root']);state.mkdir(parents=True,exist_ok=True)
    with (state/'capacity.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        completed=sum(json.loads(p.read_text()).get('model_bytes',0) for p in (state/'samples').glob('*.json')
                      if json.loads(p.read_text()).get('status')=='success')
        if completed>=c['model_target_bytes']:return False
    return True

def worker(c,row,executor='local'):
    from .worker_resources import apply as apply_worker_resources
    c=apply_worker_resources(c,executor)
    from .io_env import apply as apply_io_environment
    apply_io_environment(c)
    from .metrics import snapshot
    state=Path(c['state_root']);sample=row['sample_id'];stats=state/'samples';stats.mkdir(parents=True,exist_ok=True)
    with (stats/(sample+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        record=stats/(sample+'.json');root=sample_root(c,row)
        mode=completion_mode(c,row)
        from .archive_receipt import matches as archived_match
        if archived_match(c,row,root) is not None:return 0
        from .sequence_mode import payload_matches
        if processing_complete(root,sample,mode) and not payload_matches(c,row,root):
            raise ValueError('기존 완료본의 WGS payload가 다름. 별도 output_root/state_root를 사용해야 함: '+sample)
        if record.exists():
            previous=json.loads(record.read_text())
            if previous.get('status')=='success':
                if processing_complete(root,sample,mode):return 0
                # 잘못된 성공 상태를 모델 바이트 목표 합계에서도 제외한다.
                previous.update(status='failed',model_bytes=0,error='cleanup 완료 보고서 또는 정리 후 PASS 계약 불일치')
                atomic_json(record,previous)
        if not claim(c,sample):
            atomic_json(record,{'sample':sample,'status':'deferred','reason':'10TB 모델 입력 목표 도달'});return 75
        atomic_json(state/'worker_resources'/(sample+'.json'),c.get('_worker_resource_policy',dict(cpu_allocation_mode=c.get('cpu_allocation_mode','dynamic'),sample_cpus=c['sample_cpus'])))
        started=time.time();info={'sample':sample,'status':'running','started':started,'model_bytes':0}
        info.update(cpu_allocation_mode=c.get('cpu_allocation_mode','dynamic'),sample_cpus=c['sample_cpus'])
        atomic_json(record,info)
        lease_root=Path(c['data_root'])/'references';lease_root.mkdir(parents=True,exist_ok=True)
        from integrate_wgs_rna.gpu_compare_resources import sample_slot
        from integrate_wgs_rna.sample_cpu import sample_slot as cpu_sample_slot
        with (lease_root/'active_samples.lock').open('a') as lease, sample_slot(c,sample), cpu_sample_slot(c,sample) as fixed_cpu:
            fcntl.flock(lease,fcntl.LOCK_SH)
            try:
                resources=sample_downloads(c,row,sample_resources(ROOT,row))
                refs=json.loads((state/'references.json').read_text())
                c,resources=prepare_rna(c,row,resources,refs.get('broad'))
                index_kinds={'wgs_index','rna_index','sample_tbi','shared_tbi','illumina_panel_index'}
                for resource in sorted(resources,key=lambda r:r['kind'] in index_kinds):
                    if resource['kind'] in {'rna_alignment','rna_index'}:continue
                    bundled=bundled_assets.get(c,ROOT,resource['path']) if resource.get('shared') else None
                    if bundled is not None:
                        resource['bundled_path']=str(bundled)
                        resource['bundled_size_bytes']=bundled.stat().st_size
                    else:
                        fetch_or_index(resource,c['data_root'],state,min(4,c['sample_cpus']))
                # 이전 배치 정리와 이번 worker 시작 사이에 공유 캐시가 제거되었을 수 있다.
                required=[refs['ncbi'],refs['autosomes']]
                if rna_required(c,row):required += [refs['broad'],refs['gtf'],refs['rna_extra']]
                if row['platform']=='illumina':required += [refs['broad'],refs.get('illumina',str(lease_root/'GRCh38_full_analysis_set_plus_decoy_hla.fa'))]
                if not all(Path(p).is_file() for p in required):refs=prepare(c,[row])
                # BAM/CRAM은 모든 압축 블록을 검사한다. 기존 중간 산출물부터 재개하는
                # 샘플의 변경 없는 WGS 원본은 다시 전수 읽지 않는다.
                if int(sample_environment(c,sample).get('WGS_START_AT','01')) <= 1:
                    from .integrity_cache import validate_shared
                    for resource in resources:
                        if resource['kind'] in {'wgs_alignment','wgs_read'} and Path(resource['path']).suffix in {'.bam','.cram'}:
                            path=safe_path(c['data_root'],resource['path'])
                            key=hashlib.sha256(str(path).encode()).hexdigest()
                            validate_shared(path,state/'wgs_input_validation'/sample/(key+'.json'),c['data_root'])
                manifest,bound=bind_manifest(c,row,resources,refs)
                env=os.environ.copy()
                if executor=='slurm':
                    affinity=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else (os.cpu_count() or c['sample_cpus'])
                    pool_total=cpu_budget(c,'slurm',affinity=affinity)
                    pool_root=Path(tempfile.gettempdir())/f"mfm_cpu_pool_{os.getuid()}_{os.environ.get('SLURM_JOB_ID','job')}_{os.environ.get('SLURM_ARRAY_TASK_ID',sample)}"
                    audit_root=state/'stats'/'cpu_pool_usage';audit_root.mkdir(parents=True,exist_ok=True)
                    pool_audit=audit_root/f"slurm_{os.environ.get('SLURM_JOB_ID','job')}_{os.environ.get('SLURM_ARRAY_TASK_ID',sample)}.jsonl"
                else:
                    affinity=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else (os.cpu_count() or c['sample_cpus'])
                    pool_total=cpu_budget(c,'local',affinity=affinity)
                    pool_root=state/'cpu_pool'/c['run_id']
                    pool_audit=state/'stats'/'cpu_pool_usage.jsonl'
                env.update({'PIPELINE_MANIFEST':str(manifest),'INPUT_ROOT':c['data_root'],
                  'OUTPUT_ROOT':c['output_root'],'OUTPUT_ROOT_OVERRIDE':'1','RUN_ID':c['run_id'],
                  'MFM_REPO_ROOT':str(ROOT/'vendor/MFM-DNARNAFM'),'MFM_SAMPLE_CPU_BUDGET':str(c['sample_cpus']),
                  'MFM_CPU_BUDGET':str(pool_total),'MFM_CPU_POOL_DIR':str(pool_root),
                  'MFM_RAM_LIMIT_BYTES':str(memory_budget(c,executor)),
                  'MFM_MEMORY_PROFILE':'balanced',
                  'MFM_CPU_POOL_AUDIT':str(pool_audit),
                  'MFM_CPU_POOL_LAYOUT':'original' if ORIGINAL_LAYOUT else 'legacy',
                  'MFM_LOCAL_JOBS':str(c['local_jobs']),
                  'MFM_RESOURCE_PROFILE':str(c['resource_profile']),
                  'MFM_CPU_UTILIZATION':c.get('cpu_utilization','maximum'),
                  'MFM_CPU_POOL_PROFILE':'throughput' if c.get('cpu_utilization','maximum')=='maximum' else c['resource_profile'],
                  'MFM_ELASTIC_SAMPLE_SLOTS':str(c['resource_resolution'].get('elastic_enabled',False)).lower(),
                  'MFM_HANDOFF_METRICS':'1','REF_SEQUENCE':refs['autosomes'],'REF_DECODE':refs['ncbi'],
                  'PANEL_BASE':bound['common_vcf_dir'],
                  'PANEL_LOCAL_DIR':bound['common_vcf_dir'],
                  'ILLUMINA_REF_DECODE':bound['source_decode_reference_fasta'],'GENCODE_GTF':refs['gtf'],
                  'RNA_ADDNL_SEQNS_SOURCE':refs['broad'],'RNA_ADDNL_SEQNS_FASTA':refs['rna_extra'],
                  'CLEAN_INTERMEDIATES_ON_SUCCESS':'0' if ORIGINAL_LAYOUT else '1','CLEAN_AFTER_RNA_ONLY':'1','SEQUENCE_SOURCE':'vcf_consensus',
                  'COMPACT_VALIDATION_LEVEL':'full','STAGE12_TRACE_POLICY':'required','ALLOW_CLEANUP_WITHOUT_TRACE':'0'})
                # Baseline config allows the pool grant to grow above the old
                # per-sample default; each stage config caps these values again.
                budget=pool_total
                env.update({'STAGE04_ALIGN_THREADS':str(max(1,budget*2//3)),
                  'STAGE04_SORT_THREADS':str(max(1,budget//3)),
                  'VCF_CHROM_JOBS':str(max(1,min(6,budget//2))),
                  'VCF_THREADS_PER_CHROM':str(max(1,min(2,budget))),
                  'VCF_JOBS':str(max(1,min(8,budget))),
                  'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1',
                  'NUMEXPR_NUM_THREADS':'1'})
                explicit_threads={
                    '11':{'RNA_COLLATE_THREADS','RNA_FASTQ_THREADS'},
                    '13':{'STAR_INDEX_THREADS'},
                    '14':{'STAR_ALIGN_THREADS'},
                    '18':{'HAP_VARIANT_THREADS'},
                }
                if ORIGINAL_LAYOUT:
                    explicit_threads={f'{int(stage)+2:02d}':names for stage,names in explicit_threads.items()}
                if c.get('rna_full_decode_report'):env['RNA_FULL_DECODE_REPORT']=c['rna_full_decode_report']
                env['RNA_DECODE_REFERENCE']=refs['broad']
                user_env=sample_environment(c,sample)
                fixed_stages=sorted(stage for stage,names in explicit_threads.items()
                                    if names & set(user_env))
                env['MFM_DYNAMIC_RNA_FIXED_STAGES']=','.join(fixed_stages)
                # 번호별 재개는 선택적으로 허용하며 원본/FASTAs가 남은 실패 샘플에서만 가능하다.
                env.update({k:str(v) for k,v in user_env.items()})
                from .sequence_mode import runtime_environment
                env.update(runtime_environment(c,row))
                env.update(fixed_cpu)
                env['MFM_CPU_BUDGET']=str(min(pool_total,int(env.get('MFM_CPU_BUDGET',pool_total))))
                from .stage_resources import longread_memory_plan
                memory_plan=longread_memory_plan(row,c)
                env.pop('MFM_LONGREAD_MEMORY_PLAN',None)
                if memory_plan is not None:
                    env['MFM_LONGREAD_MEMORY_PLAN']=json.dumps(memory_plan)
                    atomic_json(root/'metadata/longread_memory_plan.json',memory_plan)
                if env.get('COMPACT_VALIDATION_LEVEL')!='full' or env.get('STAGE12_TRACE_POLICY')!='required' or env.get('ALLOW_CLEANUP_WITHOUT_TRACE')!='0':raise ValueError('전달본은 COMPACT_VALIDATION_LEVEL=full, STAGE12_TRACE_POLICY=required, ALLOW_CLEANUP_WITHOUT_TRACE=0이 필수입니다.')
                model_only=(ORIGINAL_LAYOUT and (env.get('MFM_MODEL_ONLY_CLEANUP')=='1'
                            or (env.get('MFM_WGS_PAYLOAD') in {'sequence_only','pileup_feature'} and env.get('CLEAN_INTERMEDIATES_ON_SUCCESS')=='1')))
                expected_cleanup='1' if model_only or not ORIGINAL_LAYOUT else '0'
                if env.get('CLEAN_INTERMEDIATES_ON_SUCCESS')!=expected_cleanup:
                    raise ValueError('모델 전용 정리는 MFM_MODEL_ONLY_CLEANUP=1과 CLEAN_INTERMEDIATES_ON_SUCCESS=1을 함께 지정해야 합니다.')
                (root/'metadata').mkdir(parents=True,exist_ok=True)
                (root/'metadata/download_inputs.json').write_text(json.dumps(resources,ensure_ascii=False,indent=2)+'\n')
                log=stats/(sample+'.run.log')
                with log.open('a') as output:subprocess.run(['bash',str(ROOT/'run_sample.sh'),sample],env=env,stdout=output,stderr=subprocess.STDOUT,check=True)
                if not processing_complete(root,sample,mode,require_cleanup=model_only):
                    raise ValueError('cleanup 완료 보고서의 schema/sample/mode 또는 정리 후 PASS 계약이 일치하지 않습니다.')
                info.update(status='success',model_bytes=model_bytes(root))
                if c['remove_owned_raw_inputs']:info['removed_raw_inputs']=remove_owned(resources,c['data_root'],state)
            except Exception as e:info.update(status='failed',error=str(e))
            finally:
                info.update(finished=time.time(),elapsed_seconds=time.time()-started)
                atomic_json(record,info)
        return 0 if info['status']=='success' else 1

def cpu_pool_summary(state,c):
    paths=[]
    local=state/'stats'/'cpu_pool_usage.jsonl'
    if local.is_file():paths.append(local)
    directory=state/'stats'/'cpu_pool_usage'
    if directory.is_dir():paths.extend(sorted(directory.glob('*.jsonl')))
    events=[]
    for path in paths:
        for line in path.read_text().splitlines():
            try:events.append(json.loads(line))
            except (ValueError,OSError):continue
    result={'event_count':len(events),'audit_files':[str(p) for p in paths],
      'resource_resolution':c.get('resource_resolution',{})}
    if not events:return result
    grants=[int(r.get('granted_cpus',0)) for r in events]
    maxima=[int(r.get('maximum_cpus',0)) for r in events]
    waits=[float(r.get('queue_wait_sec',0)) for r in events]
    elapsed=[float(r.get('elapsed_sec',0)) for r in events]
    stage_counts=Counter(f"{r.get('branch','')}:{r.get('stage','')}" for r in events)
    result.update(
      average_granted_cpus=round(sum(grants)/len(grants),3),
      max_granted_cpus=max(grants),
      average_queue_wait_sec=round(sum(waits)/len(waits),3),
      queue_wait_over_1s_fraction=round(sum(x>1 for x in waits)/len(waits),4),
      below_stage_cap_fraction=round(sum(g<m for g,m in zip(grants,maxima))/len(events),4),
      leased_cpu_seconds=round(sum(g*t for g,t in zip(grants,elapsed)),3),
      stage_event_counts=dict(sorted(stage_counts.items())))
    return result

def finish(c,rows):
    state=Path(c['state_root']);records=[]
    for r in rows:
        p=state/'samples'/(r['sample_id']+'.json')
        records.append(json.loads(p.read_text()) if p.exists() else {'sample':r['sample_id'],'status':'not_started'})
    total=sum(r.get('model_bytes',0) for r in records if r['status']=='success')
    result={'selection_count':len(rows),'success':sum(r['status']=='success' for r in records),
      'model_input_bytes':total,'target_bytes':None if c.get('selection_target') else c['model_target_bytes'],'records':records,
      'capacity_note':'선택한 클라우드 대상 전체 완료를 기준으로 한다.' if c.get('selection_target') else '시작한 작업은 완료한다. 동시 진행 수에 따라10TB를 초과할 수 있다.',
      'cpu_pool_summary':cpu_pool_summary(state,c)}
    state.mkdir(parents=True,exist_ok=True)
    (state/'batch_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    leftovers=[r for r in rows if next(x['status'] for x in records if x['sample']==r['sample_id'])!='success']
    planner.write_tsv(state/'remaining_samples.tsv',leftovers,fields=list(rows[0]))
    # 샘플 메타데이터의 모델/QC 적격성을 전역 목록으로 만든다.
    catalogs=[]
    for group,project in sorted({(r['data_group'],r['project_id']) for r in rows}):
        directory=Path(c['output_root'])/group/project
        if directory.exists():
            target=state/(group+'_'+project+'_catalog.tsv')
            process=subprocess.run([sys.executable,str(ROOT/('27_make_catalog.py' if ORIGINAL_LAYOUT else '25_make_catalog.py')),str(directory),'--out',str(target)],capture_output=True,text=True)
            catalogs.append({'path':str(target),'exit_code':process.returncode,'log':process.stdout[-2000:]+process.stderr[-2000:]})
    result['catalogs']=catalogs
    (state/'batch_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    # 실패 작업은 참조·중간 파일을 유지하여 재개할 수 있다.
    if c['cleanup_shared_fasta'] and result['success']==len(rows):
        result['shared_fasta_cleanup']=purge_shared(c)
        (state/'batch_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    return result

def purge_shared(c):
    ref=Path(c['data_root'])/'references';removed=[]
    with (ref/'active_samples.lock').open('a') as lease:
        try:fcntl.flock(lease,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return {'status':'deferred','reason':'다른 배치/샘플이 참조를 사용 중입니다.'}
        # 생성 자산과 다운로드 소유권이 확인된 자산만 삭제한다.
        for ownership in ref.glob('*.owned.json'):
            d=json.loads(ownership.read_text());p=Path(d['path'])
            if d.get('schema')=='mfm_derived_reference_v1' and p.parent==ref and p.exists():
                p.unlink();removed.append(str(p));ownership.unlink()
        for record in (Path(c['data_root'])/'.download_records').glob('*.json'):
            d=json.loads(record.read_text());p=Path(d['path'])
            if d.get('owned') and p.parent==ref and any(x in p.name for x in ['.fa','.fna','.fasta']) and p.exists():
                p.unlink();removed.append(str(p));d['status']='removed_after_batch';record.write_text(json.dumps(d)+'\n')
        for name in removed:
            for suffix in ['.fai','.gzi']:
                p=Path(name+suffix)
                if p.exists():p.unlink()
    return {'status':'complete','removed':removed}

def main():
    p=argparse.ArgumentParser(description='전달본: 선택 배치 다운로드→전처리→검증→정리')
    sub=p.add_subparsers(dest='command',required=True)
    for name in ['plan','doctor','dry-run','run','prepare','worker','finish','download','download-worker','download-finish','import-state','resume']:
        s=sub.add_parser(name)
        s.add_argument('--config');s.add_argument('--effective-config');s.add_argument('--selection')
        s.add_argument('--target',choices=CLOUD_TARGETS,help='클라우드 분담 명단 전체 선택; --sample으로 대상 안의 샘플 선택 가능')
        s.add_argument('--lab',choices=['tgil','medai']);s.add_argument('--batch');s.add_argument('--sample',action='append')
        s.add_argument('--executor',choices=['local','slurm'],default='local');s.add_argument('--task-index',type=int)
        s.add_argument('--resource-profile',choices=['auto','throughput','balanced','conservative'])
        s.add_argument('--local-jobs',type=int);s.add_argument('--sample-cpus',type=int)
        s.add_argument('--cpu-allocation-mode',choices=['dynamic','fixed_sample'])
        s.add_argument('--elastic-max-jobs',type=int)
        s.add_argument('--download-jobs',type=int);s.add_argument('--slurm-concurrency',type=int)
        s.add_argument('--slurm-cpus',type=int)
        s.add_argument('--cpu-pool-limit',type=int)
        acceleration.cli_arguments(s)
        if name=='import-state':s.add_argument('--from-dir',required=True)
        if name=='resume':
            s.add_argument('--from-stage',required=True)
            s.add_argument('--overwrite-prepared',action='store_true')
            s.add_argument('--rebuild-haplotag-partial',action='store_true')
        if name=='download':s.add_argument('--dry-run',action='store_true',dest='download_dry_run')
    a=p.parse_args()
    if a.command=='plan':
        if a.target:raise ValueError('클라우드 대상 확인은 dry-run --target을 사용하세요. plan은 기존 배정 재생성 명령입니다.')
        calibration=ROOT/'manifests/분담/calibration.json'
        d=json.loads(calibration.read_text()) if calibration.exists() else {}
        print(json.dumps(planner.build(ROOT,calibration=d),ensure_ascii=False,indent=2));return 0
    if a.target and not a.config and not a.effective_config:
        a.config=str(ROOT/'manifests/클라우드_분담/설정'/f'{a.target}_pileup.json')
    resource_overrides={name:getattr(a,name) for name in (
        'cpu_allocation_mode','resource_profile','local_jobs','sample_cpus','elastic_max_jobs','download_jobs',
        'slurm_concurrency','slurm_cpus','cpu_pool_limit')}
    c=acceleration.apply(load_config(a.effective_config or a.config,resource_overrides),acceleration.overrides(a))
    if a.effective_config:
        stored=c.get('selection_target')
        if a.target and stored!=a.target:raise ValueError('저장된 실행과 --target이 다릅니다.')
        a.target=a.target or stored
    elif a.target:
        c['selection_target']=a.target
        if a.sample and a.command!='resume':c['selection_scope_samples']=sorted(a.sample)
        else:c.pop('selection_scope_samples',None)
    else:
        c.pop('selection_target',None);c.pop('selection_scope_samples',None)
    if a.command=='doctor':
        from .local_queue import plan as queue_plan
        print(json.dumps({'config':c,'local_queue':queue_plan(c),'slurm_commands':{x:shutil.which(x) for x in ['sbatch','sinfo','scontrol']},
          'free_space':{k:shutil.disk_usage(Path(c[k]).parent if Path(c[k]).parent.exists() else ROOT).free for k in ['data_root','output_root','state_root']},
          'next':'실행 서버에서 sinfo -o "%P %c %m %l", scontrol show partition을 확인하고 경로·자원을 설정하세요.'},ensure_ascii=False,indent=2));return 0
    if not a.target and not a.selection and not a.lab and c.get('default_selection'):
        a.selection=str((ROOT/c['default_selection']).resolve())
    rows=selection(a)
    if not a.effective_config:c=batch_config(c,rows)
    if a.command=='import-state':
        from .state_import import import_states
        print(json.dumps(import_states(c,rows,a.from_dir,sys.modules[__name__]),ensure_ascii=False,indent=2));return 0
    if a.command=='resume':
        from .queue_recovery import enqueue
        print(json.dumps(enqueue(c,rows,a.sample,a.from_stage,a.overwrite_prepared,
                                 a.rebuild_haplotag_partial,atomic_json,ORIGINAL_LAYOUT),ensure_ascii=False,indent=2));return 0
    if a.command.startswith('download'):
        from . import download_only
        if a.command=='download':return download_only.run(c,rows,ROOT,a.executor,a.download_dry_run)
        if not a.effective_config or not a.selection:raise ValueError('다운로드 내부 명령은 effective-config와 selection이 필요함')
        if a.command=='download-worker':
            if a.task_index is None or not 0<=a.task_index<len(rows):raise ValueError('download array task index 범위 오류')
            return download_only.worker(c,rows[a.task_index],ROOT)
        result=download_only.finish(c,rows);print(json.dumps(result,ensure_ascii=False,indent=2));return 0 if result['success']==len(rows) else 1
    state=Path(c['state_root']);config_path=state/'effective_config.json';batch_path=state/'selected_samples.tsv'
    if a.command=='dry-run':
        d={'selected':len(rows),'sample_ids':[r['sample_id'] for r in rows], 'estimated_model_bytes':sum(int(r['estimated_model_bytes']) for r in rows),
          'config':c,'downloads_first_sample':sample_downloads(c,rows[0],sample_resources(ROOT,rows[0])),'side_effects':'다운로드·설치·제출·전처리를 실행하지 않는다.'}
        if a.target:d.update(selection_target=a.target,completion_policy='all_selected')
        if a.executor=='local':
            from .local_queue import plan as queue_plan
            d['local_queue']=queue_plan(c)
        if a.executor=='slurm':d['slurm']=slurm.submit(c,ROOT,config_path,batch_path,len(rows),True)
        print(json.dumps(d,ensure_ascii=False,indent=2));return 0
    if a.command=='run' and a.executor=='local':
        from .local_queue import run
        return run(c,rows,sys.modules[__name__])
    state.mkdir(parents=True,exist_ok=True)
    if not a.effective_config:
        config_path.write_text(json.dumps(c,ensure_ascii=False,indent=2)+'\n');planner.write_tsv(batch_path,rows)
    if a.command=='run' and a.executor=='slurm':
        print(json.dumps(slurm.submit(c,ROOT,config_path,batch_path,len(rows)),ensure_ascii=False,indent=2));return 0
    if a.command=='prepare':prepare(c,rows);return 0
    if a.command=='worker':
        if a.task_index is None or not 0<=a.task_index<len(rows):raise ValueError('array task index 범위 오류')
        return worker(c,rows[a.task_index],a.executor)
    if a.command=='finish':
        result=finish(c,rows);print(json.dumps(result,ensure_ascii=False,indent=2));return 0 if result['success']==len(rows) else 1
    raise ValueError('지원하지 않는 실행 명령: '+a.command)

if __name__=='__main__':
    try:sys.exit(main())
    except Exception as e:print('[오류] '+str(e),file=sys.stderr);sys.exit(1)
