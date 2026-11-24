"""선택 입력만 내려받음. 전처리·환경 설치·원본 정리·모델 용량 제한을 호출하지 않음."""
from __future__ import annotations
import concurrent.futures,fcntl,hashlib,json,os,re,shlex,subprocess,tempfile,time,uuid
from contextlib import contextmanager
from pathlib import Path
from . import bundled_assets,planner,slurm
from .downloads import fetch,safe_path,sample_resources
from .resume import sample_downloads
from integrate_wgs_rna.memory_policy import usable_memory,limit_worker_address_space
from integrate_wgs_rna.pool_leases import prune_leases
from .download_budget import memory_plan

INDEX_KINDS={'wgs_index','rna_index','sample_tbi','shared_tbi','illumina_panel_index'}
DONE={'success','success_with_missing_indexes'}
DOWNLOAD_RAM_BYTES=1_000_000_000

def atomic_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile('w',dir=path.parent,delete=False,encoding='utf-8') as f:
        json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n');temporary=Path(f.name)
    os.replace(temporary,path)

def configure(config):
    c=dict(config);c['operation']='download_only'
    c['download_jobs']=c.get('download_jobs',c['local_jobs'])
    if not isinstance(c['download_jobs'],int) or isinstance(c['download_jobs'],bool) or c['download_jobs']<1:
        raise ValueError('download_jobs는 양의 정수여야 함')
    c['state_root']=str(Path(c['state_root'])/'download_only'/c['selection_fingerprint'][:16])
    return c

def plan(c,rows,root):
    base,bundles=bundled_assets.inventory(c,root);unique={}
    for row in rows:
        for r in sample_downloads(c,row,sample_resources(root,row)):
            key=r['path']
            if key in unique and unique[key]['url']!=r['url']:
                raise ValueError('같은 다운로드 경로에 다른 URL이 있음: '+key)
            item=dict(r);safe_path(c['data_root'],key)
            if r.get('shared') and key in bundles:
                path=safe_path(base,key)
                if path.is_file() and not path.is_symlink():
                    item['bundled_path']=str(path);item['bundled_size_bytes']=bundles[key]['size_bytes']
                else:
                    item['asset_restore_required']=True
                    item['size_bytes']=bundles[key]['size_bytes']
            unique[key]=item
    return {'operation':'download_only','selected':len(rows),'sample_ids':[r['sample_id'] for r in rows],
      'data_root':c['data_root'],'download_state_root':c['state_root'],'resource_count':len(unique),
      'bundled_resources':sum('bundled_path' in r for r in unique.values()),
      'asset_restore_resources':sum(r.get('asset_restore_required',False) for r in unique.values()),
      'estimated_download_bytes':sum(int(r.get('size_bytes') or 0) for r in unique.values() if 'bundled_path' not in r),
      'unknown_size_files':sum(not r.get('size_bytes') for r in unique.values() if 'bundled_path' not in r),
      'model_target_applied':False,'raw_inputs_removed':False,
      'first_sample_resources':sample_downloads(c,rows[0],sample_resources(root,rows[0]))}

def initialize(c,rows):
    data=Path(c['data_root']);state=Path(c['state_root']);data.mkdir(parents=True,exist_ok=True);state.mkdir(parents=True,exist_ok=True)
    marker=data/'.mfm_handoff_store.json'
    with (data/'.mfm_handoff_store.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if marker.exists():
            if json.loads(marker.read_text()).get('schema')!='mfm_handoff_store_v1':raise ValueError('입력 저장소 소유권 마커가 다름')
        else:atomic_json(marker,{'schema':'mfm_handoff_store_v1','package':c['package_root'],'created':time.time()})
    fingerprint=state/'selection_fingerprint.txt'
    if fingerprint.exists() and fingerprint.read_text().strip()!=c['selection_fingerprint']:raise ValueError('다운로드 선택 fingerprint 불일치')
    fingerprint.write_text(c['selection_fingerprint']+'\n')
    atomic_json(state/'effective_config.json',c)
    planner.write_tsv(state/'selected_samples.tsv',rows,fields=list(rows[0]))

@contextmanager
def download_memory_lease(c,sample):
    estimate=memory_plan(c);memory=estimate['reservation_bytes']
    pool=c.get('pipeline_env',{}).get('MFM_CPU_POOL_DIR')
    if not pool:
        yield
        return
    root=Path(pool);directory=root/'leases';directory.mkdir(parents=True,exist_ok=True)
    path=directory/('download_'+sample+'_'+uuid.uuid4().hex+'.json')
    while True:
        with (root/'pool.lock').open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            prune_leases(root)
            reserved=sum(int(json.loads(p.read_text()).get('memory_bytes',0))
                         for p in directory.glob('*.json'))
            if reserved+memory<=usable_memory(root):
                atomic_json(path,{'sample':sample,'branch':'DOWNLOAD','stage':'download',
                                  'cpus':0,'memory_bytes':memory,'download_memory_plan':estimate,
                                  'manager_pid':os.getpid(),'child_pid':0,'started':time.time()})
                break
        time.sleep(.2)
    try:
        # CLI 단독과 범위 병렬 버퍼, 로컬 스레드 풀의 동시 worker 수를 함께 계상한다.
        limit_worker_address_space(memory * c.get('_download_worker_concurrency',1))
        yield
    finally:
        with (root/'pool.lock').open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            path.unlink(missing_ok=True)

def worker(c,row,root):
    from .io_env import apply as apply_io_environment
    apply_io_environment(c)
    if c.get('operation')!='download_only':raise ValueError('다운로드 전용 effective_config가 필요함')
    value=c.get('http_range_workers')
    if value is not None and (isinstance(value,bool) or not isinstance(value,int) or not 1 <= value <= 16):
        raise ValueError('http_range_workers는 1~16 정수여야 함')
    # RAM 사전 산정과 실제 전송에 같은 범위 수를 사용함.
    # pipeline_env/운영 정책도 누락 없이 새 다운로드 worker에 반영함.
    os.environ['MFM_HTTP_RANGE_WORKERS']=str(memory_plan(c)['range_workers'])
    state=Path(c['state_root']);sample=row['sample_id'];directory=state/'samples';directory.mkdir(parents=True,exist_ok=True)
    started=time.time();info={'schema':'mfm_download_only_sample_v1','sample':sample,'status':'running','started':started,'files':[],'missing_indexes':[]}
    with (directory/(sample+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        record=directory/(sample+'.json');atomic_json(record,info)
        try:
            with download_memory_lease(c,sample):
                resources=sample_downloads(c,row,sample_resources(root,row))
                for r in sorted(resources,key=lambda r:r['kind'] in INDEX_KINDS):
                    item=dict(r)
                    bundled=bundled_assets.get(c,root,r['path']) if r.get('shared') else None
                    if bundled is not None:
                        item.update(local_path=str(bundled),status='bundled',actual_bytes=bundled.stat().st_size)
                    else:
                        try:p=fetch(r,c['data_root'],state)
                        except subprocess.CalledProcessError as exc:
                            if r['kind'] not in INDEX_KINDS or exc.returncode!=22:raise
                            item.update(local_path=str(safe_path(c['data_root'],r['path'])),status='index_unavailable',error=str(exc))
                            info['missing_indexes'].append(item)
                            info['files'].append(item);atomic_json(record,info);continue
                        item.update(local_path=str(p),status='verified',actual_bytes=p.stat().st_size)
                    info['files'].append(item);atomic_json(record,info)
                info['status']='success_with_missing_indexes' if info['missing_indexes'] else 'success'
        except Exception as exc:info.update(status='failed',error=str(exc))
        finally:
            info.update(finished=time.time(),elapsed_seconds=time.time()-started)
            atomic_json(record,info);atomic_json(state/'attempts'/f'{sample}.{time.time_ns()}.json',info)
    return 0 if info['status'] in DONE else 1

def finish(c,rows):
    if c.get('operation')!='download_only':raise ValueError('다운로드 전용 effective_config가 필요함')
    state=Path(c['state_root']);records=[]
    for row in rows:
        p=state/'samples'/(row['sample_id']+'.json')
        records.append(json.loads(p.read_text()) if p.exists() else {'sample':row['sample_id'],'status':'not_started','files':[],'missing_indexes':[]})
    unique={r['local_path']:r for info in records for r in info['files'] if r['status'] in {'verified','bundled'}}
    warnings=[dict(sample=info['sample'],**r) for info in records for r in info['missing_indexes']]
    result={'schema':'mfm_download_only_batch_v1','operation':'download_only','selection_count':len(rows),
      'success':sum(r['status'] in DONE for r in records),'failed':sum(r['status']=='failed' for r in records),
      'not_started':sum(r['status']=='not_started' for r in records),'data_root':c['data_root'],'download_state_root':str(state),
      'verified_local_bytes':sum(r['actual_bytes'] for r in unique.values() if r['status']=='verified'),
      'bundled_bytes_reused':sum(r['actual_bytes'] for r in unique.values() if r['status']=='bundled'),
      'missing_index_count':len(warnings),'model_target_applied':False,'raw_inputs_removed':False,'records':records}
    state.mkdir(parents=True,exist_ok=True)
    atomic_json(state/'download_summary.json',result)
    failed={r['sample'] for r in records if r['status'] not in DONE}
    planner.write_tsv(state/'remaining_samples.tsv',[r for r in rows if r['sample_id'] in failed],fields=list(rows[0]))
    planner.write_tsv(state/'missing_indexes.tsv',warnings,fields=['sample','kind','path','url','local_path','status','error'])
    return result

def commands(c,root,count):
    root=Path(root);state=Path(c['state_root']);log=state/'slurm_logs'
    config=state/'effective_config.json';selection=state/'selected_samples.tsv'
    array=['sbatch',*slurm.flags(c,1,'2G',c['slurm'].get('download_time',c['slurm']['time'])),
      '--job-name=mfm_download','--array=0-'+str(count-1)+'%'+str(c['slurm']['concurrency']),
      '--output='+str(log/'download_%A_%a.log'),
      '--wrap='+shlex.join(['bash',str(root/'slurm/download_task.sh'),str(config),str(selection)])]
    summary=['sbatch',*slurm.flags(c,1,'2G','01:00:00'),'--job-name=mfm_download_summary',
      '--dependency=afterany:ARRAY_JOB_ID','--output='+str(log/'summary_%j.log'),
      '--wrap='+shlex.join(['bash',str(root/('run_handoff.sh' if (root/'09_personal_fasta.sh').is_file() else 'run_pipeline.sh')),'download-finish','--effective-config',str(config),'--selection',str(selection)])]
    return [array,summary]

def submit(c,root,count,dry_run=False):
    plans=commands(c,root,count)
    if dry_run:return {'commands':plans,'shell_commands':[shlex.join(x) for x in plans],'submitted':False}
    (Path(c['state_root'])/'slurm_logs').mkdir(parents=True,exist_ok=True);ids=[]
    for command in plans:
        cmd=[x.replace('ARRAY_JOB_ID',ids[0]) for x in command] if ids else command
        job=subprocess.check_output(cmd,text=True).strip().split(';')[0]
        if not re.fullmatch(r'\d+',job):raise ValueError('sbatch job ID를 읽지 못함: '+job)
        ids.append(job);atomic_json(Path(c['state_root'])/'slurm_jobs.json',{'job_ids':ids,'commands':plans})
    return {'job_ids':ids,'submitted':True}

def run(config,rows,root,executor='local',dry_run=False):
    c=configure(config);preview=plan(c,rows,root)
    if dry_run:
        preview['side_effects']='다운로드·설치·제출·전처리·cleanup을 실행하지 않음'
        if executor=='slurm':preview['slurm']=submit(c,root,len(rows),True)
        print(json.dumps(preview,ensure_ascii=False,indent=2));return 0
    initialize(c,rows)
    if executor=='slurm':
        preview['slurm']=submit(c,root,len(rows));print(json.dumps(preview,ensure_ascii=False,indent=2));return 0
    shared_worker_config=dict(c,_download_worker_concurrency=c['download_jobs'])
    with concurrent.futures.ThreadPoolExecutor(max_workers=c['download_jobs']) as pool:
        results=list(pool.map(lambda row:worker(shared_worker_config,row,root),rows))
    result=finish(c,rows);print(json.dumps(result,ensure_ascii=False,indent=2))
    return 0 if result['success']==len(rows) and all(x==0 for x in results) else 1
