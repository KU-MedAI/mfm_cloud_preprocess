"""개인 인덱스와 paired BAM 준비를 한 예약의 고정 하위 예산으로 중첩함."""
from pathlib import Path
from dataclasses import asdict,replace
from concurrent.futures import ThreadPoolExecutor
import argparse,fcntl,json,os,resource,signal,subprocess,sys,tempfile,time

from ..memory_policy import GIB,process_rss
from ..pool_leases import read_pool_records,cgroup_cpu_limit
from ..gpu_bwa_index import cache_root,ready_hint,ensure_index,fingerprint
from ..gpu_wgs import scratch_directory
from .config import PersonalRealignConfig,PreparedPairedInput
from .phases import _atomic

PATH_FIELDS=['output_root','sample_list','alignment','hap1_fasta_dir','hap2_fasta_dir','reference_fasta','bwa_mem2_index_root','summary_tsv','metadata_json']

def child_environment(cpus,memory):
    env=dict(os.environ)
    for name in ['MFM_CPU_POOL_LEASE_FILE','MFM_GPU_LOCAL_LEASE_FILE','MFM_GPU_HAP_ID','MFM_MEMORY_STAGE_BASE_BYTES']:
        env.pop(name,None)
    env.update(MFM_CPU_BUDGET=str(cpus),MFM_SAMPLE_CPU_BUDGET=str(cpus),MFM_CPU_POOL_GRANTED=str(cpus),
               MFM_RAM_LIMIT_BYTES=str(memory*10//9),MFM_GPU_COMPARE_UNRESERVED='0',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    return env

def _disk_fits(components,others,reserve):
    paths=[Path(p) for s in components.values() for p in s['write_paths']]
    if len({p.parent.stat().st_dev for p in paths})!=1:return False
    device=paths[0].parent.stat().st_dev
    claims=list(components.values())
    for row in others:
        claims.extend(row.get('hap_resources',{}).values());claims.extend(row.get('prepare_components',{}).values())
    needed=0
    for row in claims:
        files=[Path(p) for p in row.get('write_paths',[])]
        if not files or not files[0].parent.exists() or files[0].parent.stat().st_dev!=device:continue
        needed+=max(0,int(row.get('disk_peak_bytes',0))-sum(p.stat().st_size for p in files if p.is_file()))
    st=os.statvfs(paths[0].parent)
    return needed+reserve<=st.f_bavail*st.f_frsize

def commit_components(phases,components,*,admit=False):
    """진행 중 부분은 고정 grant를 유지하고 종료된 부분의 몫만 반환함."""
    cpus=1+sum(s['cpus'] for s in components.values())
    memory=4*GIB+sum(s['memory_bytes'] for s in components.values())
    if phases.path:
        root=phases.path.parent.parent
        with (root/'pool.lock').open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX);row=json.loads(phases.path.read_text())
            if row.get('realign_phase')!='gpu_bam_index_prepare':raise RuntimeError('중첩 준비 lease 단계 변경')
            if cpus>row['cpus'] or memory>row['memory_bytes']:raise ValueError('고정 준비 예약의 비인가 증설')
            if admit:
                others=[r for p,r in read_pool_records(root/'leases') if p!=phases.path]
                policy=json.loads((root/'policy.json').read_text()) if (root/'policy.json').exists() else {}
                floor=int(policy.get('adaptive_prepare_storage',{}).get('reserve_bytes',150_000_000_000))
                if not _disk_fits(components,others,floor):return False
            row.update(cpus=cpus,memory_bytes=memory,prepare_components=components)
            _atomic(phases.path,row)
    elif admit and not _disk_fits(components,[],150_000_000_000):return False
    event=dict(phase='gpu_prepare_components',time=time.time(),reserved_cpus=cpus,reserved_memory_bytes=memory,
               active_components=list(components),reservation_applied=not phases.unreserved)
    phases.events.append(event)
    return True

def _group_members(group):
    found=[]
    for path in Path('/proc').iterdir():
        if not path.name.isdigit():continue
        try:
            row=(path/'stat').read_text().rsplit(')',1)[1].split()
            if int(row[2])==group and row[0]!='Z':found.append(int(path.name))
        except (OSError,ValueError,IndexError):continue
    return found

def stop_group(group):
    try:os.killpg(group,signal.SIGTERM)
    except ProcessLookupError:return
    deadline=time.monotonic()+5
    while _group_members(group):
        if time.monotonic()>deadline:
            try:os.killpg(group,signal.SIGKILL)
            except ProcessLookupError:pass
        time.sleep(.05)

def overlap_seconds(primary,parts):
    intervals=sorted((max(primary['_work_started_at'],p['_work_started_at']),min(primary['_work_finished_at'],p['_work_finished_at'])) for p in parts)
    total=0.;end=float('-inf')
    for left,right in intervals:
        if right<=left:continue
        total+=max(0,right-max(left,end));end=max(end,right)
    return total


def prepare(config,sample,alignment,contigs,hap_dirs,realigned_root,work_dir,log_dir,required,phases,metrics):
    """지원/자원이 부족하면 None을 반환해 기존 순차 경로를 사용함."""
    total=min(phases.total,cgroup_cpu_limit() or phases.total)
    if total<7:return None
    phases.reserve('target_fasta_early',max_units=min(2,len(required)),memory_per_unit=2*GIB)
    from ..personal_realign import build_target_fasta
    targets={h:realigned_root/f'HAP{h}'/'target.fa' for h in required}
    def target(h):build_target_fasta(config,hap_dirs[h],targets[h],contigs,' '.join(contigs))
    with ThreadPoolExecutor(max_workers=min(2,len(required))) as pool:list(pool.map(target,required))
    scratch=scratch_directory(work_dir);storage=cache_root(scratch);storage.mkdir(parents=True,exist_ok=True)
    pending=[h for h in required if not ready_hint(targets[h],scratch)]
    if not pending:
        metrics.update(mode='indexes_already_cached',overlap_observed=False)
        phases.reserve('prepare_cache_hit',max_units=1);phases.finish()
        return None
    if storage.stat().st_dev!=scratch.stat().st_dev:
        phases.reserve('prepare_device_fallback',max_units=1);phases.finish();metrics.update(mode='device_fallback',overlap_observed=False)
        return None
    sizes={h:targets[h].stat().st_size for h in pending}
    index_ram={h:(sizes[h]*16+4*GIB)*11//10 for h in pending}
    base_cpus=len(pending)+1;base_memory=4*GIB+sum(index_ram.values())
    if total<base_cpus+4 or base_memory+8*GIB>phases.ceiling():
        phases.reserve('prepare_budget_fallback',max_units=1);phases.finish();metrics.update(mode='budget_fallback',overlap_observed=False)
        return None
    if phases.path:
        others=sum(int(r.get('memory_bytes',0)) for path,r in read_pool_records(phases.path.parent) if path!=phases.path)
        if base_memory+8*GIB>phases.ceiling()-others:
            phases.reserve('prepare_memory_fallback',max_units=1);phases.finish();metrics.update(mode='memory_fallback',overlap_observed=False)
            return None
    grant=phases.reserve('gpu_bam_index_prepare',max_units=total-base_cpus,min_units=4,
        memory_per_unit=2*GIB,base_memory=base_memory,base_cpus=base_cpus)
    input_memory=grant*2*GIB
    routes=max(1,(grant-1)//2) if os.environ.get('MFM_WGS_ROUTE_REQUESTED')=='auto' else min(config.route_jobs,max(1,(grant-1)//2))
    input_config=replace(config,route_jobs=routes,
        sort_threads=min(config.sort_threads,max(1,grant-3)),fastq_threads=min(config.fastq_threads,max(0,grant-3)))
    source_stamp=fingerprint(alignment)
    components={'input':dict(cpus=grant,memory_bytes=input_memory,write_paths=[str(scratch/'overlap.input.reserved')],
        disk_peak_bytes=max(350_000_000_000,Path(alignment).stat().st_size*8+32*GIB))}
    components.update({f'index_{h}':dict(cpus=1,memory_bytes=index_ram[h],write_paths=[str(storage/f'overlap.index_{h}.reserved')],disk_peak_bytes=sizes[h]*16+GIB) for h in pending})
    if not commit_components(phases,components,admit=True):
        commit_components(phases,{});phases.finish();metrics.update(mode='disk_fallback',overlap_observed=False)
        return None
    records={};active={};started=time.time();sampled=0;peak=0
    metrics.update(mode='bam_and_bwa_indexes',input_cpus=grant,input_ram_bytes=input_memory,
        initial_reserved_cpus=base_cpus+grant,initial_reserved_memory_bytes=base_memory+input_memory,
        index_cpus=len(pending),runs=records,started_at=started,status='running')
    metrics_path=log_dir/f'{sample}.prepare_overlap.json';_atomic(metrics_path,metrics)
    try:
        with tempfile.TemporaryDirectory(prefix='prepare_overlap_',dir=work_dir) as raw:
            directory=Path(raw)
            jobs={'input':dict(kind='input',config=asdict(input_config),sample=sample,alignment=str(alignment),contigs=contigs,work=str(work_dir))}
            jobs.update({f'index_{h}':dict(kind='index',target=str(targets[h]),work=str(scratch)) for h in pending})
            for name,job in jobs.items():
                p=directory/(name+'.job.json');result=directory/(name+'.result.json')
                p.write_text(json.dumps(job,default=str)+'\n')
                log=(log_dir/f'{sample}.prepare_overlap.{name}.log').open('w')
                spec=components[name]
                cmd=[sys.executable,'-B','-m','integrate_wgs_rna.realign.gpu_prepare_overlap','--job',str(p),'--result',str(result),'--memory',str(spec['memory_bytes'])]
                child=subprocess.Popen(cmd,env=child_environment(spec['cpus'],spec['memory_bytes']),stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                active[child.pid]=(name,child,result,log,time.time())
            failure=None;values={}
            while active:
                now=time.time()
                if now-sampled>=1:
                    peak=max(peak,sum(process_rss({name:pid for pid,(name,*_) in active.items()}).values()));sampled=now
                for pid,(name,child,result,log,began) in list(active.items()):
                    done,status,usage=os.wait4(pid,os.WNOHANG)
                    if not done:continue
                    child.returncode=os.waitstatus_to_exitcode(status);log.close()
                    stop_group(pid)
                    del active[pid];components.pop(name)
                    records[name]=dict(pid=pid,started_at=began,finished_at=time.time(),exit_code=child.returncode,
                        cpu_seconds=usage.ru_utime+usage.ru_stime,max_single_process_rss_bytes=usage.ru_maxrss*1024)
                    if child.returncode or not result.exists():
                        with Path(log.name).open('rb') as detail:
                            detail.seek(max(0,Path(log.name).stat().st_size-4096));tail=detail.read().decode(errors='replace')
                        if failure is None:failure=RuntimeError('동시 준비 worker 실패: '+name+'; 로그 '+str(log.name)+'; '+tail)
                    else:values[name]=json.loads(result.read_text())
                    commit_components(phases,components);_atomic(metrics_path,metrics)
                if active:time.sleep(.1)
            if failure is not None:raise failure
            if fingerprint(alignment)!=source_stamp:raise ValueError('중첩 준비 중 원본 변경')
            item=values['input'];from .prepared_cache import complete_metrics
            if not item['name_collated'] or not complete_metrics(item['input_counts'],item['stats']):raise ValueError('준비 BAM 성공 계약 불일치')
            overlap=overlap_seconds(values['input'],[r for n,r in values.items() if n!='input'])
            metrics.update(status='pass',finished_at=time.time(),elapsed_seconds=time.time()-started,
                overlap_seconds=overlap,overlap_observed=overlap>0,peak_children_tree_rss_sampled_bytes=peak,index_results={n:v for n,v in values.items() if n!='input'})
            return PreparedPairedInput(Path(item['path']),item['route_label'],item['stats'],item['input_counts'],item['name_collated'])
    except BaseException as error:
        metrics.update(status='failed',error=repr(error));raise
    finally:
        # 이 호출이 만든 그룹만 종료/회수하며 다른 샘플이나 기존 비교 계산은 유지함.
        for pid,(name,child,result,log,began) in active.items():
            stop_group(pid);child.wait();log.close()
        commit_components(phases,{})
        phases.finish();_atomic(metrics_path,metrics)

def worker(job_path,result_path,memory):
    if memory<1:raise ValueError('준비 worker의 RAM 상한은 양수여야 함')
    def interrupted(signum,frame):raise KeyboardInterrupt('준비 worker 중지')
    for signum in [signal.SIGTERM,signal.SIGINT]:signal.signal(signum,interrupted)
    _,hard=resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS,(min(memory,hard) if hard!=resource.RLIM_INFINITY else memory,hard))
    job=json.loads(Path(job_path).read_text())
    if job['kind']=='input':
        from ..personal_realign import prepare_illumina_paired_input
        fields=job['config']
        for field in PATH_FIELDS:
            if fields.get(field) is not None:fields[field]=Path(fields[field])
        began=time.time()
        prepared=prepare_illumina_paired_input(PersonalRealignConfig(**fields),job['sample'],1,Path(job['alignment']),job['contigs'],Path(job['work']),None,'all_to_both')
        result=asdict(prepared)
    elif job['kind']=='index':
        began=time.time();result={};reference=ensure_index(Path(job['target']),Path(job['work']),metrics=result);result['reference']=str(reference)
    else:raise ValueError('지원하지 않는 준비 worker 종류')
    result['_work_started_at']=began;result['_work_finished_at']=time.time()
    Path(result_path).write_text(json.dumps(result,default=str)+'\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--job',type=Path,required=True);parser.add_argument('--result',type=Path,required=True);parser.add_argument('--memory',type=int,required=True)
    args=parser.parse_args();worker(args.job,args.result,args.memory)
