"""준비 작업→샘플 array→최종 배치 보고서의 Slurm 의존성을 만든다."""
from __future__ import annotations
import json,re,shlex,subprocess
from pathlib import Path

def flags(config,cpus=None,memory=None,time=None):
    s=config['slurm'];out=['--parsable','--cpus-per-task='+str(cpus or s['cpus']),
      '--mem='+str(memory or s['memory']),'--time='+str(time or s['time'])]
    for k in ['partition','account','qos']:
        if s.get(k):out+=['--'+k+'='+s[k]]
    if s.get('constraint'):out+=['--constraint='+s['constraint']]
    return out

def commands(config,root,config_path,batch_path,count):
    root=Path(root);state=Path(config['state_root']);log=state/'slurm_logs'
    common=[str(root/('run_handoff.sh' if (root/'09_personal_fasta.sh').is_file() else 'run_pipeline.sh')),'--effective-config',str(config_path),'--selection',str(batch_path)]
    prep=['sbatch',*flags(config,1,'16G','08:00:00'),'--job-name=mfm_prepare','--output='+str(log/'prepare_%j.log'),
       '--wrap='+shlex.join(['bash',str(root/('run_handoff.sh' if (root/'09_personal_fasta.sh').is_file() else 'run_pipeline.sh')),'prepare','--effective-config',str(config_path),'--selection',str(batch_path)])]
    array=['sbatch',*flags(config),'--job-name=mfm_samples','--array=0-'+str(count-1)+'%'+str(config['slurm']['concurrency']),
      '--dependency=afterok:PREP_JOB_ID','--kill-on-invalid-dep=yes','--output='+str(log/'sample_%A_%a.log'),
      '--wrap='+shlex.join(['bash',str(root/'slurm/sample_task.sh'),str(config_path),str(batch_path)])]
    if config.get('slurm',{}).get('gpus_per_task'):
        count_gpu=config['slurm']['gpus_per_task']
        if isinstance(count_gpu,bool) or not isinstance(count_gpu,int) or count_gpu<1:raise ValueError('slurm.gpus_per_task는 양의 정수여야 함')
        array.insert(1,'--gpus-per-task='+str(count_gpu))
    finish=['sbatch',*flags(config,1,'4G','04:00:00'),'--job-name=mfm_summary','--dependency=afterany:ARRAY_JOB_ID',
       '--output='+str(log/'summary_%j.log'),'--wrap='+shlex.join(['bash',str(root/('run_handoff.sh' if (root/'09_personal_fasta.sh').is_file() else 'run_pipeline.sh')),'finish',
         '--effective-config',str(config_path),'--selection',str(batch_path)])]
    return [prep,array,finish]

def submit(config,root,config_path,batch_path,count,dry_run=False):
    plans=commands(config,root,config_path,batch_path,count)
    if dry_run:return {'commands':plans,'shell_commands':[shlex.join(c) for c in plans],'submitted':False}
    (Path(config['state_root'])/'slurm_logs').mkdir(parents=True,exist_ok=True)
    ids=[]
    for i,cmd in enumerate(plans):
        if i:cmd=[a.replace('PREP_JOB_ID',ids[0]).replace('ARRAY_JOB_ID',ids[1] if len(ids)>1 else '') for a in cmd]
        job=subprocess.check_output(cmd,text=True).strip().split(';')[0]
        if not re.fullmatch(r'\d+',job):raise ValueError('sbatch job ID를 읽지 못했습니다: '+job)
        ids.append(job)
        # 부분 제출 실패 시 이미 제출된 job ID도 보존한다.
        (Path(config['state_root'])/'slurm_jobs.json').write_text(json.dumps({'job_ids':ids,'commands':plans},indent=2)+'\n')
    return {'job_ids':ids,'submitted':True}
