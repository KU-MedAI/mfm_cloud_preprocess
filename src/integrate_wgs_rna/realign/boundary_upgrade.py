"""승인된 순차작업의 HAP1 완료를 보존하며 HAP2 직전만 재배정함."""
from pathlib import Path
import fcntl
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import tempfile
import time

from .phases import RealignPhases, _atomic, reference_bases
from .config import PersonalRealignConfig


def process_start(pid):
    raw = (Path('/proc')/str(pid)/'stat').read_text()
    return int(raw[raw.rfind(')')+2:].split()[19])


def prepare_boundary(command, report):
    lease_value = os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if not lease_value or len(command) != 2 or Path(command[0]).name != 'bash':
        return command, None
    lease_path = Path(lease_value);pool = lease_path.parent.parent
    requests_path = pool/'hap_boundary_requests.json'
    if not requests_path.exists():
        return command, None
    requests = json.loads(requests_path.read_text());row=json.loads(lease_path.read_text());sample=row.get('sample')
    request = requests.get('samples',{}).get(sample)
    script = Path(command[1]);report=Path(report)
    if (not request or request.get('status') != 'armed'
            or script.name != f'{sample}.HAP2.realign.sh'
            or report.name != f'{sample}.HAP2.resource.json'):
        return command,None
    if (row.get('realign_phase') != 'align_restore_sort_qc'
            or row.get('effective_alignment',{}).get('hap_jobs') != 1
            or row.get('manager_pid') != request['manager_pid']
            or process_start(row['manager_pid']) != request['manager_start']):
        return command,None
    if int(row.get('total_cpus',0)) < 5 or os.environ.get('MFM_WGS_ALIGN_REQUESTED','auto') != 'auto':
        return command,None
    root=Path(request['sample_root']).resolve();log=root/'logs/personal_realign'
    if script.resolve().parent != log or report.resolve().parent != log:
        raise ValueError('HAP 경계 요청과 실제 script/report 소유 경로가 다름')
    prior_dir=root/'realigned'/request.get('hap_read_mode','all_to_both')/'HAP1'
    prior=json.loads((prior_dir/f'{sample}.realign.json').read_text())
    if (prior.get('sample')!=sample or prior.get('haplotype')!='HAP1'
            or prior.get('hap_read_mode')!=request.get('hap_read_mode','all_to_both')
            or prior.get('qc_contract',{}).get('status')!='pass'):
        raise ValueError('HAP1 완료 QC 계약이 맞지 않음')
    fmt=prior.get('format');output=prior_dir/f'{sample}.personal.HAP1.{fmt}'
    indexes=[Path(str(output)+('.crai' if fmt=='cram' else '.bai')),output.with_suffix('.crai' if fmt=='cram' else '.bai')]
    if fmt not in {'bam','cram'} or not output.is_file() or not any(p.is_file() for p in indexes):
        raise ValueError('완료 HAP1 출력/index가 없음')
    subprocess.run([os.environ.get('SAMTOOLS_BIN','samtools'),'quickcheck',str(output)],check=True)
    original=script.read_text();align=re.search(r'\| "\$BWA_MEM2" mem [^\n]*? -t ([0-9]+)',original)
    minimap=re.search(r'\| (?:\(\s*)?"\$MINIMAP2" -ax (map-hifi|map-ont) [^\n]*? -t ([0-9]+)',original)
    platform='illumina';thread_group=1
    if align is None and minimap is not None:
        align=minimap;thread_group=2
        platform='hifi' if minimap.group(1)=='map-hifi' else 'ont'
    sort=re.search(r'"\$SAMTOOLS" sort -@ ([0-9]+) -m ([^ ]+)',original)
    if not align or not sort:
        raise ValueError('승인된 BWA/minimap2 생성 script 형식과 다름')
    sort_mem=sort.group(2).strip("'\"")
    lengths=reference_bases({2:root/'individual_fasta/HAP2'},[f'chr{n}' for n in range(1,23)])
    cfg=PersonalRealignConfig(cohort='legacy_boundary',output_root=root.parent,sample_list=None,
        samples=[sample],chroms=['1-22'],sequence_variant_policy='snv_indel',dataset_role='pre_training',
        threads=int(align.group(thread_group)),sort_threads=int(sort.group(1)),sort_mem=sort_mem,
        hap_jobs=1,fastq_threads=0,platform=platform)
    with (pool/'pool.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        current=json.loads(requests_path.read_text());record=current['samples'][sample]
        if record.get('status')!='armed' or record['manager_start']!=request['manager_start']:
            return command,None
        record.update(status='claimed',claimed_by=os.getpid(),claimed_at=time.time());_atomic(requests_path,current)
    # 구실행 진입부에서 누락을 auto로 인정한 경우 계획부에도 같은 요청을 전달함.
    # 일반 alignment_plan의 명시적 config 스레드 의미는 변경하지 않음.
    names = ('MFM_WGS_ALIGN_REQUESTED', 'MFM_WGS_SORT_REQUESTED', 'MFM_WGS_HAP_REQUESTED')
    previous = {name: os.environ.get(name) for name in names}
    os.environ.setdefault('MFM_WGS_ALIGN_REQUESTED', 'auto')
    os.environ.setdefault('MFM_WGS_SORT_REQUESTED', 'auto')
    os.environ['MFM_WGS_HAP_REQUESTED'] = '1'
    try:
        phase=RealignPhases();new=phase.align_config(cfg,lengths)
    finally:
        for name, value in previous.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value
    replacement=original[:align.start(thread_group)]+str(new.threads)+original[align.end(thread_group):]
    replacement,count=re.subn(r'("\$SAMTOOLS" sort -@ )[0-9]+ -m [^ ]+',
        lambda m:m.group(1)+str(new.sort_threads)+' -m '+shlex.quote(new.sort_mem),replacement,count=1)
    if count!=1:raise ValueError('HAP2 sort thread 교체 실패')
    replacement=re.sub(r'(--sort-threads )[0-9]+',lambda m:m.group(1)+str(new.sort_threads),replacement,count=1)
    if minimap is not None:
        replacement=re.sub(r'("\$SAMTOOLS" view)(?: -@ [0-9]+)? -u',
            lambda m:m.group(1)+' -@ '+str(new.input_decode_threads)+' -u',replacement)
    backup=script.with_suffix('.before_boundary.sh')
    if not backup.exists():shutil.copy2(script,backup)
    fd,temp=tempfile.mkstemp(prefix='.'+script.name,dir=script.parent);temp=Path(temp)
    try:
        with os.fdopen(fd,'w') as handle:handle.write(replacement);handle.flush();os.fsync(handle.fileno())
        subprocess.run(['bash','-n',str(temp)],check=True);os.chmod(temp,script.stat().st_mode&0o777);os.replace(temp,script)
    finally:temp.unlink(missing_ok=True)
    info={'schema':'mfm_hap_boundary_upgrade_v1','sample':sample,'hap':'HAP2',
          'platform':platform,'old_threads':cfg.threads,'threads':new.threads,'sort_threads':new.sort_threads,
          'old_sort_mem':cfg.sort_mem,'sort_mem':new.sort_mem,
          'input_decode_threads':new.input_decode_threads,
          'phase_plan':phase.events[-1],'original_script_sha256':hashlib.sha256(original.encode()).hexdigest(),
          'executed_script_sha256':hashlib.sha256(replacement.encode()).hexdigest(),'time':time.time()}
    _atomic(log/f'{sample}.HAP2.boundary_upgrade.json',info)
    with (pool/'pool.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);current=json.loads(requests_path.read_text());current['samples'][sample].update(status='applied',result=info);_atomic(requests_path,current)
    return command,info


def finish_boundary(info):
    if info is None:return
    # 완료된BWA pipeline은없으며,legacy부모의index최대5CPU와이후QC를함께보장함.
    phase=RealignPhases();granted=phase.reserve('output_qc',max_units=min(23,phase.total),
                                             min_units=min(5,phase.total),memory_per_unit=4*1024**3)
    info['post_alignment_qc_cpus']=granted
