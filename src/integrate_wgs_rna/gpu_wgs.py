"""기존 routing FASTQ를 받아 GPU 정렬하고 기존 sort/태그 복원으로 전달함."""
from pathlib import Path
import argparse,fcntl,json,os,shlex,shutil,subprocess,sys,tempfile,time
import pysam,resource
from contextlib import ExitStack
from .gpu_pool import GIB
from .parabricks_runtime import Runtime
from .gpu_fastq import split_fastq, cached_inputs
from .gpu_stage_resources import StageBudget
from .realign.align_memory import alignment_memory

def host_cpu_limit(fallback):
    from .gpu_compare_resources import cpu_limit
    requested=int(os.environ.get('MFM_CPU_BUDGET',str(fallback)))
    if os.environ.get('MFM_GPU_HOST_CPUS_PER_HAP'):requested=min(requested,int(os.environ['MFM_GPU_HOST_CPUS_PER_HAP']))
    return cpu_limit(requested)


def scratch_directory(work):
    root=os.environ.get('MFM_GPU_SCRATCH_ROOT')
    if not root:return Path(work)
    import hashlib
    directory=Path(root).expanduser().resolve()/hashlib.sha256(str(Path(work).resolve()).encode()).hexdigest()[:20]
    directory.mkdir(parents=True,exist_ok=True)
    return directory


def profile_hifi_input(fastq):
    """read를 바꾸지 않고 실제 길이에 맞춰 GPU seed/chain 버퍼를 산정함."""
    if os.environ.get('MFM_GPU_HIFI_AUTO_PROFILE','1')=='0':return None
    if os.environ.get('MFM_GPU_HIFI_AUTO_PROFILE','1')!='1':raise ValueError('MFM_GPU_HIFI_AUTO_PROFILE는0/1이어야 함')
    lengths=[]
    with pysam.FastxFile(str(fastq)) as source:
        for read in source:
            lengths.append(len(read.sequence))
            if len(lengths)==1000:break
    return profile_hifi_lengths(lengths)

def profile_hifi_lengths(lengths):
    flag=os.environ.get('MFM_GPU_HIFI_AUTO_PROFILE','1')
    if flag=='0':return None
    if flag!='1':raise ValueError('MFM_GPU_HIFI_AUTO_PROFILE는0/1이어야 함')
    if not lengths:return None
    lengths=sorted(lengths);p95=lengths[(95*len(lengths)+99)//100-1]
    estimated=max(15000,p95)
    values={'MFM_GPU_MINIMAP_ESTIMATED_READ_LENGTH':estimated,
            'MFM_GPU_MINIMAP_ESTIMATED_ANCHOR_PER_READ':max(7000,2*estimated)}
    for key,value in values.items():os.environ.setdefault(key,str(value))
    return {'sampled_records':len(lengths),'mean_read_length':sum(lengths)/len(lengths),
            'p95_read_length':p95,'max_read_length':max(lengths),
            'effective_parameters':{key:int(os.environ[key]) for key in values},
            'note':'앞1000개 길이 기반 성능 추정임. read 필터/서열/태그는 바꾸지 않으며 나머지 입력의 recovery는 별도 기록함.'}

def commands(reference,fastqs,output,platform,preset,threads,sample,num_gpus=1):
    common=['--ref',str(reference),'--in-fq',*map(str,fastqs),'--out-bam',str(output),
            '--num-gpus',str(num_gpus),'--gpusort','--gpuwrite','--read-group-sm',sample,
            '--read-group-id-prefix',sample]
    if platform=='illumina':
        result=['fq2bam',*common,'--no-markdups','--bwa-options=-Y -C',
                '--bwa-primary-cpus',str(threads),'--bwa-cpu-thread-pool','1',
                '--memory-limit',str(int(os.environ.get('MFM_GPU_HOST_MEMORY_BYTES',str(128*GIB)))//10**9)]
        tuning_keys=('MFM_GPU_CIGAR_ON_GPU','MFM_GPU_BWA_PARTITION')
        default='1' if num_gpus==4 and threads>=8 and not any(k in os.environ for k in tuning_keys) else '0'
        for key,option in [('MFM_GPU_CIGAR_ON_GPU','--cigar-on-gpu'),('MFM_GPU_BWA_PARTITION','--run-partition')]:
            value=os.environ.get(key,default)
            if value not in {'0','1'}:raise ValueError(key+'는0/1이어야 함')
            if value=='1':
                if option=='--run-partition':
                    if num_gpus<2 or num_gpus%2:raise ValueError('BWA partition은2GPU씩나눌수있는짝수장치가필요함')
                    result.extend(['--bwa-gpu-num-per-partition','2'])
                result.append(option)
        return result
    if preset not in {'map-hifi','map-ont','lr:hq'}:raise ValueError('GPU WGS preset 미지원: '+preset)
    result=['minimap2',*common,'--preset',preset,'--md','--copy-comment','--num-threads',str(threads),
            '--seeding-chaining-on-gpu','on']
    for variable,option in [('MFM_GPU_MINIMAP_ESTIMATED_READ_LENGTH','--estimated-read-length'),
                            ('MFM_GPU_MINIMAP_ESTIMATED_ANCHOR_PER_READ','--estimated-anchor-per-read'),
                            ('MFM_GPU_MINIMAP_LARGE_CHAINING_SIZE','--stage2-large-chaining-size')]:
        if variable in os.environ:
            value=int(os.environ[variable])
            if value<=0:raise ValueError(variable+'는양수여야 함')
            result.extend([option,str(value)])
    if os.environ.get('MFM_GPU_MINIMAP_VERBOSE')=='1':result.append('--verbose')
    return result

def _restore_bam_python(output, read_group, sample, duplicate_mode=None, destination='-',output_mode='wb0'):
    """RG와 운반용 OD를 한번에 복원하고 SAM 텍스트 변환 없이 BAM을 전달함."""
    fields=read_group.replace('\\t','\t').split('\t')
    rg=dict(field.split(':',1) for field in fields if ':' in field)
    with pysam.AlignmentFile(str(output),'rb') as source:
        header=source.header.to_dict();header['RG']=[rg]
        if header.get('HD',{}).get('SO')!='coordinate':raise ValueError('Parabricks 정렬 BAM의 좌표 정렬 보증이 없음')
        with pysam.AlignmentFile(str(destination),output_mode,header=header) as dest:
            for read in source:
                read.set_tag('RG',rg['ID'])
                if duplicate_mode is not None and read.has_tag('OD'):
                    if int(read.get_tag('OD'))==1 and duplicate_mode=='preserve_source':
                        read.is_duplicate=True
                    read.set_tag('OD',None)
                dest.write(read)

def restore_bam(output, read_group, sample, duplicate_mode=None, destination='-',output_mode='wb0'):
    mode=os.environ.get('MFM_GPU_TAG_RESTORE','auto')
    if mode not in {'auto','native','python'}:raise ValueError('태그 복원은 auto|native|python이어야 함')
    started=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN);own=resource.getrusage(resource.RUSAGE_SELF)
    def metrics(backend,decode=0,tool=None):
        after=resource.getrusage(resource.RUSAGE_CHILDREN);current=resource.getrusage(resource.RUSAGE_SELF)
        return {'backend':backend,'decode_threads':decode,'native_tool':str(tool) if tool else None,
                'wall_seconds_including_downstream_write':time.monotonic()-started,
                'cpu_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime+current.ru_utime+current.ru_stime-own.ru_utime-own.ru_stime,
                'scope':'BAM 해독·RG/OD 복원·pipe 출력 구간임. downstream CRAM 쓰기 역압을 wall에 포함함.'}
    if mode=='python':
        _restore_bam_python(output,read_group,sample,duplicate_mode,destination,output_mode);return metrics('python')
    from .record_stats_parallel import native_tool,NativeUnavailable
    try:tool=native_tool('restore_tags')
    except NativeUnavailable:
        if mode=='native':raise
        _restore_bam_python(output,read_group,sample,duplicate_mode,destination,output_mode);return metrics('python_fallback')
    rg=dict(field.split(':',1) for field in read_group.replace('\\t','\t').split('\t') if ':' in field)
    with pysam.AlignmentFile(str(output),'rb') as source:header=source.header.to_dict()
    if header.get('HD',{}).get('SO')!='coordinate':raise ValueError('Parabricks 정렬 BAM의 좌표 정렬 보증이 없음')
    header['RG']=[rg]
    with tempfile.TemporaryDirectory(prefix='mfm_restore_') as tmp:
        path=Path(tmp)/'header.sam';path.write_text(str(pysam.AlignmentHeader.from_dict(header)))
        decode=max(0,int(os.environ.get('MFM_GPU_RESTORE_DECODE_THREADS','0')))
        subprocess.run([str(tool),str(output),str(destination),str(path),rg['ID'],duplicate_mode or 'unset',output_mode,str(decode)],check=True)
    return metrics('native_htslib',decode,tool)


def bwa_reference(target,work,budget=None,*,metrics=None):
    from .gpu_bwa_index import ensure_index
    return ensure_index(target,work,budget,metrics)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ['reference','platform','preset','sample','read-group','work','metadata']:p.add_argument('--'+arg,required=True)
    p.add_argument('--threads',required=True,type=int)
    p.add_argument('--out-cram')
    p.add_argument('--input-source')
    p.add_argument('--input-cache-key')
    p.add_argument('--duplicate-mode',choices=['preserve_source','none'])
    a=p.parse_args(argv)
    from .stage_gpu_config import apply_stage_environment
    apply_stage_environment(os.environ)
    from .gpu_compare_resources import enabled,cpu_limit
    if enabled():
        a.threads=host_cpu_limit(68)
        os.environ['MFM_GPU_STAGE_CPUS']=str(a.threads)
    os.environ['MFM_GPU_HOST_CPU_LIMIT']=str(a.threads);runtime=Runtime();runtime.validate()
    if not os.environ.get('MFM_GPU_HAP_ID') and os.environ.get('MFM_CPU_POOL_LEASE_FILE'):
        # 이미 생성된 순차 shell도 입력 EOF 뒤에는 안전하게 장치 대기 예약을 반환한다.
        import re
        row=json.loads(Path(os.environ['MFM_CPU_POOL_LEASE_FILE']).read_text())
        match=re.fullmatch(re.escape(a.sample)+r'\.HAP([12])\.gpu\.json',Path(a.metadata).name)
        if match and row.get('effective_alignment',{}).get('hap_jobs')==1:
            os.environ['MFM_GPU_HAP_ID']=match.group(1)
    os.environ.setdefault('MFM_GPU_STAGE_CPUS',str(host_cpu_limit(a.threads)))
    budget=StageBudget()
    work=Path(a.work).resolve();work.mkdir(parents=True,exist_ok=True)
    reference=Path(a.reference).resolve()
    scratch=scratch_directory(work)
    input_root=Path(os.environ.get('MFM_GPU_INPUT_CACHE_ROOT',str(scratch))).expanduser().resolve();input_root.mkdir(parents=True,exist_ok=True)
    gpu_count=int(os.environ.get('MFM_GPU_ALIGNMENT_GPUS','1'))
    if gpu_count<1:raise ValueError('정렬 GPU 개수는 양수여야 함')
    index_info={}
    def prepare_index(target,where,stage_budget=None):
        return bwa_reference(target,where,stage_budget,metrics=index_info)
    with tempfile.TemporaryDirectory(prefix='gpu_align_',dir=scratch) as tmp, ExitStack() as stack:
        directory=Path(tmp);paired=a.platform=='illumina'
        started=time.time();reused=False;stream=None;expected_input_bytes=0
        preparation=None
        if a.input_source:
            if not a.input_cache_key:raise ValueError('공통 FASTQ의 입력 지문이 없음')
            identity=json.loads(a.input_cache_key)
            expected_input_bytes=identity['alignment'][1] or 0
            input_mode=os.environ.get('MFM_GPU_INPUT_MODE','auto')
            if input_mode not in {'auto','cache','stream'}:raise ValueError('GPU 입력은 auto|cache|stream이어야 함')
            overlap=os.environ.get('MFM_GPU_PREPARE_OVERLAP','1')
            if overlap not in {'0','1'}:raise ValueError('MFM_GPU_PREPARE_OVERLAP은0/1이어야 함')
            index_ready=None
            if paired and overlap=='1':
                from .gpu_bwa_index import lookup
                budget.phase('bwa_index_lookup',1,4*GIB)
                index_ready=lookup(reference,scratch,index_info)
                if index_ready is not None:reference=index_ready
            if paired and index_ready is None and overlap=='1' and input_mode in {'auto','cache','stream'}:
                from .gpu_prepare import prepare_short
                from .gpu_bwa_index import cache_root
                prepared=prepare_short(reference,scratch,a.input_source,a.input_cache_key,input_root,budget,
                    int(os.environ.get('MFM_GPU_FASTQ_PRODUCER_CPUS','3')),expected_input_bytes,prepare_index,
                    index_root=cache_root(scratch))
                if prepared:
                    reference,fastqs,count,reused,preparation=prepared
                    input_mode='prepared'
            if paired and preparation is None and index_ready is None:
                reference=prepare_index(reference,scratch,budget)
            if input_mode=='auto':
                producer=int(os.environ.get('MFM_GPU_FASTQ_PRODUCER_CPUS','3'))
                from .gpu_stream_inputs import stream_cpu_budget
                required=stream_cpu_budget(budget.cpus,paired,producer)
                input_mode='stream' if required<=budget.cpus-2 else 'cache'
            if input_mode not in {'cache','stream','prepared'}:raise ValueError('GPU 입력은 auto|cache|stream이어야 함')
            if input_mode=='stream':
                from .gpu_stream_inputs import StreamInputs
                stream=stack.enter_context(StreamInputs(a.input_source,a.input_cache_key,input_root,paired,budget,
                    int(os.environ.get('MFM_GPU_FASTQ_PRODUCER_CPUS','3')),expected_input_bytes))
                fastqs=stream.paths;reused=stream.reused;count=stream.result[1] if reused else None
            elif input_mode=='cache':
                fastqs,count,reused=cached_inputs(a.input_source,a.input_cache_key,input_root,paired,budget,
                    int(os.environ.get('MFM_GPU_FASTQ_PRODUCER_CPUS','3')),expected_input_bytes)
        else:
            if paired:reference=prepare_index(reference,scratch,budget)
            fastqs=[directory/f'input_R{i+1}.fq.gz' for i in range(2 if paired else 1)]
            count=split_fastq(sys.stdin.buffer,fastqs,paired,cpus=a.threads)
        input_seconds=time.time()-started
        input_profile=None
        if a.platform=='hifi' and a.preset=='map-hifi':
            input_profile=profile_hifi_lengths(stream.lengths) if stream and not reused else profile_hifi_input(fastqs[0])
        output=directory/'gpu.bam';native_cpus=max(1,a.threads-1-(stream.cpus if stream and not reused else 0))
        if native_cpus<1:raise ValueError('입력 생산자와 GPU 보조 CPU 예산 초과')
        budget.write_paths=[output]
        if os.environ.get('MFM_GPU_FINAL_PART'):budget.write_paths.append(Path(os.environ['MFM_GPU_FINAL_PART']))
        budget.disk_peak_bytes=(expected_input_bytes*4 if stream and not reused else sum(p.stat().st_size for p in fastqs)*3)+32*GIB
        with budget.gpu(int(os.environ.get('MFM_GPU_ALIGNMENT_VRAM_BYTES',str(24*GIB))),purpose='wgs_alignment',count=gpu_count,
                        producer_cpus=stream.cpus if stream and not reused else 0) as gpu:
            native_cpus=max(1,budget.cpus-1-(stream.cpus if stream and not reused else 0))
            if a.platform=='illumina':
                explicit=any(k in os.environ for k in ['MFM_GPU_BWA_PARTITION','MFM_GPU_CIGAR_ON_GPU'])
                if not explicit:
                    default='1' if gpu_count==4 and native_cpus>=8 else '0'
                    os.environ.update(MFM_GPU_BWA_PARTITION=default,MFM_GPU_CIGAR_ON_GPU=default)
                partition=os.environ.get('MFM_GPU_BWA_PARTITION','0')=='1'
                pools=gpu_count//2 if partition else 1
                if native_cpus<pools+(pools if partition else 1):raise ValueError('GPU BWA primary/pool CPU 예산 부족')
                native_cpus=native_cpus-pools
                if partition:native_cpus=max(pools,native_cpus//pools*pools)
            command=commands(reference,fastqs,output,a.platform,a.preset,native_cpus,a.sample,gpu_count)
            command+=['--tmp-dir',str(directory)]
            info=runtime.run(command,gpu,[reference,*fastqs,output,directory],log=Path(a.metadata).with_suffix('.log'),check=stream.check if stream else None,cpu_limit=budget.cpus)
            budget.output()  # GPU 반환 전에 CPU 후처리 몫을 원자 교환함.
        if stream:cache_paths,count,reused=stream.finish()
        if not output.is_file() or output.stat().st_size==0:raise RuntimeError('GPU 정렬 BAM 미생성')
        legacy_paired=paired and not a.input_source
        if a.out_cram:
            from .gpu_cram import encode
            output_restore=encode(output,reference,a.out_cram,a.read_group,a.sample,a.duplicate_mode,cpus=budget.cpus,work=scratch)
        else:
            output_restore=restore_bam(output,a.read_group,a.sample,a.duplicate_mode,output_mode='wh' if legacy_paired else 'wb0')
        info.update(schema='mfm_gpu_wgs_alignment_v1',input_records=count,platform=a.platform,preset=a.preset,
                    routing='upstream 기존 routing FASTQ 사용함',duplicate_policy='upstream OD + 기존 restore_duplicate_flags 사용함' if paired else '기존 long-read 정책 사용함',
                    supplementary_clipping='Parabricks preset 기본값; supplementary는 기존 모델 입력에서 제외함',
                    no_bqsr=True,no_new_duplicate_marking=True,input_seconds=input_seconds,input_preparation_overlap=preparation,bwa_index_cache=index_info if paired else None,
                    input_cache_receipt=str((cache_paths if stream else fastqs)[0].parent/'complete.json'),input_reused=reused,input_streamed=bool(stream and not reused),native_cpu_threads=native_cpus,input_producer_cpu_budget=stream.cpus if stream and not reused else 0,gpu_input_performance_profile=input_profile,host_resource_phases=budget.events,scratch_directory=str(scratch),
                    local_intermediate_cram=os.environ.get('MFM_GPU_LOCAL_CRAM'),
                    output_restore=output_restore,output_stream='cram' if a.out_cram else ('sam_legacy' if legacy_paired else 'uncompressed_bam'))
        path=Path(a.metadata);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(info,ensure_ascii=False,indent=2)+'\n')
    return 0

def write_script(config,sample,hap,alignment,target,output,contigs,work,log,prepared):
    from .personal_realign import q,routing_fastq_commands,illumina_fastq_source_command
    from .realign.policy import sam_rg_platform
    paired=prepared is not None
    if paired:
        source=illumina_fastq_source_command(prepared,log,sample,hap,config.fastq_threads)
    else:
        source='''(
  grp() {
    local code=$1; shift
    "$SAMTOOLS" view -u "${VIEW_REF[@]}" "$@" | "$SAMTOOLS" fastq -T HP,PS -
    local s=("${PIPESTATUS[@]}")
    [[ "${s[0]}" == 0 && "${s[1]}" == 0 ]] || exit "$code"
  }
'''+routing_fastq_commands(config,hap,' '.join(q(c) for c in contigs))
        if config.include_unmapped:
            source+='''  "$SAMTOOLS" view -u "${VIEW_REF[@]}" -F 0x900 "$INPUT" '*' | "$SAMTOOLS" fastq -T HP,PS - | awk 'NR%4==1{print $0 "\\tXU:i:1"; next} {print}'
  s=("${PIPESTATUS[@]}")
  [[ "${s[0]}" == 0 && "${s[1]}" == 0 && "${s[2]}" == 0 ]] || exit 93
'''
        source+=')'
        source=source.replace('"$SAMTOOLS" view -u',f'"$SAMTOOLS" view -@ {getattr(config,"input_decode_threads",0)} -u')
    rg=f'@RG\\tID:{sample}\\tSM:{sample}\\tPL:{sam_rg_platform(config.platform)}\\tDS:personal_HAP{hap};backend=parabricks'
    command=[sys.executable,'-B','-m','integrate_wgs_rna.gpu_wgs','--reference',str(target),
             '--platform','illumina' if paired else config.platform,'--preset',config.preset,
             '--sample',sample,'--read-group',rg,'--work',str(work),
             '--metadata',str(log/f'{sample}.HAP{hap}.gpu.json'),'--threads',str(config.threads)]
    # all_to_both는 개인 참조만 다르고 입력 read/태그는 같으므로 한번 준비한다.
    from .realign.policy import effective_hap_read_mode
    identity={}
    for label,path in [('alignment',prepared.path if paired else alignment),('reference',config.reference_fasta)]:
        if path:
            stat=Path(path).stat() if Path(path).exists() else None
            identity[label]=[str(Path(path).resolve()),stat.st_size if stat else None,stat.st_mtime_ns if stat else None]
    identity.update(schema='gpu_fastq_contract_v2',contigs=contigs,mode=effective_hap_read_mode(config),
                    read_routing_policy=config.read_routing_policy,unknown_read_policy=config.unknown_read_policy,
                    include_unmapped=config.include_unmapped,paired=paired,
                    name_collated=getattr(prepared,'name_collated',False),duplicate_mode=config.illumina_duplicate_mode)
    if config.read_routing_policy!='all_to_both':identity['hap']=hap
    key=json.dumps(identity,sort_keys=True)
    source_path=log/f'{sample}.HAP{hap}.input.sh'
    source_path.write_text('#!/usr/bin/env bash\nset -euo pipefail\n'+
        f'SAMTOOLS={q(config.samtools_bin)}\nINPUT={q(prepared.path if paired else alignment)}\n'+
        f"VIEW_REF=({'-T '+q(config.reference_fasta) if config.reference_fasta else ''})\n"+
        "UNTAGGED_EXPR='[HP]==0 || ![HP]'\n"+f'COLLATE_PREFIX={q(work/f"collate.{sample}.HAP{hap}")}\n'+source+'\n')
    command+=['--input-source',str(source_path),'--input-cache-key',key]
    if paired:command+=['--duplicate-mode',config.illumina_duplicate_mode]
    fai=Path(str(target)+'.fai')
    bases=sum(int(line.split('\t')[1]) for line in fai.read_text().splitlines()) if fai.exists() else max(1,Path(target).stat().st_size if Path(target).exists() else 3_000_000_000)
    memory=alignment_memory(config,{hap:bases},1,config.threads,config.sort_threads,
                            dict(os.environ,MFM_WGS_BACKEND='parabricks'))['reserved_memory_bytes']
    part=output.with_suffix(output.suffix+'.part')
    local_cram=None
    if os.environ.get('MFM_GPU_SCRATCH_ROOT'):
        local_cram=scratch_directory(work)/f'{sample}.HAP{hap}.{time.time_ns()}.cram'
        part=local_cram.with_suffix('.cram.part')
    publish=f'mv {q(part)} {q(output)}'
    if local_cram is not None:
        # 새 중간 파일만 로컬에 생성하며 기존 출력·원본을 이동하거나 덮어쓰지 않는다.
        publish=f'mv {q(part)} {q(local_cram)}\nln -s {q(local_cram)} {q(output)}'
    from .gpu_compare_resources import enabled,cpu_limit
    output_cpus=host_cpu_limit(config.threads+config.sort_threads+5);cram_threads=config.sort_threads;restore_decode=0
    if enabled():
        output_cpus=host_cpu_limit(68)
        restore_decode=min(8,max(0,(output_cpus-3)//4))
        cram_threads=max(0,output_cpus-restore_decode-3)
    # CRAM 인코딩은 sort 메모리가 아닌 read/quality·컨테이너 버퍼로 산정함.
    output_memory=(4*GIB+cram_threads*512*(1<<20)+restore_decode*64*(1<<20))*11//10
    command+=['--out-cram',str(part)]
    text=f'''#!/usr/bin/env bash
set -euo pipefail
SAMTOOLS={q(config.samtools_bin)}
INPUT={q(prepared.path if paired else alignment)}
VIEW_REF=({'-T '+q(config.reference_fasta) if config.reference_fasta else ''})
UNTAGGED_EXPR='[HP]==0 || ![HP]'
COLLATE_PREFIX={q(work/f'collate.{sample}.HAP{hap}')}
export MFM_GPU_HAP_ID={hap}
export MFM_GPU_FINAL_PART={q(part)}
export MFM_GPU_LOCAL_CRAM={q(local_cram) if local_cram is not None else ""}
export MFM_GPU_STAGE_CPUS={output_cpus}
export MFM_GPU_STAGE_MEMORY_BYTES={memory}
export MFM_GPU_OUTPUT_CPUS={output_cpus}
export MFM_GPU_RESTORE_DECODE_THREADS={restore_decode}
export MFM_GPU_OUTPUT_MEMORY_BYTES={output_memory}
export MFM_GPU_FASTQ_PRODUCER_CPUS={getattr(config,'input_decode_threads',0)+3}
{shlex.join(command)}
{publish}
'''
    path=log/f'{sample}.HAP{hap}.realign.sh';path.write_text(text);path.chmod(0o755);return path

if __name__=='__main__':raise SystemExit(main())
