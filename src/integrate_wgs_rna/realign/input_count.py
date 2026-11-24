"""인덱스별전수입력집계를분할하되placed/unplaced unmapped와HP정책을보존함."""
from pathlib import Path
import json,os,multiprocessing,subprocess
import pysam
from concurrent.futures import ProcessPoolExecutor
from ..memory_policy import limit_worker_address_space,GIB,memory_limit
from ..input_parallel import cpu_budget
from .output_scan import _output_scan_reader

KEYS=('hp1','hp2','untagged','all_primary_mapped','unmapped')

def partition_counts(config,alignment,regions,selected,invalid_hp_as_untagged):
    result=dict.fromkeys((*KEYS,'_raw_unmapped','_raw_mapped'),0)
    with _output_scan_reader(str(alignment),str(config.reference_fasta) if config.reference_fasta else None,regions,0) as source:
        for region in regions:
            if isinstance(region,tuple):
                chrom,start,end=region;reads=source.fetch(chrom,start,end)
            else:start=end=None;reads=source.fetch(region)
            for read in reads:
                if start is not None and not start<=read.reference_start<end:continue
                if read.is_unmapped:result['_raw_unmapped']+=1
                elif read.reference_name in selected:result['_raw_mapped']+=1
                if read.is_secondary or read.is_supplementary:continue
                if read.is_unmapped:result['unmapped']+=1;continue
                if read.reference_name not in selected:continue
                result['all_primary_mapped']+=1
                if not read.has_tag('HP'):result['untagged']+=1;continue
                try:hp=int(read.get_tag('HP'))
                except (TypeError,ValueError):
                    if invalid_hp_as_untagged:raise
                    continue
                if hp==1:result['hp1']+=1
                elif hp==2:result['hp2']+=1
                elif hp==0 or invalid_hp_as_untagged:result['untagged']+=1
    return result

def count_indexed(config,alignment,contigs,invalid_hp_as_untagged):
    with pysam.AlignmentFile(str(alignment),reference_filename=str(config.reference_fasta) if config.reference_fasta else None) as source:
        names=list(source.references);lengths=dict(zip(source.references,source.lengths))
        statistics=source.get_index_statistics();expected_mapped=sum(r.mapped for r in statistics if r.contig in set(contigs));expected_unmapped=sum(r.unmapped for r in statistics)+source.nocoordinate
    chosen=set(contigs);major=sorted([c for c in names if c in chosen],key=lambda c:lengths[c],reverse=True)
    other=[c for c in names if c not in chosen]
    chunk=int(os.environ.get('MFM_REALIGN_INPUT_CHUNK_BP',16*1024**2))
    if chunk<1:raise ValueError('MFM_REALIGN_INPUT_CHUNK_BP는양수여야함')
    regions=[[(c,start,min(lengths[c],start+chunk))] for c in major for start in range(0,lengths[c],chunk)]+([other] if other else [])+[['*']]
    tasks=[(config,str(alignment),r,chosen,invalid_hp_as_untagged) for r in regions]
    workers=min(len(tasks),cpu_budget());ram=memory_limit()*9//10
    for key,factor in [('SLURM_MEM_PER_NODE',1),('SLURM_MEM_PER_CPU',workers)]:
        if int(os.environ.get(key,'0'))>0:ram=min(ram,int(os.environ[key])*1024**2*factor*9//10)
    lease=os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    if lease:
        row=json.loads(Path(lease).read_text());workers=min(len(tasks),int(row['cpus']));ram=min(ram,int(row.get('memory_bytes',0)))
    if ram<6*GIB:raise ValueError('indexed입력집계최소RAM6GiB부족')
    workers=min(workers,max(1,(ram-4*GIB)//(2*GIB)))
    result=dict.fromkeys((*KEYS,'_raw_unmapped','_raw_mapped'),0)
    # caller의단계예약을유지함. 단계재배분이꺼진경우도후속정렬예약을잃지않음.
    with ProcessPoolExecutor(max_workers=max(1,workers),mp_context=multiprocessing.get_context('spawn'),initializer=limit_worker_address_space,initargs=(2*GIB,)) as pool:
        futures=[pool.submit(partition_counts,*task) for task in tasks]
        for future in futures:
            counts=future.result()
            for key in result:result[key]+=counts[key]
    if Path(alignment).suffix.lower()=='.bam' and (result['_raw_unmapped']!=expected_unmapped or result['_raw_mapped']!=expected_mapped):
        # 좌표없는특수배치/일치하지않는index의누락을허용하지않음.
        from .qc import _count_input_reads_one_pass
        print('[INPUT-COUNT] index coverage 불일치로원전수집계수행함',flush=True)
        return _count_input_reads_one_pass(config,Path(alignment),contigs,invalid_hp_as_untagged=invalid_hp_as_untagged,_force_serial=True)
    if Path(alignment).suffix.lower()=='.cram':
        # CRAI는BAM식record통계가없음. 모든좌표의primary unmapped만독립전수확인함.
        command=[getattr(config,'samtools_bin','samtools'),'view','-c','-f','4','-F','0x900','-@',str(max(0,min(4,workers-1)))]
        if config.reference_fasta:command+=['-T',str(config.reference_fasta)]
        command.append(str(alignment));result['unmapped']=int(subprocess.check_output(command,text=True))
    result.pop('_raw_mapped');result.pop('_raw_unmapped')
    print('[INPUT-COUNT] '+json.dumps({'backend':'indexed_intervals','tasks':len(tasks),'workers':workers,'chunk_bp':chunk,'reserved_memory_bytes':ram,'counts':result}),flush=True)
    return result
