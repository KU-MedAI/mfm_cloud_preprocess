"""정상 예약과 승인된 비교가 같은 단계별 CPU 설정을 사용함."""
def apply(environment,branch,stage,cpus):
    if branch.upper()=='RNA' and environment.get('MFM_CPU_ALLOCATION_MODE')=='fixed_sample':
        legacy=environment.get('MFM_CPU_POOL_LAYOUT')=='legacy'
        number=str(stage).zfill(2)
        names={('13' if legacy else '15'):('STAR_INDEX_THREADS',),
               ('14' if legacy else '16'):('STAR_ALIGN_THREADS',),
               ('11' if legacy else '13'):('RNA_COLLATE_THREADS','RNA_FASTQ_THREADS'),
               ('17' if legacy else '19'):('RNA_BAM_INDEX_THREADS',)}
        for name in names.get(number,()):environment[name]=str(cpus)
    if branch.upper()!='WGS':return
    if str(stage).zfill(2)=='06':environment['MFM_HAPLOTAG_CPU_BUDGET']=str(cpus)
    if str(stage).zfill(2)=='03':
        environment.update(VCF_CHROM_JOBS=str(max(1,min(22,cpus//3))),
                           VCF_THREADS_PER_CHROM=str(min(2,max(0,cpus-1))),
                           VCF_JOBS=str(max(1,min(22,cpus))))
