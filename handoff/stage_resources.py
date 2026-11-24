"""public VCF 기반 HiFi의 도구·buffer별 보수적 투입 RAM 예측임."""

GIB = 1024**3


def longread_memory_plan(row, config=None):
    config = config or {}
    env = dict(config.get('pipeline_env', {}))
    env.update(config.get('sample_pipeline_env', {}).get(row.get('sample_id',''), {}))
    mode = env.get('MFM_LONGREAD_MEMORY_MODEL', 'stages_v1')
    if mode not in {'stages_v1', 'legacy'}:
        raise ValueError('MFM_LONGREAD_MEMORY_MODEL은 stages_v1|legacy여야 함')
    sequence = env.get('SEQUENCE_SOURCE', row.get('sequence_source', 'vcf_consensus'))
    if (mode == 'legacy' or env.get('MFM_MEMORY_PROFILE','balanced') != 'balanced' or row.get('platform') != 'hifi' or row.get('read_type') != 'long_read'
            or row.get('vcf_source') != 'public' or row.get('phase_source') != 'provided'
            or sequence != 'vcf_consensus'):
        return None
    # 설치 전 doctor/dry-run에서 pysam을 요구하지 않도록 실제 RAM 계산 때 읽음.
    from integrate_wgs_rna.public_vcf_parallel import entry_memory
    cpus = config.get('sample_cpus', 72)
    cpus = 72 if cpus == 'auto' else max(1, int(cpus))
    # GRCh38/decoy 전체 reference를 포함하는 보수적 상한. 개인 HAP는 상염색체만 사용함.
    reference = max(4_000_000_000, int(env.get('MFM_LONGREAD_REFERENCE_BASES', 4_000_000_000)))
    sort_threads = min(24, max(1, cpus//3))
    prepare = reference*16 + 40*GIB + cpus*(GIB//2) + 2*GIB*(sort_threads+1)
    prepare_estimate=(prepare*11+9)//10
    if row.get('wgs_input_format') == 'bam':
        # aligned BAM의기존transcode는64GiB로예약함. 정렬버퍼최대8×2GiB와native집계28GiB가서로다른단계로실행됨.
        prepare_estimate=64*GIB
    align = 2*(reference*16 + 28*GIB + 10*GIB) + cpus*(GIB//4)
    stages = {'prepare_alignment': prepare_estimate,
              'public_vcf_prepare': entry_memory(cpus),
              'haplotag': 128*GIB,
              'consensus': max(128*GIB,4*GIB*(cpus//2+1)),
              'consensus_entry': 4*GIB,
              'personal_minimap_sort_qc': (align*11+9)//10,
              'compact_entry': 4*GIB}
    adaptive_realign=env.get('MFM_REALIGN_PHASE_RESOURCES','auto')!='0'
    if adaptive_realign:stages['personal_realign_entry']=4*GIB
    # compact의 실제 입력통계별worker예약은 실행시 별도 global pool에서 동적으로 결정함.
    floor = 256*GIB
    if row.get('data_group') == 'wgs_rna_pair':
        from integrate_wgs_rna.rna_memory import estimate_star_memory
        stages['rna_star_index'] = estimate_star_memory('genomeGenerate',cpus,reference_bases=reference)['memory_bytes']
        stages['rna_star_align'] = estimate_star_memory('alignReads',cpus,reference_bases=reference)['memory_bytes']
        stages['rna_star_and_sort'] = max(stages['rna_star_index'],stages['rna_star_align'])
    peak = max(floor, *stages.values())
    peak = (peak+GIB-1)//GIB*GIB
    return {'schema': 'mfm_longread_stage_memory_v1', 'mode': mode,
            'reference_upper_bases': reference, 'cpu_upper': cpus,
            'stage_estimates_bytes': stages, 'admission_bytes': peak,
            'input_is_streamed': True, 'assembly_and_self_call_included': False,
            'adaptive_compact_reserved_separately': True,
            'adaptive_realign':adaptive_realign,
            'stage_admission':env.get('MFM_STAGE_ADMISSION','1')=='1'}


def starting_memory(row,config=None):
    """시작단계만예약하고다음단계는공유풀에서다시확보함."""
    from .sequence_mode import selected,needs_wgs
    if selected(config or {},row) and not needs_wgs(config or {},row):return 20*GIB
    plan=longread_memory_plan(row,config)
    if plan is None:return None
    config=config or {};env=dict(config.get('pipeline_env',{}));env.update(config.get('sample_pipeline_env',{}).get(row.get('sample_id',''),{}))
    if env.get('MFM_STAGE_ADMISSION','1') not in {'0','1'}:raise ValueError('MFM_STAGE_ADMISSION은0|1이어야함')
    if env.get('MFM_STAGE_ADMISSION','1')=='0':return plan['admission_bytes']
    stage=str(env.get('WGS_START_AT','01')).zfill(2)
    labels={'01':'prepare_alignment','02':None,'03':'public_vcf_prepare','06':'haplotag',
            '09':'consensus_entry','10':'personal_realign_entry' if plan['adaptive_realign'] else 'personal_minimap_sort_qc','11':'windows','12':'compact_entry'}
    label=labels.get(stage)
    if stage=='02':return 64*GIB
    if stage=='11':return 8*GIB
    if stage=='23':return 64*GIB
    return plan['stage_estimates_bytes'].get(label,plan['admission_bytes'])
