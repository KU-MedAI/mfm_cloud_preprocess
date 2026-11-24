"""short-read 서열전용 출력과 RNA의 WGS 임시근거를 명시적으로 구분함."""
from pathlib import Path

LAYOUT='sequence_4bit_acgtnp__windowed'
DEFERRED_SCOPE='dna_sequence_only_rna_deferred'

def output_mode(config,row):
    from .resume import sample_environment
    env=sample_environment(config,row.get('sample_id',''))
    key='MFM_LONG_READ_OUTPUT_MODE' if row.get('read_type')=='long_read' else 'MFM_SHORT_READ_OUTPUT_MODE'
    value=env.get(key,'pileup_feature')
    if value not in {'full','sequence_only','pileup_feature'}:raise ValueError(f'{key}의 출력 모드 오류: {value}')
    return value

def selected(config,row):
    return output_mode(config,row)=='sequence_only'

def defer_rna(config,row):
    from .resume import sample_environment
    value=sample_environment(config,row.get('sample_id','')).get('MFM_SEQUENCE_DEFER_RNA','0')
    if value not in {'0','1'}:raise ValueError('MFM_SEQUENCE_DEFER_RNA는 0/1임')
    return value=='1' and output_mode(config,row) in {'sequence_only','pileup_feature'} and bool(row.get('rna_path'))

def rna_required(config,row):
    return bool(row.get('rna_path')) and not defer_rna(config,row)

def completion_mode(config,row):
    return 'wgs' if defer_rna(config,row) else ('rna' if row.get('data_group')=='wgs_rna_pair' else 'wgs')

def needs_wgs(config,row):
    from .resume import sample_environment
    from rna_signal import output_mode
    return not selected(config,row) or (rna_required(config,row) and output_mode(sample_environment(config,row['sample_id']))!='coverage_only')

def validate_selection(config,rows):
    from .resume import sample_environment
    for row in rows:
        from rna_signal import output_mode as rna_output_mode
        mode=rna_output_mode(sample_environment(config,row['sample_id']))
        if mode=='coverage_only' and defer_rna(config,row):raise ValueError('RNA coverage_only와RNA보류는동시에설정할수없음')
        if output_mode(config,row)=='pileup_feature' and sample_environment(config,row['sample_id']).get('SEQUENCE_SOURCE','vcf_consensus')!='vcf_consensus':raise ValueError('pileup_feature는 검증된 consensus chain이 필요함')
        if not selected(config,row):continue
        if sample_environment(config,row['sample_id']).get('SEQUENCE_SOURCE','vcf_consensus')!='vcf_consensus':raise ValueError('sequence_only는vcf_consensus만지원함')
        public_long=row.get('read_type')=='long_read' and row.get('vcf_source')=='public'
        if row.get('phase_source')!='provided' or not (public_long or (row.get('platform')=='illumina' and row.get('vcf_source')=='illumina_panel')):
            raise ValueError('서열전용 short는 제공된 phased Illumina panel이 필요함')
        if rna_required(config,row) and mode!='coverage_only' and sample_environment(config,row['sample_id']).get('MFM_SEQUENCE_RNA_WGS_EVIDENCE')!='temporary':
            raise ValueError('RNA ASE/CNV 유지에는 pair의 WGS 임시처리가 필요함. MFM_SEQUENCE_RNA_WGS_EVIDENCE=temporary 정책 확인이 필요함')

def runtime_environment(config,row):
    mode=output_mode(config,row)
    if mode=='full':return {}
    from integrate_wgs_rna.pileup_features.schema import LAYOUT as PILEUP_LAYOUT
    layout=LAYOUT if mode=='sequence_only' else PILEUP_LAYOUT
    return dict(MFM_WGS_PAYLOAD=mode,COMPACT_LAYOUT_NAME=layout,
                MFM_SEQUENCE_DEFER_RNA='1' if defer_rna(config,row) else '0')

def is_sequence_package(root):
    import json
    manifests=list((Path(root)/'compact').glob('*/manifest.json'))
    return len(manifests)==1 and json.loads(manifests[0].read_text()).get('wgs_payload')=='sequence_only'

def uses_bundled_only(config,row,root):
    """전송 없는 준비의 입장만 판정함. 실제 worker는 기존 SHA 검증을 수행함."""
    if not selected(config,row) or needs_wgs(config,row):return False
    from .bundled_assets import inventory
    from .downloads import sample_resources,safe_path
    from .resume import sample_downloads
    base,records=inventory(config,root)
    resources=sample_downloads(config,row,sample_resources(root,row))
    for resource in resources:
        record=records.get(resource['path'])
        if resource.get('shared') and record is not None:
            path=safe_path(base,resource['path'])
            if not path.is_symlink() and path.is_file() and path.stat().st_size==record['size_bytes']:continue
        # 이미 검증된 개인 VCF도 전송 없이 worker의 기존 검증을 재사용함.
        import hashlib,json
        path=safe_path(config['data_root'],resource['path'])
        if not path.is_file() or path.is_symlink():return False
        receipt=Path(config['data_root'])/'.download_records'/(hashlib.sha256(str(path).encode()).hexdigest()+'.json')
        try:
            cached=json.loads(receipt.read_text());stat=path.stat()
        except (OSError,ValueError):return False
        signature=dict(size=stat.st_size,mtime_ns=stat.st_mtime_ns,ctime_ns=stat.st_ctime_ns)
        if cached.get('status') not in {'verified','generated_index'} or cached.get('url')!=resource.get('url') or cached.get('signature')!=signature:return False
    return bool(resources)

def payload_matches(config,row,root):
    import json
    path=Path(root)/'qc/final_package_validation.json'
    report=json.loads(path.read_text()) if path.is_file() else {}
    if defer_rna(config,row) and report.get('completion_scope')!=('dna_pileup_feature_only_rna_deferred' if output_mode(config,row)=='pileup_feature' else DEFERRED_SCOPE):return False
    from .resume import sample_environment
    from rna_signal import output_mode as rna_output_mode
    if rna_required(config,row) and report.get('rna_payload','full')!=rna_output_mode(sample_environment(config,row['sample_id'])):return False
    observed=report.get('wgs_payload','full')
    requested=output_mode(config,row)
    if requested=='pileup_feature':
        from integrate_wgs_rna.pileup_features.schema import FEATURE_SET
        if report.get('pileup_feature_set')!=FEATURE_SET:return False
    return observed==requested
