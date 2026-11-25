# 검토용 원문 발췌임. 단독 실행 모듈이 아니며 import/호출 문맥은 전체 전달본을 사용함.

# 원본: handoff/sequence_mode.py:7-13
def output_mode(config,row):
    from .resume import sample_environment
    env=sample_environment(config,row.get('sample_id',''))
    key='MFM_LONG_READ_OUTPUT_MODE' if row.get('read_type')=='long_read' else 'MFM_SHORT_READ_OUTPUT_MODE'
    value=env.get(key,'pileup_feature')
    if value not in {'full','sequence_only','pileup_feature'}:raise ValueError(f'{key}의 출력 모드 오류: {value}')
    return value

# 원본: handoff/sequence_mode.py:50-56
def runtime_environment(config,row):
    mode=output_mode(config,row)
    if mode=='full':return {}
    from integrate_wgs_rna.pileup_features.schema import LAYOUT as PILEUP_LAYOUT
    layout=LAYOUT if mode=='sequence_only' else PILEUP_LAYOUT
    return dict(MFM_WGS_PAYLOAD=mode,COMPACT_LAYOUT_NAME=layout,
                MFM_SEQUENCE_DEFER_RNA='1' if defer_rna(config,row) else '0')

# 원본: handoff/sequence_mode.py:88-101
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
