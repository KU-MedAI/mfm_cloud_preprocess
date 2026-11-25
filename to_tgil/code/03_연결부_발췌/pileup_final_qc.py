# 검토용 원문 발췌임. 단독 실행 모듈이 아니며 import/호출 문맥은 전체 전달본을 사용함.

# 원본: scripts/validate_sample_package.py:199-208
if payload=='pileup_feature':
    from integrate_wgs_rna.pileup_features.schema import SCHEMA, CHANNELS, FEATURE_SET, LEGACY_SCHEMA, LEGACY_CHANNELS
    if manifest.get('schema')==SCHEMA and manifest.get('channels')==list(CHANNELS) and manifest.get('feature_set')==FEATURE_SET:
        a.pileup_feature_set=FEATURE_SET
    elif manifest.get('schema')==LEGACY_SCHEMA and manifest.get('channels')==list(LEGACY_CHANNELS):
        a.pileup_feature_set='legacy_reference12_v1'
    else:raise ValueError('알 수 없는 pileup 채널 계약')
    metadata=read_json(a.root/'metadata/sample.json')
    if manifest.get('schema')==SCHEMA and (metadata.get('pileup_feature_set')!=FEATURE_SET or metadata.get('pileup_feature_schema')!=SCHEMA):
        raise ValueError('pileup metadata/실제 채널 계약 불일치')

# 원본: scripts/validate_sample_package.py:232-234
if payload=='pileup_feature':
    required=['contig_metadata.json','windows.tsv','hap_sequence.zarr','pileup.zarr','filter_qc.json']
    if sequence_results[index][1] is not None:raise ValueError(sequence_results[index][1])
