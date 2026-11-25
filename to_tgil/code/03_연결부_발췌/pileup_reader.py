# 검토용 원문 발췌임. 단독 실행 모듈이 아니며 import/호출 문맥은 전체 전달본을 사용함.

# 원본: src/integrate_wgs_rna/window_reader.py:50-57
if self.pileup:
    from .pileup_features.schema import SCHEMA, LEGACY_SCHEMA, LEGACY_CHANNELS, resolve_channels
    self.pileup_schema=manifest.get('schema')
    self.pileup_feature_set=manifest.get('feature_set','legacy_reference12_v1')
    expected=(list(resolve_channels(manifest.get('channels'),self.pileup_feature_set))
              if self.pileup_schema==SCHEMA else list(LEGACY_CHANNELS) if self.pileup_schema==LEGACY_SCHEMA else None)
    if expected is None or manifest.get('channels')!=expected:raise ValueError('pileup manifest schema/채널 불일치')
    self.pileup_channels=expected

# 원본: src/integrate_wgs_rna/window_reader.py:92-95
if self.pileup:
    if depth.attrs.get('schema')!=self.pileup_schema or depth.attrs.get('channels')!=self.pileup_channels:raise ValueError('pileup schema/채널 불일치')
    from .pileup_features.schema import SCHEMA
    if self.pileup_schema==SCHEMA and depth.attrs.get('feature_set')!=self.pileup_feature_set:raise ValueError('pileup feature set 불일치')

# 원본: src/integrate_wgs_rna/window_reader.py:178-197
if self.pileup:
    g=entry['pileup']
    features=np.ascontiguousarray(g['features'][a:b])
    kind=self._chunks.read(g['mapping_kind'],a,b)
    quality=self._chunks.read(g['baseq_depth'],a,b)
    if features.shape!=(padded,len(self.pileup_channels)) or features.dtype!=np.dtype('float32') or not np.isfinite(features).all():raise ValueError('pileup feature shape/dtype/값 오류')
    result.update(pileup_features=features,pileup_channels=list(self.pileup_channels),pileup_feature_set=self.pileup_feature_set,
                  pileup_reference_mapped_mask=kind==1,pileup_insertion_mask=kind==2,
                  pileup_baseq_valid_mask=quality>0,
                  supplementary_depth=self._chunks.read(g['supplementary_depth'],a,b),
                  max_deletion_length=self._chunks.read(g['max_deletion_length'],a,b),
                  wgs_depth_coordinate_system='reference_projected_to_hap')
    from .pileup_features.schema import SCHEMA
    if self.pileup_schema==SCHEMA:
        for name in ('hp1_depth','hp2_depth','hp1_supplementary_depth','hp2_supplementary_depth'):
            result[name]=self._chunks.read(g[name],a,b)
        result.update(pileup_reference_base_code=self._chunks.read(g['reference_base_code'],a,b),
                      pileup_baf_valid_mask=self._chunks.read(g['baf_valid'],a,b)>0,
                      pileup_baf_site_mask=self._chunks.read(g['baf_site'],a,b)>0,
                      pileup_baf_heterozygous_mask=self._chunks.read(g['baf_heterozygous'],a,b)>0)
