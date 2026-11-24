"""과거 기본30과 공통 BAF·non-ref·strand 보완의 저장 계약임."""
import numpy as np

LEGACY_SCHEMA = 'mfm_reference_pileup_projected_v1'
LEGACY_LAYOUT = 'hap_sequence_reference_pileup_v1'
LEGACY_CHANNELS = (
    'log1p_depth', 'mean_mapq', 'mean_baseq', 'strand_forward_ratio',
    'read_A_ratio', 'read_C_ratio', 'read_G_ratio', 'read_T_ratio', 'read_N_ratio',
    'insertion_ratio', 'deletion_after_event_ratio', 'deletion_after_max_len',
)
HAPLOTYPE_CHANNELS = (
    'log1p_depth', 'read_A_ratio', 'read_C_ratio', 'read_G_ratio', 'read_T_ratio',
    'read_N_ratio', 'insertion_ratio', 'deletion_after_event_ratio', 'deletion_after_max_len',
)
ALLELE_CHANNELS = ('vcf_ref_allele_fraction', 'baf', 'minor_allele_fraction')
EXTRA_CHANNELS = ('read_supports_variant_ratio', 'base_diff_from_ref_ratio', 'strand_reverse_ratio')
CHANNELS = (LEGACY_CHANNELS + tuple('hp1_'+name for name in HAPLOTYPE_CHANNELS)
            + tuple('hp2_'+name for name in HAPLOTYPE_CHANNELS) + ALLELE_CHANNELS + EXTRA_CHANNELS)
SCHEMA = 'mfm_reference_pileup_projected_v2'
LAYOUT = 'hap_sequence_reference_pileup_v2'
FEATURE_SET = 'reference_hp_baf_v2_36ch'

COUNT_FIELDS = (
    'A', 'C', 'G', 'T', 'N', 'baseq_n', 'forward', 'deletion', 'insertion',
    'deletion_event', 'max_deletion', 'supplementary',
)
SUM_FIELDS = ('mapq', 'baseq')
ARRAYS = {'features': np.dtype('float32'), 'mapping_kind': np.dtype('uint8'),
          'reference_base_code': np.dtype('uint8'), 'baf_site': np.dtype('uint8'),
          'baf_heterozygous': np.dtype('uint8')}
for prefix in ('', 'hp1_', 'hp2_'):
    ARRAYS.update({prefix + name: np.dtype('uint32') for name in
                   ('depth', 'supplementary_depth', 'max_deletion_length')})
ARRAYS['baseq_depth'] = np.dtype('uint32')
ARRAYS['baf_valid'] = np.dtype('uint8')

INSERTION_CHANNELS = tuple(i for i, name in enumerate(CHANNELS)
                          if name.removeprefix('hp1_').removeprefix('hp2_') == 'insertion_ratio')

POLICY = dict(reference_basis='reference_aligned_projected_to_hap', baq=False,
              min_base_quality=0, mate_overlap='count_each_alignment', depth_cap=None,
              insertion='anchor_insertion_events_only', deletion_length='observed_read_cigar',
              missing_quality='zero_with_baseq_depth',
              segment_policy='primary_features_supplementary_separate',
              haplotype_basis='input_bam_HP_tag_1_or_2; untagged_in_shared_only',
              baf='biallelic_ACGT_SNV_ALT_count/(REF_count+ALT_count)',
              baf_missing='zero_with_explicit_site_and_group_valid_masks',
              non_reference='ACGT_reference_only; reference_base_code_marks_ambiguous_sites',
              feature_set=FEATURE_SET)

assert len(LEGACY_CHANNELS) == 12 and len(HAPLOTYPE_CHANNELS) == 9 and len(CHANNELS) == 36
