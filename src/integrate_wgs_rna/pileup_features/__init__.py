"""Reference 관측을 기존 개인 HAP window에 대응시키는 출력 모드임."""

MODE = 'pileup_feature'
LAYOUT = 'hap_sequence_reference_pileup_v1'
CHANNELS = (
    'log1p_depth', 'mean_mapq', 'mean_baseq', 'strand_forward_ratio',
    'read_A_ratio', 'read_C_ratio', 'read_G_ratio', 'read_T_ratio', 'read_N_ratio',
    'insertion_ratio', 'deletion_after_event_ratio', 'deletion_after_max_len',
)
SCHEMA = 'mfm_reference_pileup_projected_v1'
