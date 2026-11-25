# 검토용 원문 발췌임. 단독 실행 모듈이 아니며 import/호출 문맥은 전체 전달본을 사용함.

# 원본: config_v2.sh:535-539
sample_output_mode() {
  local mode=${MFM_SHORT_READ_OUTPUT_MODE:-pileup_feature}
  [[ "$(sample_read_type "$1")" != long_read ]] || mode=${MFM_LONG_READ_OUTPUT_MODE:-pileup_feature}
  printf '%s\n' "$mode"
}

# 원본: config_v2.sh:540-540
sample_pileup_feature() { [[ "$(sample_output_mode "$1")" == pileup_feature ]]; }

# 원본: config_v2.sh:541-545
sample_hap_realign_required() {
  # 실측에서만 HAP 재정렬 비용을 추가함. reference pileup 입력은 변경하지 않음.
  if [[ "${MFM_PILEUP_BENCHMARK_REALIGN:-0}" == 1 ]] && sample_pileup_feature "$1"; then return 0; fi
  [[ "$(sample_output_mode "$1")" == full ]] || { [[ -n "$(sample_rna "$1")" && "$(sample_rna "$1")" != - ]] && ! sample_rna_deferred "$1" && [[ "${MFM_RNA_OUTPUT_MODE:-full}" != coverage_only ]]; }
}
