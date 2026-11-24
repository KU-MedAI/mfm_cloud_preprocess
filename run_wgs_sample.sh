#!/usr/bin/env bash
# WGS sample 1개를 Stage 01–12, 최종 검증과 선택적 cleanup까지 실행한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./lib/fixed_cpu_entry.sh
pipeline_fixed_cpu_entry "${1:?sample이 필요함}" bash "$(pwd)/run_wgs_sample.sh" "$@"
source ./config_v2.sh
source ./lib/state.sh
source ./lib/stage_runner.sh

S=${1:?사용법: $0 <sample>}
apply_sample_payload "$S"
WGS_START_AT=${WGS_START_AT:-01}
WGS_ALLOWED_STAGES=(01 02 03 04 05 06 07 08 09 10 11 12 23 24)
pipeline_validate_stage_range WGS "$WGS_START_AT" "${WGS_STOP_AFTER:-}" \
  "${WGS_ALLOWED_STAGES[@]}"

T0=$(date +%s)
require_manifest_sample_contract "$S"
PLATFORM=$(normalize_platform "$(sample_platform "$S")")
SEQ_SOURCE=$(sequence_source_for_sample "$S")
SEQ_VCF_SOURCE=$(longread_sequence_vcf_source_for_sample "$S")
HAPLOTAG_VCF_SOURCE=$(longread_haplotag_vcf_source_for_sample "$S")
log "═══ WGS 시작: $S platform=$PLATFORM sequence_source=$SEQ_SOURCE vcf_policy=$VCF_SOURCE_POLICY seq_vcf=$SEQ_VCF_SOURCE haplotag_vcf=$HAPLOTAG_VCF_SOURCE project=$(sample_project_id "$S") ═══"

RUN_STATUS_DIR=${RUN_STATUS_DIR:-$(sample_status_root "$S")}
mkdir -p "$RUN_STATUS_DIR"
PIPELINE_SAMPLE_ID=$S
PIPELINE_STAGE_LABEL=WGS
PIPELINE_START_AT=$WGS_START_AT
PIPELINE_STOP_AFTER=${WGS_STOP_AFTER:-}
STATE_FILE=$RUN_STATUS_DIR/$S.state.tsv
STAGE_FILE=$RUN_STATUS_DIR/$S.stages.tsv
CURRENT_STAGE=init
LAST_SUCCESS_STAGE=none
pipeline_install_sample_trap

if sample_wgs_reads_required "$S"; then
  pipeline_run_stage 01 01_prepare_analysis_bam.sh
  pipeline_run_stage 02 02_preflight.sh
else
  pipeline_skip_stage 01 "sequence_only: WGS read 입력 처리 제외"
  pipeline_skip_stage 02 "sequence_only: read QC는 비적용, 개인서열/VCF QC는 유지"
fi
pipeline_run_stage 03 03_prepare_vcf.sh
if longread_needs_clair3 "$S"; then
  pipeline_run_stage 04 04_call_longread_clair3.sh
  pipeline_run_stage 05 05_clair3_qc.sh
elif is_long_read_sample "$S"; then
  pipeline_skip_stage 04 "검증된 public VCF 사용으로 Clair3 생략"
  pipeline_skip_stage 05 "public VCF는 Clair3 QC 대상 아님"
else
  pipeline_skip_stage 04 "short-read는 Clair3 long-read self-calling 대상 아님"
  pipeline_skip_stage 05 "short-read는 Clair3 QC 대상 아님"
fi
if sample_wgs_reads_required "$S"; then
  pipeline_run_stage 06 06_haplotag.sh
else
  pipeline_skip_stage 06 "sequence_only: WGS haplotag 불필요"
fi
if is_long_read_sample "$S" && [[ "$SEQ_SOURCE" == denovo ]]; then
  pipeline_run_stage 07 07_denovo_assembly.sh
  pipeline_run_stage 08 08_anchor_assembly.sh
else
  pipeline_skip_stage 07 "SEQUENCE_SOURCE가 denovo인 long-read만 실행"
  pipeline_skip_stage 08 "de novo assembly anchoring 대상 아님"
fi
pipeline_run_stage 09 09_personal_fasta.sh
# 모델 입력 window는 read가 개인 HAP에 재정렬된 뒤 정의한다.
if sample_hap_realign_required "$S"; then
  pipeline_run_stage 10 10_realign_wgs.sh
else
  pipeline_skip_stage 10 "$(sample_output_mode "$S"): 개인 HAP WGS 재정렬 비적용"
fi
pipeline_run_stage 11 11_windows.sh
pipeline_run_stage 12 12_compact_zarr.sh

if [[ "${WGS_DEFER_FINALIZE:-${WGS_DEFER_CLEANUP:-0}}" != 1 ]]; then
  pipeline_run_stage 23 23_validate_sample_package.sh wgs
  if [[ "${CLEAN_INTERMEDIATES_ON_SUCCESS:-1}" == 1 ]]; then
    pipeline_run_stage 24 24_cleanup_intermediates.sh wgs
  else
    pipeline_skip_stage 24 "CLEAN_INTERMEDIATES_ON_SUCCESS=0"
  fi
else
  pipeline_skip_stage 23 "WGS-RNA pair는 RNA 완료 뒤 full package를 검증"
  pipeline_skip_stage 24 "WGS-RNA pair는 RNA 완료 뒤 cleanup"
fi

pipeline_complete_sample
log "═══ WGS 완료: $S / $(( ($(date +%s)-T0)/60 ))분 ═══"
