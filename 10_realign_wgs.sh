#!/usr/bin/env bash
# 10 — WGS read를 HAP1/HAP2 personal FASTA에 재정렬.
#
# 입력은 HP/PS-tagged BAM이다. long-read production은 검증된 public VCF+WhatsHap
# 산출물을 사용하며, Clair3+WhatsHap은 명시적 재현/비교 실행에만 사용한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh
export WGS_AUTO_THREADS

S=${1:?사용법: $0 <sample>}
ROOT=$(sample_final_root "$S")
COHORT=$(sample_cohort "$S")
PLATFORM=$(normalize_platform "$(sample_platform "$S")")
SEQUENCE_SOURCE_RESOLVED=$(sequence_source_for_sample "$S")
# v3 기본 경로는 VCF consensus이다. assembly-direct 선택 기능이 있는
# 배포본과도 호환되도록 helper가 정의된 경우에만 assembly 상태를 조회한다.
ASSEMBLY_PHASE_STATUS=not_applicable
if declare -F sample_assembly_phase_label_status >/dev/null; then
  ASSEMBLY_PHASE_STATUS=$(sample_assembly_phase_label_status "$S")
fi
DECODE_REF=$(sample_decode_reference "$S")
BAM=$(sample_haplotagged_alignment "$S")
require_file "$BAM"
have_index "$BAM" || die "missing haplotagged BAM index: $BAM"
REQUIRE_CHAIN=0
[[ "$SEQUENCE_VARIANT_POLICY" == snv_indel && "$SEQUENCE_SOURCE_RESOLVED" != assembly_direct ]] && REQUIRE_CHAIN=1
require_hap_contig_assets "$ROOT/$S/individual_fasta" "$CHROMS" "$REQUIRE_CHAIN"
BWA_INDEX_ARGS=()
ASSEMBLY_PHASE_FLAG=
if [[ "$SEQUENCE_SOURCE_RESOLVED" == assembly_direct && "$ASSEMBLY_PHASE_STATUS" == validated ]]; then
  ASSEMBLY_PHASE_FLAG=--assembly-phase-labels-validated
fi
if [[ "$STAGE11_BWA_INDEX_MODE" == shared ]]; then
  BWA_INDEX_ARGS=(--bwa-mem2-index-root "$STAGE22_BWA_INDEX_ROOT")
fi

log "10 personal realignment: sample=$S platform=$PLATFORM memory_profile=$MFM_MEMORY_PROFILE_RESOLVED hap_read_mode=$HAP_READ_MODE illumina_aligner=$STAGE22_ILLUMINA_ALIGNER duplicate_mode=$STAGE22_ILLUMINA_DUPLICATE_MODE bwa_index_mode=$STAGE11_BWA_INDEX_MODE route_jobs=$WGS_ROUTE_JOBS hap_jobs=$WGS_HAP_JOBS align_threads=$WGS_REALIGN_THREADS sort_threads=$WGS_SORT_THREADS sort_mem=$WGS_SORT_MEM source=$BAM"
log "중간 CRAM 무손실 압축 level=$MFM_WGS_CRAM_LEVEL"
run_cli create-personal-realigned-alignments \
  --cohort "$COHORT" \
  --output-root "$ROOT" \
  --realigned-subdir "$REALIGNED_SUBDIR" \
  --sample "$S" \
  --chroms "$CHROMS" \
  --sequence-variant-policy "$SEQUENCE_VARIANT_POLICY" \
  --sequence-source "$SEQUENCE_SOURCE_RESOLVED" \
  $ASSEMBLY_PHASE_FLAG \
  --dataset-role "$(sample_dataset_role "$S")" \
  --alignment "$BAM" \
  --reference-fasta "$DECODE_REF" \
  --platform "$PLATFORM" \
  --preset auto \
  --hap-read-mode "$HAP_READ_MODE" \
  --format "$WGS_REALIGN_FORMAT" \
  --threads "$WGS_REALIGN_THREADS" \
  --sort-threads "$WGS_SORT_THREADS" \
  --sort-mem "$WGS_SORT_MEM" \
  --jobs "$WGS_REALIGN_JOBS" \
  --hap-jobs "$WGS_HAP_JOBS" \
  --route-jobs "$WGS_ROUTE_JOBS" \
  --min-mapped-frac "$WGS_MIN_MAPPED_FRAC" \
  --ratio-qc-policy "$QC_RATIO_POLICY" \
  --samtools-bin "$SAMTOOLS_BIN" \
  --minimap2-bin "$MINIMAP2_BIN" \
  --bwa-mem2-bin "$BWA_MEM2_BIN" \
  --illumina-aligner "$STAGE22_ILLUMINA_ALIGNER" \
  "${BWA_INDEX_ARGS[@]}" \
  --illumina-duplicate-mode "$STAGE22_ILLUMINA_DUPLICATE_MODE"
log "10 완료: $S"
