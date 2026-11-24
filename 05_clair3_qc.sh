#!/usr/bin/env bash
# 05 — long-read Clair3 output QC.
#
# 입력: stage 04에서 생성한 Clair3 VCF + haplotagged BAM.
# 출력: <final_root>/<sample>/qc/clair3_output_qc.{json,tsv}
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
if ! is_long_read_sample "$S"; then
  log "05 skip: short-read는 Clair3 QC 대상 아님 sample=$S platform=$(sample_platform "$S")"
  exit 0
fi
if ! longread_needs_clair3 "$S"; then
  log "05 skip: Clair3를 사용하지 않는 sample=$S source=$(sample_vcf_source "$S")"
  exit 0
fi

ROOT=$(sample_final_root "$S")
COHORT=$(sample_cohort "$S")
PLATFORM=$(normalize_platform "$(sample_platform "$S")")
MODEL_PATH=$(clair3_model_path_for_sample "$S")
[[ -d "$MODEL_PATH" ]] || die "missing Clair3 model path for sample=$S platform=$PLATFORM: $MODEL_PATH"
mkdir -p "$ROOT/$S/qc"

SUMMARY_TSV="$ROOT/$S/qc/clair3_output_qc.tsv"
METADATA_JSON="$ROOT/$S/qc/clair3_output_qc.metadata.json"
log "05 Clair3 QC: sample=$S platform=$PLATFORM model=$MODEL_PATH"
run_cli qc-longread-clair3 \
  --cohort "$COHORT" \
  --output-root "$ROOT" \
  --sample "$S" \
  --platform "$PLATFORM" \
  --model-path "$MODEL_PATH" \
  --reference-fasta "$REF_DECODE" \
  --samtools-bin "$CLAIR3_SAMTOOLS_BIN" \
  --inspect-reads "$CLAIR3_QC_INSPECT_READS" \
  --min-mq "$CLAIR3_MIN_MQ" \
  --min-coverage "$CLAIR3_MIN_COVERAGE" \
  --qual "$CLAIR3_QUAL" \
  --min-hp-fraction "$CLAIR3_QC_MIN_HP_FRACTION" \
  --summary-tsv "$SUMMARY_TSV" \
  --metadata-json "$METADATA_JSON"
log "05 완료: $S summary=$SUMMARY_TSV"
