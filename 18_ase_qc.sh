#!/usr/bin/env bash
# 구번호 Stage 18 직접 실행기. 현재 runner는 Stage 23 `23_ase_qc.py`와 variant-aware QC 계약을 쓴다.
# 18 — ASE 깊이·재현성·CNV·multimapper 민감도 QC를 수행한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
require_manifest_sample_contract "$S"
F=$(v2_final "$S")
W=$(v2_work "$S")
MTAG=$(sample_rna_evidence_label "$S")
AC=$F/rna_personal/$MTAG/allelic_counts.tsv
CNV_RD=$F/rna_personal/$MTAG/gene_rd.tsv
mkdir -p "$W/logs"
require_file "$AC"
require_file "${AC%.tsv}.mmap_adjusted.tsv"

CNV_ARGS=()
if [[ "$RNA_CNV_FILTER" == 1 ]]; then
  require_file "$CNV_RD"
  CNV_ARGS=(--cnv-rd "$CNV_RD")
fi

"$PYTHON" ./18_ase_qc.py --counts "$AC" \
  --mmap-adjusted "${AC%.tsv}.mmap_adjusted.tsv" \
  ${CNV_ARGS[@]+"${CNV_ARGS[@]}"} \
  --min-cnv-rd "${RNA_CNV_MIN_RD:-0.5}" --max-cnv-rd "${RNA_CNV_MAX_RD:-1.5}" \
  --min-total-count "${ASE_MIN_TOTAL_COUNT:-20}" \
  --min-split-count "${ASE_MIN_SPLIT_COUNT:-5}" \
  --max-split-ratio-delta "${ASE_MAX_SPLIT_RATIO_DELTA:-0.25}" \
  --max-fdr "${ASE_MAX_FDR:-0.05}" \
  --min-informative-fraction "${ASE_MIN_INFORMATIVE_FRACTION:-0.01}" \
  --min-eligible-genes "${ASE_MIN_ELIGIBLE_GENES:-100}" \
  --min-split-correlation "${ASE_MIN_SPLIT_CORRELATION:-0.30}" \
  2>&1 | tee -a "$W/logs/18.log"
[[ -s "$(dirname "$AC")/ase_qc.json" ]] || die "$S ASE QC JSON 없음"
[[ -s "$(dirname "$AC")/ase_training_eligibility.tsv" ]] || die "$S ASE 학습 적격성 TSV 없음"
log "18 완료: $S"
