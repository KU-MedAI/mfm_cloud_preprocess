#!/usr/bin/env bash
# 구번호 Stage 17 직접 실행기. 현재 RNA runner의 CNV 단계는 Stage 22 `22_cnv_filter.py`다.
# 17 — 개인 HAP WGS depth로 gene별 CNV 의심 구간을 표시한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
require_manifest_sample_contract "$S"
if [[ "$RNA_CNV_FILTER" != 1 ]]; then
  log "17 skip: RNA_CNV_FILTER=0 (CNV 필터 없이 진행)"
  exit 0
fi

F=$(v2_final "$S")
W=$(v2_work "$S")
MTAG=$(sample_rna_evidence_label "$S")
D1=$F/rna_personal/$MTAG/HAP1
D2=$F/rna_personal/$MTAG/HAP2
CNV_RD=$F/rna_personal/$MTAG/gene_rd.tsv
W1=$(sample_hap_alignment "$S" 1)
W2=$(sample_hap_alignment "$S" 2)
R1=$(sample_hap_target "$S" 1)
R2=$(sample_hap_target "$S" 2)
mkdir -p "$W/logs"

[[ -s "$W1" ]] || die "$S CNV용 HAP1 WGS 정렬 없음: $W1"
[[ -s "$W2" ]] || die "$S CNV용 HAP2 WGS 정렬 없음: $W2"
[[ -s "$R1" ]] || die "$S CNV용 HAP1 target.fa 없음: $R1"
[[ -s "$R2" ]] || die "$S CNV용 HAP2 target.fa 없음: $R2"
require_file "$D1/exons.bygene.tsv"
require_file "$D2/exons.bygene.tsv"

"$PYTHON" ./17_cnv_filter.py \
  --exons "$D1/exons.bygene.tsv" --bam "$W1" --reference "$R1" \
  --exons "$D2/exons.bygene.tsv" --bam "$W2" --reference "$R2" \
  --min-rd "${RNA_CNV_MIN_RD:-0.5}" --max-rd "${RNA_CNV_MAX_RD:-1.5}" \
  --samtools "$SAMTOOLS" --out "$CNV_RD" 2>&1 | tee -a "$W/logs/17.log"
[[ -s "$CNV_RD" ]] || die "$S CNV rd 파일 없음: $CNV_RD"
log "17 완료: $S"
