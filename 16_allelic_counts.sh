#!/usr/bin/env bash
# 구번호 Stage 16 직접 실행기. 현재 RNA runner는 Stage 21의 `21_allelic_counts.py`를 호출하며,
# Stage 20에서 검증된 HAP variant를 입력으로 요구한다. 이 파일은 과거 분석 재현용으로 보존한다.
# 16 — HAP1/HAP2의 fragment AS를 비교해 gene별 대립 카운트를 만든다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
require_manifest_sample_contract "$S"
F=$(v2_final "$S")
W=$(v2_work "$S")
MTAG=$(sample_rna_evidence_label "$S")
D1=$F/rna_personal/$MTAG/HAP1
D2=$F/rna_personal/$MTAG/HAP2
AC=$F/rna_personal/$MTAG/allelic_counts.tsv
mkdir -p "$W/logs"
for path in "$D1/$S.rna.personal.HAP1.bam" "$D2/$S.rna.personal.HAP2.bam" \
            "$D1/exons.bygene.tsv" "$D2/exons.bygene.tsv"; do
  require_file "$path"
done

ASE_ARGS=()
[[ "${RNA_ASE_INCLUDE_DUPLICATES:-0}" == 1 ]] && ASE_ARGS+=(--include-duplicates)
"$PYTHON" ./16_allelic_counts.py \
  --bam1 "$D1/$S.rna.personal.HAP1.bam" --bam2 "$D2/$S.rna.personal.HAP2.bam" \
  --exons1 "$D1/exons.bygene.tsv" --exons2 "$D2/exons.bygene.tsv" \
  ${ASE_ARGS[@]+"${ASE_ARGS[@]}"} --out "$AC" 2>&1 | tee -a "$W/logs/16.log"

[[ -s "$AC" ]] || die "$S ASE 결과 파일 없음"
awk 'END { exit (NR > 1 ? 0 : 1) }' "$AC" || die "$S ASE 결과가 header-only입니다"
for sidecar in "${AC%.tsv}.split_half.tsv" "${AC%.tsv}.both_haps.tsv" \
               "${AC%.tsv}.mmap_adjusted.tsv" "${AC%.tsv}.meta.json"; do
  [[ -s "$sidecar" ]] || die "$S ASE sidecar 없음: $sidecar"
done
log "16 완료: $S"
