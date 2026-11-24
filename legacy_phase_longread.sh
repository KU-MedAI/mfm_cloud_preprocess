#!/usr/bin/env bash
# legacy — 리드 기반 phasing (long read 샘플 전용).
#
#   ./legacy_phase_longread.sh <sample>
#
#  새 규칙(preprocessing_rule_260825.md): *"read 마다의 haplotyping 은 long read bam 가지고
#  새로 한다 — 기존 vcf 기반 haplotagging 말고 진짜로 bam 의 read 기반 haplotyping"*
#
#  변이 자체는 **Illumina 패널 VCF** 를 쓴다(short read 가 variant calling 정확도가 높다).
#  바뀌는 것은 **phase 를 어디서 얻는가** 다: 패널의 통계 phasing → **리드가 실제로 잇는 phase**.
#
#  ⚠ 리드 phase 는 **블록 단위**다. 이를 어떻게 쓸지 두 갈래를 **모두 산출**한다(jhleee 지시):
#     anchored.vcf.gz — 패널 골격 + 블록 내부만 리드로 정제(블록 방향은 패널 다수결로 앵커링)
#                       → 염색체 전체 H1/H2 정의 가능
#     readonly.vcf.gz — 리드 phase 만. 블록 간 방향은 임의이므로 block_table.tsv 를 함께 낸다
#                       → 블록 안에서만 유효. 리드가 못 잇는 het 은 제외(임의 배정 금지)
#     어느 쪽을 하류가 쓸지는 PHASE_MODE 로 고른다(config_v2.sh 의 v2_phase_vcf).
#
#  입력  : work/<S>/panel.het.vcf.gz  (02 산출 tag VCF = het + biallelic)
#          long read CRAM/BAM (매니페스트 3열)
#  산출물: v2 work/<S>/phase/{read_phased.vcf.gz, anchored.vcf.gz, anchor_stats.tsv, compare.tsv}
#          → 이후 21(haplotag)·03(개인 서열)이 anchored.vcf.gz 를 쓴다
set -uo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
v2_platform "$S" >/dev/null || die "$S 가 매니페스트에 없다: $V2_MANIFEST"
v2_is_long "$S" || die "$S 는 long read 샘플이 아니다(platform=$(v2_platform "$S")). 20 은 건너뛴다"
WGS=$(v2_wgs "$S"); [[ -s "$WGS" ]] || die "$S long read 파일 없음: $WGS"

W=$(v2_work "$S"); OUT=$W/phase; LOG=$W/logs
mkdir -p "$OUT" "$LOG"
PANEL=$W/panel.het.vcf.gz
[[ -s "$PANEL.tbi" ]] || die "패널 het VCF 없음: $PANEL (03_prepare_vcf.sh 먼저)"
STAMP=$OUT/phase.done.json
[[ -s "$STAMP" ]] && { log "[skip] $S 이미 phasing 완료"; exit 0; }

# CRAM 디코딩 reference — 규칙 1: long read CRAM 은 Vienna 로 압축돼 있다
REF_DECODE=$VIENNA_SRC
log "$S: 리드 기반 phasing (platform=$(v2_platform "$S")) / 염색체 $PHASE_CHR_JOBS 병렬"

# ── 염색체별 whatshap phase ──
#  ⚠ 염색체 단위로 쪼개는 이유: whatshap phase 는 전장 1회 실행이 매우 느리고 실패 시 전부
#    잃는다. 염색체별로 나누면 병렬화되고 재시도가 싸다.
#  --ignore-read-groups: 우리 CRAM 의 @RG SM 이 VCF 샘플명과 다를 수 있다(단일샘플 VCF 라 안전).
one_chr() {
  local c=$1
  local o=$OUT/read_phased.$c.vcf.gz
  [[ -s "$o.tbi" ]] && { echo "  [skip] $c"; return 0; }
  "$BCFTOOLS" view -r "$c" -Oz -o "$OUT/panel.$c.vcf.gz" "$PANEL" 2>/dev/null
  "$BCFTOOLS" index -t "$OUT/panel.$c.vcf.gz"
  "$WHATSHAP" phase --reference "$REF_DECODE" --ignore-read-groups \
      --chromosome "$c" --indels \
      -o "$o.part" "$OUT/panel.$c.vcf.gz" "$WGS" > "$LOG/phase.$c.log" 2>&1 \
    || { echo "  [FAIL] $c ($(tail -1 "$LOG/phase.$c.log"))"; rm -f "$o.part"; return 1; }
  mv "$o.part" "$o"; "$BCFTOOLS" index -t "$o"
  local nb
  nb=$(grep -oE "Found [0-9]+ .*blocks" "$LOG/phase.$c.log" | tail -1)
  echo "  [ok] $c ${nb:-}"
}
export -f one_chr; export OUT LOG PANEL WGS REF_DECODE BCFTOOLS WHATSHAP

parallel_map "$PHASE_CHR_JOBS" one_chr $AUTOSOMES 2>&1 | tee "$LOG/20_phase.log"
grep -q '^\s*\[FAIL' "$LOG/20_phase.log" && die "$S phasing 실패 — $LOG 확인"

# ── 합치기 ──
"$BCFTOOLS" concat -Oz -o "$OUT/read_phased.vcf.gz" \
    $(for c in $AUTOSOMES; do echo "$OUT/read_phased.$c.vcf.gz"; done) 2>/dev/null
"$BCFTOOLS" index -t "$OUT/read_phased.vcf.gz"

# ── 모드 A: 블록을 패널 방향에 앵커링 → 염색체 전체 일관 H1/H2 ──
"$PYTHON" ./anchor_blocks.py --mode anchored \
    --panel "$PANEL" --read "$OUT/read_phased.vcf.gz" \
    --sample "$S" --out "$OUT/anchored.vcf.gz" --stats "$OUT/anchor_stats.tsv" \
    --min-conf "$PHASE_ANCHOR_CONF" --min-het "$PHASE_ANCHOR_MIN_HET" \
    2>&1 | tee "$LOG/anchor.log"
"$BCFTOOLS" index -t "$OUT/anchored.vcf.gz"

# ── 모드 B: 리드 phase 만 (앵커링 없음) + 블록 경계 표 ──
"$PYTHON" ./anchor_blocks.py --mode readonly \
    --panel "$PANEL" --read "$OUT/read_phased.vcf.gz" \
    --sample "$S" --out "$OUT/readonly.vcf.gz" --block-table "$OUT/block_table.tsv" \
    2>&1 | tee "$LOG/readonly.log"
"$BCFTOOLS" index -t "$OUT/readonly.vcf.gz"

# ── 패널 vs 리드 phasing 정확도 비교 (switch error·블록 통계) ──
#  이 값이 "리드 기반 phasing 을 할 가치가 있는가" 의 직접 증거다.
"$WHATSHAP" compare --sample "$S" --names panel,readbased \
    --tsv-pairwise "$OUT/compare.tsv" "$PANEL" "$OUT/read_phased.vcf.gz" \
    > "$LOG/compare.log" 2>&1 || warn "$S compare 실패(비교만 실패, phasing 은 유효)"

n_het=$("$BCFTOOLS" index -n "$PANEL")
n_blk=$(( $(wc -l < "$OUT/anchor_stats.tsv") - 1 ))
n_corr=$(grep -oE "교정한 het: [0-9,]+" "$LOG/anchor.log" | grep -oE "[0-9,]+" | tr -d ,)
cat > "$STAMP" <<JSON
{ "schema": "mfm_v2_phase_v1", "sample": "$S", "platform": "$(v2_platform "$S")",
  "panel_het": $n_het, "read_blocks": $n_blk, "corrected_het": ${n_corr:-0},
  "anchor": { "min_conf": $PHASE_ANCHOR_CONF, "min_het": $PHASE_ANCHOR_MIN_HET },
  "policy": "변이는 Illumina 패널, phase 는 long read. 두 모드를 모두 산출한다",
  "outputs": { "anchored_vcf": "$OUT/anchored.vcf.gz",
               "readonly_vcf": "$OUT/readonly.vcf.gz",
               "block_table": "$OUT/block_table.tsv",
               "anchor_stats": "$OUT/anchor_stats.tsv",
               "compare": "$OUT/compare.tsv" } }
JSON
log "legacy phase 완료: $S / 패널 het $n_het / 리드 블록 $n_blk / 교정 ${n_corr:-0}"
