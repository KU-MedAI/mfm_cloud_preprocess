#!/usr/bin/env bash
# 11 — HAP1/HAP2 window index 생성.
#
# 기본 window size는 1M(1,048,576 bp)이고 시작점은 deterministic random offset이다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
ROOT=$(sample_final_root "$S")
COHORT=$(sample_cohort "$S")
SEQUENCE_SOURCE_RESOLVED=$(sequence_source_for_sample "$S")
FIRST_CHROM=$(first_chrom "$CHROMS")
require_file "$ROOT/$S/individual_fasta/HAP1/$FIRST_CHROM.fa"
require_file "$ROOT/$S/individual_fasta/HAP2/$FIRST_CHROM.fa"

# Stage 11은 sequence-only window 계산이지만, 파이프라인 계약상 Stage 10의 두 HAP
# 재정렬 산출물·sample/HAP/output 계약을 확인한 뒤에만 실행한다. mapped fraction은
# record_only에서 값만 확인하고 strict에서만 하한을 강제한다. 이 guard가 없으면 START_AT=11
# 또는 남은 오래된 FASTA 때문에 재정렬 없이 window/compact가 만들어질 수 있다.
if sample_hap_realign_required "$S"; then
for h in 1 2; do
  REALIGNED=$(sample_hap_alignment "$S" "$h")
  REALIGN_QC="$(dirname "$REALIGNED")/$S.realign.json"
  require_file "$REALIGNED"
  have_index "$REALIGNED" || die "11 requires Stage 10 alignment index: $REALIGNED"
  require_file "$REALIGN_QC"
  "$PYTHON" - "$REALIGN_QC" "$S" "$h" "$REALIGNED" "$QC_RATIO_POLICY" "$WGS_MIN_MAPPED_FRAC" <<'PY'
import json, pathlib, sys
path, sample, hap, alignment, policy, minimum = sys.argv[1:]
data = json.load(open(path))
assert data.get("sample") == sample, f"sample mismatch: {data.get('sample')} != {sample}"
assert data.get("haplotype") == f"HAP{hap}", f"hap mismatch: {data.get('haplotype')}"
assert pathlib.Path(data.get("output", "")).resolve() == pathlib.Path(alignment).resolve(), "output path mismatch"
fraction = float(data.get("routed", {}).get("mapped_frac", -1))
assert 0.0 <= fraction <= 1.0, f"invalid mapped fraction: {fraction}"
if policy == "strict":
    assert fraction >= float(minimum), f"mapped fraction {fraction} < {minimum}"
PY
done
elif sample_pileup_feature "$S"; then
  export PYTHONPATH="$PWD:$PWD/src:$PWD/vendor/MFM-DNARNAFM:${PYTHONPATH:-}"
  "$PYTHON" -m integrate_wgs_rna.pileup_features.export --check-inputs \
    --root "$ROOT/$S" --alignment "$(sample_haplotagged_alignment "$S")" \
    --reference "$(sample_sequence_reference "$S")" --alignment-reference "$(sample_decode_reference "$S")" \
    --read-type "$(sample_read_type "$S")" --contigs "$(chroms_space "$CHROMS")"
fi

args=(
  create-window-index
  --cohort "$COHORT"
  --output-root "$ROOT"
  --sample "$S"
  --chroms "$CHROMS"
  --sequence-variant-policy "$SEQUENCE_VARIANT_POLICY"
  --sequence-source "$SEQUENCE_SOURCE_RESOLVED"
  --window-size "$WINDOW_SIZE"
  --seed "$WINDOW_SEED"
  --max-blacklist-frac "$WINDOW_MAX_BLACKLIST_FRAC"
  --jobs "$WINDOW_JOBS"
)
if [[ "$WINDOW_DROP_N" == 0 ]]; then
  args+=(--no-drop-n)
fi
if [[ "$SEQUENCE_SOURCE_RESOLVED" == "assembly_direct" ]]; then
  args+=(--disable-reference-blacklist)
  warn "11 assembly_direct: GRCh38 multi-chain 투영 미구현으로 reference blacklist를 적용하지 않음"
elif [[ -n "$BLACKLIST_BED" ]]; then
  args+=(--blacklist-bed "$BLACKLIST_BED")
fi
if [[ "$FIXED_WINDOW_START" == 1 ]]; then
  args+=(--fixed-window-start)
elif [[ -n "$OFFSET_UNIT" ]]; then
  args+=(--offset-unit "$OFFSET_UNIT")
fi

log "11 window index 생성: sample=$S window=$WINDOW_SIZE seed=$WINDOW_SEED"
run_cli "${args[@]}"
log "11 완료: $S"
