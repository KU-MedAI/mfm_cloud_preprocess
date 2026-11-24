#!/usr/bin/env bash
# 25 — WGS–RNA full package를 검증하고 Stage 26 cleanup 허가 JSON을 만든다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample> [wgs|rna]}
MODE=${2:-wgs}
apply_sample_payload "$S"
[[ "$MODE" == wgs || "$MODE" == rna ]] || die "mode는 wgs|rna: $MODE"
ROOT="$(sample_final_root "$S")/$S"
WORK=$(sample_work_root "$S")
mkdir -p "$ROOT/logs/stage_logs" "$ROOT/qc" "$ROOT/metadata"
# 검증 시점에 log를 final package로 먼저 복제한다. 실패 sample도 진단 log를 보존한다.
if [[ -d "$WORK/logs" ]]; then
  cp -a "$WORK/logs/." "$ROOT/logs/stage_logs/"
fi
ARGS=()
[[ "${ALLOW_CLEANUP_WITHOUT_TRACE:-0}" == 1 ]] && ARGS+=(--allow-without-trace)
RNA_LABEL=
[[ "$MODE" == rna ]] && RNA_LABEL=$(sample_rna_evidence_label "$S")
"$PYTHON" ./scripts/validate_sample_package.py \
  --sample "$S" --root "$ROOT" --mode "$MODE" --layout "$COMPACT_LAYOUT_NAME" --wgs-payload "${MFM_WGS_PAYLOAD:-full}" --rna-payload "${MFM_RNA_OUTPUT_MODE:-full}" \
  --chroms "$(chroms_space "$CHROMS")" --hap-read-mode "$HAP_READ_MODE" \
  --rna-label "$RNA_LABEL" --trace-policy "$STAGE12_TRACE_POLICY" \
  --output "$ROOT/qc/final_package_validation.json" \
  --completeness-tsv "$ROOT/metadata/completeness.tsv" "${ARGS[@]}"
