#!/usr/bin/env bash
# 24 — Stage 23 PASS sample만 기존의 보수적 cleanup 구현으로 전달한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample> [wgs|rna]}
MODE=${2:-wgs}
ROOT="$(sample_final_root "$S")/$S"
REPORT=$ROOT/qc/final_package_validation.json
[[ -s "$REPORT" ]] || die "24 cleanup 거부: Stage 23 보고서 없음: $REPORT"
"$PYTHON" - "$REPORT" "$S" "$MODE" <<'PY'
import json,sys
p,sample,mode=sys.argv[1:]
d=json.load(open(p))
assert d.get("schema")=="mfm_final_package_validation_v1", "validation schema mismatch"
assert d.get("sample")==sample and d.get("mode")==mode, "validation sample/mode mismatch"
assert d.get("status")=="pass" and not d.get("failed_checks"), "validation did not pass"
PY
if [[ "${MFM_MODEL_ONLY_CLEANUP:-0}" == 1 ]] || [[ "$(sample_output_mode "$S")" != full ]]; then
  ARGS=()
  [[ "${CLEAN_INTERMEDIATES_DRY_RUN:-0}" == 1 ]] && ARGS+=(--dry-run)
  export PYTHONPATH="$PWD:$PWD/src:$PWD/vendor/MFM-DNARNAFM:${PYTHONPATH:-}"
  "$PYTHON" -m handoff.cleanup --root "$ROOT" --sample "$S" --mode "$MODE" \
    --work "$(sample_work_root "$S")" --output-root "$OUTPUT_ROOT" "${ARGS[@]}"
else
  ./20_cleanup_intermediates.sh "$S" "$MODE"
fi
