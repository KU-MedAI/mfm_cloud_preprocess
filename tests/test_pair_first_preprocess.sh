#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT="$ROOT_DIR/share_package/scripts/05_preprocess_phase.sh"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

ASSIGN="$TMP_DIR/phase_assignments.tsv"
OUT_ROOT="$TMP_DIR/out"
TRACE="$TMP_DIR/trace.txt"
export TRACE

printf '%s\n' \
  $'sample_id\tphase\tselected_platform\trecommended_training_use' \
  $'PAIR01\t01\tIllumina_30x\tpost_training_wgs_rna_pair' \
  $'WGS01\t01\tIllumina_30x\tpre_training_wgs_only' \
  $'PAIR02\t02\tPacBio_HiFi\tpost_training_wgs_rna_pair' \
  $'WGS02\t02\tPacBio_HiFi\tpre_training_wgs_only' > "$ASSIGN"

MODE=pair_first \
ASSIGN_TSV="$ASSIGN" \
OUT_ROOT="$OUT_ROOT" \
RUN_PREPROCESS=1 \
PAIR_WGS_CMD='printf "%s:%s\\n" "$PRIORITY_STAGE" "$PHASE" >> "$TRACE"' \
PAIR_RNA_CMD='printf "%s:%s\\n" "$PRIORITY_STAGE" "$PHASE" >> "$TRACE"' \
WGS_ONLY_CMD='printf "%s:%s\\n" "$PRIORITY_STAGE" "$PHASE" >> "$TRACE"' \
bash "$SCRIPT" all

expected=$'pair_wgs:01\npair_wgs:02\npair_rna:01\npair_rna:02\nwgs_only:01\nwgs_only:02'
actual="$(cat "$TRACE")"
[[ "$actual" == "$expected" ]] || {
  echo "[FAIL] pair-first order mismatch" >&2
  printf 'expected:\n%s\nactual:\n%s\n' "$expected" "$actual" >&2
  exit 1
}
[[ -s "$OUT_ROOT/phase_01/execution_order.pair_first.platform_all.tsv" ]]
[[ -s "$OUT_ROOT/phase_02/execution_order.pair_first.platform_all.tsv" ]]
echo "[PASS] pair-first execution order"
