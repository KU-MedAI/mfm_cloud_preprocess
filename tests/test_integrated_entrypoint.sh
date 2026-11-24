#!/usr/bin/env bash
set -euo pipefail
TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "$TEST_DIR/.." && pwd)"

export HAP_READ_MODE=hap_specific
source "$PIPELINE_DIR/config_v2.sh"
sample_final_root() { printf '/tmp/final\n'; }

expected=/tmp/final/S1/realigned/hap_specific/HAP1/S1.personal.HAP1.cram
actual=$(sample_hap_alignment S1 1)
[[ "$actual" == "$expected" ]] || { echo "[FAIL] mode path: $actual" >&2; exit 1; }
[[ "$COMPACT_LAYOUT_NAME" == "full_shared_hap_4bit_acgtnp_depth__hap_specific" ]]   || { echo "[FAIL] compact layout: $COMPACT_LAYOUT_NAME" >&2; exit 1; }

grep -F -- './run_pipeline.sh sample <sample>' "$PIPELINE_DIR/run_pipeline.sh" >/dev/null
grep -F -- './run_pipeline.sh sample-list <sample_list>' "$PIPELINE_DIR/run_pipeline.sh" >/dev/null
grep -F -- './run_wgs_sample.sh "$S"' "$PIPELINE_DIR/run_sample.sh" >/dev/null
grep -F -- './run_rna_sample.sh "$S"' "$PIPELINE_DIR/run_sample.sh" >/dev/null
grep -F -- 'sample) SCRIPT_REL=run_sample.sh' "$PIPELINE_DIR/runner_local.sh" >/dev/null
grep -F -- 'validate-stage70-trace' "$PIPELINE_DIR/12_compact_zarr.sh" >/dev/null
grep -F -- '--workers "$TRACE_JOBS"' "$PIPELINE_DIR/12_compact_zarr.sh" >/dev/null
grep -F -- 'assert_source_trace_ready' "$PIPELINE_DIR/20_cleanup_intermediates.sh" >/dev/null

echo "[PASS] integrated one-command and HAP-mode path contract"
