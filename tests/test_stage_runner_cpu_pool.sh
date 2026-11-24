#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/lib" "$TMP/scripts" "$TMP/status" "$TMP/final"
cp "$ROOT_DIR/run_wgs_sample.sh" "$ROOT_DIR/lib/state.sh" "$ROOT_DIR/lib/stage_runner.sh" "$TMP/"
cp "$ROOT_DIR/lib/state.sh" "$ROOT_DIR/lib/stage_runner.sh" "$TMP/lib/"
if [[ -f "$ROOT_DIR/lib/fixed_cpu_entry.sh" ]]; then
  cp "$ROOT_DIR/lib/fixed_cpu_entry.sh" "$TMP/lib/"
fi
cp "$ROOT_DIR/scripts/cpu_pool.py" "$TMP/scripts/"
chmod +x "$TMP/run_wgs_sample.sh"
cat > "$TMP/config_v2.sh" <<'EOF'
die() { echo "$*" >&2; exit 1; }
log() { echo "$*" >&2; }
require_manifest_sample_contract() { :; }
apply_sample_payload() { :; }
sample_wgs_reads_required() { return 0; }
normalize_platform() { echo illumina; }
sample_platform() { echo illumina; }
sequence_source_for_sample() { echo vcf_consensus; }
longread_sequence_vcf_source_for_sample() { echo panel; }
longread_haplotag_vcf_source_for_sample() { echo panel; }
sample_project_id() { echo test; }
sample_status_root() { echo "$TEST_STATUS"; }
sample_final_root() { echo "$TEST_FINAL"; }
longread_needs_clair3() { return 1; }
is_long_read_sample() { return 1; }
VCF_SOURCE_POLICY=provided
EOF
cat > "$TMP/01_prepare_analysis_bam.sh" <<'EOF'
#!/usr/bin/env bash
printf 'grant=%s\n' "$MFM_SAMPLE_CPU_BUDGET"
EOF
chmod +x "$TMP/01_prepare_analysis_bam.sh"
export TEST_STATUS="$TMP/status" TEST_FINAL="$TMP/final"
export PYTHONPATH="$ROOT_DIR/src:$ROOT_DIR${PYTHONPATH:+:$PYTHONPATH}"
export MFM_CPU_ALLOCATION_MODE=dynamic
export MFM_CPU_BUDGET=6 MFM_CPU_POOL_DIR="$TMP/pool" MFM_CPU_POOL_AUDIT="$TMP/audit.jsonl"
(
  cd "$TMP"
  WGS_START_AT=01 WGS_STOP_AFTER=01 RUN_STATUS_DIR="$TMP/status" bash ./run_wgs_sample.sh SAMPLE
)
grep -F 'grant=6/6' "$TMP/final/SAMPLE/logs/stage_logs/WGS/01."*.log >/dev/null
grep -F 'grant=6' "$TMP/final/SAMPLE/logs/stage_logs/WGS/01."*.log >/dev/null
python3 - "$TMP/audit.jsonl" <<'PY'
import json,sys
row=json.loads(open(sys.argv[1],encoding='utf-8').read())
assert row['sample']=='SAMPLE' and row['stage']=='01' and row['granted_cpus']==6
PY
echo '[PASS] stage runner CPU lease/env/audit integration'
