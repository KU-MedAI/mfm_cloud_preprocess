#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP=$(mktemp -d "$ROOT_DIR/tests/.tmp_common_runner.XXXXXX")
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/lib" "$TMP/status"
cp "$ROOT_DIR/run_wgs_sample.sh" "$TMP/run_wgs_sample.sh"
cp "$ROOT_DIR/lib/state.sh" "$ROOT_DIR/lib/stage_runner.sh" "$TMP/lib/"
chmod +x "$TMP/run_wgs_sample.sh"

cat > "$TMP/config_v2.sh" <<'EOF'
die() { echo "[FATAL] $*" >&2; exit 1; }
log() { echo "$*" >&2; }
require_manifest_sample_contract() { :; }
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

for stage in 01_prepare_analysis_bam 02_preflight 03_prepare_vcf \
             04_call_longread_clair3 05_clair3_qc 06_haplotag \
             07_denovo_assembly 08_anchor_assembly 09_personal_fasta \
             10_realign_wgs 11_windows 12_compact_zarr \
             23_validate_sample_package 24_cleanup_intermediates; do
  cat > "$TMP/$stage.sh" <<'EOF'
#!/usr/bin/env bash
name=$(basename "$0")
printf '%s:%s\n' "$name" "$*" >> "$TEST_TRACE"
[[ "${TEST_FAIL_STAGE:-}" != "${name%%_*}" ]] || exit 23
EOF
  chmod +x "$TMP/$stage.sh"
done

export TEST_STATUS="$TMP/status"
export TEST_FINAL="$TMP/final"
export TEST_TRACE="$TMP/trace"
(
  cd "$TMP"
  CLEAN_INTERMEDIATES_ON_SUCCESS=1 ./run_wgs_sample.sh SAMPLE
)
expected=$'01_prepare_analysis_bam.sh:SAMPLE\n02_preflight.sh:SAMPLE\n03_prepare_vcf.sh:SAMPLE\n06_haplotag.sh:SAMPLE\n09_personal_fasta.sh:SAMPLE\n10_realign_wgs.sh:SAMPLE\n11_windows.sh:SAMPLE\n12_compact_zarr.sh:SAMPLE\n23_validate_sample_package.sh:SAMPLE wgs\n24_cleanup_intermediates.sh:SAMPLE wgs'
[[ "$(cat "$TEST_TRACE")" == "$expected" ]] || {
  echo "[FAIL] WGS 기본 단계 순서/cleanup mode 불일치" >&2
  exit 1
}
grep -F $'SAMPLE\t04\tskipped\tshort-read는 Clair3 long-read self-calling 대상 아님' \
  "$TEST_STATUS/SAMPLE.stages.tsv" >/dev/null
grep -F $'SAMPLE\t07\tskipped\tSEQUENCE_SOURCE가 denovo인 long-read만 실행' \
  "$TEST_STATUS/SAMPLE.stages.tsv" >/dev/null
RUNTIME="$TEST_FINAL/SAMPLE/metadata/stage_runtime.tsv"
[[ "$(head -1 "$RUNTIME")" == $'run_id\tsample_id\tbranch\tstage\tstage_name\tstatus\tstarted_at\tfinished_at\telapsed_seconds\telapsed_minutes\texit_code\treason' ]]
awk -F'\t' 'NR>1 && $4=="01" {found=($6=="success" && $9 ~ /^[0-9]+$/ && $10 ~ /^[0-9]+\.[0-9]{3}$/)} END {exit !found}' "$RUNTIME"
awk -F'\t' 'NR>1 && $4=="04" {found=($6=="skipped" && $9=="0" && $11=="0")} END {exit !found}' "$RUNTIME"

: > "$TEST_TRACE"
(
  cd "$TMP"
  WGS_START_AT=11 WGS_STOP_AFTER=12 ./run_wgs_sample.sh SAMPLE
)
expected=$'11_windows.sh:SAMPLE\n12_compact_zarr.sh:SAMPLE'
[[ "$(cat "$TEST_TRACE")" == "$expected" ]] || {
  echo "[FAIL] WGS start/stop 순서 불일치" >&2
  exit 1
}
awk -F'\t' 'NR==2 {exit !($2=="success" && $3=="12" && $4=="12" && $5=="0")}' \
  "$TEST_STATUS/SAMPLE.state.tsv"

: > "$TEST_TRACE"
if (
  cd "$TMP"
  TEST_FAIL_STAGE=10 WGS_START_AT=09 ./run_wgs_sample.sh SAMPLE >/dev/null 2>&1
); then
  echo "[FAIL] stage 실패가 성공으로 처리됨" >&2
  exit 1
fi
awk -F'\t' 'NR==2 {exit !($2=="failed" && $3=="10" && $4=="09" && $5=="23")}' \
  "$TEST_STATUS/SAMPLE.state.tsv"
[[ "$(tail -1 "$TEST_TRACE")" == '10_realign_wgs.sh:SAMPLE' ]] || {
  echo "[FAIL] 실패 뒤 다음 stage가 실행됨" >&2
  exit 1
}

awk -F'\t' 'NR>1 && $4=="10" && $6=="failed" && $11=="23" {found=1} END {exit !found}' "$RUNTIME"
[[ "$(grep -c '^run_id' "$RUNTIME")" == 1 ]]

PROVENANCE="$TMP/provenance.jsonl"
RUNTIME_FROM_PROVENANCE="$TMP/stage_runtime_from_provenance.tsv"
cat > "$PROVENANCE" <<'EOF'
{"event":"start","timestamp":"2026-09-29T00:00:00+00:00","sample_id":"SAMPLE","branch":"wgs","stage":"01","stage_name":"01_prepare_analysis_bam.sh","run_id":"test"}
{"event":"finish","timestamp":"2026-09-29T00:01:30+00:00","sample_id":"SAMPLE","branch":"wgs","stage":"01","stage_name":"01_prepare_analysis_bam.sh","run_id":"test","status":"success","exit_code":0,"elapsed_seconds":90}
{"event":"skip","timestamp":"2026-09-29T00:01:31+00:00","sample_id":"SAMPLE","branch":"wgs","stage":"04","stage_name":"conditional_skip","run_id":"test","status":"skipped","reason":"public VCF"}
EOF
python3 "$ROOT_DIR/scripts/summarize_stage_runtime.py" --input "$PROVENANCE" --output "$RUNTIME_FROM_PROVENANCE"
awk -F'\t' 'NR==2 {exit !($1=="test" && $4=="01" && $6=="success" && $9=="90" && $10=="1.500")}' "$RUNTIME_FROM_PROVENANCE"
awk -F'\t' 'NR==3 {exit !($4=="04" && $6=="skipped" && $9=="0" && $11=="0")}' "$RUNTIME_FROM_PROVENANCE"

echo "[PASS] common WGS stage/state/start-stop/failure/runtime contract"
