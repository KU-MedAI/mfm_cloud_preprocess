#!/usr/bin/env bash
# 직접 실행도 기본 공유 CPU 풀을 사용하고 OMP의1thread 환경에 축소되지 않는지 검사한다.
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_TMP=$(mktemp -d)
trap 'rm -rf "$TEST_TMP"' EXIT
mkdir -p "$TEST_TMP/lib" "$TEST_TMP/scripts" "$TEST_TMP/status" "$TEST_TMP/final"
cp "$ROOT_DIR/run_wgs_sample.sh" "$TEST_TMP/"
cp "$ROOT_DIR/lib/state.sh" "$ROOT_DIR/lib/stage_runner.sh" "$TEST_TMP/lib/"
cp "$ROOT_DIR/scripts/cpu_pool.py" "$TEST_TMP/scripts/"
cat > "$TEST_TMP/config_v2.sh" <<'EOF'
die() { echo "$*" >&2; exit 1; }
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
cat > "$TEST_TMP/01_prepare_analysis_bam.sh" <<'EOF'
#!/usr/bin/env bash
printf 'actual_grant=%s\n' "$MFM_SAMPLE_CPU_BUDGET"
EOF
chmod +x "$TEST_TMP/01_prepare_analysis_bam.sh"
export TEST_STATUS="$TEST_TMP/status" TEST_FINAL="$TEST_TMP/final" OUTPUT_ROOT="$TEST_TMP/output"
export WGS_START_AT=01 WGS_STOP_AFTER=01 RUN_STATUS_DIR="$TEST_TMP/status" OMP_NUM_THREADS=1
unset MFM_CPU_BUDGET MFM_CPU_POOL_DIR MFM_CPU_POOL_AUDIT MFM_CPU_POOL_PROFILE MFM_CPU_UTILIZATION
unset MFM_CPU_POOL_LAYOUT SLURM_JOB_ID SLURM_CPUS_PER_TASK RUN_ID
(
  cd "$TEST_TMP"
  python3 - <<'PY'
import os
allowed = set(sorted(os.sched_getaffinity(0))[:4])
os.sched_setaffinity(0, allowed)
os.execvp('bash', ['bash', './run_wgs_sample.sh', 'SAMPLE'])
PY
)
python3 - "$OUTPUT_ROOT/_resource/default/cpu_pool/usage.jsonl" <<'PY'
import json,sys
row=json.loads(open(sys.argv[1],encoding='utf-8').read())
assert row['granted_cpus']==row['total_cpus']
assert 1 <= row['total_cpus'] <= 4
assert row['minimum_cpus']==min(4,row['total_cpus'])
PY
echo '[PASS] 직접 runner 기본 공유 풀/전체 affinity/OMP 환경 독립'
