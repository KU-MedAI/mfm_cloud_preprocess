#!/usr/bin/env bash
# compact 자동 작업 수의 소수 예산 회귀와 명시적 설정 보존 검증.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
for budget in 6 11 18 36 61 72; do
  env -i PATH="$PATH" MFM_SAMPLE_CPU_BUDGET="$budget" MFM_CPU_UTILIZATION=maximum bash -c '
    source ./config_v2.sh
    expected=$MFM_SAMPLE_CPU_BUDGET
    (( expected > 22 )) && expected=22
    [[ "$COMPACT_JOBS" == "$expected" && "$MFM_BLOSC_THREADS" == 1 && "$MFM_TRACE_BLOSC_THREADS" == 1 ]]
  '
done
env -i PATH="$PATH" MFM_SAMPLE_CPU_BUDGET=61 MFM_CPU_UTILIZATION=maximum COMPACT_JOBS=3 bash -c '
  source ./config_v2.sh
  [[ "$COMPACT_JOBS" == 3 ]]
'
echo 'PASS: 소수 CPU 배정과 명시적 compact 작업 수'
