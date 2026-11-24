#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

fail() { echo "FAIL: $*" >&2; exit 1; }

env -i PATH="$PATH" MFM_CPU_UTILIZATION=bounded MFM_SAMPLE_CPU_BUDGET=36 MFM_HAP_PARALLEL_ALLOWED=0 bash -c '
  source ./config_v2.sh
  [[ "$DENOVO_THREADS" == 36 ]]
  [[ "$ASSEMBLY_ANCHOR_THREADS" == 36 ]]
  [[ "$WGS_HAP_JOBS" == 1 ]]
  [[ "$WGS_REALIGN_THREADS" == 32 ]]
  [[ "$WGS_SORT_THREADS" == 4 ]]
  [[ "$COMPACT_JOBS" == 6 ]]
  [[ "$TRACE_JOBS" == 32 ]]
  [[ "$MFM_BLOSC_THREADS" == 6 ]]
  [[ "$MFM_TRACE_BLOSC_THREADS" == 1 ]]
' || fail "WGS 직렬 HAP과 compact 자동 CPU 예산"

env -i PATH="$PATH" MFM_CPU_UTILIZATION=bounded MFM_SAMPLE_CPU_BUDGET=36 bash -c '
  source ./config_v2.sh
  [[ "$WGS_HAP_JOBS" == 2 ]]
  [[ "$WGS_REALIGN_THREADS" == 14 ]]
  [[ "$WGS_SORT_THREADS" == 4 ]]
  [[ $(((WGS_REALIGN_THREADS + WGS_SORT_THREADS) * WGS_HAP_JOBS)) == 36 ]]
' || fail "WGS HAP 자동 병렬 합산 CPU 예산"

env -i PATH="$PATH" MFM_CPU_UTILIZATION=bounded MFM_SAMPLE_CPU_BUDGET=36 WGS_HAP_JOBS=1 bash -c '
  source ./config_v2.sh
  [[ "$WGS_HAP_JOBS" == 1 ]]
  [[ $(((WGS_REALIGN_THREADS + WGS_SORT_THREADS) * WGS_HAP_JOBS)) == 36 ]]
' || fail "명시적 직렬 HAP 설정"

env -i PATH="$PATH" MFM_CPU_UTILIZATION=bounded MFM_SAMPLE_CPU_BUDGET=6 bash -c '
  source ./config_v2.sh
  [[ "$STAR_ALIGN_THREADS" == 6 ]]
  [[ "$HAP_VARIANT_THREADS" == 6 ]]
  [[ "$STAGE04_ALIGN_THREADS" == 6 ]]
  [[ "$WGS_HAP_JOBS" == 1 ]]
  [[ "$WGS_REALIGN_THREADS" == 5 ]]
  [[ "$WGS_SORT_THREADS" == 1 ]]
  [[ "$COMPACT_JOBS" == 1 ]]
  [[ "$TRACE_JOBS" == 2 ]]
  [[ "$MFM_BLOSC_THREADS" == 6 ]]
  [[ "$MFM_TRACE_BLOSC_THREADS" == 3 ]]
' || fail "RNA/고동시성 sample thread cap"

env -i PATH="$PATH" MFM_CPU_UTILIZATION=bounded bash -c '
  source ./config_v2.sh
  [[ "$WGS_HAP_JOBS" == 2 && "$WGS_ROUTE_JOBS" == 8 && "$WGS_REALIGN_THREADS" == 6 ]]
  export MFM_SAMPLE_CPU_BUDGET=18
  source ./config_v2.sh
  [[ "$WGS_ROUTE_JOBS" == 16 ]]
  [[ "$WGS_HAP_JOBS" == 2 ]]
  [[ "$WGS_REALIGN_THREADS" == 7 ]]
  [[ "$WGS_SORT_THREADS" == 2 ]]
  [[ "$COMPACT_JOBS" == 3 ]]
  [[ "$TRACE_JOBS" == 14 ]]
  [[ "$MFM_BLOSC_THREADS" == 6 ]]
  [[ "$MFM_TRACE_BLOSC_THREADS" == 1 ]]
  [[ $(((WGS_REALIGN_THREADS + WGS_SORT_THREADS) * WGS_HAP_JOBS)) == 18 ]]
' || fail "다중 entrypoint source 뒤 child CPU 예산 재계산"

# 소수 CPU 배정과 실제 염색체 병렬 계산: 약수 조건으로 1개로 축소하지 않는다.
for budget in 6 11 18 36 61 72; do
  env -i PATH="$PATH" MFM_SAMPLE_CPU_BUDGET="$budget" MFM_CPU_UTILIZATION=maximum bash -c '
    source ./config_v2.sh
    expected=$MFM_SAMPLE_CPU_BUDGET
    (( expected > 22 )) && expected=22
    [[ "$COMPACT_JOBS" == "$expected" ]]
    [[ "$MFM_BLOSC_THREADS" == 1 ]]
    [[ "$MFM_TRACE_BLOSC_THREADS" == 1 ]]
  ' || fail "최대 CPU 모드 소수 배정 $budget"
done

env -i PATH="$PATH" MFM_SAMPLE_CPU_BUDGET=61 MFM_CPU_UTILIZATION=maximum COMPACT_JOBS=3 bash -c '
  source ./config_v2.sh
  [[ "$COMPACT_JOBS" == 3 ]]
' || fail "명시적 compact 작업 수 유지"

grep -F 'LOCAL_JOBS=${LOCAL_JOBS:-auto}' runner_local.sh >/dev/null || fail "동적 기본 local jobs"
grep -F 'MFM_CPU_BUDGET를 현재 CPU affinity로 제한' runner_local.sh >/dev/null || fail "CPU affinity oversubscription 차단"
# HAP 동시수와 RAM·디스크 조건은 test_hap_gpu_pipeline/test_align_resources에서 실제 계획으로 검사한다.
grep -F 'MFM_CPU_POOL_DIR=${MFM_CPU_POOL_DIR:-$RUN_STATUS_DIR/cpu_pool}' runner_local.sh >/dev/null || fail "공용 단계 CPU pool"
grep -F 'sample_cpu_budget=dynamic-stage-lease' runner_local.sh >/dev/null || fail "단계별 동적 CPU 임대"
grep -F 'ALIGN_SAMPLE_JOBS=${WGS_ALIGN_SAMPLE_JOBS:-4}' run_wgs_phased.sh >/dev/null || fail "Stage 11 단계형 sample jobs"
grep -F 'COMPACT_SAMPLE_JOBS=${WGS_COMPACT_SAMPLE_JOBS:-2}' run_wgs_phased.sh >/dev/null || fail "Stage 12 단계형 sample jobs"
grep -F 'wgs-phased-list)' run_pipeline.sh >/dev/null || fail "단계형 WGS 진입점"
! grep -F 'LOCAL_JOBS=${LOCAL_JOBS:-32}' run_pipeline.sh >/dev/null || fail "과거 32 jobs override 잔존"
echo "PASS: 72코어 자동 CPU 예산"
