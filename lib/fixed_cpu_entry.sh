#!/usr/bin/env bash
# 직접 샘플 실행도 내장 worker와 같은 CPU 소유권을 사용함.
pipeline_fixed_cpu_entry() {
  local sample=$1
  shift
  [[ ${MFM_CPU_ALLOCATION_MODE:-dynamic} == fixed_sample ]] || return 0
  if [[ -n ${MFM_FIXED_SAMPLE_TOKEN:-} ]]; then
    [[ "$sample" == "${MFM_FIXED_SAMPLE_ID:-}" ]] || { echo '[FATAL] 상속한 CPU 소유권과 sample이 다름' >&2; return 1; }
    return 0
  fi
  local python=${PYTHON_BIN:-${PYTHON:-python}}
  local root
  root=$(pwd)
  : "${MFM_SAMPLE_CPU_BUDGET:?fixed_sample에는 MFM_SAMPLE_CPU_BUDGET이 필요함}"
  : "${MFM_CPU_BUDGET:?fixed_sample에는 전체 MFM_CPU_BUDGET이 필요함}"
  : "${MFM_CPU_POOL_DIR:?fixed_sample에는 공유 MFM_CPU_POOL_DIR이 필요함}"
  export PYTHONPATH="$root/src:$root${PYTHONPATH:+:$PYTHONPATH}"
  exec "$python" -m integrate_wgs_rna.sample_cpu --sample "$sample" \
    --sample-cpus "$MFM_SAMPLE_CPU_BUDGET" --total-cpus "$MFM_CPU_BUDGET" \
    --pool-dir "$MFM_CPU_POOL_DIR" -- "$@"
}
