#!/usr/bin/env bash
# 일반 실행과 실물 직렬/병렬 비교를 같은 stage barrier로 수행함.
# 두 HAP가 같은 총 CPU 예산을 나눔. 첫 HAP의 전체 풀 선점을 방지함.
pipeline_rna_index_group() {
  local stage=$1 callback=$2 jobs=$3 budget=${MFM_SAMPLE_CPU_BUDGET:-${MFM_CPU_BUDGET:-1}}
  if [[ -n ${MFM_CPU_POOL_STAGE_MAX:-} ]] && (( MFM_CPU_POOL_STAGE_MAX < budget )); then budget=$MFM_CPU_POOL_STAGE_MAX; fi
  local cap=$((budget / jobs))
  (( cap >= 1 )) || cap=1
  # 함수 호출 한정 환경으로 전달해 전역 상태·errexit를 유지하고 자동 원복함.
  MFM_CPU_POOL_STAGE_MAX=$cap pipeline_run_hap_group "$stage" "$callback" "$jobs"
}

pipeline_run_rna_index_stage() {
  local stage=$1 callback=$2 saved_stop=${PIPELINE_STOP_AFTER:-} started serial_seconds parallel_seconds
  if [[ ${RNA_STAR_INDEX_BENCHMARK:-0} != 1 ]] || (( 10#$stage < 10#$PIPELINE_START_AT )); then
    pipeline_rna_index_group "$stage" "$callback" "$RNA_STAR_INDEX_HAP_JOBS"
    return
  fi
  PIPELINE_STOP_AFTER=
  local report=$F/qc/rna_index_benchmark.json
  local helper=$PIPELINE_DIR/src/integrate_wgs_rna/rna_index_benchmark.py
  RNA_INDEX_FORCE_ARGS=(--force)
  "$PYTHON" "$helper" prepare --work-dir "$W/rna" --report "$report" --seconds 0
  started=$SECONDS
  pipeline_rna_index_group "$stage" "$callback" 1
  serial_seconds=$((SECONDS-started))
  "$PYTHON" "$helper" snapshot --work-dir "$W/rna" --report "$report" --seconds "$serial_seconds"
  RNA_INDEX_FORCE_ARGS=(--force)
  "$PYTHON" "$helper" prepare --work-dir "$W/rna" --report "$report" --seconds 0
  started=$SECONDS
  pipeline_rna_index_group "$stage" "$callback" 2
  parallel_seconds=$((SECONDS-started))
  RNA_INDEX_FORCE_ARGS=()
  "$PYTHON" "$helper" compare --work-dir "$W/rna" --report "$report" --seconds "$parallel_seconds"
  PIPELINE_STOP_AFTER=$saved_stop
  pipeline_finish_if_stop "$stage"
}
