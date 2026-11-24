#!/usr/bin/env bash

# 각 단계별 통계는 시간 TSV/provenance와 함께 누적한다. 실패한 관측도 이유를 남긴다.
pipeline_collect_statistics() {
  [[ "${MFM_HANDOFF_METRICS:-0}" == 1 ]] || return 0
  local stage=$1 root bam= ref= personal_stage=08 args=()
  [[ ! -f "$PIPELINE_DIR/09_personal_fasta.sh" ]] || personal_stage=10
  root="$(sample_final_root "$PIPELINE_SAMPLE_ID")/$PIPELINE_SAMPLE_ID"
  if [[ "$PIPELINE_STAGE_LABEL" == WGS && "$stage" == 01 ]]; then
    bam=$(sample_alignment "$PIPELINE_SAMPLE_ID")
    ref=$(sample_alignment_reference "$PIPELINE_SAMPLE_ID")
  elif [[ "$PIPELINE_STAGE_LABEL" == WGS && "$stage" == 06 ]]; then
    bam=$(find "$root/haplotagged" -maxdepth 1 -type f -name '*.bam' -print -quit)
  elif [[ "$PIPELINE_STAGE_LABEL" == WGS && "$stage" == "$personal_stage" ]]; then
    bam=$(sample_hap_alignment "$PIPELINE_SAMPLE_ID" 1)
    ref=$(sample_hap_target "$PIPELINE_SAMPLE_ID" 1)
  fi
  [[ -z "$bam" ]] || args+=(--bam "$bam")
  [[ -z "$ref" ]] || args+=(--reference "$ref")
  PYTHONPATH="$PIPELINE_DIR:${PYTHONPATH:-}" "$PYTHON" -m handoff.metrics --root "$root" \
    --sample "$PIPELINE_SAMPLE_ID" --stage "$stage" --branch "$PIPELINE_STAGE_LABEL" "${args[@]}" || \
    warn "통계 수집 실패: sample=$PIPELINE_SAMPLE_ID stage=$stage; 단계 로그를 확인하세요."
}

# WGS/RNA 번호 단계의 start/stop, 상태 기록, provenance 실행을 공통화한다.

pipeline_provenance_path() {
  printf '%s/%s/metadata/execution_provenance.jsonl\n' "$(sample_final_root "$PIPELINE_SAMPLE_ID")" "$PIPELINE_SAMPLE_ID"
}

# execution_provenance.jsonl은 기계 판독용 이력이다. 이 TSV는 sample별 단계 시간의
# 사람이 읽기 쉬운 요약이며, 재실행은 별도 행으로 보존한다.
pipeline_stage_runtime_path() {
  printf '%s/%s/metadata/stage_runtime.tsv\n' "$(sample_final_root "$PIPELINE_SAMPLE_ID")" "$PIPELINE_SAMPLE_ID"
}

pipeline_runtime_field() {
  printf '%s' "${1:-}" | tr '\t\r\n' '   '
}

pipeline_record_stage_runtime() {
  local stage=$1 stage_name=$2 status=$3 started_at=$4 finished_at=$5
  local elapsed_seconds=$6 exit_code=$7 reason=${8:-}
  local path elapsed_minutes
  path=$(pipeline_stage_runtime_path)
  mkdir -p "$(dirname "$path")"
  elapsed_minutes=$(awk -v seconds="$elapsed_seconds" 'BEGIN {printf "%.3f", seconds / 60}')

  # WGS와 RNA는 한 sample 안에서 직렬이지만, 재시작/외부 후처리와의 경합에도
  # header 중복이나 행 interleave가 없도록 file lock을 사용한다.
  {
    flock -x 9
    if [[ ! -s "$path" ]]; then
      printf 'run_id\tsample_id\tbranch\tstage\tstage_name\tstatus\tstarted_at\tfinished_at\telapsed_seconds\telapsed_minutes\texit_code\treason\n' >&9
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$(pipeline_runtime_field "${RUN_ID:-}")" \
      "$(pipeline_runtime_field "$PIPELINE_SAMPLE_ID")" \
      "$(pipeline_runtime_field "$PIPELINE_STAGE_LABEL")" \
      "$(pipeline_runtime_field "$stage")" \
      "$(pipeline_runtime_field "$stage_name")" \
      "$(pipeline_runtime_field "$status")" \
      "$(pipeline_runtime_field "$started_at")" \
      "$(pipeline_runtime_field "$finished_at")" \
      "$(pipeline_runtime_field "$elapsed_seconds")" \
      "$(pipeline_runtime_field "$elapsed_minutes")" \
      "$(pipeline_runtime_field "$exit_code")" \
      "$(pipeline_runtime_field "$reason")" >&9
    flock -u 9
  } 9>>"$path"
}

pipeline_command_json() {
  "${PYTHON:-python3}" -c 'import json,sys; print(json.dumps(sys.argv[1:], ensure_ascii=False))' "$@"
}

# 최종 Stage 23까지 기다리지 않고 stdout/stderr를 sample package에 기록한다.
# mktemp는 재시도의 원본 로그를 덮어쓰지 않으며 tee의 실패도 성공으로 처리하지 않는다.
pipeline_new_stage_log() {
  local stage=$1 directory
  directory="$(sample_final_root "$PIPELINE_SAMPLE_ID")/$PIPELINE_SAMPLE_ID/logs/stage_logs/$PIPELINE_STAGE_LABEL"
  mkdir -p "$directory"
  mktemp "$directory/$stage.XXXXXXXX.log"
}

pipeline_execute_logged() {
  local path=$1
  shift
  local -a statuses
  "$@" 2>&1 | tee -a "$path"
  statuses=("${PIPESTATUS[@]}")
  [[ "${statuses[0]}" == 0 ]] || return "${statuses[0]}"
  return "${statuses[1]}"
}

# 단계 시작 시 CPU grant를 받고 종료 시 반환한다. 이미 실행 중인 작업을
# 중간에 선점하지 않으며 다음 단계/대기 sample이 반환 코어를 사용할 수 있다.
pipeline_execute_stage_logged() {
  local path=$1 stage=$2 name=$3
  shift 3
  if [[ -f ./scripts/cpu_pool.py ]]; then
    # CLI를 거치지 않은 직접 sample 실행도 같은 run의 공유 풀을 사용한다.
    local available cpu_python=${PYTHON:-python3}
    available=$("$cpu_python" -c 'import os; n=len(os.sched_getaffinity(0)) if hasattr(os,"sched_getaffinity") else (os.cpu_count() or 1); s=os.environ.get("SLURM_CPUS_PER_TASK"); print(min(n,int(s)) if os.environ.get("SLURM_JOB_ID") and s else n)')
    MFM_CPU_BUDGET=${MFM_CPU_BUDGET:-$available}
    (( MFM_CPU_BUDGET <= available )) || MFM_CPU_BUDGET=$available
    MFM_CPU_UTILIZATION=${MFM_CPU_UTILIZATION:-maximum}
    MFM_CPU_POOL_DIR=${MFM_CPU_POOL_DIR:-${OUTPUT_ROOT:-$(sample_final_root "$PIPELINE_SAMPLE_ID")}/_resource/${RUN_ID:-default}/cpu_pool}
    MFM_CPU_POOL_AUDIT=${MFM_CPU_POOL_AUDIT:-$MFM_CPU_POOL_DIR/usage.jsonl}
    if [[ -z "${MFM_CPU_POOL_LAYOUT:-}" ]]; then
      if [[ -f ./09_personal_fasta.sh ]]; then MFM_CPU_POOL_LAYOUT=original; else MFM_CPU_POOL_LAYOUT=legacy; fi
    fi
    export MFM_CPU_BUDGET MFM_CPU_UTILIZATION MFM_CPU_POOL_DIR MFM_CPU_POOL_AUDIT MFM_CPU_POOL_LAYOUT
    local -a pool_args
    pool_args=(--pool-dir "$MFM_CPU_POOL_DIR" --total-cpus "${MFM_CPU_BUDGET:-$(nproc)}"
      --branch "$PIPELINE_STAGE_LABEL" --stage "$stage" --sample "$PIPELINE_SAMPLE_ID")
    [[ -z "${MFM_CPU_POOL_AUDIT:-}" ]] || pool_args+=(--audit "$MFM_CPU_POOL_AUDIT")
    pipeline_execute_logged "$path" "${PYTHON:-python3}" ./scripts/cpu_pool.py \
      "${pool_args[@]}" -- "$@"
  else
    pipeline_execute_logged "$path" "$@"
  fi
}

pipeline_record_failure() {
  local stage=$1 name=$2 command_json=$3 rc=$4 path=$5
  [[ -f ./scripts/record_stage_failure.py ]] || return 0
  "${PYTHON:-python3}" ./scripts/record_stage_failure.py \
    --output "$(sample_final_root "$PIPELINE_SAMPLE_ID")/$PIPELINE_SAMPLE_ID/metadata/stage_failures.jsonl" \
    --sample "$PIPELINE_SAMPLE_ID" --branch "$PIPELINE_STAGE_LABEL" \
    --stage "$stage" --name "$name" --command-json "$command_json" \
    --run-id "${RUN_ID:-}" --exit-code "$rc" --log "$path"
}

pipeline_write_provenance() {
  local event=$1 stage=$2 name=$3 command_json=$4
  shift 4
  # 최소 단위 shell 테스트는 helper 없이 runner만 복사한다. 실제 배포본에는
  # helper가 필수이며 Stage 23이 provenance 존재를 다시 검증한다.
  [[ -f ./scripts/stage_provenance.py ]] || return 0
  local args
  args=("$PYTHON" ./scripts/stage_provenance.py "$event")
  args+=(--output "$(pipeline_provenance_path)" --sample "$PIPELINE_SAMPLE_ID")
  args+=(--branch "$PIPELINE_STAGE_LABEL" --stage "$stage" --name "$name")
  args+=(--run-id "${RUN_ID:-}" --command-json "$command_json")
  # 버전 조회는 외부 process를 여러 번 띄우므로 sample당 첫 start event에만 한다.
  # 이후 event는 정확한 command/시간/status만 추가해 대규모 cohort의 고정 overhead를 줄인다.
  if [[ "$event" == start && "${PIPELINE_TOOL_VERSIONS_RECORDED:-0}" != 1 ]]; then
    args+=(--tool "samtools=${SAMTOOLS_BIN:-${SAMTOOLS:-}}")
    args+=(--tool "bcftools=${BCFTOOLS_BIN:-${BCFTOOLS:-}}")
    args+=(--tool "minimap2=${MINIMAP2_BIN:-${MINIMAP2:-}}")
    args+=(--tool "bwa-mem2=${BWA_MEM2_BIN:-}")
    args+=(--tool "STAR=${STAR_BIN:-}" --tool "Picard=${PICARD_BIN:-}")
    PIPELINE_TOOL_VERSIONS_RECORDED=1
  fi
  "${args[@]}" "$@"
}

pipeline_validate_stage_range() {
  local label=$1 start=$2 stop=$3
  shift 3
  local allowed stage start_ok=0 stop_ok=0
  allowed=("$@")
  for stage in "${allowed[@]}"; do
    [[ "$start" == "$stage" ]] && start_ok=1
    [[ -n "$stop" && "$stop" == "$stage" ]] && stop_ok=1
  done
  (( start_ok == 1 )) || die "${label}_START_AT 오류: $start"
  if [[ -n "$stop" ]]; then
    (( stop_ok == 1 )) || die "${label}_STOP_AFTER 오류: $stop"
    (( 10#$stop >= 10#$start )) || die "${label}_STOP_AFTER가 ${label}_START_AT보다 앞섬: start=$start stop=$stop"
  fi
}

pipeline_finish_if_stop() {
  local stage=$1
  if [[ -n "${PIPELINE_STOP_AFTER:-}" && "$PIPELINE_STOP_AFTER" == "$stage" ]]; then
    trap - EXIT
    pipeline_write_sample_state success 0
    log "${PIPELINE_STAGE_LABEL}_STOP_AFTER=$stage 요청으로 중단"
    exit 0
  fi
}

pipeline_run_stage() {
  local stage=$1 script=$2
  shift 2
  if (( 10#$stage < 10#$PIPELINE_START_AT )); then
    log "── $stage skip: ${PIPELINE_STAGE_LABEL}_START_AT=$PIPELINE_START_AT ──"
    pipeline_record_stage "$stage" skipped "${PIPELINE_STAGE_LABEL}_START_AT=$PIPELINE_START_AT"
    return 0
  fi
  CURRENT_STAGE=$stage
  pipeline_write_sample_state running 0
  pipeline_record_stage "$stage" running
  log "── $stage $script ──"
  local started started_at finished finished_at command_json rc=0 status stage_log reason=""
  started=$(date +%s)
  started_at=$(date -Is)
  command_json=$(pipeline_command_json "./$script" "$PIPELINE_SAMPLE_ID" "$@")
  stage_log=$(pipeline_new_stage_log "$stage")
  pipeline_write_provenance start "$stage" "$script" "$command_json" --log "$stage_log"
  pipeline_execute_stage_logged "$stage_log" "$stage" "$script" "./$script" "$PIPELINE_SAMPLE_ID" "$@" || rc=$?
  status=success
  if [[ "$rc" != 0 ]]; then
    status=failed
    reason="exit=$rc; log=$stage_log"
    pipeline_record_failure "$stage" "$script" "$command_json" "$rc" "$stage_log"
  fi
  finished=$(date +%s)
  finished_at=$(date -Is)
  pipeline_write_provenance finish "$stage" "$script" "$command_json" --status "$status" --exit-code "$rc" --elapsed-seconds "$(( finished - started ))" --log "$stage_log" --reason "$reason"
  pipeline_record_stage_runtime "$stage" "$script" "$status" "$started_at" "$finished_at" "$(( finished - started ))" "$rc" "$reason"
  [[ "$rc" == 0 ]] || return "$rc"
  pipeline_collect_statistics "$stage"
  LAST_SUCCESS_STAGE=$stage
  pipeline_record_stage "$stage" success
  pipeline_write_sample_state running 0
  pipeline_finish_if_stop "$stage"
}

pipeline_run_command() {
  local stage=$1 name=$2
  shift 2
  if (( 10#$stage < 10#$PIPELINE_START_AT )); then
    log "── $stage skip: ${PIPELINE_STAGE_LABEL}_START_AT=$PIPELINE_START_AT ──"
    pipeline_record_stage "$stage" skipped "${PIPELINE_STAGE_LABEL}_START_AT=$PIPELINE_START_AT"
    return 0
  fi
  if [[ "${PIPELINE_GROUP_CHILD:-0}" != 1 ]]; then
    CURRENT_STAGE=$stage
    pipeline_write_sample_state running 0
    pipeline_record_stage "$stage" running
  fi
  local started started_at finished finished_at command_json rc=0 status stage_log reason=""
  started=$(date +%s)
  started_at=$(date -Is)
  command_json=$(pipeline_command_json "$@")
  stage_log=$(pipeline_new_stage_log "$stage")
  pipeline_write_provenance start "$stage" "$name" "$command_json" --log "$stage_log"
  log "── $stage $name ──"
  pipeline_execute_stage_logged "$stage_log" "$stage" "$name" "$@" || rc=$?
  status=success
  if [[ "$rc" != 0 ]]; then
    status=failed
    reason="exit=$rc; log=$stage_log"
    pipeline_record_failure "$stage" "$name" "$command_json" "$rc" "$stage_log"
  fi
  finished=$(date +%s)
  finished_at=$(date -Is)
  pipeline_write_provenance finish "$stage" "$name" "$command_json" --status "$status" --exit-code "$rc" --elapsed-seconds "$(( finished - started ))" --log "$stage_log" --reason "$reason"
  pipeline_record_stage_runtime "$stage" "$name" "$status" "$started_at" "$finished_at" "$(( finished - started ))" "$rc" "$reason"
  [[ "$rc" == 0 ]] || return "$rc"
  [[ "${PIPELINE_GROUP_CHILD:-0}" != 1 ]] || return 0
  pipeline_collect_statistics "$stage"
  LAST_SUCCESS_STAGE=$stage
  pipeline_record_stage "$stage" success
  pipeline_write_sample_state running 0
  pipeline_finish_if_stop "$stage"
}

pipeline_skip_stage() {
  local stage=$1 reason=$2 stop_on_skip=${3:-0}
  local skipped_at
  if (( 10#$stage < 10#$PIPELINE_START_AT )); then
    log "── $stage skip: ${PIPELINE_STAGE_LABEL}_START_AT=$PIPELINE_START_AT ──"
    pipeline_record_stage "$stage" skipped "${PIPELINE_STAGE_LABEL}_START_AT=$PIPELINE_START_AT"
    return 0
  fi
  CURRENT_STAGE=$stage
  pipeline_write_sample_state running 0
  log "── $stage skip: $reason ──"
  skipped_at=$(date -Is)
  pipeline_record_stage "$stage" skipped "$reason"
  pipeline_write_provenance skip "$stage" conditional_skip "$(pipeline_command_json)" --status skipped --reason "$reason"
  pipeline_record_stage_runtime "$stage" conditional_skip skipped "$skipped_at" "$skipped_at" 0 0 "$reason"
  if [[ "$stop_on_skip" == 1 ]]; then
    pipeline_finish_if_stop "$stage"
  fi
}
