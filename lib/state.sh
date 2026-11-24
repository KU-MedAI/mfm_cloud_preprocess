#!/usr/bin/env bash
# sample runner 상태 파일을 원자적으로 기록하는 공통 함수.
# 호출자는 PIPELINE_SAMPLE_ID/STATE_FILE/CURRENT_STAGE/LAST_SUCCESS_STAGE를 설정한다.

pipeline_write_sample_state() {
  local status=$1 exit_code=${2:-0}
  local tmp=$STATE_FILE.part.$$
  printf 'sample_id\tstatus\tcurrent_stage\tlast_success_stage\texit_code\tupdated_at\n%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$PIPELINE_SAMPLE_ID" "$status" "$CURRENT_STAGE" "$LAST_SUCCESS_STAGE" \
    "$exit_code" "$(date -Is)" > "$tmp"
  mv "$tmp" "$STATE_FILE"
}

pipeline_record_stage() {
  local stage=$1 status=$2 reason=${3:--}
  [[ -n "${STAGE_FILE:-}" ]] || return 0
  if [[ ! -s "$STAGE_FILE" ]]; then
    printf 'sample_id\tstage\tstatus\treason\tupdated_at\n' > "$STAGE_FILE"
  fi
  # reason은 TSV 한 cell이어야 한다. 호출자가 외부 오류 문자열을 넘겨도 열 수를 보존한다.
  reason=${reason//$'\t'/ }
  reason=${reason//$'\n'/ }
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$PIPELINE_SAMPLE_ID" "$stage" "$status" "$reason" "$(date -Is)" >> "$STAGE_FILE"
}

pipeline_sample_exit_trap() {
  local rc=$?
  if [[ "$rc" -ne 0 ]]; then
    pipeline_record_stage "$CURRENT_STAGE" failed "exit_code=$rc"
    pipeline_write_sample_state failed "$rc"
  fi
  exit "$rc"
}

pipeline_install_sample_trap() {
  trap pipeline_sample_exit_trap EXIT
  pipeline_write_sample_state running 0
}

pipeline_complete_sample() {
  CURRENT_STAGE=complete
  trap - EXIT
  pipeline_write_sample_state success 0
}

# 통합 WGS→RNA runner는 공개 TSV 열이 달라 별도 writer를 쓰되 원자적 교체 규칙은 공유한다.
pipeline_write_integrated_state() {
  local status=$1 exit_code=${2:-0}
  local tmp=$INTEGRATED_STATE_FILE.part.$$
  printf 'sample_id\tstatus\tcurrent_branch\texit_code\thas_rna\tupdated_at\n%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$PIPELINE_SAMPLE_ID" "$status" "$CURRENT_BRANCH" "$exit_code" \
    "$HAS_RNA" "$(date -Is)" > "$tmp"
  mv "$tmp" "$INTEGRATED_STATE_FILE"
}

pipeline_integrated_exit_trap() {
  local rc=$?
  if [[ "$rc" -ne 0 ]]; then
    pipeline_write_integrated_state failed "$rc"
  fi
  exit "$rc"
}

pipeline_install_integrated_trap() {
  trap pipeline_integrated_exit_trap EXIT
  pipeline_write_integrated_state running 0
}

pipeline_complete_integrated() {
  CURRENT_BRANCH=complete
  trap - EXIT
  pipeline_write_integrated_state success 0
}
