#!/usr/bin/env bash
# 부모가 stage 상태를 확정하고 두 HAP의 명령/원장은 독립적으로 기록함.
pipeline_group_cancel() {
  [[ ${#PIPELINE_GROUP_PIDS[@]} -gt 0 ]] || return 0
  "${PYTHON:-python3}" - "${PIPELINE_GROUP_PIDS[@]}" <<'PY'
import os,signal,sys
roots={int(x) for x in sys.argv[1:]};parents={}
for name in os.listdir('/proc'):
    if not name.isdigit():continue
    try:
        text=open('/proc/'+name+'/stat').read();parts=text[text.rfind(')')+2:].split();parents[int(name)]=int(parts[1])
    except (OSError,ValueError,IndexError):continue
members=set(roots)
while True:
    updated=members|{pid for pid,parent in parents.items() if parent in members}
    if updated==members:break
    members=updated
managers=set();leaves=set()
for pid in members-roots:
    try:args=open('/proc/'+str(pid)+'/cmdline','rb').read().split(b'\0')
    except OSError:continue
    if any(arg.endswith(b'/cpu_pool.py') for arg in args):managers.add(pid)
    if pid not in parents.values() and args and os.path.basename(args[0])!=b'tee':leaves.add(pid)
# 관리자가 자식을 wait하여 lease와zombie를정리하도록 관리자부터 종료함.
for pid in (managers if managers else leaves)|roots:
    try:os.kill(pid,signal.SIGTERM)
    except ProcessLookupError:pass
PY
}

pipeline_group_child_signal() {
  local rc=$1
  trap '' INT TERM
  local -a PIPELINE_GROUP_PIDS=("$BASHPID")
  pipeline_group_cancel
  wait || :
  exit "$rc"
}

pipeline_group_signal() {
  local rc=$1 pid
  trap '' INT TERM
  pipeline_group_cancel
  for pid in "${PIPELINE_GROUP_PIDS[@]}"; do wait "$pid" 2>/dev/null || :; done
  exit "$rc"
}

pipeline_run_hap_group() {
  local stage=$1 callback=$2 jobs=${3:-2} h pid rc=0 part_rc=0
  [[ "$jobs" == 1 || "$jobs" == 2 ]] || { echo 'HAP 동시작업은 1 또는 2여야 함' >&2; return 2; }
  if (( 10#$stage < 10#$PIPELINE_START_AT )); then
    pipeline_record_stage "$stage" skipped "${PIPELINE_STAGE_LABEL}_START_AT=$PIPELINE_START_AT"
    return 0
  fi
  CURRENT_STAGE=$stage
  pipeline_write_sample_state running 0
  pipeline_record_stage "$stage" running
  local old_int old_term
  old_int=$(trap -p INT); old_term=$(trap -p TERM)
  local -a PIPELINE_GROUP_PIDS=()
  trap 'pipeline_group_signal 130' INT
  trap 'pipeline_group_signal 143' TERM
  for h in 1 2; do
    (
      trap - EXIT
      trap 'pipeline_group_child_signal 130' INT
      trap 'pipeline_group_child_signal 143' TERM
      PIPELINE_GROUP_CHILD=1
      PIPELINE_STOP_AFTER=
      "$callback" "$h"
    ) &
    PIPELINE_GROUP_PIDS+=("$!")
    if [[ "$jobs" == 1 ]]; then
      part_rc=0; wait "$!" || part_rc=$?
      [[ "$part_rc" == 0 ]] || { rc=$part_rc; break; }
      PIPELINE_GROUP_PIDS=()
    fi
  done
  # 실패한 HAP가 있어도 이미 실행 중인 형제의 유효한 완료 결과는 보존함.
  for pid in "${PIPELINE_GROUP_PIDS[@]}"; do
    part_rc=0; wait "$pid" || part_rc=$?
    [[ "$part_rc" == 0 || "$rc" != 0 ]] || rc=$part_rc
  done
  trap - INT TERM
  [[ -z "$old_int" ]] || eval "$old_int"
  [[ -z "$old_term" ]] || eval "$old_term"
  if [[ "$rc" != 0 ]]; then
    pipeline_record_stage "$stage" failed "HAP stage group exit=$rc"
    return "$rc"
  fi
  pipeline_collect_statistics "$stage"
  LAST_SUCCESS_STAGE=$stage
  pipeline_record_stage "$stage" success
  pipeline_write_sample_state running 0
  pipeline_finish_if_stop "$stage"
}
