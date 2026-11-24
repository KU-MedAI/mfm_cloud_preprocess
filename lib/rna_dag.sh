#!/usr/bin/env bash
# 내장 runner의 독립 분기만 겹침. 명시한 START/STOP은 번호순 경로를 사용함.
rna_dag_enabled() {
  local mode=${RNA_PIPELINE_OVERLAP:-auto}
  [[ "$mode" == auto || "$mode" == 0 || "$mode" == 1 ]] || { echo 'RNA_PIPELINE_OVERLAP는 auto/0/1임' >&2; return 2; }
  [[ "$mode" != 0 && -z "${PIPELINE_STOP_AFTER:-}" && "${RNA_STAR_INDEX_BENCHMARK:-0}" != 1 ]] || return 1
  (( 10#$PIPELINE_START_AT <= 14 )) || return 1
  [[ "$RNA_STAR_INDEX_HAP_JOBS" == 2 && "$RNA_MARKDUP_HAP_JOBS" == 2 && "$RNA_ALIGN_HAP_JOBS" == 2 && "$RNA_COVERAGE_HAP_JOBS" == 2 ]] || return 1
  RNA_DAG_CPU=$("$PYTHON" -c 'from integrate_wgs_rna.gpu_compare_resources import cpu_limit;import os;print(cpu_limit(int(os.environ.get("MFM_CPU_BUDGET",os.environ.get("MFM_SAMPLE_CPU_BUDGET","1")))))') || return
  (( RNA_DAG_CPU >= 8 ))
}

rna_dag_wait() {
  local name=$1
  until [[ -f "$RNA_DAG_DIR/$name" ]]; do
    [[ ! -f "$RNA_DAG_DIR/failed" ]] || return 1
    sleep 0.1
  done
}

rna_dag_cap() {
  local role=$1 evidence=$((RNA_DAG_CPU/4)) reserved=0 haps=0 cap
  if [[ "$role" == variants || "$role" == cnv || "$role" == ase ]]; then cap=$evidence
  else
    [[ -f "$RNA_DAG_DIR/variants" ]] || reserved=$evidence
    if [[ -f "$RNA_DAG_DIR/marked1" || -f "$RNA_DAG_DIR/marked2" ]]; then
      [[ -f "$RNA_DAG_DIR/ase" ]] || reserved=$evidence
    fi
    [[ -f "$RNA_DAG_DIR/cnv" ]] || reserved=$((reserved+evidence))
    [[ -f "$RNA_DAG_DIR/hap1" ]] || haps=$((haps+1))
    [[ -f "$RNA_DAG_DIR/hap2" ]] || haps=$((haps+1))
    (( haps>0 )) || haps=1
    cap=$(((RNA_DAG_CPU-reserved)/haps))
  fi
  (( cap>0 )) || cap=1
  if [[ -n "${MFM_CPU_POOL_STAGE_MAX:-}" ]] && (( MFM_CPU_POOL_STAGE_MAX < cap )); then cap=$MFM_CPU_POOL_STAGE_MAX; fi
  printf '%s\n' "$cap"
}

rna_dag_step() {
  local role=$1 stage=$2 callback=$3 cap dynamic=$RNA_DAG_CPU rc=0
  shift 3
  cap=$(rna_dag_cap "$role") || return
  CURRENT_STAGE=$stage
  if [[ -n ${MFM_CPU_POOL_STAGE_MAX:-} ]] && (( MFM_CPU_POOL_STAGE_MAX < dynamic )); then dynamic=$MFM_CPU_POOL_STAGE_MAX; fi
  if [[ ${MFM_CPU_UTILIZATION:-maximum} != maximum || ${MFM_CPU_POOL_PROFILE:-throughput} != throughput ]]; then dynamic=$cap; fi
  case "$role" in variants|cnv|ase) ;; *) dynamic=$cap ;; esac
  printf '%s\t%s\t%s\t%s\t%s\n' "$role" "$stage" "$cap" "$(date +%s.%N)" "$dynamic" >> "$RNA_DAG_DIR/$role.execution.tsv"
  case "$role" in
    variants|cnv|ase) MFM_RNA_DAG_ROLE=$role MFM_RNA_DAG_DYNAMIC_CPU_MAX=$dynamic MFM_CPU_POOL_STAGE_MAX=$cap "$callback" "$@" || rc=$? ;;
    *) MFM_CPU_POOL_STAGE_MAX=$cap "$callback" "$@" || rc=$? ;;
  esac
  printf '%s\t%s\t%s\t%s\n' "$role" "$stage" "$(date +%s.%N)" "$rc" >> "$RNA_DAG_DIR/$role.finish.tsv"
  if [[ "$rc" != 0 ]]; then touch "$RNA_DAG_DIR/failed"; return "$rc"; fi
}

rna_dag_role() {
  local role=$1 h=${2:-} rc=0
  trap 'rna_role_rc=$?; if [[ "$rna_role_rc" != 0 ]]; then printf "%s\t%s\t%s\n" "$role" "$CURRENT_STAGE" "$rna_role_rc" >> "$RNA_DAG_DIR/failure.tsv"; touch "$RNA_DAG_DIR/failed"; fi' EXIT
  trap 'pipeline_group_child_signal 130' INT
  trap 'pipeline_group_child_signal 143' TERM
  PIPELINE_GROUP_CHILD=1
  PIPELINE_STOP_AFTER=
  STATE_FILE=$RNA_DAG_DIR/$role.state.tsv
  STAGE_FILE=$RNA_DAG_DIR/$role.stages.tsv
  case "$role" in
    variants) rna_dag_step variants 20 rna_variants || return; touch "$RNA_DAG_DIR/variants" ;;
    cnv) rna_dag_step cnv 22 rna_cnv || return; touch "$RNA_DAG_DIR/cnv" ;;
    ase)
      rna_dag_wait variants || return
      rna_dag_wait marked1 || return
      rna_dag_wait marked2 || return
      rna_dag_step ase 21 rna_allelic || return; touch "$RNA_DAG_DIR/ase" ;;
    hap*)
      # 짧은 근거 계산 뒤 index를 시작해, 고정 샘플 CPU가 4+4에 묶이지 않게 함.
      if [[ ${MFM_CPU_ALLOCATION_MODE:-dynamic} == fixed_sample ]]; then
        rna_dag_wait variants || return
        rna_dag_wait cnv || return
      fi
      rna_dag_step "$role" 15 rna_star_index_hap "$h" || return
      rna_dag_step "$role" 16 rna_star_align_hap "$h" || return
      rna_dag_step "$role" 17 rna_markdup_hap "$h" || return
      touch "$RNA_DAG_DIR/marked$h"
      rna_dag_step "$role" 18 rna_coverage_hap "$h" || return
      rna_dag_step "$role" 19 rna_qc_hap "$h" || return
      touch "$RNA_DAG_DIR/hap$h" ;;
  esac
}

pipeline_rna_dag() {
  RNA_DAG_DIR=$(mktemp -d "$W/rna/dag.XXXXXXXX") || return
  local old_int old_term pid rc=0 part_rc role
  old_int=$(trap -p INT);old_term=$(trap -p TERM)
  local -a PIPELINE_GROUP_PIDS=()
  trap 'pipeline_group_signal 130' INT
  trap 'pipeline_group_signal 143' TERM
  CURRENT_STAGE=15;pipeline_write_sample_state running 0
  for role in variants cnv hap1 hap2 ase; do
    (rna_dag_role "$role" "${role#hap}") &
    PIPELINE_GROUP_PIDS+=("$!")
  done
  for pid in "${PIPELINE_GROUP_PIDS[@]}"; do
    part_rc=0;wait "$pid" || part_rc=$?
    [[ "$part_rc" == 0 || "$rc" != 0 ]] || rc=$part_rc
  done
  trap - INT TERM
  [[ -z "$old_int" ]] || eval "$old_int"
  [[ -z "$old_term" ]] || eval "$old_term"
  "$PYTHON" - "$RNA_DAG_DIR" "$F/qc/rna_dag_execution.json" "$RNA_DAG_CPU" "$rc" <<'PY'
from pathlib import Path
import json,sys
folder,out,cpu,rc=sys.argv[1:];rows=[];finished={}
for path in Path(folder).glob('*.finish.tsv'):
    for line in path.read_text().splitlines():
        role,stage,ended,code=line.split('\t');finished[role,stage]=dict(finished_at=float(ended),exit_code=int(code))
for path in Path(folder).glob('*.execution.tsv'):
    for line in path.read_text().splitlines():
        values=line.split('\t');role,stage,cap,started=values[:4];rows.append(dict(role=role,stage=stage,cpu_cap=int(cap),elastic_cpu_upper=int(values[4]) if len(values)>4 else int(cap),started_at=float(started),**finished.get((role,stage),{})))
target=Path(out);target.parent.mkdir(parents=True,exist_ok=True)
target.write_text(json.dumps(dict(schema='mfm_rna_dag_v1',status='pass' if rc=='0' else 'failed',cpu_pool=int(cpu),phases=sorted(rows,key=lambda x:x['started_at']),resource_note='같은 공유풀을 사용함. 예약 예외도 분기 상한을 적용함. cpu_cap은 입장상한이며 elastic_cpu_upper까지 확대는 native·공유lease·CPU/RAM적합성 확인시만 허용함. 비예약/대여/고정thread는 유지하며 실제예약은 CPU audit로 구분함.'),ensure_ascii=False,indent=2)+'\n')
PY
  if [[ "$rc" != 0 ]]; then
    if [[ -s "$RNA_DAG_DIR/failure.tsv" ]]; then local failed_role failed_rc; IFS=$'\t' read -r failed_role CURRENT_STAGE failed_rc < "$RNA_DAG_DIR/failure.tsv"; fi
    return "$rc"
  fi
  for role in 15 16 17 18; do pipeline_collect_statistics "$role";pipeline_record_stage "$role" success; done
  LAST_SUCCESS_STAGE=18;CURRENT_STAGE=18;pipeline_write_sample_state running 0
}

pipeline_rna_dag_commit() {
  local stage
  for stage in 20 21 22; do
    pipeline_collect_statistics "$stage"
    if [[ "$stage" == 22 && "${RNA_CNV_FILTER:-1}" != 1 ]]; then pipeline_record_stage "$stage" skipped 'RNA_CNV_FILTER=0'
    else pipeline_record_stage "$stage" success; fi
  done
  LAST_SUCCESS_STAGE=22;CURRENT_STAGE=22;pipeline_write_sample_state running 0
}
