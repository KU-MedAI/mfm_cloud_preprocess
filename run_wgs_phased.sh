#!/usr/bin/env bash
# WGS 목록을 준비(01-10)·개인정렬(11)·compact(12-20) 큐로 분리한다.
# 단계별 sample 병렬도를 달리해 단일코어 준비 구간의 CPU 유휴를 줄인다.
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
source ./config_v2.sh

LIST=${1:?사용법: $0 <sample_list> [sample_count]}
SAMPLE_COUNT=${2:-all}
[[ -s "$LIST" ]] || die "sample list가 없거나 비어 있음: $LIST"
LIST=$(readlink -m "$LIST")
case "$SAMPLE_COUNT" in
  all|ALL|max|MAX) SAMPLE_COUNT=all ;;
  *) [[ "$SAMPLE_COUNT" =~ ^[0-9]+$ && "$SAMPLE_COUNT" -ge 1 ]] || die "sample_count 오류: $SAMPLE_COUNT" ;;
esac

BASE_RUN_ID=${RUN_ID:-$(date +%y%m%d)}
PHASED_RUN_STATUS_ROOT=${PHASED_RUN_STATUS_ROOT:-$OUTPUT_ROOT/run-status/$BASE_RUN_ID/wgs_phased}
PREP_SAMPLE_JOBS=${WGS_PREP_SAMPLE_JOBS:-4}
ALIGN_SAMPLE_JOBS=${WGS_ALIGN_SAMPLE_JOBS:-4}
COMPACT_SAMPLE_JOBS=${WGS_COMPACT_SAMPLE_JOBS:-2}
mkdir -p "$PHASED_RUN_STATUS_ROOT"

SELECTED="$PHASED_RUN_STATUS_ROOT/selected_samples.txt"
if [[ "$SAMPLE_COUNT" == all ]]; then
  awk 'NF && $1 !~ /^#/ {print $1}' "$LIST" > "$SELECTED.part"
else
  awk -v n="$SAMPLE_COUNT" 'NF && $1 !~ /^#/ {if (++k <= n) print $1}' "$LIST" > "$SELECTED.part"
fi
mv "$SELECTED.part" "$SELECTED"
[[ -s "$SELECTED" ]] || die "실행할 sample이 없음: $LIST"

FAILED_ALL="$PHASED_RUN_STATUS_ROOT/failed_samples.txt"
: > "$FAILED_ALL"

run_phase() {
  local name=$1 start=$2 stop=$3 jobs=$4 input=$5 output=$6
  : > "$output"
  [[ -s "$input" ]] || { warn "phased WGS $name: 입력 성공 sample이 없어 생략"; return 0; }
  local phase_status="$PHASED_RUN_STATUS_ROOT/$name"
  local status_tsv="$phase_status/status.tsv"
  log "phased WGS: phase=$name start=$start stop=$stop jobs=$jobs input=$input"
  set +e
  RUN_ID="${BASE_RUN_ID}_${name}" \
  RUN_STATUS_DIR="$phase_status" RUN_LOG_DIR="$phase_status/logs" \
  WGS_START_AT="$start" WGS_STOP_AFTER="$stop" LOCAL_JOBS="$jobs" \
    ./run_pipeline.sh wgs-list "$input" all
  local rc=$?
  set -e
  if [[ ! -s "$status_tsv" ]]; then
    warn "phased WGS $name: status.tsv 없음; 입력 전체를 실패 처리"
    awk 'NF && $1 !~ /^#/ {print $1}' "$input" >> "$FAILED_ALL"
    return 0
  fi
  awk -F '\t' 'NR>1 && $2=="success" {print $1}' "$status_tsv" > "$output.part"
  mv "$output.part" "$output"
  awk -F '\t' 'NR>1 && $2=="failed" {print $1}' "$status_tsv" >> "$FAILED_ALL"
  [[ "$rc" == 0 ]] || warn "phased WGS $name: 일부 sample 실패 rc=$rc; 성공 sample만 다음 단계로 진행"
}

PREP_OK="$PHASED_RUN_STATUS_ROOT/01_10.success.txt"
ALIGN_OK="$PHASED_RUN_STATUS_ROOT/11.success.txt"
FINAL_OK="$PHASED_RUN_STATUS_ROOT/12_20.success.txt"
run_phase 01_10 01 10 "$PREP_SAMPLE_JOBS" "$SELECTED" "$PREP_OK"
run_phase 11 11 11 "$ALIGN_SAMPLE_JOBS" "$PREP_OK" "$ALIGN_OK"
run_phase 12_20 12 20 "$COMPACT_SAMPLE_JOBS" "$ALIGN_OK" "$FINAL_OK"

sort -u "$FAILED_ALL" -o "$FAILED_ALL"
selected_count=$(wc -l < "$SELECTED")
success_count=$(wc -l < "$FINAL_OK")
failed_count=$(wc -l < "$FAILED_ALL")
printf 'selected\tsuccess\tfailed\tprep_jobs\talign_jobs\tcompact_jobs\n%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$selected_count" "$success_count" "$failed_count" \
  "$PREP_SAMPLE_JOBS" "$ALIGN_SAMPLE_JOBS" "$COMPACT_SAMPLE_JOBS" \
  > "$PHASED_RUN_STATUS_ROOT/summary.tsv"
log "phased WGS 완료: selected=$selected_count success=$success_count failed=$failed_count root=$PHASED_RUN_STATUS_ROOT"
[[ "$failed_count" == 0 && "$success_count" == "$selected_count" ]]
