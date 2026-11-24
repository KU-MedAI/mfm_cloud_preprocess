#!/usr/bin/env bash
# local runner for process_code형 integrated v2 pipeline.
#
# 사용법:
#   ./runner_local.sh sample <sample_list.txt> [sample_count]
#   ./runner_local.sh wgs <sample_list.txt> [sample_count]
#   ./runner_local.sh rna <sample_list.txt> [sample_count]
#   ./runner_local.sh clair3-qc <sample_list.txt> [sample_count]
#
# LOCAL_JOBS=auto는 CPU/RAM/임시공간 한도에서 동시 sample 수를 정한다.
set -euo pipefail
cd "$(dirname "$0")"
CODE_SOURCE_DIR=$(pwd -P)
source ./config_v2.sh
mkdir -p "$OUTPUT_ROOT"

MODE=${1:?사용법: $0 sample|wgs|rna|clair3-qc <sample_list.txt> [sample_count]}
LIST=${2:-$SAMPLE_LIST}
SAMPLE_COUNT=${3:-${SAMPLE_COUNT:-all}}
[[ -n "$LIST" && -s "$LIST" ]] || die "sample list가 없습니다: $LIST"
# CPU 예산은 runner가 아니라 각 단계 시작 시 공용 pool에서 동적으로 임대한다.
MFM_CPU_BUDGET=${MFM_CPU_BUDGET:-$(nproc)}
LOCAL_JOBS=${LOCAL_JOBS:-auto}
RUNNER_MAX_RETRIES=${RUNNER_MAX_RETRIES:-0}
case "$LOCAL_JOBS" in
  all|ALL|max|MAX)
    LOCAL_JOBS=$(awk 'NF && $1 !~ /^#/ {n++} END {print n+0}' "$LIST")
    ;;
  auto|AUTO) ;;
  *) ;;
esac
if [[ "$LOCAL_JOBS" != auto && "$LOCAL_JOBS" != AUTO ]]; then
  [[ "$LOCAL_JOBS" =~ ^[0-9]+$ && "$LOCAL_JOBS" -ge 1 ]] || die "LOCAL_JOBS는 양의 정수, auto 또는 all|max여야 합니다: $LOCAL_JOBS"
fi
[[ "$MFM_CPU_BUDGET" =~ ^[0-9]+$ && "$MFM_CPU_BUDGET" -ge 1 ]] || die "MFM_CPU_BUDGET는 양의 정수여야 합니다: $MFM_CPU_BUDGET"
CPU_AFFINITY_BUDGET=$(nproc)
if (( MFM_CPU_BUDGET > CPU_AFFINITY_BUDGET )); then
  warn "MFM_CPU_BUDGET를 현재 CPU affinity로 제한: requested=$MFM_CPU_BUDGET available=$CPU_AFFINITY_BUDGET"
  MFM_CPU_BUDGET=$CPU_AFFINITY_BUDGET
fi
case "$MODE" in
  rna)
    SAMPLE_MEMORY_GIB=${RNA_SAMPLE_MEMORY_GIB:-48}
    SAMPLE_TMP_GIB=${RNA_SAMPLE_TMP_GIB:-80}
    MIN_CPU_PER_SAMPLE=${RNA_MIN_CPU_PER_SAMPLE:-6}
    ;;
  sample)
    SAMPLE_MEMORY_GIB=${SAMPLE_MEMORY_GIB:-64}
    SAMPLE_TMP_GIB=${SAMPLE_TMP_GIB:-250}
    MIN_CPU_PER_SAMPLE=${MIN_CPU_PER_SAMPLE:-8}
    ;;
  wgs|clair3-qc)
    SAMPLE_MEMORY_GIB=${WGS_SAMPLE_MEMORY_GIB:-64}
    # HG00281 all_to_both collate 실측 peak가 250GiB에 근접해 400GiB로 보수화한다.
    SAMPLE_TMP_GIB=${WGS_SAMPLE_TMP_GIB:-400}
    MIN_CPU_PER_SAMPLE=${WGS_MIN_CPU_PER_SAMPLE:-8}
    ;;
esac
if [[ -n "${SAMPLE_MEMORY_GIB:-}" ]]; then
  mem_available_gib=$(awk '/MemAvailable/{print int($2/1024/1024)}' /proc/meminfo)
  RESOURCE_DISK_ROOT=${MFM_SCRATCH_ROOT:-$DATASET_ROOT}
  if [[ -n "${MFM_SCRATCH_ROOT:-}" ]]; then
    mkdir -p "$MFM_SCRATCH_ROOT" "$MFM_STAGE11_SCRATCH_ROOT" "$MFM_STAGE12_SCRATCH_ROOT"
  fi
  disk_available_gib=$(df -BG --output=avail "$RESOURCE_DISK_ROOT" | tail -1 | tr -dc '0-9')
  by_mem=$((mem_available_gib / SAMPLE_MEMORY_GIB))
  by_disk=$((disk_available_gib / SAMPLE_TMP_GIB))
  by_cpu=$((MFM_CPU_BUDGET / MIN_CPU_PER_SAMPLE))
  (( by_mem >= 1 )) || by_mem=1
  (( by_disk >= 1 )) || by_disk=1
  (( by_cpu >= 1 )) || by_cpu=1
  resource_jobs=$by_mem
  (( by_disk < resource_jobs )) && resource_jobs=$by_disk
  (( by_cpu < resource_jobs )) && resource_jobs=$by_cpu
  if [[ "$LOCAL_JOBS" == auto || "$LOCAL_JOBS" == AUTO ]]; then
    LOCAL_JOBS=$resource_jobs
    warn "LOCAL_JOBS=auto: CPU/RSS/disk 예산으로 동시 sample 수 산정: jobs=$LOCAL_JOBS cpu=$MFM_CPU_BUDGET min_cpu_per_sample=$MIN_CPU_PER_SAMPLE mem=${mem_available_gib}GiB disk=${disk_available_gib}GiB"
  fi
  if (( LOCAL_JOBS > resource_jobs )); then
    warn "LOCAL_JOBS를 CPU/RSS/disk 예산으로 제한: requested=$LOCAL_JOBS effective=$resource_jobs mem=${mem_available_gib}GiB disk=${disk_available_gib}GiB"
    LOCAL_JOBS=$resource_jobs
  fi
fi
# worker 초기화/안전 정책용 예상 share. 각 실제 단계는 pool grant로 덮어쓴다.
MFM_SAMPLE_CPU_BUDGET=$((MFM_CPU_BUDGET / LOCAL_JOBS))
(( MFM_SAMPLE_CPU_BUDGET >= 1 )) || MFM_SAMPLE_CPU_BUDGET=1
# 샘플 수로 미리 금지하지 않고 실제 단계의 공유 CPU/RAM 예약으로 결정한다.
MFM_HAP_PARALLEL_ALLOWED=${MFM_HAP_PARALLEL_ALLOWED:-1}
[[ "$MFM_HAP_PARALLEL_ALLOWED" == 0 || "$MFM_HAP_PARALLEL_ALLOWED" == 1 ]] || die "MFM_HAP_PARALLEL_ALLOWED는 0 또는 1이어야 합니다."
export MFM_CPU_BUDGET MFM_SAMPLE_CPU_BUDGET MFM_HAP_PARALLEL_ALLOWED
[[ "$RUNNER_MAX_RETRIES" =~ ^[0-9]+$ ]] || die "RUNNER_MAX_RETRIES는 0 이상의 정수여야 합니다: $RUNNER_MAX_RETRIES"
case "$SAMPLE_COUNT" in
  all|ALL|max|MAX) SAMPLE_COUNT=all ;;
  *) [[ "$SAMPLE_COUNT" =~ ^[0-9]+$ && "$SAMPLE_COUNT" -ge 1 ]] || die "sample_count는 양의 정수 또는 all이어야 합니다: $SAMPLE_COUNT" ;;
esac

case "$MODE" in
  sample) SCRIPT_REL=run_sample.sh ;;
  wgs) SCRIPT_REL=run_wgs_sample.sh ;;
  rna) SCRIPT_REL=run_rna_sample.sh ;;
  clair3-qc) SCRIPT_REL=05_clair3_qc.sh ;;
  *) die "unknown mode: $MODE" ;;
esac

log "local runner: mode=$MODE jobs=$LOCAL_JOBS cpu_budget=$MFM_CPU_BUDGET sample_cpu_budget=dynamic-stage-lease hap_parallel_allowed=$MFM_HAP_PARALLEL_ALLOWED sample_count=$SAMPLE_COUNT max_retries=$RUNNER_MAX_RETRIES disk_root=${RESOURCE_DISK_ROOT:-$DATASET_ROOT} list=$LIST"
RUN_STATUS_DIR=${RUN_STATUS_DIR:-$(runner_status_root "$MODE")}
RUN_LOG_DIR=${RUN_LOG_DIR:-$RUN_STATUS_DIR/logs}
RUN_STATUS_TSV=${RUN_STATUS_TSV:-$RUN_STATUS_DIR/status.tsv}
RETRY_LIST=${RETRY_LIST:-$RUN_STATUS_DIR/retry_samples.txt}
mkdir -p "$RUN_STATUS_DIR/results" "$RUN_LOG_DIR"
MFM_CPU_POOL_DIR=${MFM_CPU_POOL_DIR:-$RUN_STATUS_DIR/cpu_pool}
MFM_CPU_POOL_AUDIT=${MFM_CPU_POOL_AUDIT:-$RUN_STATUS_DIR/cpu_pool_usage.jsonl}
export MFM_CPU_POOL_DIR MFM_CPU_POOL_AUDIT

# 실행 시작 시 코드만 run-status 아래에 고정한다. output이 v2 내부에 있어도
# output/.git/cache/patch 임시 파일은 fingerprint와 snapshot에서 제외한다.
SNAPSHOT_OUTPUT_ROOT=$(readlink -m "$OUTPUT_ROOT")
snapshot_files() {
  find "$CODE_SOURCE_DIR" \
    \( -path "$SNAPSHOT_OUTPUT_ROOT" -o -path "$CODE_SOURCE_DIR/.git" -o \
       -name __pycache__ -o -name .pytest_cache \) -prune -o \
    -type f ! -name '*.orig' ! -name '*.rej' ! -name '*.pyc' -print0
}
SNAPSHOT_FINGERPRINT=$(snapshot_files | sort -z | xargs -0 sha256sum | sha256sum | awk '{print $1}')
RUN_CODE_SNAPSHOT_DIR=${RUN_CODE_SNAPSHOT_DIR:-$RUN_STATUS_DIR/code_snapshot_${SNAPSHOT_FINGERPRINT:0:12}}
if [[ ! -x "$RUN_CODE_SNAPSHOT_DIR/$SCRIPT_REL" ]]; then
  SNAPSHOT_TMP="$RUN_STATUS_DIR/.code_snapshot_${SNAPSHOT_FINGERPRINT:0:12}.part.$$"
  [[ ! -e "$SNAPSHOT_TMP" ]] || die "code snapshot temporary path exists: $SNAPSHOT_TMP"
  mkdir -p "$SNAPSHOT_TMP"
  TAR_EXCLUDES=(--exclude='./output' --exclude='./.git' --exclude='*/__pycache__' --exclude='*.pyc' --exclude='*.orig' --exclude='*.rej')
  if [[ "$SNAPSHOT_OUTPUT_ROOT" == "$CODE_SOURCE_DIR"/* ]]; then
    TAR_EXCLUDES+=(--exclude="./${SNAPSHOT_OUTPUT_ROOT#"$CODE_SOURCE_DIR"/}")
  fi
  tar -C "$CODE_SOURCE_DIR" "${TAR_EXCLUDES[@]}" -cf - . | tar -C "$SNAPSHOT_TMP" -xf -
  (
    cd "$SNAPSHOT_TMP"
    find . -type f -not -path '*/__pycache__/*' -not -name 'CODE_SNAPSHOT.sha256' -print0 | sort -z | xargs -0 sha256sum
  ) > "$SNAPSHOT_TMP/CODE_SNAPSHOT.sha256"
  mv "$SNAPSHOT_TMP" "$RUN_CODE_SNAPSHOT_DIR"
fi
[[ -x "$RUN_CODE_SNAPSHOT_DIR/$SCRIPT_REL" ]] || die "missing code snapshot target: $RUN_CODE_SNAPSHOT_DIR/$SCRIPT_REL"
SCRIPT="$RUN_CODE_SNAPSHOT_DIR/$SCRIPT_REL"
log "local runner code snapshot: $RUN_CODE_SNAPSHOT_DIR"
export RUN_STATUS_DIR RUN_LOG_DIR RUNNER_TARGET_SCRIPT="$SCRIPT"
SELECTED_LIST="$RUN_STATUS_DIR/selected_samples.txt"
if [[ "$SAMPLE_COUNT" == all ]]; then
  awk 'NF && $1 !~ /^#/ {print $1}' "$LIST" > "$SELECTED_LIST.part"
else
  awk -v n="$SAMPLE_COUNT" 'NF && $1 !~ /^#/ {if (++k <= n) print $1}' "$LIST" > "$SELECTED_LIST.part"
fi
mv "$SELECTED_LIST.part" "$SELECTED_LIST"
[[ -s "$SELECTED_LIST" ]] || die "실행할 sample이 없습니다: $LIST"

attempt=0
CURRENT_LIST="$SELECTED_LIST"
while :; do
  export RUNNER_ATTEMPT="$attempt"
  log "local runner pass: attempt=$attempt list=$CURRENT_LIST"
  xargs -r -n1 -P "$LOCAL_JOBS" ./scripts/run_one_local.sh < "$CURRENT_LIST"

  {
    printf 'sample_id\tstatus\texit_code\telapsed_sec\tattempt\tlog_path\tmax_rss_kib\tuser_cpu_sec\tsystem_cpu_sec\tresource_path\n'
    while IFS= read -r selected_sample; do
      result_file="$RUN_STATUS_DIR/results/$selected_sample.result.tsv"
      [[ -s "$result_file" ]] && cat "$result_file"
    done < "$SELECTED_LIST"
  } > "$RUN_STATUS_TSV.part"
  mv "$RUN_STATUS_TSV.part" "$RUN_STATUS_TSV"
  awk -F'\t' 'NR>1 && $2=="failed"{print $1}' "$RUN_STATUS_TSV" > "$RETRY_LIST.part"
  mv "$RETRY_LIST.part" "$RETRY_LIST"
  failed=$(awk -F'\t' 'NR>1 && $2=="failed"{n++}END{print n+0}' "$RUN_STATUS_TSV")
  [[ "$failed" == 0 ]] && break
  [[ "$attempt" -ge "$RUNNER_MAX_RETRIES" ]] && break
  cp "$RETRY_LIST" "$RUN_STATUS_DIR/retry_attempt_$((attempt + 1)).txt"
  CURRENT_LIST="$RUN_STATUS_DIR/retry_attempt_$((attempt + 1)).txt"
  attempt=$((attempt + 1))
done

log "local runner 완료: status=$RUN_STATUS_TSV retry=$RETRY_LIST failed=$failed"
[[ "$failed" == 0 ]]
