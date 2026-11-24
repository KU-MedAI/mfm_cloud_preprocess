#!/usr/bin/env bash
# 계층화 shard 하나를 다운로드·검증하고 pair 우선으로 전처리한다.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARD_ID=${1:?사용법: $0 shard_001 [--data-root DIR] [--output-root DIR] [--shard-root DIR] [--mode all|download|verify|preprocess] [--dry-run]}
shift
SHARD_ROOT=${SHARD_ROOT:-$SCRIPT_DIR/share_package/sharded_delivery}
DATA_ROOT=${DATA_ROOT:-$SHARD_ROOT/data}
OUTPUT_ROOT=${OUTPUT_ROOT:-$SHARD_ROOT/output}
MODE=all
DRY_RUN=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-root) DATA_ROOT=${2:?}; shift 2 ;;
    --output-root) OUTPUT_ROOT=${2:?}; shift 2 ;;
    --shard-root) SHARD_ROOT=${2:?}; shift 2 ;;
    --mode) MODE=${2:?}; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    *) echo "[FATAL] 알 수 없는 옵션: $1" >&2; exit 2 ;;
  esac
done
case "$MODE" in all|download|verify|preprocess) ;; *) echo "[FATAL] mode 오류: $MODE" >&2; exit 2 ;; esac
SHARD_DIR="$SHARD_ROOT/shards/$SHARD_ID"
MANIFEST="$SHARD_DIR/shard_manifest.tsv"
PAIR_LIST="$SHARD_DIR/wgs_rna_pair.txt"
WGS_LIST="$SHARD_DIR/wgs_only.txt"
DOWNLOAD_PLAN="$SHARD_DIR/download_plan.tsv"
SHARED_PLAN="$SHARD_ROOT/shared/download_plan.tsv"
[[ -s "$MANIFEST" ]] || { echo "[FATAL] shard manifest 없음: $MANIFEST" >&2; exit 1; }
[[ -s "$DOWNLOAD_PLAN" ]] || { echo "[FATAL] download plan 없음: $DOWNLOAD_PLAN" >&2; exit 1; }
mkdir -p "$DATA_ROOT" "$OUTPUT_ROOT" "$SHARD_DIR/qc"
DATA_ROOT=$(readlink -m "$DATA_ROOT")
OUTPUT_ROOT=$(readlink -m "$OUTPUT_ROOT")
export PIPELINE_MANIFEST="$MANIFEST" INPUT_ROOT="$DATA_ROOT" OUTPUT_ROOT MFM_CPU_BUDGET=${MFM_CPU_BUDGET:-72}

DOWNLOAD_JOBS=${DOWNLOAD_JOBS:-2}
[[ "$DOWNLOAD_JOBS" =~ ^[0-9]+$ && "$DOWNLOAD_JOBS" -ge 1 ]] || { echo "[FATAL] DOWNLOAD_JOBS 오류: $DOWNLOAD_JOBS" >&2; exit 2; }
download_args=(--plan "$DOWNLOAD_PLAN" --data-root "$DATA_ROOT" --jobs "$DOWNLOAD_JOBS")
[[ "$DRY_RUN" == 0 ]] || download_args+=(--dry-run)
if [[ "$MODE" == all || "$MODE" == download ]]; then
  if [[ -s "$SHARED_PLAN" ]]; then
    shared_args=(--plan "$SHARED_PLAN" --data-root "$DATA_ROOT" --jobs "$DOWNLOAD_JOBS")
    [[ "$DRY_RUN" == 0 ]] || shared_args+=(--dry-run)
    "$SCRIPT_DIR/scripts/download_shard.py" "${shared_args[@]}"
  fi
  "$SCRIPT_DIR/scripts/download_shard.py" "${download_args[@]}"
fi
if [[ "$MODE" == verify ]]; then
  if [[ -s "$SHARED_PLAN" ]]; then
    "$SCRIPT_DIR/scripts/download_shard.py" --plan "$SHARED_PLAN" --data-root "$DATA_ROOT" --verify-only --jobs "$DOWNLOAD_JOBS"
  fi
  "$SCRIPT_DIR/scripts/download_shard.py" --plan "$DOWNLOAD_PLAN" --data-root "$DATA_ROOT" --verify-only --jobs "$DOWNLOAD_JOBS"
  exit $?
fi
[[ "$MODE" == all || "$MODE" == preprocess ]] || exit 0
if [[ "$DRY_RUN" == 1 ]]; then
  echo "[DRY-RUN] pair WGS: PIPELINE_MANIFEST=$MANIFEST INPUT_ROOT=$DATA_ROOT OUTPUT_ROOT=$OUTPUT_ROOT ./run_pipeline.sh wgs-phased-list $PAIR_LIST all"
  echo "[DRY-RUN] pair RNA: PIPELINE_MANIFEST=$MANIFEST INPUT_ROOT=$DATA_ROOT OUTPUT_ROOT=$OUTPUT_ROOT LOCAL_JOBS=${LOCAL_JOBS_RNA:-12} ./run_pipeline.sh rna-list $PAIR_LIST all"
  echo "[DRY-RUN] WGS-only: PIPELINE_MANIFEST=$MANIFEST INPUT_ROOT=$DATA_ROOT OUTPUT_ROOT=$OUTPUT_ROOT ./run_pipeline.sh wgs-phased-list $WGS_LIST all"
  exit 0
fi

failures=0
run_list() {
  local mode=$1 list=$2 jobs=$3
  [[ -s "$list" ]] || return 0
  echo "[RUN] shard=$SHARD_ID mode=$mode jobs=$jobs list=$list"
  set +e
  local list_mode="$mode-list"
  if [[ "$mode" == wgs && "${WGS_PHASED_SCHEDULER:-1}" == 1 ]]; then
    list_mode=wgs-phased-list
  fi
  LOCAL_JOBS="$jobs" "$SCRIPT_DIR/run_pipeline.sh" "$list_mode" "$list" all
  local rc=$?
  set -e
  if [[ "$rc" -ne 0 ]]; then
    failures=$((failures + 1))
    echo "[WARN] 일부 sample 실패: shard=$SHARD_ID mode=$mode rc=$rc; 다음 단계/샘플 계속" >&2
  fi
}
# WGS–RNA pair를 먼저 완성한다. WGS는 기본적으로 단계형 4/4/2 sample
# 병렬도를 사용하되 runner가 CPU/RSS/scratch 여유에 맞춰 자동으로 낮춘다.
run_list wgs "$PAIR_LIST" "${LOCAL_JOBS_WGS:-2}"
run_list rna "$PAIR_LIST" "${LOCAL_JOBS_RNA:-12}"
run_list wgs "$WGS_LIST" "${LOCAL_JOBS_WGS:-2}"
printf '{"shard_id":"%s","status":"%s","failed_batches":%d,"completed_at":"%s"}\n'   "$SHARD_ID" "$([[ "$failures" == 0 ]] && echo pass || echo partial_fail)" "$failures" "$(date -Is)"   > "$SHARD_DIR/qc/run_status.json.part"
mv "$SHARD_DIR/qc/run_status.json.part" "$SHARD_DIR/qc/run_status.json"
[[ "$failures" == 0 ]]
