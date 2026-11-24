#!/usr/bin/env bash
# 현재 shard 계산과 다음 shard 다운로드/checksum을 중첩하는 실행기.
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARD_LIST=${1:?사용법: $0 shard_ids.txt [--data-root DIR] [--output-root DIR] [--shard-root DIR]}
shift
SHARD_ROOT=${SHARD_ROOT:-$SCRIPT_DIR/share_package/sharded_delivery}
DATA_ROOT=${DATA_ROOT:-$SHARD_ROOT/data}
OUTPUT_ROOT=${OUTPUT_ROOT:-$SHARD_ROOT/output}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-root) DATA_ROOT=${2:?}; shift 2 ;;
    --output-root) OUTPUT_ROOT=${2:?}; shift 2 ;;
    --shard-root) SHARD_ROOT=${2:?}; shift 2 ;;
    *) echo "[FATAL] 알 수 없는 옵션: $1" >&2; exit 2 ;;
  esac
done
[[ -s "$SHARD_LIST" ]] || { echo "[FATAL] shard list 없음: $SHARD_LIST" >&2; exit 1; }
mapfile -t SHARDS < <(awk 'NF && $1 !~ /^#/ {print $1}' "$SHARD_LIST")
[[ ${#SHARDS[@]} -gt 0 ]] || { echo "[FATAL] 실행할 shard가 없음" >&2; exit 1; }
mkdir -p "$SHARD_ROOT/queue_logs" "$DATA_ROOT" "$OUTPUT_ROOT"
STATUS="$SHARD_ROOT/queue_logs/queue_status.tsv"
printf 'shard_id\tdownload\tverify\tpreprocess\n' > "$STATUS.part"

download_shard() {
  local shard=$1 log=$2
  "$SCRIPT_DIR/run_shard.sh" "$shard" --shard-root "$SHARD_ROOT" \
    --data-root "$DATA_ROOT" --output-root "$OUTPUT_ROOT" --mode download >"$log" 2>&1
}

declare -A DOWNLOAD_STATE
first=${SHARDS[0]}
first_log="$SHARD_ROOT/queue_logs/$first.download.log"
echo "[QUEUE] 최초 다운로드: $first"
if download_shard "$first" "$first_log"; then
  DOWNLOAD_STATE["$first"]=pass
else
  DOWNLOAD_STATE["$first"]=fail
  echo "[WARN] 최초 다운로드 일부 실패: $first" >&2
fi

failures=0
for ((i=0; i<${#SHARDS[@]}; i++)); do
  current=${SHARDS[$i]}
  download_status=${DOWNLOAD_STATE[$current]:-fail}
  next_pid=
  next=
  if (( i + 1 < ${#SHARDS[@]} )); then
    next=${SHARDS[$((i + 1))]}
    next_log="$SHARD_ROOT/queue_logs/$next.download.log"
    echo "[QUEUE] 계산=$current / 선다운로드=$next"
    download_shard "$next" "$next_log" &
    next_pid=$!
  else
    echo "[QUEUE] 계산=$current / 마지막 shard"
  fi

  verify_status=pass
  preprocess_status=pass
  "$SCRIPT_DIR/run_shard.sh" "$current" --shard-root "$SHARD_ROOT" \
    --data-root "$DATA_ROOT" --output-root "$OUTPUT_ROOT" --mode verify \
    >"$SHARD_ROOT/queue_logs/$current.verify.log" 2>&1 || verify_status=fail
  if [[ "$verify_status" == pass ]]; then
    "$SCRIPT_DIR/run_shard.sh" "$current" --shard-root "$SHARD_ROOT" \
      --data-root "$DATA_ROOT" --output-root "$OUTPUT_ROOT" --mode preprocess \
      >"$SHARD_ROOT/queue_logs/$current.preprocess.log" 2>&1 || preprocess_status=partial_fail
  else
    preprocess_status=skipped
  fi
  if [[ -n "$next_pid" ]]; then
    if wait "$next_pid"; then
      DOWNLOAD_STATE["$next"]=pass
    else
      DOWNLOAD_STATE["$next"]=fail
      echo "[WARN] 다음 shard 다운로드 일부 실패: $next" >&2
    fi
  fi
  printf '%s\t%s\t%s\t%s\n' "$current" "$download_status" "$verify_status" "$preprocess_status" >> "$STATUS.part"
  if [[ "$download_status" != pass || "$verify_status" != pass || "$preprocess_status" != pass ]]; then
    failures=$((failures + 1))
  fi
done
mv "$STATUS.part" "$STATUS"
echo "[QUEUE] 완료: shards=${#SHARDS[@]} failed=$failures status=$STATUS"
[[ "$failures" == 0 ]]
