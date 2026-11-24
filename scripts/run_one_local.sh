#!/usr/bin/env bash
# runner_local.sh가 호출하는 sample별 오류 격리 wrapper.
set -u

sample=${1:?sample을 지정하세요}
target=${RUNNER_TARGET_SCRIPT:?RUNNER_TARGET_SCRIPT를 지정하세요}
status_dir=${RUN_STATUS_DIR:?RUN_STATUS_DIR를 지정하세요}
log_dir=${RUN_LOG_DIR:?RUN_LOG_DIR를 지정하세요}
attempt=${RUNNER_ATTEMPT:-0}
mkdir -p "$status_dir/results" "$log_dir"
log_path="$log_dir/$sample.attempt_${attempt}.log"
resource_path="$log_dir/$sample.attempt_${attempt}.resource.json"
resource_runner="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/resource_runner.py"
python_bin=${PYTHON_BIN:-$(command -v python3 || command -v python || true)}
if [[ -z "$python_bin" || ! -x "$python_bin" ]]; then
  echo "[FATAL] PYTHON_BIN is not executable and no python3/python fallback was found" >&2
  exit 2
fi
result="$status_dir/results/$sample.result.tsv"
started=$(date +%s)

set +e
"$python_bin" "$resource_runner" \
  --output "$resource_path" -- "$target" "$sample" >"$log_path" 2>&1
rc=$?
set -e
ended=$(date +%s)
if [[ "$rc" == 0 ]]; then status=success; else status=failed; fi
read -r max_rss_kib user_cpu_sec system_cpu_sec < <(
  "$python_bin" -c 'import json,sys; d=json.load(open(sys.argv[1])); print(d["max_rss_kib"], d["user_cpu_sec"], d["system_cpu_sec"])' "$resource_path"
)
tmp="$result.part.$$"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$sample" "$status" "$rc" "$((ended-started))" "$attempt" "$log_path" \
  "${max_rss_kib:-}" "${user_cpu_sec:-}" "${system_cpu_sec:-}" "$resource_path" > "$tmp"
mv "$tmp" "$result"

# 한 sample 실패가 xargs와 다른 sample 실행을 중단시키지 않게 wrapper 자체는 성공으로 끝낸다.
exit 0
