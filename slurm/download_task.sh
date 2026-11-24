#!/usr/bin/env bash
# 선택 manifest의 한 sample 입력만 내려받고 전처리 없이 종료함.
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CONFIG=${1:?다운로드 effective_config.json 경로가 필요함}
SELECTION=${2:?선택 manifest 경로가 필요함}
TASK=${SLURM_ARRAY_TASK_ID:?Slurm array task ID가 필요함}
exec bash "$ROOT/run_handoff.sh" download-worker --effective-config "$CONFIG" --selection "$SELECTION" --task-index "$TASK"
