#!/usr/bin/env bash
# array index는 선택 manifest의 0-based 행이다. 샘플마다 독립 실행/상태 파일을 만든다.
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CONFIG=${1:?effective_config.json 경로가 필요합니다.}
SELECTION=${2:?선택 manifest 경로가 필요합니다.}
TASK=${SLURM_ARRAY_TASK_ID:?Slurm array task ID가 필요합니다.}
exec bash "$ROOT/run_handoff.sh" worker --effective-config "$CONFIG" --selection "$SELECTION" --task-index "$TASK" --executor slurm
