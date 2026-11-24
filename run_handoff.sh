#!/usr/bin/env bash
# original 번호와 Slurm 호출의 호환 진입점임. 실제 처리는 공용 run_pipeline.sh에 위임함.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$ROOT/run_pipeline.sh" "$@"
