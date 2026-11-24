#!/usr/bin/env bash
# 리드별 목록을 기존 배치 단위로 실행함. run은 jj_mfm 환경 확인/설치를 기존 진입점에 위임함.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export PYTHONPATH="$PWD:$PWD/src:${PYTHONPATH:-}"
if [[ -n "${MFM_CORE_ENV:-}" && -x "$MFM_CORE_ENV/bin/python" ]]; then
  export PYTHON_BIN="$MFM_CORE_ENV/bin/python"
fi
exec "${PYTHON_BIN:-python3}" -B -m handoff.read_groups "$@"
