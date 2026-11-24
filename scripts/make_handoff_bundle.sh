#!/usr/bin/env bash
# 검증한 클라우드 코드 inventory만 내보냄. 실행 데이터/캐시/설치 바이너리는 제외함.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN=${PYTHON_BIN:-$(command -v python3 || command -v python)}
exec "$PYTHON_BIN" "$ROOT/tools/check_cloud_package.py" --export-to "${1:?사용법: make_handoff_bundle.sh 새_출력_경로}"
