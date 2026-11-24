#!/usr/bin/env bash
# 호환 wrapper: window 생성의 production 번호는 11이다.
set -euo pipefail
cd "$(dirname "$0")"
echo "[WARN] 10_windows.sh는 deprecated입니다. 11_windows.sh를 실행합니다." >&2
exec ./11_windows.sh "$@"
