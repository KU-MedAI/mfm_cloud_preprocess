#!/usr/bin/env bash
# 호환 wrapper: personal realignment의 production 번호는 10이다.
set -euo pipefail
cd "$(dirname "$0")"
echo "[WARN] 11_realign_wgs.sh는 deprecated입니다. 10_realign_wgs.sh를 실행합니다." >&2
exec ./10_realign_wgs.sh "$@"
