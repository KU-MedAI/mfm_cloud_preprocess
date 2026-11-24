#!/usr/bin/env bash
# 구번호 Stage 19 직접 실행기. 현재 RNA runner의 metadata 단계는 Stage 24다.
# 19 — RNA 산출물 경로와 완료 여부를 sample.json에 반영한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
require_manifest_sample_contract "$S"
F=$(v2_final "$S")
W=$(v2_work "$S")
mkdir -p "$W/logs"
"$PYTHON" ./19_stamp_rna_metadata.py "$F" 2>&1 | tee -a "$W/logs/19.log"
log "19 완료: $S"
