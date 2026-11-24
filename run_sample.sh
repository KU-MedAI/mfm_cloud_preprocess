#!/usr/bin/env bash
# manifest의 RNA 입력 유무에 따라 WGS-only 또는 WGS→RNA를 한 번에 실행한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./lib/fixed_cpu_entry.sh
pipeline_fixed_cpu_entry "${1:?sample이 필요함}" bash "$(pwd)/run_sample.sh" "$@"
source ./config_v2.sh
source ./lib/state.sh

S=${1:?사용법: $0 <sample>}
apply_sample_payload "$S"
require_manifest_sample_contract "$S"
RNA=$(sample_rna "$S" 2>/dev/null || true)
STATUS_ROOT=${RUN_STATUS_DIR:-$(sample_status_root "$S")}
mkdir -p "$STATUS_ROOT"
PIPELINE_SAMPLE_ID=$S
INTEGRATED_STATE_FILE=$STATUS_ROOT/$S.integrated.state.tsv
CURRENT_BRANCH=wgs
HAS_RNA=0
if [[ -n "$RNA" && "$RNA" != "-" ]] && ! sample_rna_deferred "$S"; then HAS_RNA=1; fi
pipeline_install_integrated_trap

if [[ "$HAS_RNA" == 1 ]]; then
  log "통합 실행: $S WGS–RNA pair; WGS 01–12 후 RNA 13–26, 검증 25, cleanup 26 실행"
  # pair는 RNA 산출물까지 포함한 package를 한 번에 검증한 뒤 정리한다.
  WGS_DEFER_FINALIZE=1 ./run_wgs_sample.sh "$S"
  CURRENT_BRANCH=rna
  pipeline_write_integrated_state running 0
  ./run_rna_sample.sh "$S"
else
  if sample_rna_deferred "$S"; then log "RNA 보류: $S DNA 서열/window/QC만 실행하며 WGS/RNA read 처리는 생략함"; fi
  log "통합 실행: $S WGS-only; WGS 01–12 실행"
  ./run_wgs_sample.sh "$S"
fi

pipeline_complete_integrated
log "통합 실행 완료: sample=$S has_rna=$HAS_RNA"
