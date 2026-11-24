#!/usr/bin/env bash
# 테스트 전체 실행.
#
#   ./tests/run_all.sh
#
# 이 러너가 없던 동안 테스트는 "각자 알아서 돌리는" 상태였다. python 테스트 중 일부만
# PYTHONPATH=src 를 요구해서(src 패키지를 import 하는 것들) 처음 돌리는 사람이 반드시 밟았다.
# 여기서 한 번에 맞춰 준다. pytest 는 필요 없다 — 전부 standalone 이다.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
PIPELINE_DIR=$(pwd)

# 이 프로젝트의 검증 기준 환경은 jj_mfm이다. 호출 셸이 base conda를
# 활성화한 상태여도 그 CONDA_PREFIX를 암묵적으로 채택하지 않는다.
if [[ -z "${MFM_CORE_ENV:-}" && -x /home/aibio_1/miniconda3/envs/jj_mfm/bin/python ]]; then
  MFM_CORE_ENV=/home/aibio_1/miniconda3/envs/jj_mfm
fi
MFM_CORE_ENV=${MFM_CORE_ENV:-${CONDA_PREFIX:-}}
PYTHON_BIN=${PYTHON_BIN:-${MFM_CORE_ENV:+$MFM_CORE_ENV/bin/python}}
export MFM_CORE_ENV PYTHON_BIN
[[ -x "$PYTHON_BIN" ]] || PYTHON_BIN=$(command -v python || command -v python3 || true)
[[ -x "$PYTHON_BIN" ]] || { echo "python을 찾을 수 없습니다. PYTHON_BIN을 지정하세요." >&2; exit 2; }
echo "python: $PYTHON_BIN"

if [[ -z "${MFM_REPO_ROOT:-}" ]]; then
  for candidate in     "$PIPELINE_DIR/vendor/MFM-DNARNAFM"     "$(readlink -m "$PIPELINE_DIR/../../../pre-training/jungjin/MFM-DNARNAFM")"
  do
    if [[ -d "$candidate/pre_training/preprocess_raw_readset" ]]; then MFM_REPO_ROOT=$candidate; break; fi
  done
fi
MFM_REPO_ROOT=${MFM_REPO_ROOT:-$PIPELINE_DIR/vendor/MFM-DNARNAFM}
[[ -d "$MFM_REPO_ROOT/pre_training/preprocess_raw_readset" ]] || { echo "MFM_REPO_ROOT가 올바르지 않습니다: $MFM_REPO_ROOT" >&2; exit 2; }
export PATH="$(dirname "$PYTHON_BIN"):$PATH"
# 단계 파일 일부는 v2 루트의 공통 모듈(rna_cnv_policy 등)을 직접 import한다.
export PYTHONPATH="$PIPELINE_DIR:$PIPELINE_DIR/src:$MFM_REPO_ROOT"

only=${1:-}
pass=0; fail=0; failed=()

run_one() {
  local label=$1; shift
  [[ -n "$only" && "$label" != *"$only"* ]] && return 0
  printf '  %-36s ' "$label"
  local out
  if out=$("$@" 2>&1); then
    echo "PASS"; pass=$((pass + 1))
  else
    echo "FAIL"; fail=$((fail + 1)); failed+=("$label")
    printf '%s\n' "$out" | tail -12 | sed 's/^/        /'
  fi
}

echo "python (jj_mfm; pipeline src + MFM repo PYTHONPATH)"
for f in tests/*.py; do run_one "$(basename "$f")" "$PYTHON_BIN" "$f"; done
echo "shell"
for f in tests/*.sh; do
  [[ "$(basename "$f")" == "run_all.sh" ]] && continue
  run_one "$(basename "$f")" timeout 600 bash "$f"
done

echo "─────────────────────────────────────────"
printf '%d PASS / %d FAIL\n' "$pass" "$fail"
if (( fail )); then
  printf '실패: %s\n' "${failed[*]}"
  exit 1
fi
