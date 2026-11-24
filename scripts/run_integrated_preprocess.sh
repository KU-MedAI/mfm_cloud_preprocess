#!/usr/bin/env bash
set -euo pipefail
# native/NumPy 보조 thread는 기본1로 두고 단계별 병렬화 예산을 따름.
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-${CONDA_PREFIX:+$CONDA_PREFIX/bin/python}}"
if [[ -z "${PYTHON_BIN:-}" && -n "${MFM_CORE_ENV:-}" ]]; then
  PYTHON_BIN="$MFM_CORE_ENV/bin/python"
fi
if [[ -z "${MFM_REPO_ROOT:-}" ]]; then
  for candidate in     "$PKG_ROOT/vendor/MFM-DNARNAFM"
  do
    if [[ -d "$candidate/pre_training/preprocess_raw_readset" ]]; then
      MFM_REPO_ROOT=$candidate
      break
    fi
  done
  MFM_REPO_ROOT=${MFM_REPO_ROOT:-$PKG_ROOT/vendor/MFM-DNARNAFM}
fi
[[ -x "$PYTHON_BIN" ]] || PYTHON_BIN=$(command -v python || command -v python3 || true)
[[ -n "$PYTHON_BIN" && -x "$PYTHON_BIN" ]] || { echo "[FATAL] python 실행 파일을 찾지 못했습니다" >&2; exit 1; }
"$PYTHON_BIN" -c 'import pysam, zarr, numcodecs' >/dev/null 2>&1 || {
  echo "[FATAL] core Python 모듈(pysam/zarr/numcodecs)을 불러올 수 없습니다: $PYTHON_BIN" >&2
  echo "[FATAL] MFM_CORE_ENV 또는 PYTHON_BIN을 jj_mfm 환경으로 지정하세요" >&2
  exit 1
}

# compact writer가 참조하는 pre_training.preprocess_raw_readset은 전달본의
# vendor/MFM-DNARNAFM 또는 사용자가 지정한 MFM_REPO_ROOT에서 불러온다.
if [[ -d "$MFM_REPO_ROOT/pre_training" ]]; then
  export PYTHONPATH="${PKG_ROOT}/src:${MFM_REPO_ROOT}:${PYTHONPATH:-}"
else
  export PYTHONPATH="${PKG_ROOT}/src:${PYTHONPATH:-}"
fi

exec "${PYTHON_BIN}" -m integrate_wgs_rna.cli "$@"
