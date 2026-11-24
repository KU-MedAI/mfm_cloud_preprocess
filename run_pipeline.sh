#!/usr/bin/env bash
# 한 번의 호출로 환경 확인/설치, 선택 배치 다운로드, 전처리, 최종 정리를 실행한다.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
if (($# == 0)); then set -- run; fi
MODE=${1:-help}
# native/NumPy 보조 thread는 기본1로 두고 단계별 병렬화 예산을 따름.
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}
# 설치 경로를 가정하지 않고 명시 환경 또는 현재 PATH에서 Python을 선택함.
if [[ -n "${MFM_CORE_ENV:-}" && -x "$MFM_CORE_ENV/bin/python" ]]; then
  export PYTHON_BIN="$MFM_CORE_ENV/bin/python"
elif [[ -n "${CONDA_PREFIX:-}" && -x "$CONDA_PREFIX/bin/python" ]]; then
  export PYTHON_BIN="$CONDA_PREFIX/bin/python"
else
  export PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || command -v python || true)}"
fi
if [[ "$MODE" == help || "$MODE" == --help || "$MODE" == -h ]]; then
  printf '%s\n' '사용법: bash run_pipeline.sh dry-run|run --lab tgil|medai --batch 001 [--executor local|slurm] [--config 파일]' \
    '클라우드 분담: bash run_pipeline.sh dry-run|run --target azure_poc|azure_main|aws_main|server [--sample ID]' \
    '클라우드 입력 계획: bash run_pipeline.sh download --target azure_poc --dry-run' \
    '분담 목록 재생성: bash run_pipeline.sh plan' '클러스터/경로 확인: bash run_pipeline.sh doctor' \
    '샘플 선택: bash run_pipeline.sh run --lab tgil --sample ID --config handoff_config.json' \
    '다운로드만: bash run_pipeline.sh download --lab tgil --batch 003 --executor local|slurm --config 파일 [--dry-run]' \
    '단일 다운로드: bash run_pipeline.sh download --lab medai --sample HG00281 --config 파일' \
    '자원 override: --resource-profile auto|throughput|balanced|conservative --elastic-max-jobs 6' \
    '전체 CPU 상한: --cpu-pool-limit 68 / 샘플 재개: resume --sample ID --from-stage 01' \
    '직접 지정: --local-jobs 4 --sample-cpus 16 --download-jobs 4 [--slurm-concurrency 4 --slurm-cpus 16]'
  exit 0
fi
# 과거의 직접 branch/목록/Slurm/smoke 명령은 원본의 호환 dispatcher에 보존함.
case "$MODE" in
  run|dry-run|doctor|plan|resume|import-state|download|download-worker|download-finish|prepare|worker|finish) ;;
  *) exec bash "$PWD/scripts/run_compat_pipeline.sh" "$@" ;;
esac
# 다운로드는 Python/curl만 사용하고 전처리 환경을 설치·검증하지 않음.
if [[ "$MODE" == download || "$MODE" == download-worker || "$MODE" == download-finish ]]; then
  if [[ -n "${MFM_CORE_ENV:-}" && -x "$MFM_CORE_ENV/bin/python" ]]; then
    export PYTHON_BIN="$MFM_CORE_ENV/bin/python"
    export PATH="$MFM_CORE_ENV/bin:$PATH"
  fi
  export PYTHONPATH="$PWD:$PWD/src:$PWD/vendor/MFM-DNARNAFM:${PYTHONPATH:-}"
  exec "${PYTHON_BIN:-python3}" -m handoff.cli "$@"
fi
if [[ "$MODE" == run || "$MODE" == prepare || "$MODE" == worker || "$MODE" == finish ]]; then
  if [[ "$MODE" == run ]]; then
    # lab/batch 오류는 설치나 제출 이전에 판정한다.
    export PYTHONPATH="$PWD:$PWD/src:$PWD/vendor/MFM-DNARNAFM:${PYTHONPATH:-}"
    "$PYTHON_BIN" -m handoff.cli dry-run "${@:2}" >/dev/null
  fi
  source ./handoff/environment.sh
  handoff_environment
fi
export PYTHONPATH="$PWD:$PWD/src:$PWD/vendor/MFM-DNARNAFM:${PYTHONPATH:-}"
exec "${PYTHON_BIN:-python3}" -m handoff.cli "$@"
