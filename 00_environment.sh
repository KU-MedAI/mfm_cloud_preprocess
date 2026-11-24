#!/usr/bin/env bash
# 00 — integrated v2 실행 환경 설치/검증.
#
# 사용법:
#   ./00_environment.sh                # core + Clair3 설치/보강 후 도구 검증
#   ./00_environment.sh verify         # core + Clair3 + RNA + de novo + reference 검증
#   ./00_environment.sh verify-core    # jj_mfm/core tool 검증
#   ./00_environment.sh verify-clair3  # Clair3 전용 env 검증
#   ./00_environment.sh verify-rna     # RNA 관련 tool 검증(선택)
#   ./00_environment.sh verify-denovo  # hifiasm 등 de novo assembler 검증
#   ./00_environment.sh install-core   # jj_mfm core package 설치/보강
#   ./00_environment.sh install-clair3 # mamba 없이 conda-only Clair3 env 생성
#   ./00_environment.sh install-rna    # STAR/Picard/bedtools/pyBigWig 등 RNA tool 설치/보강
#   ./00_environment.sh install-denovo # hifiasm 설치/보강
#   ./00_environment.sh install-all    # core(RNA·de novo 포함) + Clair3 설치/보강
#   ./00_environment.sh print          # 현재 config 값 출력
#
# 정책:
# - integrated wrapper/backend는 jj_mfm을 기본으로 사용한다.
# - Clair3는 dependency 충돌을 줄이기 위해 별도 env(clair3_v2)를 사용한다.
# - mamba/micromamba는 사용하지 않는다.
# - 인수 없는 실행과 install* 모드만 conda env를 수정한다.
set -euo pipefail
cd "$(dirname "$0")"
if [[ -n "${CORE_ENV_NAME:-}" && -n "${MFM_CORE_ENV_NAME:-}" && "$CORE_ENV_NAME" != "$MFM_CORE_ENV_NAME" ]]; then
  echo '[FATAL] CORE_ENV_NAME과 MFM_CORE_ENV_NAME이 서로 다릅니다.' >&2
  exit 1
fi
export MFM_CORE_ENV_NAME=${MFM_CORE_ENV_NAME:-${CORE_ENV_NAME:-jj_mfm}}
source ./config_v2.sh
export PYTHONPATH="$PIPELINE_DIR/src:$MFM_REPO_ROOT:${PYTHONPATH:-}"

MODE=${1:-setup}
CORE_ENV_NAME=$MFM_CORE_ENV_NAME
STAR_BIN=${STAR_BIN:-$(find_tool STAR)}
BEDTOOLS_BIN=${BEDTOOLS_BIN:-$(find_tool bedtools)}

conda_env_exists() {
  "$CONDA_BIN" env list | awk '{print $1}' | grep -qx "$1"
}

require_conda() {
  [[ -n "$CONDA_BIN" ]] || die 'conda가 필요합니다. Conda/Miniforge를 먼저 설치하세요.'
}

env_prefix_for_name() {
  "$CONDA_BIN" env list --json | "$PYTHON_BIN" -c '
import json, os, sys
name = sys.argv[1]
print(next((p for p in json.load(sys.stdin)["envs"] if os.path.basename(p) == name), ""))
' "$1"
}

write_versions() {
  local prefix=$1 output=$2 selected=$3
  "$CONDA_BIN" list --prefix "$prefix" --export |
    awk -F= -v selected="$selected" 'BEGIN {print "package\tversion\tbuild"}
      $1 ~ selected {print $1 "\t" $2 "\t" $3}' > "$output.tmp"
  mv "$output.tmp" "$output"
}

refresh_core_paths() {
  local prefix=$1
  [[ -x "$prefix/bin/python" ]] || die "설치된 core Python을 찾지 못했습니다: $prefix"
  MFM_CORE_ENV=$prefix
  PYTHON_BIN=$prefix/bin/python
  export MFM_CORE_ENV PYTHON_BIN PATH="$prefix/bin:$PATH"
  BCFTOOLS_BIN=$prefix/bin/bcftools
  BWA_MEM2_BIN=$prefix/bin/bwa-mem2
  MINIMAP2_BIN=$prefix/bin/minimap2
  SAMTOOLS_BIN=$prefix/bin/samtools
  WHATSHAP_BIN=$prefix/bin/whatshap
  STAR_BIN=$prefix/bin/STAR
  PICARD_BIN=$prefix/bin/picard
  BEDTOOLS_BIN=$prefix/bin/bedtools
  CURL_BIN=$prefix/bin/curl
  HIFIASM_BIN=$prefix/bin/hifiasm
}

install_core() {
  local prefix
  require_conda
  log "Install/repair core env: $CORE_ENV_NAME"
  if conda_env_exists "$CORE_ENV_NAME"; then
    "$CONDA_BIN" env update -n "$CORE_ENV_NAME" -f "$PIPELINE_DIR/environment.yml"
  else
    "$CONDA_BIN" env create -n "$CORE_ENV_NAME" -f "$PIPELINE_DIR/environment.yml"
  fi
  prefix=$(env_prefix_for_name "$CORE_ENV_NAME")
  [[ -n "$prefix" ]] || die "설치된 core env 경로를 찾지 못했습니다: $CORE_ENV_NAME"
  refresh_core_paths "$prefix"
  "$PYTHON_BIN" -m pip install --no-deps -r "$PIPELINE_DIR/requirements.txt"
  write_versions "$prefix" "$PIPELINE_DIR/TOOL_VERSIONS.tsv" \
    '^(python|pip|numpy|pandas|pysam|zarr|numcodecs|biopython|pybigwig|bcftools|samtools|minimap2|bwa-mem2|whatshap|star|picard|openjdk|bedtools|subread|curl|hifiasm|pigz|htslib|openssl|cxx-compiler|cxx_linux-64|gxx_impl_linux-64|pytest)$'
}

install_clair3() {
  local prefix
  require_conda
  log "Install Clair3 with conda-only: $CLAIR3_ENV_NAME"
  if conda_env_exists "$CLAIR3_ENV_NAME"; then
    "$CONDA_BIN" env update -n "$CLAIR3_ENV_NAME" -f "$PIPELINE_DIR/clair3_environment.yml"
  else
    "$CONDA_BIN" env create -n "$CLAIR3_ENV_NAME" -f "$PIPELINE_DIR/clair3_environment.yml"
  fi
  prefix=$(env_prefix_for_name "$CLAIR3_ENV_NAME")
  [[ -n "$prefix" ]] || die "설치된 Clair3 env 경로를 찾지 못했습니다: $CLAIR3_ENV_NAME"
  CLAIR3_ENV=$prefix
  CLAIR3_PYTHON_BIN=$prefix/bin/python
  CLAIR3_PY=$prefix/bin/run_clair3.py
  CLAIR3_WHATSHAP_BIN=$prefix/bin/whatshap
  CLAIR3_SAMTOOLS_BIN=$prefix/bin/samtools
  CLAIR3_PARALLEL_BIN=$prefix/bin/parallel
  CLAIR3_PYPY_BIN=$prefix/bin/pypy3
  CLAIR3_TABIX_BIN=$prefix/bin/tabix
  CLAIR3_MODELS_DIR=$prefix/bin/models
  CLAIR3_MODEL_PATH_ONT=$CLAIR3_MODELS_DIR/ont
  CLAIR3_MODEL_PATH_HIFI=$CLAIR3_MODELS_DIR/hifi
  CLAIR3_MODEL_PATH=${CLAIR3_MODEL_PATH_OVERRIDE:-$CLAIR3_MODEL_PATH_ONT}
  write_versions "$prefix" "$PIPELINE_DIR/TOOL_VERSIONS_CLAIR3.tsv" \
    '^(python|clair3|whatshap|samtools|parallel|pypy3|htslib)$'
}

install_denovo() {
  install_core
}

install_rna() {
  install_core
}

verify_core() {
  local miss=0
  for item in \
    "$PYTHON_BIN" \
    "$RUNNER" \
    "$BCFTOOLS_BIN" \
    "$BWA_MEM2_BIN" \
    "$MINIMAP2_BIN" \
    "$SAMTOOLS_BIN" \
    "${CURL_BIN:-$(find_tool curl)}" \
    "$WHATSHAP_BIN"
  do
    if [[ ! -x "$item" ]]; then
      warn "missing core executable: $item"
      miss=1
    fi
  done
  "$PYTHON_BIN" - <<'PYCORE' >/dev/null || miss=1
import pysam, zarr, numcodecs, numpy, pandas
assert zarr.__version__.split('.')[0] == '2', 'Zarr 2 환경이 필요합니다.'
import pre_training.preprocess_raw_readset.config
PYCORE
  if (( miss != 0 )); then
    die "core environment verification failed"
  fi
  log "core Python: $($PYTHON_BIN --version 2>&1)"
  log "bcftools: $($BCFTOOLS_BIN --version 2>/dev/null | head -1 || true)"
  log "minimap2: $($MINIMAP2_BIN --version 2>/dev/null || true)"
}

verify_clair3() {
  local miss=0
  for item in \
    "$CLAIR3_PY" \
    "$CLAIR3_PYTHON_BIN" \
    "$CLAIR3_WHATSHAP_BIN" \
    "$CLAIR3_SAMTOOLS_BIN" \
    "$CLAIR3_PARALLEL_BIN" \
    "$CLAIR3_PYPY_BIN" \
    "$CLAIR3_TABIX_BIN"
  do
    if [[ ! -x "$item" ]]; then
      warn "missing Clair3 executable: $item"
      miss=1
    fi
  done
  [[ -d "$CLAIR3_MODELS_DIR" ]] || { warn "missing Clair3 models: $CLAIR3_MODELS_DIR"; miss=1; }
  [[ -d "$CLAIR3_MODEL_PATH_ONT" ]] || { warn "missing ONT Clair3 model: $CLAIR3_MODEL_PATH_ONT"; miss=1; }
  [[ -d "$CLAIR3_MODEL_PATH_HIFI" ]] || { warn "missing HiFi Clair3 model: $CLAIR3_MODEL_PATH_HIFI"; miss=1; }
  if (( miss != 0 )); then
    die "Clair3 environment verification failed"
  fi
  log "Clair3 version: $($CLAIR3_PYTHON_BIN "$CLAIR3_PY" --version 2>/dev/null || true)"
  log "WhatsHap version: $($CLAIR3_WHATSHAP_BIN --version 2>/dev/null || true)"
  log "samtools version: $($CLAIR3_SAMTOOLS_BIN --version 2>/dev/null | head -1 || true)"
  log "ONT Clair3 model: $CLAIR3_MODEL_PATH_ONT"
  log "HiFi Clair3 model: $CLAIR3_MODEL_PATH_HIFI"
  case "$(normalize_platform "$DEFAULT_PLATFORM")" in
    ont|hifi) log "selected Clair3 model for DEFAULT_PLATFORM=$DEFAULT_PLATFORM: $(clair3_model_path_for_platform "$DEFAULT_PLATFORM")" ;;
    clr) log "PacBio CLR: Clair3 model selection skipped; production requires vcf_source=public" ;;
    *) log "selected Clair3 model: skipped because DEFAULT_PLATFORM=$DEFAULT_PLATFORM is not long-read" ;;
  esac
}

verify_denovo() {
  local miss=0
  if [[ ! -x "$HIFIASM_BIN" ]]; then
    warn "missing hifiasm executable: $HIFIASM_BIN"
    miss=1
  fi
  if (( miss != 0 )); then
    die "de novo assembly environment verification failed"
  fi
  log "hifiasm: $($HIFIASM_BIN --version 2>/dev/null || true)"
}

verify_rna() {
  local miss=0
  for item in "$STAR_BIN" "$PICARD_BIN" "$BEDTOOLS_BIN" "$SAMTOOLS_BIN"; do
    if [[ ! -x "$item" ]]; then
      warn "missing RNA executable: $item"
      miss=1
    fi
  done
  "$PYTHON_BIN" - <<'PYRNA' >/dev/null || miss=1
import pyBigWig, numpy, zarr
PYRNA
  if (( miss != 0 )); then
    die "RNA environment verification failed"
  fi
  [[ "$($STAR_BIN --version 2>/dev/null)" == 2.7.11b ]] || die 'RNA profile은 STAR 2.7.11b를 요구합니다.'
  log "STAR: $($STAR_BIN --version 2>/dev/null || true)"
  log "Picard: $($PICARD_BIN MarkDuplicates --version 2>&1 | head -1 || true)"
  log "bedtools: $($BEDTOOLS_BIN --version 2>/dev/null | head -1 || true)"
}

verify_references() {
  local miss=0
  [[ -s "$REF_DECODE" ]] || { warn "missing REF_DECODE: $REF_DECODE"; miss=1; }
  [[ -s "$REF_SEQUENCE" ]] || { warn "missing REF_SEQUENCE: $REF_SEQUENCE"; miss=1; }
  if (( miss != 0 )); then
    die "reference verification failed"
  fi
}

case "$MODE" in
  install-gpu) "$PYTHON_BIN" "$PIPELINE_DIR/tools/setup_gpu.py" ;;
  install-gpu-tools) "$PYTHON_BIN" "$PIPELINE_DIR/tools/setup_gpu.py" --tools ;;
  verify-gpu) "$PYTHON_BIN" "$PIPELINE_DIR/tools/setup_gpu.py" --verify ;;
  setup|install|install-all)
    install_core
    install_clair3
    verify_core
    verify_clair3
    verify_rna
    verify_denovo
    log '도구 설치·검증 통과. 샘플별 참조/VCF 자산은 별도 확인하세요.'
    ;;
  install-core)
    install_core
    verify_core
    ;;
  install-clair3)
    install_clair3
    verify_clair3
    ;;
  install-rna)
    install_rna
    verify_rna
    ;;
  install-denovo)
    install_denovo
    verify_denovo
    ;;
  verify)
    verify_core
    verify_clair3
    verify_denovo
    verify_rna
    verify_references
    ;;
  verify-core)
    verify_core
    ;;
  verify-clair3)
    verify_clair3
    ;;
  verify-rna)
    verify_rna
    ;;
  verify-denovo)
    verify_denovo
    ;;
  print)
    print_config
    ;;
  *)
    die "unknown mode: $MODE (use setup|install|install-core|install-clair3|install-rna|install-denovo|install-all|verify|verify-core|verify-clair3|verify-rna|verify-denovo|print)"
    ;;
esac
