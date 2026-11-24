#!/usr/bin/env bash
# jj_mfm이 없으면 conda/Miniforge와 환경을 설치한다. 기존 환경은 필수 도구를 검증한다.
handoff_environment() {
  local env_prefix=${MFM_CORE_ENV:-${CONDA_PREFIX:-}} conda_cmd=${CONDA_EXE:-} base bootstrap version asset
  if [[ -n "$env_prefix" && "$(basename "$env_prefix")" == jj_mfm && -x "$env_prefix/bin/python" ]]; then
    :
  else
    [[ -n "$conda_cmd" ]] || conda_cmd=$(command -v conda || true)
    if [[ -z "$conda_cmd" ]]; then
      base=${MFM_CONDA_ROOT:-$PWD/.conda_miniforge}
      if [[ ! -x "$base/bin/conda" ]]; then
        [[ "$(uname -s)" == Linux ]] || { echo '[오류] 자동 설치는 Linux 서버를 지원합니다.' >&2; return 1; }
        version=26.7.2-0;asset="Miniforge3-$version-Linux-$(uname -m).sh"
        bootstrap=$PWD/.bootstrap;mkdir -p "$bootstrap"
        curl --fail --location --retry 5 --output "$bootstrap/$asset" "https://github.com/conda-forge/miniforge/releases/download/$version/$asset"
        curl --fail --location --retry 5 --output "$bootstrap/$asset.sha256" "https://github.com/conda-forge/miniforge/releases/download/$version/$asset.sha256"
        (cd "$bootstrap" && sha256sum --check "$asset.sha256")
        bash "$bootstrap/$asset" -b -p "$base"
      fi
      conda_cmd=$base/bin/conda
    fi
    base=$("$conda_cmd" info --base)
    source "$base/etc/profile.d/conda.sh"
    if ! conda env list | awk '{print $1}' | grep -qx jj_mfm; then
      # 기본 solver는 conda libmamba를 사용하지만 별도 mamba 명령은 필요 없다.
      conda env create -n jj_mfm -f "$PWD/environment.yml"
    fi
    conda activate jj_mfm
    env_prefix=$CONDA_PREFIX
  fi
  export MFM_CORE_ENV=$env_prefix PYTHON_BIN=$env_prefix/bin/python
  export PATH="$env_prefix/bin:$PATH"
  export PYTHONPATH="$PWD:$PWD/src:$PWD/vendor/MFM-DNARNAFM:${PYTHONPATH:-}"
  "$PYTHON_BIN" - <<'PY'
import pysam,zarr,numcodecs,numpy,pandas
assert zarr.__version__.split('.')[0]=='2', 'Zarr 2 환경이 필요합니다.'
import pre_training.preprocess_raw_readset.config
PY
  local tool
  for tool in samtools bcftools minimap2 bwa-mem2 whatshap STAR picard bedtools curl; do
    command -v "$tool" >/dev/null || { echo "[오류] jj_mfm 필수 도구 누락: $tool. conda env update -n jj_mfm -f environment.yml 로 보강하세요." >&2; return 1; }
  done
  [[ "$(STAR --version)" == 2.7.11b ]] || { echo '[오류] RNA profile은 STAR2.7.11b를 요구합니다. environment.yml의 버전을 사용하세요.' >&2; return 1; }
}
