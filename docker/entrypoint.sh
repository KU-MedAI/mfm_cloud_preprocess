#!/usr/bin/env bash
# 환경 재설치 없이 기존 파이프라인으로 인수를 전달함.
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export HOME=${HOME:-/scratch/home}
export TMPDIR=${TMPDIR:-/scratch/tmp}
export XDG_CACHE_HOME=${XDG_CACHE_HOME:-/cache}
export MFM_NATIVE_CACHE_ROOT=${MFM_NATIVE_CACHE_ROOT:-$XDG_CACHE_HOME/native}
mkdir -p -- "$HOME" "$TMPDIR" "$XDG_CACHE_HOME" "$MFM_NATIVE_CACHE_ROOT"

# native_tool과 같은 파일 잠금으로 기존 cache를 보존하며 빌드 결과를 재사용함.
seed_root=${MFM_IMAGE_NATIVE_CACHE:-/opt/mfm-image/native}
for source in "$seed_root"/*; do
  [[ -f "$source" && "$source" != *.lock ]] || continue
  target=$MFM_NATIVE_CACHE_ROOT/${source##*/}
  (
    flock -x 9
    [[ ! -e "$target" ]] || exit 0
    temporary=$(mktemp "$MFM_NATIVE_CACHE_ROOT/.seed.XXXXXXXX")
    trap 'rm -f -- "$temporary"' EXIT
    cp -- "$source" "$temporary"
    chmod 0700 "$temporary"
    mv -- "$temporary" "$target"
  ) 9>"$target.lock"
done
if (($# == 0)); then set -- help; fi
exec bash "$ROOT/run_pipeline.sh" "$@"
