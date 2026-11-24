#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash scripts/03_verify_phase.sh <phase> [all|shard_id]

Verify downloaded CRAM/BAM files for one phase/shard.

Environment variables:
  DOWNLOAD_ROOT=downloads
  RUN_SAMTOOLS=1          Run samtools quickcheck when samtools is available.

Examples:
  bash scripts/03_verify_phase.sh 1
  bash scripts/03_verify_phase.sh 1 shard_001
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 1 ]]; then
  usage
  exit $([[ $# -lt 1 ]] && echo 1 || echo 0)
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PHASE="$(printf "%02d" "${1#0}")"
SHARD_ARG="${2:-all}"
DOWNLOAD_ROOT="${DOWNLOAD_ROOT:-$PKG_ROOT/downloads}"
RUN_SAMTOOLS="${RUN_SAMTOOLS:-1}"
PHASE_ROOT="$DOWNLOAD_ROOT/phase_$PHASE"

if [[ ! -d "$PHASE_ROOT" ]]; then
  echo "[ERROR] missing phase download dir: $PHASE_ROOT" >&2
  exit 1
fi

dirs=()
if [[ "$SHARD_ARG" == "all" ]]; then
  while IFS= read -r d; do dirs+=("$d"); done < <(find "$PHASE_ROOT" -mindepth 1 -maxdepth 1 -type d | sort)
else
  shard="${SHARD_ARG#shard_}"
  shard="$(printf "%03d" "${shard#0}")"
  dirs+=("$PHASE_ROOT/shard_$shard")
fi

status=0
for dir in "${dirs[@]}"; do
  [[ -d "$dir" ]] || { echo "[WARN] missing shard dir: $dir"; status=1; continue; }
  echo "[INFO] verifying $dir"
  echo "[INFO] file count: $(find "$dir" -type f | wc -l)"

  while IFS= read -r cram; do
    if [[ ! -f "${cram}.crai" && ! -f "${cram%.cram}.crai" ]]; then
      echo "[MISSING_INDEX] $cram" >&2
      status=1
    fi
  done < <(find "$dir" -type f -name '*.cram' | sort)

  while IFS= read -r bam; do
    if [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" || -f "${bam}.csi" || -f "${bam%.bam}.csi" ]]; then
      true
    else
      echo "[WARN_NO_BAM_INDEX] $bam" >&2
    fi
  done < <(find "$dir" -type f -name '*.bam' | sort)

  if [[ "$RUN_SAMTOOLS" == "1" ]] && command -v samtools >/dev/null 2>&1; then
    find "$dir" -type f \( -name '*.cram' -o -name '*.bam' \) -print0 \
      | xargs -0 -r -n 64 samtools quickcheck -v || status=1
  else
    echo "[INFO] samtools quickcheck skipped"
  fi
done

exit "$status"
