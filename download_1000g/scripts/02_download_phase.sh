#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash scripts/02_download_phase.sh <phase> [all|shard_id]

Download URL lists produced by 01_make_phase_url_lists.sh.

Environment variables:
  DOWNLOADER=aria2c          aria2c | wget
  DRY_RUN=0                  If 1, print what would be downloaded.
  DOWNLOAD_ROOT=downloads    Output directory under package root.
  OUT_ROOT=download_lists    URL list directory under package root.
  ARIA2_CONNECTIONS=4        aria2c -x/-s value.
  ARIA2_JOBS=8               aria2c parallel jobs.
  WGET_OPTS=-c               wget options.

Examples:
  DRY_RUN=1 bash scripts/02_download_phase.sh 1
  bash scripts/02_download_phase.sh 1 shard_001
  DOWNLOADER=wget bash scripts/02_download_phase.sh 2 all
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
OUT_ROOT="${OUT_ROOT:-$PKG_ROOT/download_lists}"
DOWNLOAD_ROOT="${DOWNLOAD_ROOT:-$PKG_ROOT/downloads}"
DOWNLOADER="${DOWNLOADER:-aria2c}"
DRY_RUN="${DRY_RUN:-0}"
PHASE_DIR="$OUT_ROOT/phase_$PHASE"
SHARD_DIR="$PHASE_DIR/shards"

if [[ ! -f "$PHASE_DIR/all.urls" ]]; then
  echo "[INFO] URL list not found; generating phase $PHASE lists"
  bash "$SCRIPT_DIR/01_make_phase_url_lists.sh" "$PHASE"
fi

lists=()
if [[ "$SHARD_ARG" == "all" ]]; then
  while IFS= read -r f; do lists+=("$f"); done < <(find "$SHARD_DIR" -name 'shard_*.urls' | sort)
else
  shard="${SHARD_ARG#shard_}"
  shard="$(printf "%03d" "${shard#0}")"
  lists+=("$SHARD_DIR/shard_$shard.urls")
fi

for list in "${lists[@]}"; do
  if [[ ! -s "$list" ]]; then
    echo "[WARN] empty or missing URL list: $list" >&2
    continue
  fi
  shard_name="$(basename "$list" .urls)"
  dest="$DOWNLOAD_ROOT/phase_$PHASE/$shard_name"
  mkdir -p "$dest"
  url_count="$(wc -l < "$list")"
  echo "[INFO] phase=$PHASE shard=$shard_name urls=$url_count dest=$dest"
  if [[ "$DRY_RUN" == "1" ]]; then
    sed -n '1,10p' "$list"
    [[ "$url_count" -gt 10 ]] && echo "[DRY_RUN] ... ($url_count URLs total)"
    continue
  fi
  case "$DOWNLOADER" in
    aria2c)
      command -v aria2c >/dev/null 2>&1 || { echo "[ERROR] aria2c not found; set DOWNLOADER=wget" >&2; exit 1; }
      aria2c -x "${ARIA2_CONNECTIONS:-4}" -s "${ARIA2_CONNECTIONS:-4}" -j "${ARIA2_JOBS:-8}" -c --auto-file-renaming=false --allow-overwrite=false -i "$list" -d "$dest"
      ;;
    wget)
      command -v wget >/dev/null 2>&1 || { echo "[ERROR] wget not found" >&2; exit 1; }
      wget ${WGET_OPTS:--c} -i "$list" -P "$dest"
      ;;
    *) echo "[ERROR] unsupported DOWNLOADER=$DOWNLOADER" >&2; exit 1 ;;
  esac
done
