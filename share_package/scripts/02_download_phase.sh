#!/usr/bin/env bash
set -euo pipefail

[[ $# -ge 1 ]] || { echo "Usage: bash scripts/02_download_phase.sh <phase>" >&2; exit 1; }
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PHASE="$(printf '%02d' "${1#0}")"
OUT_ROOT="${OUT_ROOT:-$PKG_ROOT/download_lists}"
DOWNLOAD_ROOT="${DOWNLOAD_ROOT:-$PKG_ROOT/downloads}"
DOWNLOADER="${DOWNLOADER:-aria2c}"
DRY_RUN="${DRY_RUN:-0}"
LIST="$OUT_ROOT/phase_$PHASE/all.urls"
DEST="$DOWNLOAD_ROOT/phase_$PHASE"

[[ -s "$LIST" ]] || {
  echo "[ERROR] missing URL list: $LIST" >&2
  echo "[INFO] current safe example: PLATFORM_FILTER=Illumina_30x bash scripts/01_make_phase_url_lists.sh $PHASE" >&2
  exit 1
}
mkdir -p "$DEST"
echo "[INFO] phase=$PHASE urls=$(wc -l < "$LIST") dest=$DEST"
[[ "$DRY_RUN" == "1" ]] && { sed -n '1,20p' "$LIST"; exit 0; }

case "$DOWNLOADER" in
  aria2c)
    command -v aria2c >/dev/null 2>&1 || { echo "[ERROR] aria2c not found" >&2; exit 1; }
    aria2c -x "${ARIA2_CONNECTIONS:-4}" -s "${ARIA2_CONNECTIONS:-4}" -j "${ARIA2_JOBS:-8}" -c --auto-file-renaming=false --allow-overwrite=false -i "$LIST" -d "$DEST"
    ;;
  wget)
    command -v wget >/dev/null 2>&1 || { echo "[ERROR] wget not found" >&2; exit 1; }
    wget -c -i "$LIST" -P "$DEST"
    ;;
  *) echo "[ERROR] unsupported DOWNLOADER=$DOWNLOADER" >&2; exit 1;;
esac

