#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash scripts/04_download_common_vcf.sh

Download 2022 1kGP high-coverage phased SNV/INDEL/SV panel VCF.

Environment variables:
  CHROMS="1-22"             Chromosomes to download. Examples: "1-22", "1 2 22".
  OUT_DIR=downloads/common/vcf
  DRY_RUN=0
  DOWNLOADER=aria2c         aria2c | wget
  ARIA2_CONNECTIONS=4
  ARIA2_JOBS=8
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CHROMS="${CHROMS:-1-22}"
OUT_DIR="${OUT_DIR:-$PKG_ROOT/downloads/common/vcf}"
DOWNLOADER="${DOWNLOADER:-aria2c}"
DRY_RUN="${DRY_RUN:-0}"
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
LIST_DIR="$PKG_ROOT/download_lists/common"
LIST="$LIST_DIR/1kgp2022_phased_vcf.urls"
mkdir -p "$LIST_DIR" "$OUT_DIR"
: > "$LIST"

expand_chroms() {
  local spec="$1"
  if [[ "$spec" =~ ^([0-9]+)-([0-9]+)$ ]]; then
    seq "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"
  else
    for c in $spec; do echo "$c"; done
  fi
}

while read -r chrom; do
  [[ -n "$chrom" ]] || continue
  url="$BASE_URL/1kGP_high_coverage_Illumina.chr${chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  echo "$url" >> "$LIST"
  echo "$url.tbi" >> "$LIST"
done < <(expand_chroms "$CHROMS")

sort -u "$LIST" -o "$LIST"
echo "[INFO] VCF URL list: $LIST ($(wc -l < "$LIST") URLs)"

if [[ "$DRY_RUN" == "1" ]]; then
  sed -n '1,20p' "$LIST"
  exit 0
fi

case "$DOWNLOADER" in
  aria2c)
    command -v aria2c >/dev/null 2>&1 || { echo "[ERROR] aria2c not found; set DOWNLOADER=wget" >&2; exit 1; }
    aria2c -x "${ARIA2_CONNECTIONS:-4}" -s "${ARIA2_CONNECTIONS:-4}" -j "${ARIA2_JOBS:-8}" -c --auto-file-renaming=false --allow-overwrite=false -i "$LIST" -d "$OUT_DIR"
    ;;
  wget)
    command -v wget >/dev/null 2>&1 || { echo "[ERROR] wget not found" >&2; exit 1; }
    wget -c -i "$LIST" -P "$OUT_DIR"
    ;;
  *) echo "[ERROR] unsupported DOWNLOADER=$DOWNLOADER" >&2; exit 1 ;;
esac
