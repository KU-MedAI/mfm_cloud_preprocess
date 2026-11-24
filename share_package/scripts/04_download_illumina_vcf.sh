#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CHROMS="${CHROMS:-1-22}"
OUT_DIR="${OUT_DIR:-$PKG_ROOT/downloads/common/illumina_vcf}"
DRY_RUN="${DRY_RUN:-0}"
DOWNLOADER="${DOWNLOADER:-aria2c}"
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
LIST="$PKG_ROOT/download_lists/common/illumina_1kgp2022_vcf.urls"
mkdir -p "$(dirname "$LIST")" "$OUT_DIR"
: > "$LIST"

expand_chroms() {
  if [[ "$1" =~ ^([0-9]+)-([0-9]+)$ ]]; then seq "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"; else for c in $1; do echo "$c"; done; fi
}

while read -r chrom; do
  url="$BASE_URL/1kGP_high_coverage_Illumina.chr${chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  printf '%s\n%s\n' "$url" "$url.tbi" >> "$LIST"
done < <(expand_chroms "$CHROMS")

echo "[WARN] this VCF is for selected 1KGP Illumina donors, not a PacBio project VCF"
[[ "$DRY_RUN" == "1" ]] && { cat "$LIST"; exit 0; }

case "$DOWNLOADER" in
  aria2c) aria2c -x "${ARIA2_CONNECTIONS:-4}" -s "${ARIA2_CONNECTIONS:-4}" -j "${ARIA2_JOBS:-8}" -c --auto-file-renaming=false --allow-overwrite=false -i "$LIST" -d "$OUT_DIR" ;;
  wget) wget -c -i "$LIST" -P "$OUT_DIR" ;;
  *) echo "[ERROR] unsupported DOWNLOADER=$DOWNLOADER" >&2; exit 1;;
esac

