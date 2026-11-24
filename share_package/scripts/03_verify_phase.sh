#!/usr/bin/env bash
set -euo pipefail

[[ $# -ge 1 ]] || { echo "Usage: bash scripts/03_verify_phase.sh <phase>" >&2; exit 1; }
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PHASE="$(printf '%02d' "${1#0}")"
ROOT="${DOWNLOAD_ROOT:-$PKG_ROOT/downloads}/phase_$PHASE"
LIST_ROOT="${OUT_ROOT:-$PKG_ROOT/download_lists}"
WITH_KIND="$LIST_ROOT/phase_$PHASE/urls.with_kind.tsv"
PYTHON_BIN="${PYTHON_BIN:-python}"
[[ -d "$ROOT" ]] || { echo "[ERROR] missing download directory: $ROOT" >&2; exit 1; }
[[ -s "$WITH_KIND" ]] || { echo "[ERROR] missing URL metadata list: $WITH_KIND" >&2; exit 1; }

"$PYTHON_BIN" "$SCRIPT_DIR/07_verify_frozen_downloads.py" \
  --package-root "$PKG_ROOT" --download-root "$ROOT" --url-list "$WITH_KIND"

status=0
while IFS= read -r cram; do
  [[ -f "$cram.crai" || -f "${cram%.cram}.crai" ]] || { echo "[MISSING_INDEX] $cram" >&2; status=1; }
done < <(find "$ROOT" -type f -name '*.cram' | sort)

while IFS= read -r bam; do
  [[ -f "$bam.bai" || -f "${bam%.bam}.bai" || -f "$bam.csi" ]] || echo "[WARN_NO_INDEX] $bam" >&2
done < <(find "$ROOT" -type f -name '*.bam' | sort)

missing_vcf_index=0
while IFS= read -r vcf; do
  if [[ ! -s "$vcf.tbi" && ! -s "$vcf.csi" ]]; then
    if command -v bcftools >/dev/null 2>&1; then
      echo "[INFO] 공개 index가 없는 VCF의 TBI 생성: $vcf" >&2
      bcftools index -t "$vcf" || status=1
    else
      echo "[ERROR] VCF index가 없고 bcftools도 사용할 수 없음: $vcf" >&2
      missing_vcf_index=1
    fi
  fi
done < <(find "$ROOT" -type f -name '*.vcf.gz' | sort)
(( missing_vcf_index == 0 )) || status=1

if command -v samtools >/dev/null 2>&1; then
  find "$ROOT" -type f \( -name '*.cram' -o -name '*.bam' \) -print0 | xargs -0 -r -n 64 samtools quickcheck -v || status=1
else
  echo "[WARN] samtools unavailable; quickcheck skipped" >&2
fi

echo "[INFO] files=$(find "$ROOT" -type f | wc -l) root=$ROOT"
exit "$status"

