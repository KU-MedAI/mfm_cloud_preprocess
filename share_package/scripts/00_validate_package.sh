#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

required=(
  samples_3212.tsv pacbio_294_samples.tsv illumina_2918_samples.tsv rna_pair_449.tsv
  pacbio_295_samples.tsv illumina_2917_samples.tsv
  summary_counts.tsv summary.json phase_assignments.tsv phase_summary.tsv
  manifests/illumina_30x_files.tsv manifests/rna_geuvadis_files.tsv
  manifests/pacbio_file_manifest.tsv manifests/pacbio_active_file_manifest.tsv manifests/public_vcf_sources.tsv
  manifests/public_vcf_file_manifest.tsv manifests/reference_file_manifest.tsv
  PACBIO_FREEZE_PROVENANCE.md
)

for rel in "${required[@]}"; do
  [[ -s "$PKG_ROOT/$rel" ]] || { echo "[ERROR] missing or empty: $rel" >&2; exit 1; }
done

count_rows() { awk 'END { print NR-1 }' "$1"; }
assert_count() {
  local expected="$1" actual="$2" label="$3"
  [[ "$actual" == "$expected" ]] || {
    echo "[ERROR] $label expected=$expected actual=$actual" >&2
    exit 1
  }
}

assert_count 3212 "$(count_rows "$PKG_ROOT/samples_3212.tsv")" samples_3212
assert_count 294 "$(count_rows "$PKG_ROOT/pacbio_294_samples.tsv")" pacbio_294_active
assert_count 2918 "$(count_rows "$PKG_ROOT/illumina_2918_samples.tsv")" illumina_2918_active
assert_count 449 "$(count_rows "$PKG_ROOT/rna_pair_449.tsv")" rna_pair_449
assert_count 3212 "$(count_rows "$PKG_ROOT/phase_assignments.tsv")" phase_assignments

duplicates="$(tail -n +2 "$PKG_ROOT/samples_3212.tsv" | cut -f1 | sort | uniq -d | wc -l)"
assert_count 0 "$duplicates" duplicate_sample_ids

pacbio_ready="$(awk -F '\t' '
  FNR==NR {if(FNR>1) want[$1]=1; next}
  FNR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
  $(h["sample_id"]) in want && $(h["status"])=="ready" {seen[$(h["sample_id"])]=1}
  END{for(s in seen)n++; print n+0}' "$PKG_ROOT/pacbio_294_samples.tsv" "$PKG_ROOT/manifests/pacbio_file_manifest.tsv")"
assert_count 294 "$pacbio_ready" pacbio_ready_active_donors

active_clr="$(awk -F '\t' 'NR>1 && $3=="PacBio_CLR" {n++} END{print n+0}' "$PKG_ROOT/samples_3212.tsv")"
assert_count 0 "$active_clr" active_clr_donors

if rg -n '/NHNHOME/' "$PKG_ROOT" --glob '*.md' --glob '*.tsv' --glob '*.json' >/dev/null; then
  echo "[ERROR] package contains server-local absolute path" >&2
  rg -n '/NHNHOME/' "$PKG_ROOT" --glob '*.md' --glob '*.tsv' --glob '*.json' >&2
  exit 1
fi

python_bin="${PYTHON_BIN:-python}"
"$python_bin" -m json.tool "$PKG_ROOT/summary.json" >/dev/null
"$python_bin" "$SCRIPT_DIR/06_validate_pacbio_freeze.py" --package-root "$PKG_ROOT"
bash -n "$SCRIPT_DIR"/*.sh

echo "[OK] package rows, IDs, JSON, paths and shell syntax validated"
echo "[INFO] active PacBio HiFi ready donors: $pacbio_ready/294"
echo "[INFO] historical 295-donor freeze is preserved and independently validated"
