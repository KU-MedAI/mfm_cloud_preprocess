#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash scripts/00_make_phase_assignments.sh

Create balanced phase assignments from samples_3202.tsv.

Environment variables:
  PHASES=5                         Number of logical phases.
  BALANCE_MODE=aligned_ready_fallback
                                   aligned_ready_fallback | primary_url | all_manifest_runs
  INCLUDE_RNA=1                     Add RNA BAM/available source size for strict RNA-pair samples.
  INCLUDE_BIGWIG=0                  Legacy option. BigWig is excluded by default.
  SAMPLES_TSV=samples_3202.tsv      Input TSV. Default: package samples_3202.tsv.
  OUT_TSV=phase_assignments.tsv     Output phase assignment TSV.
  SUMMARY_TSV=phase_summary.tsv     Output phase summary TSV.

Recommended default:
  PHASES=5 BALANCE_MODE=aligned_ready_fallback INCLUDE_RNA=1 INCLUDE_BIGWIG=0 \
    bash scripts/00_make_phase_assignments.sh
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
SAMPLES_TSV="${SAMPLES_TSV:-$PKG_ROOT/samples_3202.tsv}"
OUT_TSV="${OUT_TSV:-$PKG_ROOT/phase_assignments.tsv}"
SUMMARY_TSV="${SUMMARY_TSV:-$PKG_ROOT/phase_summary.tsv}"
PHASES="${PHASES:-5}"
BALANCE_MODE="${BALANCE_MODE:-aligned_ready_fallback}"
INCLUDE_RNA="${INCLUDE_RNA:-1}"
INCLUDE_BIGWIG="${INCLUDE_BIGWIG:-0}"

if [[ ! -f "$SAMPLES_TSV" ]]; then
  echo "[ERROR] missing SAMPLES_TSV: $SAMPLES_TSV" >&2
  exit 1
fi

case "$BALANCE_MODE" in
  aligned_ready_fallback|primary_url|all_manifest_runs) ;;
  *) echo "[ERROR] unsupported BALANCE_MODE: $BALANCE_MODE" >&2; exit 1 ;;
esac

if [[ "$INCLUDE_BIGWIG" == "1" ]]; then
  echo "[ERROR] recount3 BigWig is excluded from this package. Use RNA BAM sources only." >&2
  exit 1
fi

TMP_EST="$(mktemp)"
trap 'rm -f "$TMP_EST"' EXIT

awk -v mode="$BALANCE_MODE" -v include_rna="$INCLUDE_RNA" -v include_bigwig="$INCLUDE_BIGWIG" '
BEGIN { FS=OFS="\t"; GiB=1073741824 }
NR==1 {
  for (i=1; i<=NF; i++) h[$i]=i
  next
}
function val(name) { return (name in h) ? $(h[name]) : "" }
function num(x) { gsub(/[^0-9.]/, "", x); return x+0 }
function fallback(plat, source) {
  if (plat == "Illumina_30x") return 15*GiB
  if (plat == "ONT" && source == "HPRC_ONT") return 50*GiB
  if (plat == "ONT") return 6*GiB
  if (plat == "HiFi_CCS") return 40*GiB
  return 15*GiB
}
{
  sample=val("sample_id")
  use=val("recommended_training_use")
  plat=val("assigned_wgs_platform")
  source=val("wgs_download_source")
  total=num(val("wgs_total_bytes_for_assigned_platform"))
  primary=num(val("wgs_primary_bytes"))
  if (mode == "all_manifest_runs") wgs=total
  else if (mode == "primary_url") wgs=primary
  else {
    if (plat == "Illumina_30x" || source == "ONT_Vienna") wgs=(primary>0 ? primary : total)
    else wgs=fallback(plat, source)
  }
  if (wgs <= 0) wgs=fallback(plat, source)

  rna=0
  if (include_rna == "1" && use == "post_training_wgs_rna_pair") {
    if (mode == "all_manifest_runs") rna=num(val("rna_total_bytes_available"))
    else rna=num(val("rna_primary_bytes"))
    if (rna <= 0) rna=5*GiB
    if (include_bigwig == "1" && val("has_rna_bigwig") == "True") rna += int(0.1*GiB)
  }
  est=wgs+rna
  printf "%.0f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", est, sample, use, plat, source, val("rna_pair_status"), val("has_rna_bam_strict"), wgs, rna
}' "$SAMPLES_TSV" | sort -nr -k1,1 > "$TMP_EST"

awk -v phases="$PHASES" -v summary="$SUMMARY_TSV" -v mode="$BALANCE_MODE" -v include_rna="$INCLUDE_RNA" -v include_bigwig="$INCLUDE_BIGWIG" '
BEGIN {
  FS=OFS="\t"; GiB=1073741824
  print "sample_id","phase","estimated_download_bytes","estimated_download_gib","recommended_training_use","assigned_wgs_platform","wgs_download_source","rna_pair_status","has_rna_bam_strict","estimated_wgs_bytes","estimated_rna_bytes","balance_mode","include_rna","include_bigwig"
}
{
  est=$1+0; sample=$2; use=$3; plat=$4; source=$5; rna_status=$6; has_bam=$7; wgs=$8+0; rna=$9+0
  min=1
  for (i=2; i<=phases; i++) if (sum[i] < sum[min]) min=i
  sum[min]+=est; count[min]++
  if (use == "post_training_wgs_rna_pair") post[min]++
  else pre_only[min]++
  plat_count[min,plat]++
  if (est > max_sample[min]) max_sample[min]=est
  printf "%s\t%02d\t%.0f\t%.3f\t%s\t%s\t%s\t%s\t%s\t%.0f\t%.0f\t%s\t%s\t%s\n", sample, min, est, est/GiB, use, plat, source, rna_status, has_bam, wgs, rna, mode, include_rna, include_bigwig
}
END {
  print "phase","samples","estimated_download_bytes","estimated_download_tib","post_training_wgs_rna_pair","pre_training_wgs_only","HiFi_CCS","ONT","Illumina_30x","max_sample_gib","balance_mode","include_rna","include_bigwig" > summary
  for (i=1; i<=phases; i++) {
    printf "%02d\t%d\t%.0f\t%.3f\t%d\t%d\t%d\t%d\t%d\t%.3f\t%s\t%s\t%s\n", i, count[i]+0, sum[i], sum[i]/GiB/1024, post[i]+0, pre_only[i]+0, plat_count[i,"HiFi_CCS"]+0, plat_count[i,"ONT"]+0, plat_count[i,"Illumina_30x"]+0, max_sample[i]/GiB, mode, include_rna, include_bigwig >> summary
  }
}
' "$TMP_EST" > "$OUT_TSV"

echo "[OK] wrote $OUT_TSV"
echo "[OK] wrote $SUMMARY_TSV"
