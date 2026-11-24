#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash scripts/01_make_phase_url_lists.sh <phase|all>

Create URL lists for one phase or all phases.

Environment variables:
  SHARD_GIB=1500                 Max estimated raw size per shard in GiB.
  INCLUDE_WGS=1                  Include WGS URLs.
  INCLUDE_RNA=1                  Include RNA URLs for strict RNA-pair samples.
  INCLUDE_BIGWIG=0               recount3 BigWig is excluded; setting 1 is rejected.
  WGS_SOURCE_MODE=assigned       assigned | primary | all_manifest_runs | aligned_only
  RNA_SOURCE_MODE=strict_bam     strict_bam | primary | none
  OUT_ROOT=download_lists        Output directory under package root.

Examples:
  bash scripts/01_make_phase_url_lists.sh 1
  SHARD_GIB=750 INCLUDE_RNA=1 bash scripts/01_make_phase_url_lists.sh 1
  WGS_SOURCE_MODE=all_manifest_runs bash scripts/01_make_phase_url_lists.sh 2
  WGS_SOURCE_MODE=aligned_only RNA_SOURCE_MODE=strict_bam bash scripts/01_make_phase_url_lists.sh all
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 1 ]]; then
  usage
  exit $([[ $# -lt 1 ]] && echo 1 || echo 0)
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PHASE_ARG="$1"
SAMPLES_TSV="${SAMPLES_TSV:-$PKG_ROOT/samples_3202.tsv}"
ASSIGN_TSV="${ASSIGN_TSV:-$PKG_ROOT/phase_assignments.tsv}"
WGS_MANIFEST="${WGS_MANIFEST:-$PKG_ROOT/manifests/wgs_urls.tsv}"
RNA_MANIFEST="${RNA_MANIFEST:-$PKG_ROOT/manifests/rna_urls.tsv}"
OUT_ROOT="${OUT_ROOT:-$PKG_ROOT/download_lists}"
SHARD_GIB="${SHARD_GIB:-1500}"
INCLUDE_WGS="${INCLUDE_WGS:-1}"
INCLUDE_RNA="${INCLUDE_RNA:-1}"
INCLUDE_BIGWIG="${INCLUDE_BIGWIG:-0}"
WGS_SOURCE_MODE="${WGS_SOURCE_MODE:-assigned}"
RNA_SOURCE_MODE="${RNA_SOURCE_MODE:-strict_bam}"

case "$WGS_SOURCE_MODE" in assigned|primary|all_manifest_runs|aligned_only) ;; *) echo "[ERROR] bad WGS_SOURCE_MODE=$WGS_SOURCE_MODE" >&2; exit 1;; esac
case "$RNA_SOURCE_MODE" in strict_bam|primary|none) ;; *) echo "[ERROR] bad RNA_SOURCE_MODE=$RNA_SOURCE_MODE" >&2; exit 1;; esac

if [[ "$INCLUDE_BIGWIG" == "1" ]]; then
  echo "[ERROR] recount3 BigWig is excluded from this package. Use RNA BAM sources only." >&2
  exit 1
fi

if [[ ! -f "$ASSIGN_TSV" ]]; then
  echo "[INFO] phase_assignments.tsv not found; creating default assignments"
  bash "$SCRIPT_DIR/00_make_phase_assignments.sh"
fi

make_one_phase() {
  local phase_raw="$1"
  local phase
  phase="$(printf "%02d" "${phase_raw#0}")"
  local phase_dir="$OUT_ROOT/phase_$phase"
  local shard_dir="$phase_dir/shards"
  mkdir -p "$phase_dir" "$shard_dir"
  rm -f "$phase_dir"/*.urls "$phase_dir"/*.tsv "$shard_dir"/*.urls "$shard_dir"/*.samples

  awk -v phase="$phase" -v limit_gib="$SHARD_GIB" '
  BEGIN { FS=OFS="\t"; GiB=1073741824; limit=limit_gib*GiB; shard=1 }
  NR==1 { for (i=1;i<=NF;i++) h[$i]=i; print "sample_id","shard","estimated_download_bytes","recommended_training_use","assigned_wgs_platform","wgs_download_source"; next }
  function val(name) { return (name in h) ? $(h[name]) : "" }
  val("phase") == phase {
    est=val("estimated_download_bytes")+0
    if (current > 0 && current + est > limit) { shard++; current=0 }
    current += est
    printf "%s\t%03d\t%.0f\t%s\t%s\t%s\n", val("sample_id"), shard, est, val("recommended_training_use"), val("assigned_wgs_platform"), val("wgs_download_source")
  }' "$ASSIGN_TSV" > "$phase_dir/phase_samples.tsv"

  awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1 > (dir "/shard_" $2 ".samples")}' dir="$shard_dir" "$phase_dir/phase_samples.tsv"

  local url_with_shard="$phase_dir/urls.with_shard.tsv"
  : > "$url_with_shard"

  if [[ "$INCLUDE_WGS" == "1" ]]; then
    if [[ "$WGS_SOURCE_MODE" == "primary" || "$WGS_SOURCE_MODE" == "assigned" || "$WGS_SOURCE_MODE" == "aligned_only" ]]; then
      awk -v mode="$WGS_SOURCE_MODE" 'BEGIN{FS=OFS="\t"}
        FNR==NR { if (FNR>1) shard[$1]=$2; next }
        FNR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
        function val(name){ return (name in h) ? $(h[name]) : "" }
        (val("sample_id") in shard) {
          source=val("wgs_download_source"); plat=val("assigned_wgs_platform")
          if (mode=="aligned_only" && !(source=="ONT_Vienna" || plat=="Illumina_30x")) next
          url=val("wgs_primary_url"); idx=val("wgs_index_url_if_predictable")
          if (url != "") print shard[val("sample_id")],"wgs",url
          if (idx != "") print shard[val("sample_id")],"wgs",idx
        }' "$phase_dir/phase_samples.tsv" "$SAMPLES_TSV" >> "$url_with_shard"
    fi

    if [[ "$WGS_SOURCE_MODE" == "all_manifest_runs" ]]; then
      awk 'BEGIN{FS=OFS="\t"}
        FNR==NR { if (FNR>1) shard[$1]=$2; next }
        FILENAME==samples && FNR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
        FILENAME==samples {
          sid=$(h["sample_id"]); if (sid in shard) { src[sid]=$(h["wgs_download_source"]); primary[sid]=$(h["wgs_primary_url"]); idx[sid]=$(h["wgs_index_url_if_predictable"]) }
          next
        }
        FILENAME==manifest && FNR==1 { for(i=1;i<=NF;i++) mh[$i]=i; next }
        FILENAME==manifest {
          sid=$(mh["sample"]); s=$(mh["source"]); url=$(mh["url"])
          if (s == "HGSVC_HiFi" && (url ~ /Kinnex/ || url ~ /\.flnc\.bam$/)) next
          if (sid in shard && s==src[sid] && url != "") {
            print shard[sid],"wgs",url
            if (url ~ /\.cram$/) print shard[sid],"wgs",url ".crai"
          }
        }
        END {
          for (sid in shard) {
            if (src[sid] == "ENA_PRJEB31736_1kGP_30x" && primary[sid] != "") {
              print shard[sid],"wgs",primary[sid]
              if (idx[sid] != "") print shard[sid],"wgs",idx[sid]
            }
          }
        }' samples="$SAMPLES_TSV" manifest="$WGS_MANIFEST" "$phase_dir/phase_samples.tsv" "$SAMPLES_TSV" "$WGS_MANIFEST" >> "$url_with_shard"
    fi
  fi

  if [[ "$INCLUDE_RNA" == "1" && "$RNA_SOURCE_MODE" != "none" ]]; then
    if [[ "$RNA_SOURCE_MODE" == "primary" ]]; then
      awk 'BEGIN{FS=OFS="\t"}
        FNR==NR { if (FNR>1) shard[$1]=$2; next }
        FNR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
        function val(name){ return (name in h) ? $(h[name]) : "" }
        (val("sample_id") in shard && val("recommended_training_use") == "post_training_wgs_rna_pair") {
          url=val("rna_primary_url")
          if (url != "") {
            print shard[val("sample_id")],"rna_bam",url
            if (url ~ /\.bam$/) print shard[val("sample_id")],"rna_bam",url ".bai"
          }
        }' "$phase_dir/phase_samples.tsv" "$SAMPLES_TSV" >> "$url_with_shard"
    else
      awk -v mode="$RNA_SOURCE_MODE" 'BEGIN{FS=OFS="\t"}
        FNR==NR { if (FNR>1) { shard[$1]=$2; use[$1]=$4 } next }
        FILENAME==manifest && FNR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
        FILENAME==manifest {
          sid=$(h["sample"]); source=$(h["source"]); url=$(h["url"])
          if (!(sid in shard) || use[sid] != "post_training_wgs_rna_pair" || url == "") next
          if (mode=="strict_bam" && source != "Geuvadis_TOPMed") next
          kind=(source=="Geuvadis_TOPMed" ? "rna_bam" : "rna_other")
          print shard[sid],kind,url
          if (source=="Geuvadis_TOPMed" && url ~ /\.bam$/) print shard[sid],kind,url ".bai"
        }' manifest="$RNA_MANIFEST" "$phase_dir/phase_samples.tsv" "$RNA_MANIFEST" >> "$url_with_shard"
    fi
  fi

  if [[ "$INCLUDE_BIGWIG" == "1" ]]; then
    awk 'BEGIN{FS=OFS="\t"}
      FNR==NR { if (FNR>1) shard[$1]=$2; next }
      FNR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
      function val(name){ return (name in h) ? $(h[name]) : "" }
      (val("sample_id") in shard && val("rna_bigwig_url_validated_subset") != "") {
        print shard[val("sample_id")],"rna_bigwig",val("rna_bigwig_url_validated_subset")
      }' "$phase_dir/phase_samples.tsv" "$SAMPLES_TSV" >> "$url_with_shard"
  fi

  sort -u "$url_with_shard" -o "$url_with_shard"
  awk -v phase_dir="$phase_dir" -v shard_dir="$shard_dir" 'BEGIN{FS=OFS="\t"}
    NF>=3 {
      print $3 >> (phase_dir "/all.urls")
      print $3 >> (phase_dir "/" $2 ".urls")
      print $3 >> (shard_dir "/shard_" $1 ".urls")
    }' "$url_with_shard"

  for f in "$phase_dir"/*.urls "$shard_dir"/*.urls; do
    [[ -e "$f" ]] || continue
    sort -u "$f" -o "$f"
  done

  awk 'BEGIN{FS=OFS="\t"} NR>1 {c[$2]++} END{for (s in c) print s,c[s]}' "$phase_dir/phase_samples.tsv" | sort -k1,1 > "$phase_dir/shard_summary.tsv.tmp"
  { printf "shard\tsamples\n"; cat "$phase_dir/shard_summary.tsv.tmp"; } > "$phase_dir/shard_summary.tsv"
  rm -f "$phase_dir/shard_summary.tsv.tmp"
  echo "[OK] phase $phase URL lists: $phase_dir"
  echo "     samples: $(($(wc -l < "$phase_dir/phase_samples.tsv") - 1)), urls: $(wc -l < "$phase_dir/all.urls" 2>/dev/null || echo 0), shards: $(find "$shard_dir" -name 'shard_*.samples' | wc -l)"
}

if [[ "$PHASE_ARG" == "all" ]]; then
  awk 'BEGIN{FS="\t"} NR>1{print $2}' "$ASSIGN_TSV" | sort -u | while read -r p; do make_one_phase "$p"; done
else
  make_one_phase "$PHASE_ARG"
fi
