#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash scripts/01_make_phase_url_lists.sh <phase|all>

Environment:
  PLATFORM_FILTER=all          all | Illumina_30x | PacBio | PacBio_HiFi | PacBio_CLR
  INCLUDE_RNA=1               include strict GEUVADIS BAM/BAI
  OUT_ROOT=download_lists

PacBio read·VCF·reference URL은 동결 manifest에서 함께 생성된다.
USAGE
}

[[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && { usage; exit 0; }
[[ $# -ge 1 ]] || { usage >&2; exit 1; }

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PHASE_ARG="$1"
PLATFORM_FILTER="${PLATFORM_FILTER:-all}"
INCLUDE_RNA="${INCLUDE_RNA:-1}"
OUT_ROOT="${OUT_ROOT:-$PKG_ROOT/download_lists}"
ASSIGN="$PKG_ROOT/phase_assignments.tsv"
ILLUMINA="$PKG_ROOT/manifests/illumina_30x_files.tsv"
PACBIO="$PKG_ROOT/manifests/pacbio_file_manifest.tsv"
PUBLIC_VCF="$PKG_ROOT/manifests/public_vcf_file_manifest.tsv"
REFERENCE="$PKG_ROOT/manifests/reference_file_manifest.tsv"
RNA="$PKG_ROOT/manifests/rna_geuvadis_files.tsv"

case "$PLATFORM_FILTER" in all|Illumina_30x|PacBio|PacBio_HiFi|PacBio_CLR) ;; *) echo "[ERROR] bad PLATFORM_FILTER" >&2; exit 1;; esac
case "$INCLUDE_RNA" in 0|1) ;; *) echo "[ERROR] INCLUDE_RNA must be 0 or 1" >&2; exit 1;; esac

make_phase() {
  local raw="$1" phase out selected
  phase="$(printf '%02d' "${raw#0}")"
  out="$OUT_ROOT/phase_$phase"
  mkdir -p "$out"
  selected="$out/samples.tsv"

  awk -F '\t' -v OFS='\t' -v p="$phase" -v platform="$PLATFORM_FILTER" '
    function selected(value) { return platform=="all" || value==platform || (platform=="PacBio" && value ~ /^PacBio_/) }
    NR==1 {for(i=1;i<=NF;i++) h[$i]=i; print "sample_id","selected_platform","selected_project"; next}
    $(h["phase"])==p && selected($(h["selected_platform"])) {
      print $(h["sample_id"]),$(h["selected_platform"]),$(h["selected_project"])
    }' "$ASSIGN" > "$selected"

  local pacbio_selected
  pacbio_selected="$(awk -F '\t' 'NR>1 && $2 ~ /^PacBio_/ {n++} END{print n+0}' "$selected")"
  if (( pacbio_selected > 0 )); then
    local pacbio_ready
    pacbio_ready="$(awk -F '\t' '
      FNR==NR {if(FNR>1 && $2 ~ /^PacBio_/) want[$1]=1; next}
      FNR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
      $(h["sample_id"]) in want && $(h["status"])=="ready" && $(h["file_url"])!="" {seen[$(h["sample_id"])]=1}
      END{for(s in seen)n++; print n+0}' "$selected" "$PACBIO")"
    if (( pacbio_ready < pacbio_selected )); then
      echo "[ERROR] phase=$phase PacBio selected=$pacbio_selected donors_with_ready_files=$pacbio_ready" >&2
      echo "[ERROR] freeze run-level PacBio URL/accession/size/MD5 before generating mixed/PacBio lists" >&2
      exit 2
    fi
  fi

  local with_kind="$out/urls.with_kind.tsv"
  printf 'sample_id\tkind\turl\n' > "$with_kind"

  awk -F '\t' -v OFS='\t' '
    FNR==NR {if(FNR>1 && $2=="Illumina_30x") want[$1]=1; next}
    FNR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
    $(h["sample_id"]) in want {
      if($(h["file_url"])!="") print $(h["sample_id"]),"wgs_cram",$(h["file_url"])
      if($(h["index_url"])!="") print $(h["sample_id"]),"wgs_index",$(h["index_url"])
    }' "$selected" "$ILLUMINA" >> "$with_kind"

  if (( pacbio_selected > 0 )); then
    awk -F '\t' -v OFS='\t' '
      FNR==NR {if(FNR>1 && $2 ~ /^PacBio_/) want[$1]=1; next}
      FNR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
      $(h["sample_id"]) in want && $(h["status"])=="ready" {
        if($(h["file_url"])!="") print $(h["sample_id"]),"wgs_pacbio",$(h["file_url"])
        if($(h["index_url"])!="") print $(h["sample_id"]),"wgs_index",$(h["index_url"])
      }' "$selected" "$PACBIO" >> "$with_kind"
  fi

  if (( pacbio_selected > 0 )); then
    awk -F '\t' -v OFS='\t' '
      FNR==NR {if(FNR>1 && $2 ~ /^PacBio_/) {want[$1]=1; projects[$3]=1}; next}
      FNR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
      $(h["project"]) in projects && ($(h["sample_id"])=="" || $(h["sample_id"]) in want) && $(h["status"])=="ready" {
        kind=$(h["file_type"]); if(kind=="vcf")kind="public_vcf"; else if(kind=="tbi")kind="public_vcf_index"; else if(kind=="tsv")kind="phase_tags";
        label=$(h["sample_id"]); if(label=="")label="__" $(h["project"]);
        print label,kind,$(h["file_url"])
      }' "$selected" "$PUBLIC_VCF" >> "$with_kind"
    awk -F '\t' -v OFS='\t' '
      NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
      $(h["status"])=="ready" && ($(h["resource_type"])=="fasta_gz" || $(h["resource_type"])=="fai") {
        print "__REFERENCE__","reference_" $(h["resource_type"]),$(h["file_url"])
      }' "$REFERENCE" >> "$with_kind"
  fi

  if [[ "$INCLUDE_RNA" == "1" ]]; then
    awk -F '\t' -v OFS='\t' '
      FNR==NR {if(FNR>1) want[$1]=1; next}
      FNR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
      $(h["sample_id"]) in want {
        if($(h["file_url"])!="") print $(h["sample_id"]),"rna_bam",$(h["file_url"])
        if($(h["index_url"])!="") print $(h["sample_id"]),"rna_index",$(h["index_url"])
      }' "$selected" "$RNA" >> "$with_kind"
  fi

  tail -n +2 "$with_kind" | cut -f3 | sort -u > "$out/all.urls"
  if ! awk -F/ '{base=$NF; if(base in seen && seen[base]!=$0){print "[ERROR] URL basename collision: " base > "/dev/stderr"; bad=1} seen[base]=$0} END{exit bad}' "$out/all.urls"; then
    exit 3
  fi
  echo "[OK] phase=$phase samples=$(($(wc -l < "$selected")-1)) urls=$(wc -l < "$out/all.urls")"
}

if [[ "$PHASE_ARG" == "all" ]]; then
  awk -F '\t' 'NR>1{print $2}' "$ASSIGN" | sort -u | while read -r phase; do make_phase "$phase"; done
else
  make_phase "$PHASE_ARG"
fi

