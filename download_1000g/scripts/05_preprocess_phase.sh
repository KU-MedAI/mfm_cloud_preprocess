#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash scripts/05_preprocess_phase.sh <phase|all> [all|shard_id]

Create sample lists for preprocessing, optionally run a user-provided preprocessing command.
This script does not download anything. It is the no-download preprocessing entry point.

Environment variables:
  MODE=all                    all | pre | post | wgs_only | rna_pair
  PLATFORM=all                all | Illumina_30x | ONT | HiFi_CCS
  RUN_PREPROCESS=0            If 1, execute PREPROCESS_CMD.
  PREPROCESS_CMD=...          Shell command template. Available env vars:
                              SAMPLE_LIST, SAMPLES_SPACE, PHASE, SHARD, MODE, PLATFORM, PACKAGE_ROOT
  PROCESS_CODE_DIR=           Optional directory to cd into before running PREPROCESS_CMD.
  OUT_ROOT=preprocess_lists   Output sample list directory.

Examples:
  bash scripts/05_preprocess_phase.sh 1
  MODE=post PLATFORM=ONT bash scripts/05_preprocess_phase.sh 1
  RUN_PREPROCESS=1 PROCESS_CODE_DIR=../process_code \
    PREPROCESS_CMD='SAMPLES="$SAMPLES_SPACE" ./run_cohort.sh igsr_ont_wgs_only' \
    bash scripts/05_preprocess_phase.sh 1 shard_001
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 1 ]]; then
  usage
  exit $([[ $# -lt 1 ]] && echo 1 || echo 0)
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PHASE_ARG="$1"
SHARD_ARG="${2:-all}"
ASSIGN_TSV="${ASSIGN_TSV:-$PKG_ROOT/phase_assignments.tsv}"
OUT_ROOT="${OUT_ROOT:-$PKG_ROOT/preprocess_lists}"
MODE="${MODE:-all}"
PLATFORM="${PLATFORM:-all}"
RUN_PREPROCESS="${RUN_PREPROCESS:-0}"
PREPROCESS_CMD="${PREPROCESS_CMD:-}"
PROCESS_CODE_DIR="${PROCESS_CODE_DIR:-}"

case "$MODE" in all|pre|post|wgs_only|rna_pair) ;; *) echo "[ERROR] bad MODE=$MODE" >&2; exit 1;; esac
case "$PLATFORM" in all|Illumina_30x|ONT|HiFi_CCS) ;; *) echo "[ERROR] bad PLATFORM=$PLATFORM" >&2; exit 1;; esac

if [[ ! -f "$ASSIGN_TSV" ]]; then
  echo "[INFO] phase_assignments.tsv not found; creating default assignments"
  bash "$SCRIPT_DIR/00_make_phase_assignments.sh"
fi

make_list() {
  local phase_raw="$1"
  local phase="$(printf "%02d" "${phase_raw#0}")"
  local out_dir="$OUT_ROOT/phase_$phase"
  mkdir -p "$out_dir"
  local sample_list="$out_dir/samples.mode_${MODE}.platform_${PLATFORM}.shard_${SHARD_ARG}.txt"

  awk -v phase="$phase" -v shard_arg="$SHARD_ARG" -v mode="$MODE" -v platform="$PLATFORM" '
  BEGIN { FS=OFS="\t" }
  NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
  function val(name){ return (name in h) ? $(h[name]) : "" }
  function include_mode(use) {
    if (mode=="all") return 1
    if (mode=="post" || mode=="rna_pair") return use=="post_training_wgs_rna_pair"
    if (mode=="pre" || mode=="wgs_only") return use ~ /^pre_training/
    return 0
  }
  val("phase") == phase {
    if (platform != "all" && val("assigned_wgs_platform") != platform) next
    if (!include_mode(val("recommended_training_use"))) next
    print val("sample_id")
  }' "$ASSIGN_TSV" > "$sample_list"

  if [[ "$SHARD_ARG" != "all" ]]; then
    # If a shard-specific phase sample file exists, intersect with it.
    local shard="${SHARD_ARG#shard_}"
    shard="$(printf "%03d" "${shard#0}")"
    local shard_samples="$PKG_ROOT/download_lists/phase_$phase/shards/shard_$shard.samples"
    if [[ -f "$shard_samples" ]]; then
      grep -Fxf "$shard_samples" "$sample_list" > "$sample_list.tmp" || true
      mv "$sample_list.tmp" "$sample_list"
    else
      echo "[WARN] shard sample file not found; using full phase list: $shard_samples" >&2
    fi
  fi

  echo "[OK] sample list: $sample_list ($(wc -l < "$sample_list") samples)"

  if [[ "$RUN_PREPROCESS" == "1" ]]; then
    if [[ -z "$PREPROCESS_CMD" ]]; then
      echo "[ERROR] RUN_PREPROCESS=1 requires PREPROCESS_CMD" >&2
      exit 1
    fi
    export SAMPLE_LIST="$sample_list"
    export SAMPLES_SPACE
    SAMPLES_SPACE="$(paste -sd ' ' "$sample_list")"
    export PHASE="$phase" SHARD="$SHARD_ARG" MODE PLATFORM PACKAGE_ROOT="$PKG_ROOT"
    if [[ -n "$PROCESS_CODE_DIR" ]]; then
      cd "$PROCESS_CODE_DIR"
    fi
    echo "[INFO] running PREPROCESS_CMD"
    bash -lc "$PREPROCESS_CMD"
  fi
}

if [[ "$PHASE_ARG" == "all" ]]; then
  awk 'BEGIN{FS="\t"} NR>1{print $2}' "$ASSIGN_TSV" | sort -u | while read -r p; do make_list "$p"; done
else
  make_list "$PHASE_ARG"
fi
