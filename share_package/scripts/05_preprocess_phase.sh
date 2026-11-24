#!/usr/bin/env bash
set -euo pipefail

[[ $# -ge 1 ]] || {
  echo "Usage: bash scripts/05_preprocess_phase.sh <phase|all>" >&2
  exit 1
}

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PHASE_ARG="$1"
MODE="${MODE:-all}"
PLATFORM="${PLATFORM:-all}"
OUT_ROOT="${OUT_ROOT:-$PKG_ROOT/preprocess_lists}"
ASSIGN_TSV="${ASSIGN_TSV:-$PKG_ROOT/phase_assignments.tsv}"
RUN_PREPROCESS="${RUN_PREPROCESS:-0}"
PREPROCESS_CMD="${PREPROCESS_CMD:-}"
PAIR_WGS_CMD="${PAIR_WGS_CMD:-}"
PAIR_RNA_CMD="${PAIR_RNA_CMD:-}"
WGS_ONLY_CMD="${WGS_ONLY_CMD:-}"
PROCESS_CODE_DIR="${PROCESS_CODE_DIR:-}"

case "$MODE" in
  all|wgs_only|rna_pair|pair_first) ;;
  *) echo "[ERROR] MODE must be all, wgs_only, rna_pair or pair_first" >&2; exit 1 ;;
esac
case "$PLATFORM" in
  all|Illumina_30x|PacBio_HiFi) ;;
  *) echo "[ERROR] bad PLATFORM=$PLATFORM" >&2; exit 1 ;;
esac
case "$RUN_PREPROCESS" in
  0|1) ;;
  *) echo "[ERROR] RUN_PREPROCESS must be 0 or 1" >&2; exit 1 ;;
esac
[[ -f "$ASSIGN_TSV" ]] || { echo "[ERROR] assignment TSV not found: $ASSIGN_TSV" >&2; exit 1; }

normalise_phase() {
  printf "%02d" "${1#0}"
}

build_list() {
  local phase="$1"
  local list_mode="$2"
  local out_dir out
  out_dir="$OUT_ROOT/phase_$phase"
  mkdir -p "$out_dir"
  out="$out_dir/samples.mode_$list_mode.platform_$PLATFORM.txt"

  awk -F '\t' -v p="$phase" -v mode="$list_mode" -v platform="$PLATFORM" '
    NR==1 {
      for (i=1; i<=NF; i++) header[$i]=i
      required[1]="sample_id"
      required[2]="phase"
      required[3]="selected_platform"
      required[4]="recommended_training_use"
      for (i=1; i<=4; i++) {
        if (!(required[i] in header)) {
          print "[ERROR] assignment TSV missing column: " required[i] > "/dev/stderr"
          exit 2
        }
      }
      next
    }
    $(header["phase"]) == p {
      if (platform != "all" && $(header["selected_platform"]) != platform) next
      if (mode == "wgs_only" && $(header["recommended_training_use"]) != "pre_training_wgs_only") next
      if (mode == "rna_pair" && $(header["recommended_training_use"]) != "post_training_wgs_rna_pair") next
      print $(header["sample_id"])
    }
  ' "$ASSIGN_TSV" > "$out"

  echo "[OK] $out ($(wc -l < "$out") samples)" >&2
  printf "%s\n" "$out"
}

write_order() {
  local phase="$1"
  local pair_list="$2"
  local wgs_list="$3"
  local out="$OUT_ROOT/phase_$phase/execution_order.pair_first.platform_$PLATFORM.tsv"
  {
    printf "phase\tstep\tmode\tsample_list\tsample_count\n"
    printf "%s\tpair_wgs\trna_pair\t%s\t%s\n" "$phase" "$pair_list" "$(wc -l < "$pair_list")"
    printf "%s\tpair_rna\trna_pair\t%s\t%s\n" "$phase" "$pair_list" "$(wc -l < "$pair_list")"
    printf "%s\twgs_only\twgs_only\t%s\t%s\n" "$phase" "$wgs_list" "$(wc -l < "$wgs_list")"
  } > "$out"
  echo "[OK] $out" >&2
}

run_stage() {
  local stage="$1"
  local list_mode="$2"
  local phase="$3"
  local sample_list="$4"
  local command="$5"

  if [[ ! -s "$sample_list" ]]; then
    echo "[INFO] empty list; stage skipped: $stage phase=$phase" >&2
    return 0
  fi
  [[ -n "$command" ]] || {
    echo "[ERROR] RUN_PREPROCESS=1 requires a command for stage=$stage" >&2
    return 1
  }

  (
    export SAMPLE_LIST="$sample_list"
    export SAMPLES_SPACE
    SAMPLES_SPACE="$(paste -sd ' ' "$sample_list")"
    export PHASE="$phase"
    export MODE="$list_mode"
    export PLATFORM
    export PRIORITY_STAGE="$stage"
    export PACKAGE_ROOT="$PKG_ROOT"
    echo "[INFO] running stage=$stage phase=$phase samples=$(wc -l < "$sample_list")" >&2
    if [[ -n "$PROCESS_CODE_DIR" ]]; then
      cd "$PROCESS_CODE_DIR"
    fi
    bash -lc "$command"
  )
}

phase_values() {
  if [[ "$PHASE_ARG" == "all" ]]; then
    awk -F '\t' 'NR>1 {print $2}' "$ASSIGN_TSV" | sort -u
  else
    normalise_phase "$PHASE_ARG"
  fi
}

mapfile -t PHASES < <(phase_values)
[[ "${#PHASES[@]}" -gt 0 ]] || { echo "[ERROR] no phase found" >&2; exit 1; }

if [[ "$MODE" != "pair_first" ]]; then
  for phase in "${PHASES[@]}"; do
    list="$(build_list "$phase" "$MODE")"
    if [[ "$RUN_PREPROCESS" == "1" ]]; then
      run_stage "standard" "$MODE" "$phase" "$list" "$PREPROCESS_CMD"
    fi
  done
  exit 0
fi

if [[ "$RUN_PREPROCESS" == "1" ]]; then
  [[ -n "$PAIR_WGS_CMD" ]] || { echo "[ERROR] PAIR_WGS_CMD is required for MODE=pair_first" >&2; exit 1; }
  [[ -n "$PAIR_RNA_CMD" ]] || { echo "[ERROR] PAIR_RNA_CMD is required for MODE=pair_first" >&2; exit 1; }
  [[ -n "$WGS_ONLY_CMD" ]] || { echo "[ERROR] WGS_ONLY_CMD is required for MODE=pair_first" >&2; exit 1; }
fi

PAIR_LISTS=()
WGS_LISTS=()
for phase in "${PHASES[@]}"; do
  pair_list="$(build_list "$phase" "rna_pair")"
  wgs_list="$(build_list "$phase" "wgs_only")"
  write_order "$phase" "$pair_list" "$wgs_list"
  PAIR_LISTS+=("$pair_list")
  WGS_LISTS+=("$wgs_list")
done

if [[ "$RUN_PREPROCESS" == "0" ]]; then
  echo "[INFO] pair-first lists only; set RUN_PREPROCESS=1 with three stage commands to execute" >&2
  exit 0
fi

for i in "${!PHASES[@]}"; do
  run_stage "pair_wgs" "rna_pair" "${PHASES[$i]}" "${PAIR_LISTS[$i]}" "$PAIR_WGS_CMD"
done
for i in "${!PHASES[@]}"; do
  run_stage "pair_rna" "rna_pair" "${PHASES[$i]}" "${PAIR_LISTS[$i]}" "$PAIR_RNA_CMD"
done
for i in "${!PHASES[@]}"; do
  run_stage "wgs_only" "wgs_only" "${PHASES[$i]}" "${WGS_LISTS[$i]}" "$WGS_ONLY_CMD"
done
