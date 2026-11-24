#!/usr/bin/env bash
# 08 — de novo assembly contig anchoring QC.
#
# assembly contig를 reference에 align해 contig가 어느 reference chromosome에
# 대응되는지 QC한다. 아직 chromosome-scale scaffolding/promote는 수행하지 않는다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
is_long_read_sample "$S" || die "$S is not a long-read sample: platform=$(sample_platform "$S")"
ROOT=$(sample_final_root "$S")
COHORT=$(sample_cohort "$S")
LABEL=$(denovo_label "$S")
ASSEMBLY_DIR="$ROOT/$S/denovo_assembly/$LABEL"
[[ -d "$ASSEMBLY_DIR" ]] || die "missing de novo assembly dir: $ASSEMBLY_DIR"

args=(
  anchor-denovo-assembly
  --cohort "$COHORT"
  --output-root "$ROOT"
  --sample "$S"
  --assembly-dir "$ASSEMBLY_DIR"
  --label "$LABEL"
  --reference-fasta "$ASSEMBLY_ANCHOR_REF"
  --minimap2-bin "$MINIMAP2_BIN"
  --preset "$ASSEMBLY_ANCHOR_PRESET"
  --threads "$ASSEMBLY_ANCHOR_THREADS"
  --min-mapq "$ASSEMBLY_MIN_ANCHOR_MAPQ"
  --min-contig-len "$DENOVO_MIN_CONTIG_LEN"
)
if [[ "${ASSEMBLY_ANCHOR_OVERWRITE:-0}" == 1 ]]; then
  args+=(--overwrite)
fi
if [[ "${ASSEMBLY_ANCHOR_DRY_RUN:-0}" == 1 ]]; then
  args+=(--dry-run)
fi

log "08 assembly anchoring QC: sample=$S label=$LABEL ref=$ASSEMBLY_ANCHOR_REF"
run_cli "${args[@]}"
log "08 완료: $S"
