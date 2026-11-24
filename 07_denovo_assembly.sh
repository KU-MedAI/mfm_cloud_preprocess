#!/usr/bin/env bash
# 07 — long-read de novo assembly benchmark/production-prep.
#
# 현재 구현은 hifiasm --ont 우선 테스트용이다. 이 단계는 assembly FASTA를
# 생성/검증하지만, 임의 assembly contig를 HAP1/HAP2 chromosome FASTA로
# 조용히 승격하지 않는다. 26 anchoring/scaffolding QC 이후에만 40단계로
# 연결해야 한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
is_long_read_sample "$S" || die "$S is not a long-read sample: platform=$(sample_platform "$S")"
ALIGN=$(sample_alignment "$S")
[[ -n "$ALIGN" ]] || die "alignment path is empty for $S"
require_file "$ALIGN"
ROOT=$(sample_final_root "$S")
COHORT=$(sample_cohort "$S")
PLATFORM=$(normalize_platform "$(sample_platform "$S")")
LABEL=$(denovo_label "$S")
mkdir -p "$ROOT/$S/logs"

args=(
  run-denovo-assembly
  --cohort "$COHORT"
  --output-root "$ROOT"
  --reference-fasta "$REF_DECODE"
  --platform "$PLATFORM"
  --alignment "$ALIGN"
  --sample "$S"
  --contigs "${DENOVO_CONTIGS:-$CHROMS}"
  --label "$LABEL"
  --assembler "$DENOVO_ASSEMBLER"
  --hifiasm-bin "$HIFIASM_BIN"
  --samtools-bin "$SAMTOOLS_BIN"
  --threads "$DENOVO_THREADS"
  --fastq-threads "$DENOVO_FASTQ_THREADS"
  --exclude-flags "$DENOVO_EXCLUDE_FLAGS"
  --min-contig-len "$DENOVO_MIN_CONTIG_LEN"
)

if [[ "$DENOVO_KEEP_FASTQ" == 1 ]]; then
  args+=(--keep-fastq)
else
  args+=(--no-keep-fastq)
fi
if [[ "$DENOVO_KEEP_GFA" == 1 ]]; then
  args+=(--keep-gfa)
else
  args+=(--no-keep-gfa)
fi
if [[ "$DENOVO_KEEP_BIN" == 1 ]]; then
  args+=(--keep-bin)
else
  args+=(--no-keep-bin)
fi
if [[ "${DENOVO_OVERWRITE:-0}" == 1 ]]; then
  args+=(--overwrite)
fi
if [[ "${DENOVO_DRY_RUN:-0}" == 1 ]]; then
  args+=(--dry-run)
fi

log "07 de novo assembly: sample=$S platform=$PLATFORM assembler=$DENOVO_ASSEMBLER label=$LABEL contigs=${DENOVO_CONTIGS:-$CHROMS}"
run_cli "${args[@]}"
log "07 완료: $S"
