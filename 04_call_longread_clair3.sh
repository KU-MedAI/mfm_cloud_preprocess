#!/usr/bin/env bash
# 04 — long-read Clair3 self-calling + WhatsHap haplotagging.
#
# 입력: ONT/HiFi BAM/CRAM + reference.
# 출력: Clair3 VCF + HP/PS-tagged BAM.
# 현재 production 기본 경로가 아니다. 과거 self-call 결과 재현/비교를 위해
# vcf_source=clair3와 ENABLE_CLAIR3_SELF_CALL=1을 함께 명시한 sample에서만 실행한다.
# FASTA와 haplotagging은 같은 long-read BAM/CRAM에서 Clair3로 call/phasing한 VCF를 사용한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
is_long_read_sample "$S" || die "$S is not a long-read sample: platform=$(sample_platform "$S")"
[[ "$ENABLE_CLAIR3_SELF_CALL" == 1 ]] || die "04 Clair3 self-call은 기본 비활성임: sample=$S"
[[ "$(sample_vcf_source "$S")" == clair3 ]] || die "04 Clair3는 vcf_source=clair3 sample만 허용함: sample=$S"
ALIGN=$(sample_alignment "$S")
[[ -n "$ALIGN" ]] || die "alignment path is empty for $S"
require_file "$ALIGN"
ROOT=$(sample_final_root "$S")
COHORT=$(sample_cohort "$S")
PERSAMPLE=$(sample_persample_vcf_dir "$S")
PLATFORM=$(normalize_platform "$(sample_platform "$S")")
MODEL_PATH=$(clair3_model_path_for_sample "$S")
[[ -d "$MODEL_PATH" ]] || die "missing Clair3 model path for sample=$S platform=$PLATFORM: $MODEL_PATH"
SEQ_VCF_SOURCE=$(longread_sequence_vcf_source_for_sample "$S")
HAPLOTAG_VCF_SOURCE=$(longread_haplotag_vcf_source_for_sample "$S")
mkdir -p "$ROOT/$S/logs" "$PERSAMPLE"

args=(
  call-longread-clair3
  --cohort "$COHORT"
  --output-root "$ROOT"
  --reference-fasta "$REF_DECODE"
  --model-path "$MODEL_PATH"
  --platform "$PLATFORM"
  --alignment "$ALIGN"
  --sample "$S"
  --threads "$CLAIR3_THREADS"
  --output-threads "$CLAIR3_OUTPUT_THREADS"
  --jobs "$CLAIR3_JOBS"
  --chunk-size "$CLAIR3_CHUNK_SIZE"
  --min-mq "$CLAIR3_MIN_MQ"
  --min-coverage "$CLAIR3_MIN_COVERAGE"
  --qual "$CLAIR3_QUAL"
  --clair3-bin "$CLAIR3_PY"
  --python-bin "$CLAIR3_PYTHON_BIN"
  --whatshap-bin "$CLAIR3_WHATSHAP_BIN"
  --samtools-bin "$CLAIR3_SAMTOOLS_BIN"
  --pypy-bin "$CLAIR3_PYPY_BIN"
  --parallel-bin "$CLAIR3_PARALLEL_BIN"
  --tabix-bin "$CLAIR3_TABIX_BIN"
)

if [[ -n "$CLAIR3_BED_FN" ]]; then
  args+=(--bed-fn "$CLAIR3_BED_FN")
else
  args+=(--contigs "$CHROMS")
fi
if [[ "$SEQ_VCF_SOURCE" == "clair3" ]]; then
  args+=(--export-sequence-vcf --sequence-vcf-dir "$PERSAMPLE" --sequence-vcf-export-mode "$LONGREAD_SEQ_VCF_EXPORT_MODE")
fi
if [[ "$LONGREAD_ALLOW_UNPHASED_SEQUENCE_VCF" == 1 ]]; then
  args+=(--allow-unphased-sequence-vcf)
fi
if [[ "$CLAIR3_ALLOW_EMPTY_HAPLOTAG" == 1 ]]; then
  args+=(--allow-empty-haplotag)
fi
if [[ "$CLAIR3_ALLOW_PILEUP_FALLBACK" == 1 ]]; then
  args+=(--allow-pileup-vcf-fallback)
fi
if [[ "${CLAIR3_OVERWRITE:-0}" == 1 ]]; then
  args+=(--overwrite)
fi
if [[ "$CLAIR3_DRY_RUN" == 1 ]]; then
  args+=(--dry-run --overwrite)
fi

log "04 Clair3 long-read call/haplotag: sample=$S platform=$PLATFORM seq_vcf=$SEQ_VCF_SOURCE haplotag_vcf=$HAPLOTAG_VCF_SOURCE model=$MODEL_PATH"
run_cli "${args[@]}"
if [[ "$SEQ_VCF_SOURCE" == clair3 && "$CLAIR3_DRY_RUN" != 1 ]]; then
  # Clair3 VCF는 실행 단계 02 시점에는 존재하지 않는다. 생성 직후 실제 REF allele을
  # 개인 FASTA reference와 대조해 stage 09 이전에 잘못된 sample을 중단한다.
  SEQUENCE_REF=$(sample_sequence_reference "$S")
  SEQ_VCF="$PERSAMPLE/$S.seq.vcf.gz"
  [[ -s "$SEQ_VCF" && -s "$SEQ_VCF.tbi" ]] || die "Clair3 sequence VCF 또는 index 없음: $SEQ_VCF"
  AUTOSOME_REGIONS=${AUTOSOMES// /,}
  VCF_REF_LOG="$ROOT/$S/logs/20.clair3_vcf_reference.log"
  if "$BCFTOOLS_BIN" view -r "$AUTOSOME_REGIONS" -Ou "$SEQ_VCF" 2>"$VCF_REF_LOG" \
      | "$BCFTOOLS_BIN" norm --check-ref e --fasta-ref "$SEQUENCE_REF" -Ou - 2>>"$VCF_REF_LOG" >/dev/null; then
    log "04 Clair3 VCF REF 계약 통과: sample=$S autosomes=$AUTOSOME_REGIONS"
  else
    detail=$(tr '\n' ' ' < "$VCF_REF_LOG" | cut -c1-500)
    die "Clair3 VCF REF 계약 실패: sample=$S ${detail:-bcftools norm --check-ref e 실패}"
  fi
fi
log "04 완료: $S"
