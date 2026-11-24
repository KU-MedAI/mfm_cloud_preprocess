#!/usr/bin/env bash
# 03 — VCF 준비.
#
# 역할:
# 1) 선택적으로 cohort common VCF를 만든다(CREATE_COMMON_VCF=1일 때).
# 2) sample별 seq/tag VCF를 만든다.
#
# short-read는 Illumina 1000G phased panel에서 seq/tag VCF를 만든다.
# long-read clair3는 stage 04에서 만들고, public은 여기서 원본을 검증·정규화한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
COHORT=$(sample_cohort "$S")
COMMON=$(sample_common_vcf_dir "$S")
PERSAMPLE=$(sample_persample_vcf_dir "$S")
CONCAT=${VCF_CONCAT_DIR:-$(sample_work_root "$S")/vcf_concat}
VCF_SOURCE=$(sample_vcf_source "$S")
mkdir -p "$(sample_work_root "$S")/logs" "$COMMON" "$PERSAMPLE" "$CONCAT"

if is_long_read_sample "$S"; then
  if [[ "$VCF_SOURCE" == clair3 ]]; then
    log "03 skip: long-read Clair3 VCF는 stage 04에서 생성함(sample=$S)"
    exit 0
  fi
  [[ "$VCF_SOURCE" == public ]] || die "지원하지 않는 long-read VCF source: sample=$S source=$VCF_SOURCE"
  mapfile -t PUBLIC_VCFS < <(sample_public_vcfs "$S")
  [[ "${#PUBLIC_VCFS[@]}" -gt 0 ]] || die "public_vcf_path가 비어 있음: sample=$S"
  PUBLIC_SAMPLE=$(sample_public_vcf_sample_id "$S")
  PUBLIC_REF=$(sample_vcf_reference "$S")
  PUBLIC_CALLSET=$(sample_vcf_callset "$S")
  PHASE_SOURCE_RESOLVED=$(sample_phase_source "$S")
  args=(
    prepare-public-vcf
    --sample "$S"
    --source-sample "$PUBLIC_SAMPLE"
    --reference-fasta "$PUBLIC_REF"
    --output-dir "$PERSAMPLE"
    --callset "$PUBLIC_CALLSET"
    --phase-source "$PHASE_SOURCE_RESOLVED"
    --chroms "$CHROMS"
    --sequence-variant-policy "$SEQUENCE_VARIANT_POLICY"
    --bcftools-bin "$BCFTOOLS_BIN"
  )
  for PUBLIC_VCF in "${PUBLIC_VCFS[@]}"; do args+=(--source-vcf "$PUBLIC_VCF"); done
  [[ "${PUBLIC_VCF_OVERWRITE:-0}" == 1 ]] && args+=(--overwrite)
  log "03 public VCF 준비: sample=$S source_sample=$PUBLIC_SAMPLE callset=$PUBLIC_CALLSET"
  run_cli "${args[@]}"
  require_file "$PERSAMPLE/$S.seq.vcf.gz"
  require_file "$PERSAMPLE/$S.seq.vcf.gz.tbi"
  require_file "$PERSAMPLE/$S.tag.vcf.gz"
  require_file "$PERSAMPLE/$S.tag.vcf.gz.tbi"
  require_file "$PERSAMPLE/$S.seq.vcf.source.json"
  log "03 완료: $S public VCF"
  exit 0
fi
[[ "$VCF_SOURCE" == illumina_panel ]] || die "short-read는 vcf_source=illumina_panel만 지원함: sample=$S source=$VCF_SOURCE"

if [[ "$CREATE_COMMON_VCF" == 1 ]]; then
  [[ -n "$SAMPLE_LIST" && -s "$SAMPLE_LIST" ]] || die "CREATE_COMMON_VCF=1 requires SAMPLE_LIST"
  log "03 common VCF 생성: cohort=$COHORT out=$COMMON"
  run_cli extract-cohort-common-vcf \
    --cohort "$COHORT" \
    --sample-list "$SAMPLE_LIST" \
    --output-dir "$COMMON" \
    --chroms "$CHROMS" \
    --panel-base "$PANEL_LOCAL_DIR" \
    --panel-file-template "$PANEL_FILE_TEMPLATE" \
    --filter-status PASS,. \
    --min-ac 1 \
    --jobs "$VCF_CHROM_JOBS" \
    --threads-per-chrom "$VCF_THREADS_PER_CHROM" \
    --bcftools-bin "$BCFTOOLS_BIN"
fi

log "03 per-sample VCF 생성: sample=$S out=$PERSAMPLE"
vcf_args=(
  extract-per-sample-vcfs
  --extraction-mode "$VCF_EXTRACTION_MODE"
  --cohort "$COHORT"
  --common-vcf-dir "$COMMON"
  --output-dir "$PERSAMPLE"
  --sample "$S"
  --chroms "$CHROMS"
  --sequence-variant-policy "$SEQUENCE_VARIANT_POLICY"
  --concat-dir "$CONCAT"
  --bcftools-bin "$BCFTOOLS_BIN"
  --jobs "$VCF_JOBS"
)
if [[ "${SAMPLE_VCF_FAST_VALIDATE:-0}" == 1 ]]; then
  vcf_args+=(--fast-validate)
fi
run_cli "${vcf_args[@]}"

require_file "$PERSAMPLE/$S.seq.vcf.gz"
require_file "$PERSAMPLE/$S.tag.vcf.gz"
log "03 완료: $S"
