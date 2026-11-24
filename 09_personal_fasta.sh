#!/usr/bin/env bash
# 09 — HAP1/HAP2 personal FASTA + chain 생성.
#
# 변이 집합은 sample별 seq VCF를 사용한다. v2 기본은 SNV+indel이다.
# short-read: 03_prepare_vcf.sh가 Illumina panel에서 만든 seq VCF 사용.
# long-read production: 03_prepare_vcf.sh가 검증한 제공 public VCF 사용.
# Clair3 seq VCF는 명시적 재현/비교 실행에서만 사용한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
SEQUENCE_SOURCE_RESOLVED=$(sequence_source_for_sample "$S")
ROOT=$(sample_final_root "$S")
COHORT=$(sample_cohort "$S")
PERSAMPLE=$(sample_persample_vcf_dir "$S")
VCF_SOURCE=$(sample_vcf_source "$S")
SEQUENCE_REF=$(sample_sequence_reference "$S")
FIRST_CHROM=$(first_chrom "$CHROMS")
if is_long_read_sample "$S" && [[ "$SEQUENCE_SOURCE_RESOLVED" == "denovo" ]]; then
  LABEL=$(denovo_label "$S")
  DENOVO_DIR="$ROOT/$S/denovo_assembly/$LABEL"
  require_file "$DENOVO_DIR/denovo_assembly_qc.json"
  if ! grep -q '"hap_resolved_ready": true' "$DENOVO_DIR/denovo_assembly_qc.json"; then
    die "de novo assembly exists but HAP1/HAP2 FASTA is not ready. Run/inspect 26 anchoring and decide scaffolding/promotion policy: $DENOVO_DIR"
  fi
  die "de novo HAP FASTA promotion into individual_fasta/HAP1,HAP2 is intentionally not automatic yet: $DENOVO_DIR"
fi
require_file "$PERSAMPLE/$S.seq.vcf.gz"
if is_long_read_sample "$S" && [[ "$VCF_SOURCE" == clair3 ]]; then
  require_file "$PERSAMPLE/$S.seq.vcf.source.json"
  log "09 long-read FASTA VCF source: Clair3 self-called VCF export($PERSAMPLE/$S.seq.vcf.gz)"
elif is_long_read_sample "$S"; then
  require_file "$PERSAMPLE/$S.seq.vcf.source.json"
  log "09 long-read FASTA VCF source: $VCF_SOURCE reference=$SEQUENCE_REF"
fi
log "09 personal FASTA 생성: sample=$S root=$ROOT"
run_cli create-personal-fastas \
  --cohort "$COHORT" \
  --seq-vcf-dir "$PERSAMPLE" \
  --output-root "$ROOT" \
  --reference-fasta "$SEQUENCE_REF" \
  --sample "$S" \
  --chroms "$CHROMS" \
  --sequence-variant-policy "$SEQUENCE_VARIANT_POLICY" \
  --bcftools-bin "$BCFTOOLS_BIN" \
  --samtools-bin "$SAMTOOLS_BIN" \
  --jobs 1 \
  --sanitize-invalid-bases
require_file "$ROOT/$S/individual_fasta/HAP1/$FIRST_CHROM.fa"
require_file "$ROOT/$S/individual_fasta/HAP2/$FIRST_CHROM.fa"
# Stage 09이 만든 영구 자산의 hash를 이 시점에 고정한다. final validator나
# cleanup 직전에 manifest를 새로 만들면 그 사이의 변조를 정상 상태로 오인할 수 있다.
"$PYTHON" ./tools/personal_assets.py manifest \
  --individual-fasta-dir "$ROOT/$S/individual_fasta" \
  --contigs "$(chroms_space "$CHROMS")" \
  --reference-fasta "$SEQUENCE_REF"
log "09 완료: $S"
