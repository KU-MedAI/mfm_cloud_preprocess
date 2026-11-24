#!/usr/bin/env bash
# 06 — haplotagged BAM 준비.
#
# long-read: Clair3 self-call 또는 검증된 public phased VCF를 사용한다.
# short-read: Illumina VCF 기반 WhatsHap haplotagging을 실행한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
ROOT=$(sample_final_root "$S")
COHORT=$(sample_cohort "$S")
PLATFORM=$(normalize_platform "$(sample_platform "$S")")
# canonical BAM을 새로 정렬했으면 analysis reference, 기존 정렬을 보존했으면
# source decode reference를 사용한다. 전체 BAM contig·unmapped read를 유지하며
# reference 범위는 haplotag_contract.py가 검증한다.
DECODE_REF=$(sample_alignment_reference "$S")
SOURCE_SAMPLE_ID=$(sample_source_sample_id "$S")
BAM=$(sample_haplotagged_alignment "$S")
VCF_SOURCE=$(sample_vcf_source "$S")
HAP_DIR=$(dirname "$BAM")
SUMMARY_TSV="$HAP_DIR/$S.haplotag.summary.tsv"
METADATA_JSON="$HAP_DIR/$S.haplotag.metadata.json"

# WhatsHap 산출물은 BAM·index·summary·metadata가 함께 있어야 재사용 가능하다.
# 일부만 남은 경우에는 기본적으로 중단한다. 사용자가 명시적으로
# HAPLOTAG_REBUILD_PARTIAL=1을 준 경우에만 남은 파일을 timestamped backup으로
# 옮긴 뒤 새 산출물을 만든다. 원본을 덮어 쓰거나 조용히 삭제하지 않는다.
HAPLOTAG_APPEND_OVERWRITE=0
prepare_haplotag_output() {
  local candidate backup stamp
  local -a outputs=(
    "$BAM" "$BAM.bai" "$BAM.csi" "${BAM%.bam}.bai" "${BAM%.cram}.crai"
    "$SUMMARY_TSV" "$METADATA_JSON"
    "$BAM.part" "$HAP_DIR/$S.taglist.tsv.gz" "$HAP_DIR/$S.haplotag_qc.json"
  )
  if [[ -s "$BAM" ]] && have_index "$BAM" \
      && [[ -s "$SUMMARY_TSV" && -s "$METADATA_JSON" \
            && -s "$HAP_DIR/$S.taglist.tsv.gz" && -s "$HAP_DIR/$S.haplotag_qc.json" ]]; then
    if [[ "${MFM_HAPLOTAG_KEEP_SHARDS:-0}" != 1 ]]; then
      PYTHONPATH="$PIPELINE_DIR/src${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON_BIN" -m integrate_wgs_rna.haplotag_cleanup \
        --qc-json "$HAP_DIR/$S.haplotag_qc.json" --samtools "$SAMTOOLS_BIN" \
        > /dev/null
    fi
    log "[skip] 06 complete haplotag outputs exist: $BAM"
    return 1
  fi
  for candidate in "${outputs[@]}"; do
    if [[ -e "$candidate" ]]; then
      if [[ "${HAPLOTAG_REBUILD_PARTIAL:-0}" == preserve ]]; then
        [[ "${MFM_HAPLOTAG_CPU_BUDGET:-1}" -gt 1 ]] || die "preserve 재개는 병렬 haplotag에만 허용됩니다"
        log "06 기존 부분 결과를 유지하고 병렬 완성 후보를 생성함: $HAP_DIR"
        HAPLOTAG_APPEND_OVERWRITE=1
        break
      fi
      if [[ "${HAPLOTAG_REBUILD_PARTIAL:-0}" != 1 ]]; then
        die "06 partial haplotag output detected: $candidate; inspect it or rerun with HAPLOTAG_REBUILD_PARTIAL=1 (files will be backed up, never deleted)"
      fi
      stamp=$(date '+%Y%m%dT%H%M%S%z')
      backup="$HAP_DIR/partial_backup/$S.$stamp"
      mkdir -p "$backup"
      for candidate in "${outputs[@]}"; do
        [[ -e "$candidate" ]] || continue
        mv -- "$candidate" "$backup/"
      done
      printf '%s\n' "reason=explicit HAPLOTAG_REBUILD_PARTIAL=1" > "$backup/README.txt"
      log "06 partial haplotag outputs backed up: $backup"
      HAPLOTAG_APPEND_OVERWRITE=1
      break
    fi
  done
  return 0
}

if is_long_read_sample "$S"; then
  if [[ "$VCF_SOURCE" == clair3 ]]; then
    [[ "$ENABLE_CLAIR3_SELF_CALL" == 1 ]] || die "06 Clair3 fallback은 기본 비활성임: sample=$S"
    if [[ -s "$BAM" ]] && have_index "$BAM"; then
      log "[skip] 21 long-read Clair3 haplotagged BAM exists: $BAM"
    else
      log "06 long-read Clair3 wrapper가 haplotag까지 수행함. 04 실행: $S"
      ./04_call_longread_clair3.sh "$S"
    fi
  elif [[ "$VCF_SOURCE" == public ]]; then
    ALIGN=$(sample_alignment "$S")
    PERSAMPLE=$(sample_persample_vcf_dir "$S")
    RAW_ROOT=$(sample_raw_root "$S")
    [[ -n "$RAW_ROOT" ]] || RAW_ROOT=$(dirname "$ALIGN")
    require_file "$ALIGN"
    require_file "$PERSAMPLE/$S.tag.vcf.gz"
    require_file "$PERSAMPLE/$S.tag.vcf.gz.tbi"
    require_file "$PERSAMPLE/$S.seq.vcf.source.json"
    if ! prepare_haplotag_output; then
      exit 0
    fi
    args=(
      create-haplotagged-alignments
      --cohort "$COHORT"
      --raw-root "$RAW_ROOT"
      --alignment "$ALIGN"
      --tag-vcf-dir "$PERSAMPLE"
      --output-root "$ROOT"
      --reference-fasta "$DECODE_REF"
      --platform "$PLATFORM"
      --sample "$S"
      --source-sample-id "$SOURCE_SAMPLE_ID"
      --output-prefix "{sample}.public_vcf_longread.haplotagged"
      --whatshap-bin "$WHATSHAP_BIN"
      --samtools-bin "$SAMTOOLS_BIN"
      --jobs 1
      --output-threads "$CLAIR3_OUTPUT_THREADS"
      --tag-supplementary
      --allow-longread-vcf-haplotag
      --summary-tsv "$SUMMARY_TSV"
      --metadata-json "$METADATA_JSON"
    )
    # Legacy PUBLIC_VCF_HAPLOTAG_OVERWRITE does not overwrite in place. A partial
    # output must first be moved by the explicit safe-resume switch above.
    [[ "$HAPLOTAG_APPEND_OVERWRITE" == 1 ]] && args+=(--overwrite)
    log "06 public VCF WhatsHap: sample=$S callset=$(sample_vcf_callset "$S")"
    run_cli "${args[@]}"
  else
    die "지원하지 않는 long-read VCF source: sample=$S source=$VCF_SOURCE"
  fi
  require_file "$BAM"
  have_index "$BAM" || die "missing BAM index: $BAM"
  log "06 완료: $S long-read source=$VCF_SOURCE haplotagged BAM 확인"
  exit 0
fi

ALIGN=$(sample_alignment "$S")
require_file "$ALIGN"
# 명시적 BAM을 사용하므로 입력폴더가 없으면 BAM의 부모폴더를 기본값으로 쓴다.
RAW_ROOT=$(sample_raw_root "$S")
[[ -n "$RAW_ROOT" ]] || RAW_ROOT=$(dirname "$ALIGN")
PERSAMPLE=$(sample_persample_vcf_dir "$S")
require_file "$PERSAMPLE/$S.tag.vcf.gz"
if ! prepare_haplotag_output; then
  exit 0
fi
log "06 short-read WhatsHap haplotagging: sample=$S"
args=(
  create-haplotagged-alignments \
  --cohort "$COHORT" \
  --raw-root "$RAW_ROOT" \
  --alignment "$ALIGN" \
  --tag-vcf-dir "$PERSAMPLE" \
  --output-root "$ROOT" \
  --reference-fasta "$DECODE_REF" \
  --platform illumina \
  --sample "$S" \
  --source-sample-id "$SOURCE_SAMPLE_ID" \
  --alignment-pattern "${ALIGNMENT_PATTERN:-{sample}/wgs/{sample}.hg38.cram}" \
  --output-prefix "{sample}.illumina.hg38.haplotagged" \
  --whatshap-bin "$WHATSHAP_BIN" \
  --samtools-bin "$SAMTOOLS_BIN" \
  --jobs 1 \
  --output-threads "$CLAIR3_OUTPUT_THREADS" \
  --tag-supplementary \
  --summary-tsv "$SUMMARY_TSV" \
  --metadata-json "$METADATA_JSON"
)
[[ "$HAPLOTAG_APPEND_OVERWRITE" == 1 ]] && args+=(--overwrite)
run_cli "${args[@]}"
log "06 완료: $S"
