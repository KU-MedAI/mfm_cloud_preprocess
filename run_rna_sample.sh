#!/usr/bin/env bash
# RNA sample 1개를 Stage 13–26 순서로 실행한다. WGS Stage 01–12가 선행되어야 한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./lib/fixed_cpu_entry.sh
pipeline_fixed_cpu_entry "${1:?sample이 필요함}" bash "$(pwd)/run_rna_sample.sh" "$@"
source ./config_v2.sh
# 실물 RNA STAR host RSS 약150GB를 수용하며 명시값과 전체/Slurm 예산은 유지함.
if [[ ${MFM_RNA_BACKEND:-cpu} == parabricks ]]; then
  export MFM_GPU_HOST_MEMORY_BYTES=${MFM_GPU_HOST_MEMORY_BYTES:-274877906944}
fi
source ./lib/state.sh
source ./lib/stage_runner.sh
source ./lib/hap_stage_group.sh
source ./lib/rna_index_stage.sh
source ./lib/rna_dag.sh

S=${1:?사용법: $0 <sample>}
apply_sample_payload "$S"
RNA_START_AT=${RNA_START_AT:-13}
RNA_ALLOWED_STAGES=(13 14 15 16 17 18 19 20 21 22 23 24 25 26)
pipeline_validate_stage_range RNA "$RNA_START_AT" "${RNA_STOP_AFTER:-}" "${RNA_ALLOWED_STAGES[@]}"

v2_require
require_manifest_sample_contract "$S"
[[ "$(sequence_source_for_sample "$S")" != assembly_direct ]] || \
  die "RNA assembly_direct는 multi-chain GTF lift 검증 전까지 지원하지 않음: sample=$S"

F=$(v2_final "$S")
W=$(v2_work "$S")
MTAG=$(sample_rna_evidence_label "$S")
R=$F/rna_personal/$MTAG
RNA=$(v2_rna "$S")
[[ "$RNA" != "-" && -s "$RNA" ]] || die "$S RNA alignment 없음: $RNA"
mkdir -p "$W/logs" "$W/rna" "$R" "$F/qc"
require_hap_contig_assets "$F/individual_fasta" "$CHROMS" 1
RNA_CONTIGS=$(chroms_space "$CHROMS")
[[ -n "$RNA_CONTIGS" ]] || die "$S RNA contig 목록이 비어 있음: CHROMS=$CHROMS"
SEQ_REF_FAI="$(sample_sequence_reference "$S").fai"
[[ -s "$SEQ_REF_FAI" ]] || die "$S HAP FASTA reference .fai 없음: $SEQ_REF_FAI"
RNA_DECODE_REF=$(named_path_or_default "$S" rna_reference_fasta "${RNA_DECODE_REFERENCE:-}")
SM_EXPECTED=$S
PUBLIC_ID=$(sample_public_vcf_sample_id "$S" 2>/dev/null || true)
[[ -n "$PUBLIC_ID" && "$PUBLIC_ID" != "$S" ]] && SM_EXPECTED="$S,$PUBLIC_ID"

EXTRA_FAI=-
[[ "$RNA_INDEX_ADDNL_SEQNS" == 1 ]] && EXTRA_FAI="$RNA_ADDNL_SEQNS_SOURCE.fai"
if [[ "$RNA_INDEX_ADDNL_SEQNS" == 1 ]]; then
  # 이 asset은 RNA profile fingerprint와 STAR index 양쪽에 쓰인다. Stage 02가
  # 먼저 검사하지만, RNA-only 재개에도 깨진/다른 reference를 쓰지 않게 방어한다.
  require_file "$RNA_ADDNL_SEQNS_SOURCE"
  require_file "$EXTRA_FAI"
  require_file "$RNA_ADDNL_SEQNS_FASTA"
fi
PROFILE_FP=$("$PYTHON" ./rna_profile.py "$GENCODE_GTF" "$EXTRA_FAI" "$RNA_INDEX_ADDNL_SEQNS") \
  || die "$S RNA profile fingerprint 실패"
ALIGN_FP=$("$PYTHON" ./rna_profile.py "$GENCODE_GTF" "$EXTRA_FAI" "$RNA_INDEX_ADDNL_SEQNS" align) \
  || die "$S RNA alignment fingerprint 실패"
CACHE_TAG="gtex_${GENCODE_TAG}_${ALIGN_FP:0:8}"
FP_ARGS=(--extra-contig-fai "${EXTRA_FAI/#-/}")
ADDNL_ARGS=()
[[ "$RNA_INDEX_ADDNL_SEQNS" == 1 ]] && \
  ADDNL_ARGS=(--addnl-fasta "$RNA_ADDNL_SEQNS_FASTA" --addnl-source "$RNA_ADDNL_SEQNS_SOURCE")
FORCE_ARGS=()
[[ "${RNA_FORCE:-0}" == 1 ]] && FORCE_ARGS=(--force)

RUN_STATUS_DIR=${RUN_STATUS_DIR:-$(sample_status_root "$S")}
mkdir -p "$RUN_STATUS_DIR"
PIPELINE_SAMPLE_ID=$S
PIPELINE_STAGE_LABEL=RNA
PIPELINE_START_AT=$RNA_START_AT
PIPELINE_STOP_AFTER=${RNA_STOP_AFTER:-}
STATE_FILE=$RUN_STATUS_DIR/$S.rna.state.tsv
STAGE_FILE=$RUN_STATUS_DIR/$S.rna.stages.tsv
CURRENT_STAGE=init
LAST_SUCCESS_STAGE=none
pipeline_install_sample_trap
T0=$(date +%s)
log "═══ RNA 시작: $S phase=$MTAG annotation=$GENCODE_TAG profile=${PROFILE_FP:0:8} align=${ALIGN_FP:0:8} ═══"

# 같은 번호에서 HAP1/HAP2 명령을 모두 provenance에 남긴 뒤에만 STOP_AFTER를 평가한다.
pipeline_run_part() {
  local saved=${PIPELINE_STOP_AFTER:-}
  PIPELINE_STOP_AFTER=
  pipeline_run_command "$@"
  PIPELINE_STOP_AFTER=$saved
}

SM_MATCH_ARGS=(--require-sm-match "$RNA_REQUIRE_SM_MATCH")
pipeline_run_command 13 13_rna_to_fastq.py \
  "$PYTHON" ./13_rna_to_fastq.py --sample "$S" --rna "$RNA" \
  --out-dir "$W/rna" --log-dir "$W/logs" --decode-reference "${RNA_DECODE_REF:-}" \
  --sm-expected "$SM_EXPECTED" "${SM_MATCH_ARGS[@]}" \
  --min-paired-fraction "$RNA_MIN_PAIRED_FRACTION" --samtools "$SAMTOOLS" \
  --ratio-qc-policy "$QC_RATIO_POLICY" \
  --preflight-json "$F/qc/rna_preflight.json" \
  --gtf "$GENCODE_GTF" --contigs "$RNA_CONTIGS" --hap-read-mode "$HAP_READ_MODE" \
  --collate-threads "$RNA_COLLATE_THREADS" --fastq-threads "$RNA_FASTQ_THREADS" \
  "${FORCE_ARGS[@]}"
FASTQ_QC=$W/rna/fastq_qc.json

for h in 1 2; do
  TD=$R/HAP$h
  HD=$F/individual_fasta/HAP$h
  mkdir -p "$TD"
  REUSE_ARGS=(--reuse)
  [[ "${RNA_FORCE:-0}" == 1 ]] && REUSE_ARGS=()
  if [[ "$h" == 1 ]]; then RUNNER=pipeline_run_part; else RUNNER=pipeline_run_command; fi
  "$RUNNER" 14 "14_lift_annotation.py:HAP$h" \
    "$PYTHON" ./14_lift_annotation.py --chain-dir "$HD" --gtf "$GENCODE_GTF" \
    --partner-chain-dir "$F/individual_fasta/HAP$((3-h))" --reference-fai "$SEQ_REF_FAI" \
    --out-gtf "$TD/annotation.gtf" --out-exons "$TD/exons.bygene.tsv" \
    --qc-json "$TD/gtf_lift_qc.json" --contigs "$RNA_CONTIGS" "${REUSE_ARGS[@]}"
done

D1=$R/HAP1
D2=$R/HAP2
PF_QC=$F/individual_fasta/personal_fasta_qc.json
[[ -s "$PF_QC" ]] || die "$S personal FASTA QC 없음: $PF_QC"
HAP_VCF=$("$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))["refok_vcf"])' "$PF_QC") || die "$S HAP FASTA 생성 VCF 없음"
HV=$R/hap_variants.tsv
AC=$R/allelic_counts.tsv
CNV_RD=$R/gene_rd.tsv
CNV_ARGS=()
[[ "${RNA_CNV_FILTER:-1}" != 1 ]] || CNV_ARGS=(--cnv-rd "$CNV_RD")
RNA_INDEX_FORCE_ARGS=()
rna_star_index_hap() {
  local h=$1 TD HD
  TD=$R/HAP$h
  HD=$F/individual_fasta/HAP$h
  pipeline_run_command 15 "15_star_index.py:HAP$h" \
    "$PYTHON" ./15_star_index.py --sample "$S" --hap "$h" --hap-fasta-dir "$HD" \
    --contigs "$RNA_CONTIGS" --gtf "$TD/annotation.gtf" --fastq-qc "$FASTQ_QC" \
    --work-dir "$W/rna" --log-dir "$W/logs" --cache-tag "$CACHE_TAG" \
    --sjdb-overhang "${RNA_SJDB_OVERHANG:-${STAR_SJDB_OVERHANG:-auto}}" "${ADDNL_ARGS[@]}" \
    --threads "$STAR_INDEX_THREADS" --star "$STAR_BIN" --samtools "$SAMTOOLS" "${RNA_INDEX_FORCE_ARGS[@]}"
}

rna_star_align_hap() {
  local h=$1 TD
  TD=$R/HAP$h
    pipeline_run_command 16 "16_star_align.py:HAP$h" \
      "$PYTHON" ./16_star_align.py --sample "$S" --hap "$h" \
      --index-json "$W/rna/hap$h.index.json" --fastq-qc "$FASTQ_QC" \
      --out-dir "$TD" --work-dir "$W/rna" --log-dir "$W/logs" \
      --threads "$STAR_ALIGN_THREADS" --star "$STAR_BIN"
}

rna_markdup_hap() {
  local h=$1 TD
  TD=$R/HAP$h
  pipeline_run_command 17 "17_mark_duplicates.py:HAP$h" \
    "$PYTHON" ./17_mark_duplicates.py --sample "$S" --hap "$h" \
    --align-json "$W/rna/hap$h.align.json" --out-dir "$TD" --work-dir "$W/rna" \
    --log-dir "$W/logs" --picard "$PICARD_BIN" --heap-gb "$RNA_MARKDUP_HEAP_GB" \
    --samtools "$SAMTOOLS" "${FORCE_ARGS[@]}"
}

rna_coverage_hap() {
  local h=$1 TD MD_BAM TFAI
  TD=$R/HAP$h
    MD_BAM=$("$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))["marked_bam"])' "$W/rna/hap$h.markdup.json") || return
    TFAI=$("$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))["target_fai"])' "$W/rna/hap$h.index.json") || return
    # pipeline1_update와 같은 HAP 전장 1-bp coverage를 contig별 worker로 계산한다.
    pipeline_run_command 18 "18_coverage_zarr.py:HAP$h" \
      "$PYTHON" ./18_coverage_zarr.py --bam "$MD_BAM" --fai "$TFAI" --out "$TD/coverage.zarr" \
      --windows-tsv "$F/windows/windows.tsv" --window-out "$TD/coverage_windows.zarr" \
      --contigs "$RNA_CONTIGS" --chunk "$RNA_ZARR_CHUNK" --power "$AG_NORM_POWER" \
      --clip "$AG_NORM_CLIP" --sample "$S" --hap "HAP$h" --samtools "$SAMTOOLS" \
      --jobs "${RNA_COVERAGE_JOBS:-auto}"
}

rna_qc_hap() {
  local h=$1 TD
  TD=$R/HAP$h
    pipeline_run_command 19 "19_rna_realign_qc.py:HAP$h" \
      "$PYTHON" ./19_rna_realign_qc.py hap --sample "$S" --hap "$h" --label "$MTAG" \
      --work-dir "$W/rna" --out-dir "$TD" --fastq-qc "$FASTQ_QC" \
      --lift-qc "$TD/gtf_lift_qc.json" --coverage "$TD/coverage.zarr" \
      --gencode-gtf "$GENCODE_GTF" --gencode-tag "$GENCODE_TAG" --star "$STAR_BIN" "${FP_ARGS[@]}"
}

rna_qc_summary() {
  pipeline_run_command 19 19_rna_realign_qc.py:summary \
    "$PYTHON" ./19_rna_realign_qc.py summary --sample "$S" --label "$MTAG" \
    --fastq-qc "$FASTQ_QC" --rna-dir "$R"
}

rna_variants() {
  pipeline_run_command 20 20_hap_variants.py \
    "$PYTHON" ./20_hap_variants.py \
    --vcf "$HAP_VCF" --wgs-bam "$(sample_haplotagged_alignment "$S")" \
    --reference "$(sample_analysis_reference "$S")" \
    --hap-dir1 "$F/individual_fasta/HAP1" --hap-dir2 "$F/individual_fasta/HAP2" \
    --gtf "$GENCODE_GTF" --contigs "$RNA_CONTIGS" \
    --min-mapq "$HAP_VARIANT_MIN_MAPQ" --min-depth "$HAP_VARIANT_MIN_DEPTH" \
    --min-alt-frac "$HAP_VARIANT_MIN_ALT_FRAC" --threads "$HAP_VARIANT_THREADS" \
    --out "$HV" --qc-json "$R/hap_variants_qc.json" || return
  [[ -s "$HV" ]] || die "$S HAP variant TSV 없음: $HV"
  [[ -s "$R/hap_variants_qc.json" ]] || die "$S HAP variant QC 없음: $R/hap_variants_qc.json"
}

rna_allelic() {
  pipeline_run_command 21 21_allelic_counts.py \
    "$PYTHON" ./21_allelic_counts.py \
    --bam1 "$D1/$S.rna.personal.HAP1.bam" --bam2 "$D2/$S.rna.personal.HAP2.bam" \
    --exons1 "$D1/exons.bygene.tsv" --exons2 "$D2/exons.bygene.tsv" \
    --hap-variants "$HV" --out "$R/allelic_counts.tsv" || return
  [[ -s "$AC" ]] || die "$S ASE 결과 파일 없음: $AC"
  for sidecar in "${AC%.tsv}.split_half.tsv" "${AC%.tsv}.both_haps.tsv" \
                 "${AC%.tsv}.mmap_adjusted.tsv" "${AC%.tsv}.meta.json"; do
    [[ -s "$sidecar" ]] || die "$S ASE sidecar 없음: $sidecar"
  done
}

rna_cnv() {
  local W1 W2
  W1=$(sample_hap_alignment "$S" 1) || return
  W2=$(sample_hap_alignment "$S" 2) || return
  if [[ "${RNA_CNV_FILTER:-1}" == 1 ]]; then
    # 파일 존재 조건으로 skip하지 않는다. 필수 입력 누락은 Stage 22 명령에서
    # 실패하여 단계 로그/provenance에 남고 해당 샘플만 중단된다.
    pipeline_run_command 22 22_cnv_filter.py \
      "$PYTHON" ./22_cnv_filter.py \
      --exons "$D1/exons.bygene.tsv" --bam "$W1" --reference "$(sample_hap_target "$S" 1)" \
      --exons "$D2/exons.bygene.tsv" --bam "$W2" --reference "$(sample_hap_target "$S" 2)" \
      --min-rd "$RNA_CNV_MIN_RD" --max-rd "$RNA_CNV_MAX_RD" \
      --samtools "$SAMTOOLS" --out "$CNV_RD" --jobs "${RNA_CNV_JOBS:-auto}" || return
    CNV_ARGS=(--cnv-rd "$CNV_RD")
  else
    pipeline_skip_stage 22 "RNA_CNV_FILTER=0 (명시적 비활성화)" 1
  fi
}

if [[ "${MFM_RNA_OUTPUT_MODE:-full}" != coverage_only ]] && rna_dag_enabled; then
  pipeline_rna_dag
  rna_qc_summary
  pipeline_rna_dag_commit
else
  pipeline_run_rna_index_stage 15 rna_star_index_hap
  pipeline_rna_index_group 16 rna_star_align_hap "$RNA_ALIGN_HAP_JOBS"
  pipeline_rna_index_group 17 rna_markdup_hap "$RNA_MARKDUP_HAP_JOBS"
  pipeline_rna_index_group 18 rna_coverage_hap "$RNA_COVERAGE_HAP_JOBS"
  saved_stop=${PIPELINE_STOP_AFTER:-};PIPELINE_STOP_AFTER=
  pipeline_run_hap_group 19 rna_qc_hap 2
  PIPELINE_STOP_AFTER=$saved_stop
  rna_qc_summary
  if [[ "${MFM_RNA_OUTPUT_MODE:-full}" != coverage_only ]]; then
    rna_variants
    rna_allelic
    rna_cnv
  else
    for stage in 20 21 22; do pipeline_skip_stage "$stage" "coverage_only: WGS변이근거·ASE·CNV 비적용" 1; done
  fi
fi

if [[ "${MFM_RNA_OUTPUT_MODE:-full}" != coverage_only ]]; then
pipeline_run_command 23 23_ase_qc.py \
  "$PYTHON" ./23_ase_qc.py --counts "$AC" --hap-variants "$HV" \
  --mmap-adjusted "${AC%.tsv}.mmap_adjusted.tsv" "${CNV_ARGS[@]}" \
  --min-cnv-rd "$RNA_CNV_MIN_RD" --max-cnv-rd "$RNA_CNV_MAX_RD" \
  --min-total-count "${ASE_MIN_TOTAL_COUNT:-20}" --min-split-count "${ASE_MIN_SPLIT_COUNT:-5}" \
  --max-split-ratio-delta "${ASE_MAX_SPLIT_RATIO_DELTA:-0.25}" --max-fdr "${ASE_MAX_FDR:-0.05}" \
  --min-informative-fraction "${ASE_MIN_INFORMATIVE_FRACTION:-0.01}" \
  --min-eligible-genes "${ASE_MIN_ELIGIBLE_GENES:-100}" \
  --min-split-correlation "${ASE_MIN_SPLIT_CORRELATION:-0.30}"

else
  pipeline_skip_stage 23 "coverage_only: ASE QC 비적용" 1
fi

pipeline_run_command 24 24_stamp_rna_metadata.py \
  "$PYTHON" ./24_stamp_rna_metadata.py "$F"
pipeline_run_command 25 25_validate_sample_package.sh \
  ./25_validate_sample_package.sh "$S" rna

if [[ "${CLEAN_INTERMEDIATES_ON_SUCCESS:-1}" == 1 ]]; then
  pipeline_run_stage 26 26_cleanup_intermediates.sh rna
else
  pipeline_skip_stage 26 "CLEAN_INTERMEDIATES_ON_SUCCESS=0" 1
fi

pipeline_complete_sample
log "═══ RNA 완료: $S / $(( ($(date +%s)-T0)/60 ))분 ═══"
