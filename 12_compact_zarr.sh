#!/usr/bin/env bash
# 12 — 최종 model input compact zarr 생성.
#
# personal_realigned mode에서 HAP별 readset과 per-base depth(uint32)를 저장한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
apply_sample_payload "$S"
if sample_sequence_only "$S"; then
  export PYTHONPATH="$PWD:$PWD/src:$PWD/vendor/MFM-DNARNAFM:${PYTHONPATH:-}"
  "$PYTHON" -m integrate_wgs_rna.sequence_only --root "$(sample_final_root "$S")/$S" --contigs "$(chroms_space "$CHROMS")"
  exit 0
fi
if sample_pileup_feature "$S"; then
  export PYTHONPATH="$PWD:$PWD/src:$PWD/vendor/MFM-DNARNAFM:${PYTHONPATH:-}"
  "$PYTHON" -m integrate_wgs_rna.pileup_features.export \
    --root "$(sample_final_root "$S")/$S" --alignment "$(sample_haplotagged_alignment "$S")" \
    --reference "$(sample_sequence_reference "$S")" --alignment-reference "$(sample_decode_reference "$S")" \
    --read-type "$(sample_read_type "$S")" --min-mapq "$MIN_MAPQ" --contigs "$(chroms_space "$CHROMS")"
  exit 0
fi
ROOT=$(sample_final_root "$S")
COHORT=$(sample_cohort "$S")
PROJECT=$(sample_project_id "$S")
READ_TYPE=$(sample_read_type "$S")
ROLE=$(sample_dataset_role "$S")
PLATFORM=$(normalize_platform "$(sample_platform "$S")")
SOURCE_SAMPLE_ID=$(sample_source_sample_id "$S")
DECODE_REF=$(sample_decode_reference "$S")
ALIGN=$(sample_haplotagged_alignment "$S")
HAP1_ALIGN=$(sample_hap_alignment "$S" 1)
HAP2_ALIGN=$(sample_hap_alignment "$S" 2)
if [[ -n "$MFM_STAGE12_SCRATCH_ROOT" ]]; then
  SCRATCH_KEY=$(printf '%s' "$ROOT/$S" | sha256sum | cut -c1-12)
  WORK="$MFM_STAGE12_SCRATCH_ROOT/$SCRATCH_KEY/$ROLE/$PROJECT/$S/compact"
else
  WORK=$(sample_work_root "$S")/compact
fi
PERSAMPLE=$(sample_persample_vcf_dir "$S")
SEQUENCE_REF=$(sample_sequence_reference "$S")
PERSONAL_FASTA_QC="$ROOT/$S/individual_fasta/personal_fasta_qc.json"

require_file "$ALIGN"
require_file "$HAP1_ALIGN"
require_file "$HAP2_ALIGN"
require_file "$ROOT/$S/windows/windows.tsv"
require_file "$PERSONAL_FASTA_QC"

args=(
  preprocess-full-shared-hap
  --sample-id "$S"
  --source-sample-id "$SOURCE_SAMPLE_ID"
  --dataset-role "$ROLE"
  --read-type "$READ_TYPE"
  --project-id "$PROJECT"
  --project-label "$COHORT"
  --platform "$PLATFORM"
  --reference-build "$REFERENCE_BUILD"
  --sample-status "$DEFAULT_SAMPLE_STATUS"
  --disease "$DEFAULT_DISEASE"
  --hap1-fasta-dir "$ROOT/$S/individual_fasta/HAP1"
  --hap1-chain-dir "$ROOT/$S/individual_fasta/HAP1"
  --hap2-fasta-dir "$ROOT/$S/individual_fasta/HAP2"
  --hap2-chain-dir "$ROOT/$S/individual_fasta/HAP2"
  --alignment "$ALIGN"
  --hap1-alignment "$HAP1_ALIGN"
  --hap2-alignment "$HAP2_ALIGN"
  --window-index-dir "$ROOT/$S/windows"
  --read-coordinate-mode personal_realigned
  --layout-name "$COMPACT_LAYOUT_NAME"
  --hap-read-mode "$HAP_READ_MODE"
  --contigs "$CHROMS"
  --output-sample-root "$ROOT/$S"
  --work-root "$WORK"
  --compact-mode direct
  --storage-policy model_input_only
  --validation-level "$COMPACT_VALIDATION_LEVEL"
  --compact-compressor "$COMPACT_COMPRESSOR"
  --compact-compression-level "$COMPACT_COMPRESSION_LEVEL"
  --jobs "$COMPACT_JOBS"
  --windows-per-chunk "$WINDOWS_PER_CHUNK"
  --window-size "$WINDOW_SIZE"
  --seed "$WINDOW_SEED"
  --min-mapq "$MIN_MAPQ"
  --reference-fasta "$DECODE_REF"
  --sequence-reference-fasta "$SEQUENCE_REF"
  --personal-fasta-qc "$PERSONAL_FASTA_QC"
  --samtools-bin "$SAMTOOLS_BIN"
  --bcftools-bin "$BCFTOOLS_BIN"
  --source-wgs "$ALIGN"
)
if [[ -s "$PERSAMPLE/$S.seq.vcf.gz" ]]; then
  args+=(--source-vcf "$PERSAMPLE/$S.seq.vcf.gz")
fi
if [[ "$FIXED_WINDOW_START" == 1 ]]; then
  args+=(--fixed-window-start)
elif [[ -n "$OFFSET_UNIT" ]]; then
  args+=(--offset-unit "$OFFSET_UNIT")
fi
if [[ -n "$COMPACT_MAX_WINDOWS" ]]; then
  args+=(--max-windows "$COMPACT_MAX_WINDOWS")
fi
if [[ "$COMPACT_ALLOW_UNPHASED_LONGREAD" == 1 ]]; then
  args+=(--allow-unphased-longread)
fi
if [[ "${COMPACT_OVERWRITE:-0}" == 1 ]]; then
  # stage 11 HAP CRAM을 재생성한 경우에는 이전 compact readset을 재사용하면 안 된다.
  # 기본값은 resume 유지이며, 재정렬 후속 실행에서만 명시적으로 활성화한다.
  args+=(--overwrite)
fi

COMPACT_ROOT="$ROOT/$S/compact/$COMPACT_LAYOUT_NAME"
TRACE_ROOT="$ROOT/$S/qc/stage12_trace/$HAP_READ_MODE"
TRACE_CONTIGS=$(expand_chroms "$CHROMS" | paste -sd, -)
log "12 compact zarr 생성: sample=$S root=$COMPACT_ROOT mode=$HAP_READ_MODE work=$WORK scratch=$([[ -n "$MFM_STAGE12_SCRATCH_ROOT" ]] && echo 1 || echo 0) trace_workers=$TRACE_JOBS trace_blosc_threads=$MFM_TRACE_BLOSC_THREADS trace_prefetch=$MFM_TRACE_PREFETCH_RECORDS sqlite_batch=$MFM_TRACE_SQLITE_BATCH sqlite_cache_mib_per_worker=$MFM_TRACE_SQLITE_CACHE_MIB"
run_cli "${args[@]}"

case "$STAGE12_TRACE_POLICY" in
  required)
    log "12 independent source trace 검증: sample=$S contigs=$TRACE_CONTIGS workers=$TRACE_JOBS"
    TRACE_SCRATCH_ARGS=()
    if [[ -n "$MFM_STAGE12_SCRATCH_ROOT" ]]; then
      TRACE_SCRATCH="$WORK/trace_sqlite"
      mkdir -p "$TRACE_SCRATCH"
      TRACE_SCRATCH_ARGS=(--scratch-dir "$TRACE_SCRATCH")
    fi
    run_cli validate-stage70-trace \
      --compact-root "$COMPACT_ROOT" \
      --contigs "$TRACE_CONTIGS" \
      --workers "$TRACE_JOBS" \
      --report-root "$TRACE_ROOT" \
      "${TRACE_SCRATCH_ARGS[@]}"
    ;;
  deferred)
    warn "12 source trace 검증 유예: cleanup은 trace 완료 전 차단됨 sample=$S"
    ;;
  off)
    warn "12 source trace 검증 비활성화: cleanup에는 ALLOW_CLEANUP_WITHOUT_TRACE=1이 별도로 필요함"
    ;;
esac
log "12 완료: $S"
