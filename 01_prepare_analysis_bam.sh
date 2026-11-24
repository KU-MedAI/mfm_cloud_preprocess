#!/usr/bin/env bash
# 01 — raw long-read uBAM/FASTQ를 analysis reference 좌표의 BAM으로 표준화한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh
export PYTHONPATH="$PIPELINE_DIR/src:$PIPELINE_DIR${PYTHONPATH:+:$PYTHONPATH}"

S=${1:?사용법: $0 <sample>}
# RNA quickcheck는 시간이 긴 WGS 정렬 전에 Stage 01에서 수행한다.
RNA_INPUT=$(sample_rna "$S" 2>/dev/null || true)
if [[ -n "$RNA_INPUT" && "$RNA_INPUT" != "-" ]]; then
  export PYTHONPATH="$PIPELINE_DIR:${PYTHONPATH:-}"
  "$PYTHON" -m handoff.input_check --rna "$RNA_INPUT" --sample "$S" --samtools "$SAMTOOLS_BIN" \
    --output "$(sample_final_root "$S")/$S/qc/stage01_rna_quickcheck.json"
fi
require_manifest_sample_contract "$S"
if [[ -z "${PIPELINE_MANIFEST:-}" ]]; then
  log "$S: manifest 없는 legacy 실행 — 실행 단계 01(내부 Stage04) 변환 없이 기존 alignment 재사용"
  exit 0
fi

WORK=$(sample_work_root "$S")
OUT=$(sample_alignment "$S")
PREP_DIR="$WORK/prepared_input"
CONTRACT="$PREP_DIR/$S.input_contract.json"
[[ "$OUT" == *.bam ]] || die "$S 실행 단계 01(내부 Stage04) canonical output must be BAM: $OUT"
QC=$(sample_stage04_qc "$S")
mkdir -p "$PREP_DIR" "$WORK/logs"

log "$S: WGS 입력 계약 동결 → $CONTRACT (sort_mem= markdup_heap=GiB)"
"$RUNNER" validate-wgs-input-manifest \
  --manifest "$PIPELINE_MANIFEST" --sample "$S" --output-json "$CONTRACT" \
  > "$WORK/logs/04_contract.log"

export MFM_PREPARE_INPUT_CONTRACT="$CONTRACT"
STORAGE_CPUS=$("$PYTHON" -m integrate_wgs_rna.prepare_phases "$CONTRACT")
if [[ -n "${MFM_CPU_POOL_LEASE_FILE:-}" && "$STORAGE_CPUS" =~ ^[1-9][0-9]*$ ]]; then
  export MFM_SAMPLE_CPU_BUDGET="$STORAGE_CPUS" MFM_CPU_POOL_GRANTED="$STORAGE_CPUS"
  source ./config_v2.sh
fi

args=(
  --contract-json "$CONTRACT"
  --output "$OUT"
  --qc-json "$QC"
  --minimap2 "$MINIMAP2_BIN"
  --bwa-mem2 "$BWA_MEM2_BIN"
  --illumina-aligner "$STAGE04_ILLUMINA_ALIGNER"
  --bwa-mem2-index-root "$BWA_MEM2_INDEX_ROOT"
  --picard "$PICARD_BIN"
  --illumina-duplicate-policy "$STAGE04_ILLUMINA_DUPLICATE_POLICY"
  --markdup-heap-gb "$STAGE04_MARKDUP_HEAP_GB"
  --samtools "$SAMTOOLS_BIN"
  --align-threads "$STAGE04_ALIGN_THREADS"
  --sort-threads "$STAGE04_SORT_THREADS"
  --sort-memory "$STAGE04_SORT_MEM"
  --extract-threads "$STAGE04_EXTRACT_THREADS"
  --min-quality-values "$STAGE04_MIN_QUAL_VALUES"
)
[[ "$STAGE04_OVERWRITE" == 1 ]] && args+=(--overwrite)
[[ "$STAGE04_KEEP_RUN_BAMS" == 1 ]] && args+=(--keep-run-bams)

log "$S: 실행 단계 01(내부 Stage04) 모든 입력 포맷을 canonical analysis BAM으로 표준화 → $OUT"
"$PYTHON" ./src/integrate_wgs_rna/prepare_alignment.py "${args[@]}" \
  2>&1 | tee "$WORK/logs/04_prepare_alignment.log"
[[ -s "$QC" ]] || die "$S 실행 단계 01(내부 Stage04) QC JSON 없음: $QC"
[[ -s "$OUT" ]] || die "$S 실행 단계 01(내부 Stage04) analysis alignment 없음: $OUT"
