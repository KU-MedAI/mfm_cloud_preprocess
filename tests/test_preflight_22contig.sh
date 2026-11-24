#!/usr/bin/env bash
# Stage 05의 22-contig M5/LN, VCF REF, runner sample 격리를 실제 작은 fixture로 검사한다.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SAMTOOLS_BIN=${SAMTOOLS_BIN:-$(command -v samtools)}
BCFTOOLS_BIN=${BCFTOOLS_BIN:-$(command -v bcftools)}
PYTHON_BIN=${PYTHON_BIN:-$(command -v python)}
TMP_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/mfm_preflight_test.XXXXXX")
if [[ "${KEEP_TEST_TMP:-0}" == 1 ]]; then
  echo "[INFO] preserving fixture: $TMP_ROOT" >&2
else
  trap 'rm -rf "$TMP_ROOT"' EXIT
fi

REF="$TMP_ROOT/ref.fa"
SEQ=$(printf 'A%.0s' {1..100})
for chrom in $(seq 1 22); do
  printf '>chr%s\n%s\n' "$chrom" "$SEQ" >> "$REF"
done
"$SAMTOOLS_BIN" faidx "$REF"

make_sam() {
  local out=$1 bad_m5=${2:-0} chrom m5
  {
    printf '@HD\tVN:1.6\tSO:coordinate\n'
    for chrom in $(seq 1 22); do
      m5=$("$SAMTOOLS_BIN" faidx "$REF" "chr$chrom" | awk 'NR>1 {printf "%s", $0}' | md5sum | awk '{print $1}')
      if [[ "$bad_m5" == 1 && "$chrom" == 22 ]]; then m5=00000000000000000000000000000000; fi
      printf '@SQ\tSN:chr%s\tLN:100\tM5:%s\n' "$chrom" "$m5"
    done
    printf '@RG\tID:rg1\tSM:S1\n'
    printf 'read1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rg1\n'
  } > "$out"
}

make_sam_without_m5() {
  local out=$1 sample=$2 chrom
  {
    printf '@HD\tVN:1.6\tSO:coordinate\n'
    for chrom in $(seq 1 22); do
      printf '@SQ\tSN:chr%s\tLN:100\n' "$chrom"
    done
    printf '@RG\tID:rg1\tSM:%s\n' "$sample"
    printf 'read1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rg1\n'
  } > "$out"
}

make_sam "$TMP_ROOT/good.sam"
make_sam "$TMP_ROOT/bad_m5.sam" 1
"$SAMTOOLS_BIN" view -b -o "$TMP_ROOT/good.bam" "$TMP_ROOT/good.sam"
"$SAMTOOLS_BIN" view -b -o "$TMP_ROOT/bad_m5.bam" "$TMP_ROOT/bad_m5.sam"
make_sam_without_m5 "$TMP_ROOT/no_m5.sam" SRC1
"$SAMTOOLS_BIN" view -b -o "$TMP_ROOT/no_m5.bam" "$TMP_ROOT/no_m5.sam"
"$SAMTOOLS_BIN" index "$TMP_ROOT/good.bam"
"$SAMTOOLS_BIN" index "$TMP_ROOT/bad_m5.bam"
"$SAMTOOLS_BIN" index "$TMP_ROOT/no_m5.bam"

make_vcf() {
  local out=$1 ref=$2 chrom
  {
    printf '##fileformat=VCFv4.2\n'
    for chrom in $(seq 1 22); do printf '##contig=<ID=chr%s,length=100>\n' "$chrom"; done
    printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
    printf 'chr1\t2\t.\t%s\tT\t.\tPASS\t.\tGT\t0|1\n' "$ref"
  } > "$out"
  "$BCFTOOLS_BIN" view -Oz -o "$out.gz" "$out"
  "$BCFTOOLS_BIN" index -t "$out.gz"
}

make_vcf "$TMP_ROOT/good.vcf" A
make_vcf "$TMP_ROOT/bad_ref.vcf" C

run_preflight() {
  local alignment=$1 vcf=$2 run_id=$3 source_sample_id=${4:-}
  DATASET_ROOT="$TMP_ROOT/data" RUN_ID="$run_id" \
  DEFAULT_PLATFORM=hifi DEFAULT_READ_TYPE=long_read DEFAULT_DATASET_ROLE=pre_training \
  DEFAULT_PROJECT_ID=fixture DEFAULT_COHORT=fixture ALIGNMENT_PATH="$alignment" \
  PREPARED_ALIGNMENT_PATH="$alignment" \
  VCF_SOURCE_POLICY=manifest LONGREAD_SEQUENCE_VCF_SOURCE=public LONGREAD_HAPLOTAG_VCF_SOURCE=public \
  PUBLIC_VCF_PATH="$vcf" PUBLIC_VCF_SAMPLE_ID=S1 VCF_REFERENCE_FASTA="$REF" \
  VCF_CALLSET=fixture PHASE_SOURCE=provided REF_DECODE="$REF" REF_SEQUENCE="$REF" \
  SOURCE_SAMPLE_ID="$source_sample_id" \
  SAMTOOLS_BIN="$SAMTOOLS_BIN" BCFTOOLS_BIN="$BCFTOOLS_BIN" PREFLIGHT_MIN_FREE_GIB=0 \
  "$PIPELINE_DIR/02_preflight.sh" S1
}

run_preflight "$TMP_ROOT/good.bam" "$TMP_ROOT/good.vcf.gz" good
GOOD_QC="$TMP_ROOT/data/_run_status/good/wgs/fixture/S1/preflight/S1.preflight.json"
"$PYTHON_BIN" - "$GOOD_QC" <<'PY'
import json
import sys

qc = json.load(open(sys.argv[1]))
checks = {row["check"]: row for row in qc["checks"]}
assert qc["status"] == "pass"
assert sum(name.startswith("cram_decode_m5_chr") for name in checks) == 22
assert all(checks[f"cram_decode_m5_chr{i}"]["status"] == "PASS" for i in range(1, 23))
assert checks["public_vcf_ref_alleles"]["status"] == "PASS"
PY

if run_preflight "$TMP_ROOT/bad_m5.bam" "$TMP_ROOT/good.vcf.gz" bad_m5; then
  echo "expected chr22 M5 mismatch to fail" >&2
  exit 1
fi
BAD_M5_QC="$TMP_ROOT/data/_run_status/bad_m5/wgs/fixture/S1/preflight/S1.preflight.json"
"$PYTHON_BIN" - "$BAD_M5_QC" <<'PY'
import json
import sys
checks = {row["check"]: row for row in json.load(open(sys.argv[1]))["checks"]}
assert checks["cram_decode_m5_chr22"]["status"] == "FAIL"
PY

if run_preflight "$TMP_ROOT/good.bam" "$TMP_ROOT/bad_ref.vcf.gz" bad_ref; then
  echo "expected VCF REF mismatch to fail" >&2
  exit 1
fi
BAD_REF_QC="$TMP_ROOT/data/_run_status/bad_ref/wgs/fixture/S1/preflight/S1.preflight.json"
"$PYTHON_BIN" - "$BAD_REF_QC" <<'PY'
import json
import sys
checks = json.load(open(sys.argv[1]))["checks"]
assert any(row["check"] == "public_vcf_ref_alleles" and row["status"] == "FAIL" for row in checks)
PY

# canonical BAM이 M5를 보존하지 않아도 Stage 01 QC의 source reference와 22개 autosome
# 염기서열이 일치하고 @RG alias가 허용되면 통과해야 한다.
STAGE04_DIR="$TMP_ROOT/data/_work/no_m5_alias/wgs/fixture/S1/prepared_input"
mkdir -p "$STAGE04_DIR"
"$PYTHON_BIN" - "$STAGE04_DIR/S1.stage04.json" "$TMP_ROOT/no_m5.bam" "$REF" <<'PY'
import json
import sys

out, alignment, reference = sys.argv[1:]
json.dump({
    "status": "converted_aligned",
    "canonical_format": "bam",
    "analysis_alignment": alignment,
    "normalization_reference": {"path": reference},
    "duplicate_marking": {},
}, open(out, "w"))
PY
run_preflight "$TMP_ROOT/no_m5.bam" "$TMP_ROOT/good.vcf.gz" no_m5_alias SRC1
NO_M5_QC="$TMP_ROOT/data/_run_status/no_m5_alias/wgs/fixture/S1/preflight/S1.preflight.json"
"$PYTHON_BIN" - "$NO_M5_QC" <<'PY'
import json
import sys

checks = {row["check"]: row for row in json.load(open(sys.argv[1]))["checks"]}
assert checks["stage04_reference_provenance"]["status"] == "PASS"
assert checks["bam_source_reference_autosomes"]["status"] == "PASS"
assert checks["cram_sample_id"]["status"] == "PASS"
assert all(checks[f"cram_decode_m5_chr{i}"]["status"] == "PASS" for i in range(1, 23))
PY

# raw uBAM/FASTQ를 analysis reference에 정렬한 BAM은 원본 source decode와
# 서열이 다르더라도 정상이다. 실제 정렬 reference의 M5 계약은 유지해야 한다.
RAW_DIR="$TMP_ROOT/data/_work/raw_reference/wgs/fixture/S1/prepared_input"
mkdir -p "$RAW_DIR"
ALT_REF="$TMP_ROOT/source_other.fa"
for chrom in $(seq 1 22); do
  printf '>chr%s\n%s\n' "$chrom" "$(printf 'C%.0s' {1..100})" >> "$ALT_REF"
done
"$SAMTOOLS_BIN" faidx "$ALT_REF"
"$PYTHON_BIN" - "$RAW_DIR/S1.stage04.json" "$TMP_ROOT/no_m5.bam" "$REF" <<'PY'
import json
import sys
out, alignment, reference = sys.argv[1:]
json.dump({"status": "created_raw_longread", "canonical_format": "bam",
           "analysis_alignment": alignment, "normalization_reference": {"path": reference},
           "alignment_policy": {"mode": "raw_long_read"}}, open(out, "w"))
PY
SOURCE_DECODE_REFERENCE_FASTA="$ALT_REF" \
  run_preflight "$TMP_ROOT/no_m5.bam" "$TMP_ROOT/good.vcf.gz" raw_reference SRC1
RAW_QC="$TMP_ROOT/data/_run_status/raw_reference/wgs/fixture/S1/preflight/S1.preflight.json"
"$PYTHON_BIN" - "$RAW_QC" <<'PY'
import json
import sys
qc = json.load(open(sys.argv[1]))
assert qc["status"] == "pass"
assert all(row["status"] == "PASS" for row in qc["checks"] if row["check"].startswith("cram_decode_m5_chr"))
PY

TARGET=$(command -v env)
RUN_STATUS_DIR="$TMP_ROOT/runner" RUN_LOG_DIR="$TMP_ROOT/runner/logs" RUNNER_TARGET_SCRIPT="$TARGET" \
  "$PIPELINE_DIR/scripts/run_one_local.sh" false
RUN_STATUS_DIR="$TMP_ROOT/runner" RUN_LOG_DIR="$TMP_ROOT/runner/logs" RUNNER_TARGET_SCRIPT="$TARGET" \
  "$PIPELINE_DIR/scripts/run_one_local.sh" true
awk -F'\t' '$1=="false" && $2=="failed" {bad=1} $1=="true" && $2=="success" {good=1} END {exit !(bad && good)}' \
  "$TMP_ROOT/runner/results/false.result.tsv" "$TMP_ROOT/runner/results/true.result.tsv"

echo "[OK] Stage 05 22-contig/VCF 계약 및 sample 오류 격리 검증 통과"
