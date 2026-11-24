#!/usr/bin/env bash
set -euo pipefail

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR=${PIPELINE_DIR:-$(cd "$TEST_DIR/.." && pwd)}
source "$PIPELINE_DIR/config_v2.sh"

fail() { echo "[FAIL] $*" >&2; exit 1; }
assert_eq() { [[ "$1" == "$2" ]] || fail "expected=[$2] actual=[$1]"; }

assert_eq "$(chroms_space 'chr1-chr3,5 chr22')" "chr1 chr2 chr3 chr5 chr22"
assert_eq "$(first_chrom 'chr22')" "chr22"

TMP=$(mktemp -d "$TEST_DIR/.tmp_rna_entrypoints.XXXXXX")
trap 'rm -rf "$TMP"' EXIT
for h in 1 2; do
  mkdir -p "$TMP/assets/HAP$h"
  for c in chr5 chr6; do
    printf '>%s\nA\n' "$c" > "$TMP/assets/HAP$h/$c.fa"
    printf '%s\t1\t6\t1\t2\n' "$c" > "$TMP/assets/HAP$h/$c.fa.fai"
    printf 'chain\n' > "$TMP/assets/HAP$h/$c.chain"
  done
done
require_hap_contig_assets "$TMP/assets" "chr5-chr6" 1
if (require_hap_contig_assets "$TMP/assets" "chr5-chr7" 1 >/dev/null 2>&1); then
  fail "missing chr7 assets were accepted"
fi

SLURM_DIR="$TMP/slurm"
mkdir -p "$SLURM_DIR"
cp "$PIPELINE_DIR/slurm_rna.sbatch" "$SLURM_DIR/slurm_rna.sbatch"
printf 'CONFIG_LOADED=1\nexport CONFIG_LOADED\n' > "$SLURM_DIR/config_v2.sh"
printf '%s\n' '#!/usr/bin/env bash' 'set -euo pipefail' '[[ "${CONFIG_LOADED:-0}" == 1 ]]' 'printf "%s\\n" "$1" > "${TEST_RNA_SELECTED:?}"' > "$SLURM_DIR/run_rna_sample.sh"
chmod +x "$SLURM_DIR/run_rna_sample.sh"
printf '# comment\n\nSAMPLE_A\n  # indented comment\nSAMPLE_B extra\n' > "$TMP/samples.txt"
TEST_RNA_SELECTED="$TMP/selected.txt" \
SLURM_ARRAY_TASK_ID=2 SLURM_CPUS_PER_TASK=2 SLURM_SUBMIT_DIR="$TMP" \
SAMPLE_LIST=samples.txt bash "$SLURM_DIR/slurm_rna.sbatch"
assert_eq "$(cat "$TMP/selected.txt")" "SAMPLE_B"

mkdir -p "$TMP/bin"
printf '%s\n' '#!/usr/bin/env bash' 'printf "%s\\n" "$@" > "${SBATCH_CAPTURE:?}"' > "$TMP/bin/sbatch"
chmod +x "$TMP/bin/sbatch"
SBATCH_CAPTURE="$TMP/sbatch.args" PATH="$TMP/bin:$PATH" \
  bash "$PIPELINE_DIR/run_pipeline.sh" slurm-rna "$TMP/samples.txt"
grep -Fx -- '--array=1-2%4' "$TMP/sbatch.args" >/dev/null || fail "slurm-rna array가 유효 sample 2개 기준이 아님"
grep -F -- "SAMPLE_LIST=$(readlink -f "$TMP/samples.txt")" "$TMP/sbatch.args" >/dev/null || fail "slurm-rna SAMPLE_LIST가 절대경로가 아님"

if grep -nE 'HAP1/chr1[.]fa|HD/chr1[.]fa' "$PIPELINE_DIR/run_rna_sample.sh"; then
  fail "RNA entrypoint에 chr1.fa 하드코딩이 남아 있음"
fi
grep -F -- '--contigs "$RNA_CONTIGS"' "$PIPELINE_DIR/run_rna_sample.sh" >/dev/null \
  || fail "dynamic contig argument missing in RNA runner"
grep -F -- 'CACHE_TAG="gtex_${GENCODE_TAG}_${ALIGN_FP:0:8}"' "$PIPELINE_DIR/run_rna_sample.sh" >/dev/null \
  || fail "STAR cache must use alignment-only fingerprint"

for required in '--reference "$(sample_hap_target "$S" 1)"' '--reference "$(sample_hap_target "$S" 2)"'; do
  grep -F -- "$required" "$PIPELINE_DIR/run_rna_sample.sh" >/dev/null \
    || fail "RNA CNV personal-reference contract missing: $required"
done

for script in 13_rna_to_fastq.py 14_lift_annotation.py 15_star_index.py 16_star_align.py \
              17_mark_duplicates.py 18_coverage_zarr.py 19_rna_realign_qc.py \
              20_hap_variants.py 21_allelic_counts.py 22_cnv_filter.py 23_ase_qc.py \
              24_stamp_rna_metadata.py 25_validate_sample_package.sh 26_cleanup_intermediates.sh; do
  [[ -x "$PIPELINE_DIR/$script" ]] || fail "RNA stage implementation missing/not executable: $script"
done
for required in \
  'pipeline_run_command 13 13_rna_to_fastq.py' \
  '14_lift_annotation.py:HAP$h' \
  '15_star_index.py:HAP$h' \
  '16_star_align.py:HAP$h' \
  '17_mark_duplicates.py:HAP$h' \
  '18_coverage_zarr.py:HAP' \
  '--window-out "$TD/coverage_windows.zarr"' \
  'pipeline_run_command 20 20_hap_variants.py' \
  'pipeline_run_command 21 21_allelic_counts.py' \
  'pipeline_run_command 22 22_cnv_filter.py' \
  'pipeline_run_command 23 23_ase_qc.py' \
  'pipeline_run_command 24 24_stamp_rna_metadata.py' \
  '25_validate_sample_package.sh' \
  '26_cleanup_intermediates.sh' \
  'RNA_START_AT=${RNA_START_AT:-13}' 'RNA_STOP_AFTER'; do
  grep -F -- "$required" "$PIPELINE_DIR/run_rna_sample.sh" >/dev/null \
    || fail "RNA numbered orchestration contract missing: $required"
done
grep -F -- '${WGS_DEFER_FINALIZE:-${WGS_DEFER_CLEANUP:-0}}' "$PIPELINE_DIR/run_wgs_sample.sh" >/dev/null \
  || fail "WGS deferred cleanup contract missing"
for lib in lib/state.sh lib/stage_runner.sh; do
  [[ -s "$PIPELINE_DIR/$lib" ]] || fail "common runner library missing: $lib"
done

for script in run_rna_sample.sh 20_cleanup_intermediates.sh; do
  grep -F -- 'sample_rna_evidence_label "$S"' "$PIPELINE_DIR/$script" >/dev/null \
    || fail "RNA evidence label helper missing: $script"
done
(
  is_long_read_sample() { return 0; }
  sample_vcf_source() { printf '%s\n' public; }
  PIPELINE_MANIFEST="$TMP/manifest.tsv"
  assert_eq "$(sample_rna_evidence_label S1)" "public"
)
(
  is_long_read_sample() { return 1; }
  assert_eq "$(sample_rna_evidence_label S2)" "panel"
)
cnv_parse=$($PYTHON "$PIPELINE_DIR/22_cnv_filter.py" \
  --exons "$TMP/missing1.tsv" --bam "$TMP/missing1.bam" --reference "$TMP/ref1.fa" \
  --exons "$TMP/missing2.tsv" --bam "$TMP/missing2.bam" --reference "$TMP/ref2.fa" \
  --out "$TMP/gene_rd.tsv" 2>&1 || true)
[[ "$cnv_parse" != *"AttributeError"* ]] || fail "Stage21 repeated --reference argparse contract"
[[ "$cnv_parse" == *"exon 파일 없음"* ]] || fail "Stage22 repeated input parsing did not reach file validation"
echo "[PASS] RNA entrypoint/slurm/partial-contig contract"
