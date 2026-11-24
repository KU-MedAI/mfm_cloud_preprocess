#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
ROOT=$(pwd -P)
TMP_ROOT=$(mktemp -d /tmp/mfm_portable_paths.XXXXXX)
export PIPELINE_MANIFEST="$ROOT/manifests/preprocess_manifest.tsv"
export INPUT_ROOT="$TMP_ROOT/input"
export OUTPUT_ROOT="$TMP_ROOT/output"
export OUTPUT_ROOT_OVERRIDE=1
export RUN_ID=portable_test
source ./config_v2.sh

column_number() {
  awk -F'\t' -v name="$1" 'NR==1{for(i=1;i<=NF;i++)if($i==name){print i;exit}}' "$PIPELINE_MANIFEST"
}
platform_col=$(column_number platform)
group_col=$(column_number data_group)
read_manifest_col=$(column_number read_files_manifest)
short_pair=$(awk -F'\t' -v p="$platform_col" -v g="$group_col" 'NR>1 && $p=="illumina" && $g=="wgs_rna_pair"{print $1;exit}' "$PIPELINE_MANIFEST")
long_public=$(awk -F'\t' -v p="$platform_col" 'NR>1 && $p!="illumina"{print $1;exit}' "$PIPELINE_MANIFEST")
multi_read=$(awk -F'\t' -v m="$read_manifest_col" 'NR>1 && $m!=""{print $1;exit}' "$PIPELINE_MANIFEST")

assert_prefix() {
  case "$1" in "$2"|"$2"/*) ;; *) echo "prefix mismatch: value=$1 prefix=$2" >&2; exit 1;; esac
}

assert_prefix "$(sample_source_alignment "$short_pair")" "$INPUT_ROOT/wgs"
assert_prefix "$(sample_rna "$short_pair")" "$INPUT_ROOT/rna"
assert_prefix "$(sample_common_vcf_dir "$short_pair")" "$INPUT_ROOT/vcf/1KGP_ILLUMINA"
assert_prefix "$(sample_final_root "$short_pair")" "$OUTPUT_ROOT/wgs_rna_pair"
assert_prefix "$(sample_work_root "$short_pair")" "$OUTPUT_ROOT/_work/portable_test/wgs_rna_pair"
assert_prefix "$(sample_status_root "$short_pair")" "$OUTPUT_ROOT/_run_status/portable_test/wgs_rna_pair"
assert_prefix "$(sample_public_vcf "$long_public")" "$INPUT_ROOT/public_vcf"

expected_read_manifest="$ROOT/manifests/read_files.tsv"
[[ "$(sample_read_files_manifest "$multi_read")" == "$expected_read_manifest" ]] || {
  echo "read_files manifest mismatch" >&2; exit 1;
}

rg -q -- "--exclude='./output'" runner_local.sh
rg -q "wgs_rna_pair/<project>" output/README.md
echo "portable path contract PASS"
