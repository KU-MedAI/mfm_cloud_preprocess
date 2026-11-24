#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

TMP_DIR=$(mktemp -d /tmp/mfm_clair_policy.XXXXXX)
trap 'rm -rf "$TMP_DIR"' EXIT

write_manifest() {
  local path=$1 platform=$2 read_type=$3 source=$4
  printf 'sample_id\tplatform\tread_type\tdataset_role\tproject_id\tcohort\talignment_path\tvcf_source\n' > "$path"
  printf 'S1\t%s\t%s\tpre_training\ttest\ttest\t/input.bam\t%s\n' "$platform" "$read_type" "$source" >> "$path"
}

run_source() {
  local manifest=$1
  shift
  env -i PATH="$PATH" PIPELINE_MANIFEST="$manifest" "$@" bash -c 'source ./config_v2.sh; sample_vcf_source S1'
}

write_manifest "$TMP_DIR/long_missing.tsv" hifi long_read ''
if run_source "$TMP_DIR/long_missing.tsv" >/dev/null 2>&1; then
  echo "FAIL: long-read missing vcf_source must fail" >&2; exit 1
fi

write_manifest "$TMP_DIR/long_public.tsv" hifi long_read public
[[ "$(run_source "$TMP_DIR/long_public.tsv")" == public ]] || { echo "FAIL: public VCF policy" >&2; exit 1; }

write_manifest "$TMP_DIR/long_clair.tsv" hifi long_read clair3
if run_source "$TMP_DIR/long_clair.tsv" >/dev/null 2>&1; then
  echo "FAIL: Clair3 without opt-in must fail" >&2; exit 1
fi
[[ "$(run_source "$TMP_DIR/long_clair.tsv" ENABLE_CLAIR3_SELF_CALL=1)" == clair3 ]] || {
  echo "FAIL: explicit Clair3 opt-in" >&2; exit 1;
}

write_manifest "$TMP_DIR/short.tsv" illumina short_read ''
[[ "$(run_source "$TMP_DIR/short.tsv")" == illumina_panel ]] || { echo "FAIL: short-read panel" >&2; exit 1; }

grep -F '[[ "$ENABLE_CLAIR3_SELF_CALL" == 1 ]]' 04_call_longread_clair3.sh >/dev/null
grep -F '[[ "$ENABLE_CLAIR3_SELF_CALL" == 1 ]]' 06_haplotag.sh >/dev/null
grep -F 'result.get("fail_count", 0)' src/integrate_wgs_rna/cli.py >/dev/null
grep -F 'result.get("missing_count", 0)' src/integrate_wgs_rna/cli.py >/dev/null
grep -F 'HAPLOTAG_REBUILD_PARTIAL=1' 06_haplotag.sh >/dev/null
grep -F 'partial_backup' 06_haplotag.sh >/dev/null
! grep -F '[[ "${PUBLIC_VCF_HAPLOTAG_OVERWRITE:-0}" == 1 ]]' 06_haplotag.sh >/dev/null
grep -F './00_environment.sh verify-core' run_pipeline.sh >/dev/null
grep -F 'if [[ "$ENABLE_CLAIR3_SELF_CALL" != 1 ]]' run_pipeline.sh >/dev/null
echo "PASS: Clair3 explicit opt-in and public-VCF fail-closed policy"
