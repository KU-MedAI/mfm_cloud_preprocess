#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
export PATH=/home/aibio_1/miniconda3/envs/jj_mfm/bin:$PATH
PYTHON_BIN=/home/aibio_1/miniconda3/envs/jj_mfm/bin/python
fail() { echo "FAIL: $*" >&2; exit 1; }

env -i PATH=/usr/bin:/bin MFM_CORE_ENV=/home/aibio_1/miniconda3/envs/jj_mfm bash -c '
  source ./config_v2.sh
  [[ "$PYTHON_BIN" == /home/aibio_1/miniconda3/envs/jj_mfm/bin/python ]]
  [[ "$(command -v samtools)" == /home/aibio_1/miniconda3/envs/jj_mfm/bin/samtools ]]
  [[ -d "$MFM_REPO_ROOT/pre_training/preprocess_raw_readset" ]]
' || fail "jj_mfm 환경 고정"

rg -F 'wgs-phased-list)' run_pipeline.sh >/dev/null || fail "phased entrypoint"
rg -F 'WGS_START_AT=${WGS_START_AT:-01}' run_wgs_sample.sh >/dev/null || fail "stage start selector"
rg -F 'MFM_STAGE11_SCRATCH_ROOT' src/integrate_wgs_rna/personal_realign.py >/dev/null || fail "stage11 scratch"
rg -F 'MFM_STAGE12_SCRATCH_ROOT' 12_compact_zarr.sh >/dev/null || fail "stage12 scratch"
bash -n run_wgs_phased.sh run_shard_queue.sh

TMP=$(mktemp -d /tmp/mfm_parallel_download.XXXXXX)
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/data/source" "$TMP/out/input"
printf 'alpha\n' > "$TMP/data/source/a.txt"
printf 'beta\n' > "$TMP/data/source/b.txt"
sha_a=$(sha256sum "$TMP/data/source/a.txt" | awk '{print $1}')
sha_b=$(sha256sum "$TMP/data/source/b.txt" | awk '{print $1}')
size_a=$(stat -c %s "$TMP/data/source/a.txt")
size_b=$(stat -c %s "$TMP/data/source/b.txt")
{
  printf 'source_url\trelative_path\tsize_bytes\tchecksum_type\tchecksum_value\n'
  printf 'file://%s\tinput/a.txt\t%s\tsha256\t%s\n' "$TMP/data/source/a.txt" "$size_a" "$sha_a"
  printf 'file://%s\tinput/b.txt\t%s\tsha256\t%s\n' "$TMP/data/source/b.txt" "$size_b" "$sha_b"
} > "$TMP/plan.tsv"
cp "$TMP/data/source/a.txt" "$TMP/out/input/a.txt"
cp "$TMP/data/source/b.txt" "$TMP/out/input/b.txt"
"$PYTHON_BIN" scripts/download_shard.py --plan "$TMP/plan.tsv" --data-root "$TMP/out" --verify-only --jobs 2 \
  | rg -F '[SUMMARY] pass=2 jobs=2' >/dev/null || fail "parallel checksum verification"
printf 'corrupt\n' > "$TMP/out/input/b.txt"
if "$PYTHON_BIN" scripts/download_shard.py --plan "$TMP/plan.tsv" --data-root "$TMP/out" --verify-only --jobs 2 >/dev/null; then
  fail "parallel verifier missed corruption"
fi

echo "PASS: 환경·phased scheduler·scratch·parallel download 계약"
