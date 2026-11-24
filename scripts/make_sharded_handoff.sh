#!/usr/bin/env bash
# 외부 서버 전달용 pipeline/shared/catalog/shards 단일 폴더를 원자적으로 조립한다.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DEST=${1:-}
SOURCE_MFM_REPO=${SOURCE_MFM_REPO:-$PIPELINE_ROOT/vendor/MFM-DNARNAFM}
PYTHON_BIN=${PYTHON_BIN:-$(command -v python3 || command -v python)}

die() { echo "[FATAL] $*" >&2; exit 1; }
[[ -n "$DEST" ]] || die "사용법: SOURCE_MFM_REPO=/path/MFM-DNARNAFM $0 /path/mfm_preprocess_share"
DEST=$(readlink -m "$DEST")
[[ ! -e "$DEST" ]] || die "대상 경로가 이미 존재합니다(덮어쓰지 않음): $DEST"
[[ -x "$PYTHON_BIN" ]] || die "jj_mfm Python 없음: $PYTHON_BIN"
[[ -d "$SOURCE_MFM_REPO/pre_training/preprocess_raw_readset" ]] || die "MFM source 없음: $SOURCE_MFM_REPO"
TMP="$DEST.part.$$"
[[ ! -e "$TMP" ]] || die "임시 경로가 이미 존재함: $TMP"
mkdir -p "$TMP"

SOURCE_MFM_REPO="$SOURCE_MFM_REPO" "$PIPELINE_ROOT/scripts/make_handoff_bundle.sh" "$TMP/pipeline"
"$PYTHON_BIN" "$TMP/pipeline/scripts/build_stratified_shards.py" --out-dir "$TMP/.shard_build"
mv "$TMP/.shard_build/catalog" "$TMP/catalog"
mv "$TMP/.shard_build/shards" "$TMP/shards"
mv "$TMP/.shard_build/shared" "$TMP/shared"
rmdir "$TMP/.shard_build"
mkdir -p "$TMP/shared/manifests" "$TMP/data" "$TMP/output"
cp -a "$PIPELINE_ROOT/share_package/manifests/." "$TMP/shared/manifests/"
cp "$PIPELINE_ROOT/SHARDED_HANDOFF.md" "$TMP/README.md"
cp "$PIPELINE_ROOT/scripts/share_run_shard.sh" "$TMP/run_shard.sh"
cp "$PIPELINE_ROOT/scripts/share_run_shard_queue.sh" "$TMP/run_shard_queue.sh"
chmod +x "$TMP/run_shard.sh" "$TMP/run_shard_queue.sh"

(
  cd "$TMP"
  find README.md run_shard.sh run_shard_queue.sh pipeline catalog shared shards -type f ! -name HANDOFF_SHA256SUMS -print0 |
    sort -z | xargs -0 sha256sum
) > "$TMP/HANDOFF_SHA256SUMS"
cat > "$TMP/HANDOFF_INFO.txt" <<EOF
schema=mfm_stratified_handoff_v1
created_at=$(date -Is)
samples=3212
shards=$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["shards"])' "$TMP/catalog/shard_summary.json")
target_output_tib=8.5
hard_cap_tib=10.0
EOF
mv "$TMP" "$DEST"
echo "[OK] sharded handoff: $DEST"
echo "[OK] 실행: $DEST/run_shard.sh shard_001 --data-root /data/mfm/input --output-root /data/mfm/output"
