#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

fail() { echo "FAIL: $*" >&2; exit 1; }
# 단독 실행도 run_all.sh와 같은 jj_mfm interpreter를 사용한다.
if [[ -z "${PYTHON_BIN:-}" && -x /home/aibio_1/miniconda3/envs/jj_mfm/bin/python ]]; then
  PYTHON_BIN=/home/aibio_1/miniconda3/envs/jj_mfm/bin/python
fi
PYTHON_BIN=${PYTHON_BIN:-$(command -v python3 || command -v python)}
[[ -x "$PYTHON_BIN" ]] || fail "Python interpreter not found"

env -i PATH="$PATH" MFM_MEMORY_AVAILABLE_GIB_OVERRIDE=2048 bash -c '
  source ./config_v2.sh
  [[ "$MFM_MEMORY_PROFILE_RESOLVED" == balanced ]]
  [[ "$MFM_RAM_LIMIT_BYTES" == 1650000000000 ]]
  [[ "$MFM_COMPACT_RESERVATION_MODEL" == input_v2 ]]
  [[ "$WGS_SORT_MEM" == 2G ]]
' || fail "기본 1.65TB 예산과 균형 메모리 버퍼"

env -i PATH="$PATH" MFM_COMPACT_RESERVATION_MODEL=legacy \
  MFM_VALIDATION_CACHE_CHUNKS=0 MFM_TRACE_DECODED_CACHE_MIB=0 bash -c '
  source ./config_v2.sh
  [[ "$MFM_COMPACT_RESERVATION_MODEL" == legacy ]]
  [[ "$MFM_VALIDATION_CACHE_CHUNKS" == 0 && "$MFM_TRACE_DECODED_CACHE_MIB" == 0 ]]
' || fail "기존 예약 모델 및 캐시 비활성 명시값 보존"

env -i PATH="$PATH" MFM_MEMORY_PROFILE=balanced MFM_MEMORY_AVAILABLE_GIB_OVERRIDE=128 bash -c '
  source ./config_v2.sh
  [[ "$MFM_MEMORY_PROFILE_RESOLVED" == balanced ]]
  [[ "$WGS_SORT_MEM" == 2G && "$STAGE04_SORT_MEM" == 2G ]]
  [[ "$MFM_TRACE_PREFETCH_RECORDS" == 65536 ]]
  [[ "$MFM_TRACE_SQLITE_BATCH" == 8192 ]]
  [[ "$MFM_VALIDATION_RECORD_CHUNK" == 262144 ]]
  [[ "$MFM_VALIDATION_CACHE_CHUNKS" == 4 ]]
' || fail "balanced profile"

env -i PATH="$PATH" MFM_MEMORY_PROFILE=high MFM_MEMORY_AVAILABLE_GIB_OVERRIDE=2048 bash -c '
  source ./config_v2.sh
  [[ "$MFM_MEMORY_PROFILE_RESOLVED" == high ]]
  [[ "$WGS_SORT_MEM" == 16G && "$STAGE04_SORT_MEM" == 16G ]]
  [[ "$MFM_TRACE_PREFETCH_RECORDS" == 1048576 ]]
  [[ "$MFM_TRACE_SQLITE_BATCH" == 262144 ]]
  [[ "$MFM_VALIDATION_RECORD_CHUNK" == 262144 ]]
  [[ "$MFM_VALIDATION_PACKED_BYTE_CHUNK" == 16777216 ]]
  [[ "$MFM_VALIDATION_CACHE_CHUNKS" == 4 ]]
  [[ "$STAGE04_MARKDUP_HEAP_GB" == 32 && "$RNA_MARKDUP_HEAP_GB" == 16 ]]
' || fail "high profile"

env -i PATH="$PATH" MFM_MEMORY_PROFILE=auto MFM_MEMORY_AVAILABLE_GIB_OVERRIDE=2048 bash -c '
  source ./config_v2.sh
  [[ "$MFM_MEMORY_PROFILE_RESOLVED" == high ]]
' || fail "auto high"

env -i PATH="$PATH" MFM_MEMORY_PROFILE=auto MFM_MEMORY_AVAILABLE_GIB_OVERRIDE=64 bash -c '
  source ./config_v2.sh
  [[ "$MFM_MEMORY_PROFILE_RESOLVED" == balanced ]]
' || fail "auto balanced"

PYTHONPATH=src MFM_TRACE_SQLITE_BATCH=17 MFM_TRACE_PREFETCH_RECORDS=33 \
  MFM_VALIDATION_RECORD_CHUNK=65 MFM_VALIDATION_PACKED_BYTE_CHUNK=129 \
  MFM_VALIDATION_CACHE_CHUNKS=3 \
  "$PYTHON_BIN" -c '
from integrate_wgs_rna import stage70_trace_validator as trace
from integrate_wgs_rna.compact import runtime
assert trace.SQLITE_MULTISET_BATCH == 17
assert trace.ZARR_PREFETCH_RECORDS == 33
assert runtime.VALIDATION_RECORD_CHUNK == 65
assert runtime.VALIDATION_PACKED_BYTE_CHUNK == 129
assert runtime.VALIDATION_CACHE_CHUNKS == 3
' || fail "Python runtime buffer environment"

if env PYTHONPATH=src MFM_TRACE_SQLITE_BATCH=0 "$PYTHON_BIN" -c \
  'import integrate_wgs_rna.stage70_trace_validator' >/dev/null 2>&1; then
  fail "invalid runtime buffer accepted"
fi

echo "PASS: 메모리 프로필과 처리 버퍼"
