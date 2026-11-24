#!/usr/bin/env bash
# 처리 결과의 schema/chunk layout은 건드리지 않고 실행 버퍼만 서버 RAM에 맞춘다.

mfm_memory_profile_init() {
  if [[ -z "${MFM_MEMORY_PROFILE_REQUESTED+x}" ]]; then
    MFM_MEMORY_PROFILE_REQUESTED=${MFM_MEMORY_PROFILE:-balanced}
  fi
  local requested available reserve_default
  requested=$(printf '%s' "$MFM_MEMORY_PROFILE_REQUESTED" | tr '[:upper:]' '[:lower:]' | tr -d ' ')
  case "$requested" in auto|balanced|high) ;; *)
    echo "[FATAL] MFM_MEMORY_PROFILE must be auto|balanced|high: $requested" >&2; return 1 ;;
  esac
  if [[ -n "${MFM_MEMORY_AVAILABLE_GIB_OVERRIDE:-}" ]]; then
    available=$MFM_MEMORY_AVAILABLE_GIB_OVERRIDE
  else
    available=$(awk '/MemAvailable/{print int($2/1024/1024)}' /proc/meminfo)
  fi
  [[ "$available" =~ ^[1-9][0-9]*$ ]] || {
    echo "[FATAL] invalid MemAvailable GiB: $available" >&2; return 1;
  }
  MFM_MEMORY_AVAILABLE_GIB=$available
  if [[ "$requested" == auto ]]; then
    if (( available >= 512 )); then MFM_MEMORY_PROFILE_RESOLVED=high; else MFM_MEMORY_PROFILE_RESOLVED=balanced; fi
  else
    MFM_MEMORY_PROFILE_RESOLVED=$requested
  fi
  MFM_MEMORY_PROFILE=$MFM_MEMORY_PROFILE_RESOLVED
  if [[ "$MFM_MEMORY_PROFILE_RESOLVED" == high ]]; then
    MFM_SORT_MEM_DEFAULT=16G
    MFM_STAGE04_MARKDUP_HEAP_GB_DEFAULT=32
    MFM_RNA_MARKDUP_HEAP_GB_DEFAULT=16
    MFM_TRACE_SQLITE_BATCH=${MFM_TRACE_SQLITE_BATCH:-262144}
    MFM_TRACE_PREFETCH_RECORDS=${MFM_TRACE_PREFETCH_RECORDS:-1048576}
    MFM_VALIDATION_RECORD_CHUNK=${MFM_VALIDATION_RECORD_CHUNK:-262144}
    MFM_VALIDATION_PACKED_BYTE_CHUNK=${MFM_VALIDATION_PACKED_BYTE_CHUNK:-16777216}
    MFM_VALIDATION_CACHE_CHUNKS=${MFM_VALIDATION_CACHE_CHUNKS:-4}
    MFM_WGS_SAMPLE_MEMORY_GIB_DEFAULT=256
    MFM_WGS_HAP_MEMORY_GIB_DEFAULT=256
    MFM_RNA_SAMPLE_MEMORY_GIB_DEFAULT=64
    reserve_default=$((available / 5)); (( reserve_default >= 256 )) || reserve_default=256
  else
    MFM_SORT_MEM_DEFAULT=2G
    MFM_STAGE04_MARKDUP_HEAP_GB_DEFAULT=16
    MFM_RNA_MARKDUP_HEAP_GB_DEFAULT=8
    MFM_TRACE_SQLITE_BATCH=${MFM_TRACE_SQLITE_BATCH:-8192}
    MFM_TRACE_PREFETCH_RECORDS=${MFM_TRACE_PREFETCH_RECORDS:-65536}
    MFM_VALIDATION_RECORD_CHUNK=${MFM_VALIDATION_RECORD_CHUNK:-262144}
    MFM_VALIDATION_PACKED_BYTE_CHUNK=${MFM_VALIDATION_PACKED_BYTE_CHUNK:-16777216}
    MFM_VALIDATION_CACHE_CHUNKS=${MFM_VALIDATION_CACHE_CHUNKS:-4}
    MFM_WGS_SAMPLE_MEMORY_GIB_DEFAULT=64
    MFM_WGS_HAP_MEMORY_GIB_DEFAULT=96
    MFM_RNA_SAMPLE_MEMORY_GIB_DEFAULT=48
    reserve_default=$((available / 5)); (( reserve_default >= 16 )) || reserve_default=16
  fi
  (( reserve_default < available )) || reserve_default=$((available / 4))
  (( reserve_default >= 1 )) || reserve_default=1
  MFM_MEMORY_RESERVE_GIB=${MFM_MEMORY_RESERVE_GIB:-$reserve_default}
  local name value
  for name in MFM_TRACE_SQLITE_BATCH MFM_TRACE_PREFETCH_RECORDS MFM_VALIDATION_RECORD_CHUNK \
    MFM_VALIDATION_PACKED_BYTE_CHUNK MFM_VALIDATION_CACHE_CHUNKS MFM_MEMORY_RESERVE_GIB; do
    value=${!name}
    if [[ "$name" == MFM_VALIDATION_CACHE_CHUNKS && "$value" == 0 ]]; then
      continue
    fi
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || {
      echo "[FATAL] $name must be a positive integer: $value" >&2; return 1;
    }
  done
  export MFM_MEMORY_PROFILE_REQUESTED MFM_MEMORY_PROFILE MFM_MEMORY_PROFILE_RESOLVED
  export MFM_MEMORY_AVAILABLE_GIB MFM_MEMORY_RESERVE_GIB MFM_SORT_MEM_DEFAULT
  export MFM_STAGE04_MARKDUP_HEAP_GB_DEFAULT MFM_RNA_MARKDUP_HEAP_GB_DEFAULT
  export MFM_TRACE_SQLITE_BATCH MFM_TRACE_PREFETCH_RECORDS MFM_VALIDATION_RECORD_CHUNK
  export MFM_VALIDATION_PACKED_BYTE_CHUNK MFM_VALIDATION_CACHE_CHUNKS
  export MFM_WGS_SAMPLE_MEMORY_GIB_DEFAULT MFM_WGS_HAP_MEMORY_GIB_DEFAULT MFM_RNA_SAMPLE_MEMORY_GIB_DEFAULT
}
