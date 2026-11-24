#!/usr/bin/env bash
# 공용 cleanup backend — WGS Stage 24 / RNA Stage 26에서 호출한다.
#
# 삭제 범위:
# - INPUT_ROOT 아래 원본 WGS/RNA/VCF/reference는 건드리지 않는다.
# - final sample root 안의 대형 중간 디렉터리와 sample work root만 제거한다.
# - compact zarr/RNA zarr/metadata/logs/qc는 보존한다.
#
# 환경변수:
#   CLEAN_INTERMEDIATES_ON_SUCCESS=1  # 1이면 cleanup 실행
#   CLEAN_INTERMEDIATES_DRY_RUN=0     # 1이면 삭제하지 않고 대상만 출력
#   CLEAN_AFTER_RNA_ONLY=1            # RNA sample은 WGS 직후 cleanup을 미루고 RNA 후 실행
#   CLEAN_SHARED_ASSETS_ON_SUCCESS=0  # 1이면 windows/haplotagged/longread_calling 과 consensus stderr 삭제
#                                     # (HAP FASTA/chain/personal_fasta_qc.json 은 항상 보존)
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample> [wgs|rna]}
MODE=${2:-wgs}
CLEANUP_STAGE=24
[[ "$MODE" == rna ]] && CLEANUP_STAGE=26
ROOT="$(sample_final_root "$S")/$S"
WORK="$(sample_work_root "$S")"
# 12_compact_zarr.sh/direct_compact_writer의 현재 layout 이름.
# 예전 산출물 호환을 위해 old layout도 fallback으로 허용한다.
LAYOUT_NAME=${COMPACT_LAYOUT_NAME:-full_shared_hap_4bit_acgtnp_depth}
LAYOUT="$ROOT/compact/$LAYOUT_NAME"
if [[ ! -d "$LAYOUT" && -d "$ROOT/compact/full_shared_hap_4bit_acgtnp" ]]; then
  LAYOUT="$ROOT/compact/full_shared_hap_4bit_acgtnp"
fi
DRY_RUN=${CLEAN_INTERMEDIATES_DRY_RUN:-0}
CLEANUP_REPORT="$ROOT/metadata/cleanup_intermediates.tsv"
# ── 영구 자산 (2026-09-28) ─────────────────────────────────────────────────────
# HAP1/HAP2 FASTA(.fa, .fa.fai), chain, personal_fasta_qc.json, permanent_assets.json 은 어떤 옵션에서도
# 지우지 않는다. SNV+indel 을 적용한 HAP 좌표는 reference 와 달라서, chain 이 있어야 GENCODE 가 바뀌어도
# GTF 만 다시 lift 하고 WGS 01–12 없이 RNA 13 이후만 재실행할 수 있다.
# Stage 09이 만든 크기·sha256·chain/reference manifest를 cleanup 전에 다시 검증하고,
# 없거나 변경되었으면 cleanup을 중단한다. 삭제 직전에 manifest를 재생성해 변조를 승인하지 않는다.
PERMANENT_DIR="$ROOT/individual_fasta"

if [[ "${CLEAN_INTERMEDIATES_ON_SUCCESS:-1}" != 1 ]]; then
  log "${CLEANUP_STAGE} cleanup skip: CLEAN_INTERMEDIATES_ON_SUCCESS!=1 sample=$S"
  exit 0
fi

normalize_bool() {
  case "${1:-}" in
    1|true|TRUE|yes|YES|y|Y) echo 1 ;;
    *) echo 0 ;;
  esac
}

expand_chroms_for_cleanup() {
  local raw=${1:-$CHROMS} part left right i
  raw=${raw//,/ }
  for part in $raw; do
    [[ -n "$part" ]] || continue
    if [[ "$part" == *-* ]]; then
      left=${part%-*}; right=${part#*-}
      left=${left#chr}; right=${right#chr}
      for ((i=left; i<=right; i++)); do echo "chr$i"; done
    else
      [[ "$part" == chr* ]] && echo "$part" || echo "chr$part"
    fi
  done
}

has_rna_input() {
  local rna
  rna=$(sample_rna "$S" 2>/dev/null || true)
  [[ -n "$rna" && "$rna" != "-" ]]
}

assert_wgs_final_ready() {
  local contig d missing=0
  [[ -d "$LAYOUT" ]] || { warn "${CLEANUP_STAGE} cleanup guard fail: compact layout 없음: $LAYOUT"; return 1; }
  while read -r contig; do
    [[ -n "$contig" ]] || continue
    d="$LAYOUT/$contig"
    for item in contig_metadata.json windows.tsv hap_sequence.zarr depth.zarr readset.HAP1.zarr readset.HAP2.zarr; do
      if [[ ! -e "$d/$item" ]]; then
        warn "${CLEANUP_STAGE} cleanup guard fail: missing $d/$item"
        missing=1
      fi
    done
  done < <(expand_chroms_for_cleanup "$CHROMS")
  [[ "$missing" == 0 ]]
}

assert_compact_manifest_ready() {
  local manifest="$LAYOUT/manifest.json" expected
  expected=$(expand_chroms_for_cleanup "$CHROMS" | paste -sd, -)
  [[ -s "$manifest" ]] || { warn "${CLEANUP_STAGE} cleanup guard fail: compact manifest 없음: $manifest"; return 1; }
  "$PYTHON" -c 'import json,sys; p,e,v=sys.argv[1:]; d=json.load(open(p)); exp=set(filter(None,e.split(","))); rec=d.get("records",{}); got={str(x.get("contig","")) for x in rec.values() if isinstance(x,dict) and x.get("status")=="pass"}; assert d.get("status")=="complete", "manifest status is not complete"; assert d.get("validation_level")==v=="full", "full validation required"; assert got==exp, f"pass contigs mismatch expected={sorted(exp)} got={sorted(got)}"' "$manifest" "$expected" "$COMPACT_VALIDATION_LEVEL" || { warn "${CLEANUP_STAGE} cleanup guard fail: compact manifest/full validation 불일치"; return 1; }
}

assert_source_trace_ready() {
  local trace_manifest="$ROOT/qc/stage12_trace/$HAP_READ_MODE/trace_manifest.json" expected newer
  if [[ "$STAGE12_TRACE_POLICY" == off && "$ALLOW_CLEANUP_WITHOUT_TRACE" == 1 ]]; then
    warn "${CLEANUP_STAGE} cleanup: source trace 없이 명시적 override로 진행"
    return 0
  fi
  [[ -s "$trace_manifest" ]] || { warn "${CLEANUP_STAGE} cleanup guard fail: source trace manifest 없음: $trace_manifest"; return 1; }
  expected=$(expand_chroms_for_cleanup "$CHROMS" | paste -sd, -)
  "$PYTHON" -c 'import json,os,sys; p,r,e=sys.argv[1:]; d=json.load(open(p)); exp=set(filter(None,e.split(","))); got={k for k,v in d.get("contigs",{}).items() if v=="pass"}; assert d.get("status")=="pass", "trace status is not pass"; assert os.path.realpath(d.get("compact_root",""))==os.path.realpath(r), "trace compact_root mismatch"; assert got==exp, f"trace contigs mismatch expected={sorted(exp)} got={sorted(got)}"' "$trace_manifest" "$LAYOUT" "$expected" || { warn "${CLEANUP_STAGE} cleanup guard fail: source trace 불일치"; return 1; }
  newer=$(find "$LAYOUT" -type f -newer "$trace_manifest" -print -quit)
  [[ -z "$newer" ]] || { warn "${CLEANUP_STAGE} cleanup guard fail: trace 이후 compact 파일 변경: $newer"; return 1; }
}

assert_cleanup_predelete_ready() {
  "$PYTHON" ./scripts/cleanup_guard.py --root "$ROOT" --sample "$S" --mode "$MODE" || {
    warn "${CLEANUP_STAGE} cleanup guard fail: 최종 검증/HAP ASE provenance 불일치"
    return 1
  }
}

assert_rna_final_ready() {
  local mtag d1 d2
  mtag=$(sample_rna_evidence_label "$S")
  d1="$ROOT/rna_personal/$mtag/HAP1"
  d2="$ROOT/rna_personal/$mtag/HAP2"
  [[ -d "$d1/coverage.zarr" ]] || { warn "${CLEANUP_STAGE} cleanup guard fail: RNA HAP1 genome coverage 없음: $d1/coverage.zarr"; return 1; }
  [[ -d "$d2/coverage.zarr" ]] || { warn "${CLEANUP_STAGE} cleanup guard fail: RNA HAP2 genome coverage 없음: $d2/coverage.zarr"; return 1; }
  [[ -d "$d1/coverage_windows.zarr" ]] || { warn "${CLEANUP_STAGE} cleanup guard fail: RNA HAP1 window coverage 없음: $d1/coverage_windows.zarr"; return 1; }
  [[ -d "$d2/coverage_windows.zarr" ]] || { warn "${CLEANUP_STAGE} cleanup guard fail: RNA HAP2 window coverage 없음: $d2/coverage_windows.zarr"; return 1; }
  [[ -s "$ROOT/rna_personal/$mtag/hap_variants.tsv" ]] || { warn "${CLEANUP_STAGE} cleanup guard fail: hap_variants 없음: $ROOT/rna_personal/$mtag/hap_variants.tsv"; return 1; }
  [[ -s "$ROOT/rna_personal/$mtag/hap_variants_qc.json" ]] || { warn "${CLEANUP_STAGE} cleanup guard fail: hap_variants QC 없음: $ROOT/rna_personal/$mtag/hap_variants_qc.json"; return 1; }
  [[ -s "$ROOT/rna_personal/$mtag/allelic_counts.tsv" ]] || { warn "${CLEANUP_STAGE} cleanup guard fail: allelic_counts 없음: $ROOT/rna_personal/$mtag/allelic_counts.tsv"; return 1; }
  [[ -s "$ROOT/rna_personal/$mtag/ase_qc.json" ]] || { warn "${CLEANUP_STAGE} cleanup guard fail: ASE QC 없음: $ROOT/rna_personal/$mtag/ase_qc.json"; return 1; }
  [[ -s "$ROOT/rna_personal/$mtag/ase_training_eligibility.tsv" ]] || { warn "${CLEANUP_STAGE} cleanup guard fail: ASE 적격성 없음: $ROOT/rna_personal/$mtag/ase_training_eligibility.tsv"; return 1; }
}

copy_work_logs() {
  if [[ -d "$WORK/logs" ]]; then
    mkdir -p "$ROOT/logs/stage_logs"
    cp -a "$WORK/logs/." "$ROOT/logs/stage_logs/" 2>/dev/null || true
  fi
}

safe_target() {
  local target=$1 safe_work_root
  safe_work_root=$(readlink -m "$OUTPUT_ROOT/_work")
  [[ -n "$target" ]] || return 1
  [[ "$target" != "/" ]] || return 1
  [[ "$target" != "$ROOT" ]] || return 1
  case "$target" in
    "$ROOT"/*) return 0 ;;
    "$safe_work_root"/*) return 0 ;;
    *) return 1 ;;
  esac
}

is_permanent_target() {
  # 영구 자산 디렉터리 자체, 그 상위, 또는 그 안의 FASTA/fai/chain/QC 는 영구 자산이다.
  local target perm
  target=$(readlink -m "$1"); perm=$(readlink -m "$PERMANENT_DIR")
  [[ "$perm/" == "$target/"* ]] && return 0
  if [[ "$target" == "$perm/"* ]]; then
    case "$target" in
      *.fa|*.fa.fai|*.chain|*/personal_fasta_qc.json|*/permanent_assets.json|*/HAP1|*/HAP2) return 0 ;;
    esac
  fi
  return 1
}

verify_permanent_manifest() {
  [[ -d "$PERMANENT_DIR" ]] || { warn "${CLEANUP_STAGE} cleanup guard fail: 영구 자산 디렉터리 없음: $PERMANENT_DIR"; return 1; }
  [[ -s "$PERMANENT_DIR/permanent_assets.json" ]] || {
    warn "${CLEANUP_STAGE} cleanup guard fail: Stage 09 영구 자산 manifest 없음"; return 1;
  }
  "$PYTHON" ./tools/personal_assets.py verify --individual-fasta-dir "$PERMANENT_DIR" || {
    warn "${CLEANUP_STAGE} cleanup guard fail: 영구 자산 manifest/hash 검증 실패"; return 1;
  }
  printf "%s\t%s\t%s\t%s\tpreserved_permanent\t%s\n" "$(date '+%F %T')" "$S" "$MODE" "$DRY_RUN" "$PERMANENT_DIR" >> "$CLEANUP_REPORT"
}

remove_individual_fasta_scratch() {
  # 영구 자산이 아닌 consensus stderr 로그만 logs 로 옮긴 뒤 지운다. FASTA/chain/QC 는 남긴다.
  local f
  [[ -d "$PERMANENT_DIR" ]] || return 0
  mkdir -p "$ROOT/logs/personal_fasta"
  while IFS= read -r -d '' f; do
    cp -a "$f" "$ROOT/logs/personal_fasta/" 2>/dev/null || true
    remove_target "$f"
  done < <(find "$PERMANENT_DIR" -type f -name '*.bcftools_consensus.stderr' -print0)
}

remove_target() {
  local target=$1
  [[ -e "$target" ]] || return 0
  if ! safe_target "$target"; then
    die "${CLEANUP_STAGE} cleanup unsafe target 거부: $target"
  fi
  if is_permanent_target "$target"; then
    die "${CLEANUP_STAGE} cleanup 영구 자산 삭제 거부: $target"
  fi
  if [[ "$(normalize_bool "$DRY_RUN")" == 1 ]]; then
    log "${CLEANUP_STAGE} cleanup dry-run remove: $target"
    printf "%s\t%s\t%s\t%s\tdry_run_remove\t%s\n" "$(date '+%F %T')" "$S" "$MODE" "$DRY_RUN" "$target" >> "$CLEANUP_REPORT"
  else
    log "${CLEANUP_STAGE} cleanup remove: $target"
    rm -rf -- "$target"
    printf "%s\t%s\t%s\t%s\tremoved\t%s\n" "$(date '+%F %T')" "$S" "$MODE" "$DRY_RUN" "$target" >> "$CLEANUP_REPORT"
  fi
}

main() {
  [[ -d "$ROOT" ]] || die "${CLEANUP_STAGE} cleanup sample final root 없음: $ROOT"

  if [[ "$MODE" == "wgs" && "${CLEAN_AFTER_RNA_ONLY:-1}" == 1 ]] && has_rna_input; then
    log "${CLEANUP_STAGE} cleanup skip: $S 는 RNA sample이므로 RNA 완료 후 cleanup 수행"
    exit 0
  fi

  # 예전 PASS JSON만 신뢰하지 않는다. 삭제 직전 현재 파일의 schema/hash/freshness를 직접 검사한다.
  assert_cleanup_predelete_ready || die "${CLEANUP_STAGE} cleanup 중단: 삭제 전 최종 계약 미통과 sample=$S"
  assert_wgs_final_ready || die "${CLEANUP_STAGE} cleanup 중단: WGS 최종 compact가 완전하지 않음 sample=$S"
  assert_compact_manifest_ready || die "${CLEANUP_STAGE} cleanup 중단: compact manifest/full 검증 미통과 sample=$S"
  assert_source_trace_ready || die "${CLEANUP_STAGE} cleanup 중단: 독립 source trace 미통과 sample=$S"
  if [[ "$MODE" == "rna" ]]; then
    assert_rna_final_ready || die "${CLEANUP_STAGE} cleanup 중단: RNA 최종 산출물이 완전하지 않음 sample=$S"
  fi

  mkdir -p "$ROOT/metadata"
  printf "timestamp\tsample\tmode\tdry_run\taction\tpath\n" > "$CLEANUP_REPORT"

  copy_work_logs
  verify_permanent_manifest || die "${CLEANUP_STAGE} cleanup 중단: 영구 자산(HAP FASTA/chain/QC) 검증 실패 sample=$S"

  # final model input/metadata/logs/qc/rna_personal 과 영구 자산(individual_fasta 의 FASTA/chain/QC)은 보존한다.
  remove_target "$ROOT/$REALIGNED_SUBDIR"
  if [[ "${CLEAN_SHARED_ASSETS_ON_SUCCESS:-0}" == 1 ]]; then
    remove_individual_fasta_scratch
    remove_target "$ROOT/windows"
    remove_target "$ROOT/haplotagged"
    remove_target "$ROOT/longread_calling"
  fi
  remove_target "$ROOT/work"
  remove_target "$WORK"

  log "${CLEANUP_STAGE} cleanup 완료: sample=$S mode=$MODE dry_run=$DRY_RUN"
}

main
