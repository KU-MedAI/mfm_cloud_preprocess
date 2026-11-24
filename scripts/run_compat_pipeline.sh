#!/usr/bin/env bash
# integrated/process_code v2 대표 실행 entrypoint.
#
# 예시:
#   ./run_pipeline.sh env verify
#   ./run_pipeline.sh wgs HG00100
#   ./run_pipeline.sh wgs-list /path/samples.txt
#   ./run_pipeline.sh slurm-wgs /path/samples.txt
#   ./run_pipeline.sh smoke
#   ./run_pipeline.sh smoke-wgs HG00100
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$SCRIPT_DIR"
# 배치 자동 실행은 공용 handoff와 같은 local/Slurm 진입점을 사용함.
if (($# == 0)); then set -- run; fi
case "${1:-}" in
  run|dry-run|doctor|plan|resume|import-state)
    exec bash "$SCRIPT_DIR/run_handoff.sh" "$@"
    ;;
esac
source ./config_v2.sh
DOWNLOAD_1000G_ROOT="${DOWNLOAD_1000G_ROOT:-$SCRIPT_DIR/download_1000g}"
CURRENT_SHARE_ROOT="${CURRENT_SHARE_ROOT:-$SCRIPT_DIR/share_package}"

usage() {
  cat <<'EOF'
사용법:
  ./run_pipeline.sh dry-run|run --target azure_poc|azure_main|aws_main|server [--sample ID] [--config FILE]
  ./run_pipeline.sh dry-run|run --lab medai|tgil --batch 001 --config FILE [--executor local|slurm]
  ./run_pipeline.sh doctor --config FILE
  ./run_pipeline.sh env <verify|verify-core|verify-clair3|verify-rna|verify-denovo|install-core|install-clair3|install-rna|install-denovo|install-all|print>
  ./run_pipeline.sh sample <sample> [--manifest FILE] [--input-root DIR] [--output-root DIR]
  ./run_pipeline.sh sample-list <sample_list> [sample_count] [portable options]
  ./run_pipeline.sh sample-list --manifest FILE --sample-list FILE [--input-root DIR] [--output-root DIR] [--sample-count N|all]
  ./run_pipeline.sh build-portable-manifest [--registry FILE] [--out-dir DIR]
  ./run_pipeline.sh wgs <sample>
  ./run_pipeline.sh wgs-list <sample_list> [sample_count]
  ./run_pipeline.sh wgs-phased-list <sample_list> [sample_count]
  ./run_pipeline.sh clair3-qc <sample>
  ./run_pipeline.sh clair3-qc-list <sample_list> [sample_count]
  ./run_pipeline.sh slurm-wgs <sample_list> [array_spec]
  ./run_pipeline.sh rna <sample>
  ./run_pipeline.sh rna-list <sample_list> [sample_count]
  ./run_pipeline.sh slurm-rna <sample_list> [array_spec]
  ./run_pipeline.sh download-plan
  ./run_pipeline.sh download-list <phase|all>
  ./run_pipeline.sh download <phase> [all|shard_id]
  ./run_pipeline.sh verify-download <phase> [all|shard_id]
  ./run_pipeline.sh download-vcf
  ./run_pipeline.sh current-validate-package
  ./run_pipeline.sh current-download-list <phase|all>
  ./run_pipeline.sh current-download <phase>
  ./run_pipeline.sh current-verify-download <phase>
  ./run_pipeline.sh current-download-vcf
  ./run_pipeline.sh current-preprocess-phase <phase|all>
  ./run_pipeline.sh legacy-download-plan
  ./run_pipeline.sh legacy-download-list <phase|all>
  ./run_pipeline.sh legacy-download <phase> [all|shard_id]
  ./run_pipeline.sh legacy-verify-download <phase> [all|shard_id]
  ./run_pipeline.sh legacy-download-vcf
  ./run_pipeline.sh legacy-preprocess-phase <phase|all> [all|shard_id]
  ./run_pipeline.sh preprocess-phase <phase|all> [all|shard_id]
  ./run_pipeline.sh smoke
  ./run_pipeline.sh smoke-wgs <sample>

주요 환경변수:
  PIPELINE_MANIFEST=/path/manifest.tsv
  INPUT_ROOT=/path/input  OUTPUT_ROOT=/path/output
  CHROMS=chr1-chr22
  LOCAL_JOBS=2                 # WGS/Clair3 기본; RNA 기본은 12
  WGS_STOP_AFTER=01|02|03|04|05|06|07|08|09|10|11|12|23|24   # 번호=파일 접두사=실행 순서
  VCF_SOURCE_POLICY=manifest         # production 기본값
  ENABLE_CLAIR3_SELF_CALL=0          # Clair3 재현 실행만 1
  # manifest에서는 vcf_source, public_vcf_path, public_vcf_sample_id,
  # vcf_reference_fasta, vcf_callset, phase_source 열을 사용한다.

1000G 다운로드/phase 환경변수:
  DOWNLOAD_1000G_ROOT=./download_1000g
  PHASES=5
  SHARD_GIB=1500
  WGS_SOURCE_MODE=assigned|primary|all_manifest_runs|aligned_only
  RNA_SOURCE_MODE=strict_bam|primary|available|none
  INCLUDE_BIGWIG=0  # recount3 BigWig은 운영 대상에서 제외한다. 1이면 에러 처리한다.

smoke:
  빠른 계약 검증이다. shell/python/core env를 검사한다. Clair3 dry-run은
  ENABLE_CLAIR3_SELF_CALL=1로 요청한 경우에만 수행한다.

smoke-wgs:
  실제 작은 WGS 실행 모드다. 기본은 CHROMS=chr22, COMPACT_MAX_WINDOWS=1이다.
  Clair3 실제 calling이 필요하면 CLAIR3_DRY_RUN=0과 CLAIR3_BED_FN을 명시한다.
EOF
}


run_1000g_download_script() {
  local script=${1:?script를 지정하세요}
  shift || true
  [[ -d "$DOWNLOAD_1000G_ROOT" ]] || die "1000G 다운로드 root가 없습니다: $DOWNLOAD_1000G_ROOT"
  [[ -x "$DOWNLOAD_1000G_ROOT/scripts/$script" ]] || die "1000G 다운로드 스크립트가 없습니다: $DOWNLOAD_1000G_ROOT/scripts/$script"
  (cd "$DOWNLOAD_1000G_ROOT" && bash "./scripts/$script" "$@")
}
run_current_share_script() {
  local script=${1:?script를 지정하세요}
  shift || true
  [[ -d "$CURRENT_SHARE_ROOT" ]] || die "current share package root가 없습니다: $CURRENT_SHARE_ROOT"
  [[ -x "$CURRENT_SHARE_ROOT/scripts/$script" ]] || die "current share package script가 없습니다: $CURRENT_SHARE_ROOT/scripts/$script"
  (cd "$CURRENT_SHARE_ROOT" && bash "./scripts/$script" "$@")
}

warn_legacy_download_alias() {
  warn "download-* 비접두 명령은 historical 1KGP 3,202 downloader입니다. current 3,212 cohort는 current-* 명령을 사용하세요."
}

run_1000g_preprocess_phase() {
  local phase=${1:?phase|all을 지정하세요}
  local shard=${2:-all}
  export PROCESS_CODE_DIR="${PROCESS_CODE_DIR:-$SCRIPT_DIR}"
  if [[ "${RUN_PREPROCESS:-0}" == "1" && -z "${PREPROCESS_CMD:-}" ]]; then
    case "${MODE:-all}" in
      pre|wgs_only)
        export PREPROCESS_CMD='./run_pipeline.sh wgs-list "$SAMPLE_LIST"'
        ;;
      post|rna_pair)
        die "post/rna_pair는 WGS+RNA 실행 순서를 명확히 해야 하므로 PREPROCESS_CMD를 직접 지정하세요"
        ;;
      *)
        die "RUN_PREPROCESS=1에서는 MODE=pre|wgs_only 또는 PREPROCESS_CMD를 명시하세요"
        ;;
    esac
  fi
  run_1000g_download_script 05_preprocess_phase.sh "$phase" "$shard"
}

run_syntax_check() {
  log "smoke: shell/sbatch syntax check"
  bash -n ./*.sh ./*.sbatch
  log "smoke: python compile check"
  "$PYTHON_BIN" -m py_compile ./src/integrate_wgs_rna/*.py
}

run_fast_smoke() {
  local sample align
  run_syntax_check
  ./00_environment.sh verify-core
  if [[ "$ENABLE_CLAIR3_SELF_CALL" != 1 ]]; then
    log "smoke 완료: syntax/python/core env 통과 (Clair3는 기본 비활성)"
    return 0
  fi
  ./00_environment.sh verify-clair3
  sample=${SMOKE_SAMPLE:-HG00096}
  align=${SMOKE_ALIGNMENT:-}
  if [[ ! -s "$align" ]]; then
    warn "smoke alignment 없음. Clair3 dry-run wrapper 검증 생략: $align"
    return 0
  fi
  MFM_SMOKE_ROOT=$(mktemp -d /tmp/mfm_pipeline_smoke.XXXXXX)
  trap 'rm -rf "${MFM_SMOKE_ROOT:-}"' EXIT
  log "smoke: Clair3 wrapper dry-run sample=$sample root=$MFM_SMOKE_ROOT"
  DATASET_ROOT="$MFM_SMOKE_ROOT" \
  RUN_ID=smoke \
  DEFAULT_COHORT=smoke_ont \
  DEFAULT_PROJECT_ID=smoke_ont \
  DEFAULT_PLATFORM=ont \
  DEFAULT_READ_TYPE=long_read \
  ENABLE_CLAIR3_SELF_CALL=1 \
  VCF_SOURCE_POLICY=platform_matched \
  LONGREAD_SEQUENCE_VCF_SOURCE=clair3 \
  ALIGNMENT_PATH="$align" \
  CHROMS=${SMOKE_CHROMS:-chr22} \
  CLAIR3_DRY_RUN=1 \
  CLAIR3_THREADS=${SMOKE_THREADS:-2} \
  CLAIR3_OUTPUT_THREADS=1 \
  ./04_call_longread_clair3.sh "$sample"
  log "smoke 완료: syntax/python/core env/명시적 Clair3 dry-run 통과"
}

run_wgs_smoke() {
  local sample=${1:-${SMOKE_SAMPLE:-HG00096}}
  local smoke_root source_alignment smoke_region smoke_bam smoke_bed smoke_common
  smoke_root=${SMOKE_ROOT:-$(mktemp -d /tmp/mfm_wgs_smoke.XXXXXX)}
  if [[ "${SMOKE_KEEP:-0}" != 1 ]]; then
    MFM_WGS_SMOKE_ROOT="$smoke_root"
    trap 'rm -rf "${MFM_WGS_SMOKE_ROOT:-}"' EXIT
  fi
  mkdir -p "$smoke_root/work/slice" "$smoke_root/work/common_vcf" "$smoke_root/work/persample_vcf"
  source_alignment=${SMOKE_SOURCE_ALIGNMENT:-}
  smoke_region=${SMOKE_REGION:-chr22:20000000-21000000}
  smoke_bam=${SMOKE_ALIGNMENT:-$smoke_root/work/slice/$sample.${smoke_region//[:\-]/_}.bam}
  smoke_bed=$smoke_root/work/slice/${sample}.${smoke_region//[:\-]/_}.bed
  smoke_common=${SMOKE_COMMON_VCF:-}

  [[ -s "$source_alignment" ]] || die "SMOKE_SOURCE_ALIGNMENT 없음: $source_alignment"
  [[ -s "$smoke_common" ]] || die "SMOKE_COMMON_VCF 없음: $smoke_common"
  if [[ ! -s "$smoke_bam" ]]; then
    log "smoke-wgs: CRAM/BAM slice 생성 source=$source_alignment region=$smoke_region"
    "$SAMTOOLS_BIN" view -b -T "$REF_DECODE" -o "$smoke_bam.part" "$source_alignment" "$smoke_region"
    mv "$smoke_bam.part" "$smoke_bam"
    "$SAMTOOLS_BIN" index -@ 2 "$smoke_bam"
  fi
  local contig range start end
  contig=${smoke_region%%:*}
  range=${smoke_region#*:}
  start=${range%-*}
  end=${range#*-}
  printf '%s\t%s\t%s\n' "$contig" "$((start-1))" "$end" > "$smoke_bed"
  ln -sf "$smoke_common" "$smoke_root/work/common_vcf/smoke_ont.$contig.vcf.gz"
  if [[ -s "$smoke_common.tbi" ]]; then
    ln -sf "$smoke_common.tbi" "$smoke_root/work/common_vcf/smoke_ont.$contig.vcf.gz.tbi"
  elif [[ -s "$smoke_common.csi" ]]; then
    ln -sf "$smoke_common.csi" "$smoke_root/work/common_vcf/smoke_ont.$contig.vcf.gz.csi"
  fi

  log "smoke-wgs: 전체 WGS pipeline 실행 sample=$sample root=$smoke_root region=$smoke_region"
  DATASET_ROOT="$smoke_root" \
  RUN_ID=smoke_wgs \
  DEFAULT_COHORT=smoke_ont \
  DEFAULT_PROJECT_ID=smoke_ont \
  DEFAULT_PROJECT_LABEL="smoke ONT" \
  DEFAULT_DATASET_ROLE=pre_training \
  DEFAULT_PLATFORM=ont \
  DEFAULT_READ_TYPE=long_read \
  ENABLE_CLAIR3_SELF_CALL=1 \
  VCF_SOURCE_POLICY=platform_matched \
  LONGREAD_SEQUENCE_VCF_SOURCE=clair3 \
  ALIGNMENT_PATH="$smoke_bam" \
  COMMON_VCF_DIR="$smoke_root/work/common_vcf" \
  PERSAMPLE_VCF_DIR="$smoke_root/work/persample_vcf" \
  VCF_CONCAT_DIR="$smoke_root/work/vcf_concat" \
  SAMPLE_VCF_FAST_VALIDATE=1 \
  CLAIR3_BED_FN="$smoke_bed" \
  CLAIR3_THREADS=${SMOKE_CLAIR3_THREADS:-4} \
  CLAIR3_OUTPUT_THREADS=${SMOKE_CLAIR3_OUTPUT_THREADS:-2} \
  CLAIR3_ALLOW_EMPTY_HAPLOTAG=1 \
  CLAIR3_ALLOW_PILEUP_FALLBACK=1 \
  CHROMS="${SMOKE_CHROMS:-$contig}" \
  WINDOW_JOBS=${WINDOW_JOBS:-2} \
  WGS_REALIGN_THREADS=${WGS_REALIGN_THREADS:-4} \
  WGS_SORT_THREADS=${WGS_SORT_THREADS:-2} \
  COMPACT_JOBS=${COMPACT_JOBS:-2} \
  COMPACT_ALLOW_UNPHASED_LONGREAD=1 \
  SEQUENCE_SOURCE=vcf_consensus \
  ./run_wgs_sample.sh "$sample"
  log "smoke-wgs 완료: root=$smoke_root"
}

configure_portable_roots() {
  local output_explicit=${1:-0}
  [[ -z "${PIPELINE_MANIFEST:-}" ]] || PIPELINE_MANIFEST=$(readlink -m "$PIPELINE_MANIFEST")
  [[ -z "${SAMPLE_LIST:-}" ]] || SAMPLE_LIST=$(readlink -m "$SAMPLE_LIST")
  [[ -z "${INPUT_ROOT:-}" ]] || INPUT_ROOT=$(readlink -m "$INPUT_ROOT")
  OUTPUT_ROOT=$(readlink -m "${OUTPUT_ROOT:-$MFM_DEFAULT_OUTPUT_ROOT}")
  if [[ "$output_explicit" == 1 ]]; then OUTPUT_ROOT_OVERRIDE=1; fi
  DATASET_ROOT=$OUTPUT_ROOT
  export PIPELINE_MANIFEST SAMPLE_LIST INPUT_ROOT OUTPUT_ROOT OUTPUT_ROOT_OVERRIDE DATASET_ROOT
}

parse_sample_command() {
  SAMPLE_ARG=
  local output_explicit=0
  while (($#)); do
    case "$1" in
      --sample) SAMPLE_ARG=${2:?--sample 값이 필요함}; shift 2 ;;
      --manifest) PIPELINE_MANIFEST=${2:?--manifest 값이 필요함}; shift 2 ;;
      --input-root) INPUT_ROOT=${2:?--input-root 값이 필요함}; shift 2 ;;
      --output-root) OUTPUT_ROOT=${2:?--output-root 값이 필요함}; output_explicit=1; shift 2 ;;
      --) shift; break ;;
      -*) die "알 수 없는 sample 옵션: $1" ;;
      *) [[ -z "$SAMPLE_ARG" ]] || die "sample은 하나만 지정: $1"; SAMPLE_ARG=$1; shift ;;
    esac
  done
  [[ -n "$SAMPLE_ARG" ]] || die "sample을 지정하세요"
  [[ -n "${PIPELINE_MANIFEST:-}" && -s "$PIPELINE_MANIFEST" ]] || die "--manifest 파일이 필요합니다: ${PIPELINE_MANIFEST:-}"
  configure_portable_roots "$output_explicit"
}

parse_sample_list_command() {
  LIST_ARG=
  COUNT_ARG=${SAMPLE_COUNT:-all}
  local output_explicit=0 positional=0
  while (($#)); do
    case "$1" in
      --manifest) PIPELINE_MANIFEST=${2:?--manifest 값이 필요함}; shift 2 ;;
      --sample-list) LIST_ARG=${2:?--sample-list 값이 필요함}; shift 2 ;;
      --sample-count) COUNT_ARG=${2:?--sample-count 값이 필요함}; shift 2 ;;
      --input-root) INPUT_ROOT=${2:?--input-root 값이 필요함}; shift 2 ;;
      --output-root) OUTPUT_ROOT=${2:?--output-root 값이 필요함}; output_explicit=1; shift 2 ;;
      --) shift; break ;;
      -*) die "알 수 없는 sample-list 옵션: $1" ;;
      *)
        if (( positional == 0 )); then LIST_ARG=$1; positional=1
        elif (( positional == 1 )); then COUNT_ARG=$1; positional=2
        else die "불필요한 positional 인자: $1"
        fi
        shift
        ;;
    esac
  done
  [[ -n "$LIST_ARG" && -s "$LIST_ARG" ]] || die "--sample-list 파일이 필요합니다: $LIST_ARG"
  [[ -n "${PIPELINE_MANIFEST:-}" && -s "$PIPELINE_MANIFEST" ]] || die "--manifest 파일이 필요합니다: ${PIPELINE_MANIFEST:-}"
  SAMPLE_LIST=$LIST_ARG
  configure_portable_roots "$output_explicit"
}

mode=${1:-help}
shift || true
case "$mode" in
  env)
    ./00_environment.sh "${1:-verify}"
    ;;
  sample)
    parse_sample_command "$@"
    ./run_sample.sh "$SAMPLE_ARG"
    ;;
  sample-list)
    parse_sample_list_command "$@"
    ./runner_local.sh sample "$SAMPLE_LIST" "$COUNT_ARG"
    ;;
  build-portable-manifest)
    exec "$PYTHON_BIN" "$SCRIPT_DIR/scripts/build_portable_manifest.py" "$@"
    ;;
  wgs)
    sample=${1:?sample을 지정하세요}
    ./run_wgs_sample.sh "$sample"
    ;;
  wgs-list)
    list=${1:?sample_list를 지정하세요}
    count=${2:-${SAMPLE_COUNT:-all}}
    ./runner_local.sh wgs "$list" "$count"
    ;;
  wgs-phased-list)
    list=${1:?sample_list를 지정하세요}
    count=${2:-${SAMPLE_COUNT:-all}}
    ./run_wgs_phased.sh "$list" "$count"
    ;;
  clair3-qc)
    sample=${1:?sample을 지정하세요}
    ./05_clair3_qc.sh "$sample"
    ;;
  clair3-qc-list)
    list=${1:?sample_list를 지정하세요}
    count=${2:-${SAMPLE_COUNT:-all}}
    ./runner_local.sh clair3-qc "$list" "$count"
    ;;
  slurm-wgs)
    list=${1:?sample_list를 지정하세요}
    array=${2:-1-$(wc -l < "$list")%4}
    sbatch --array="$array" --export=ALL,SAMPLE_LIST="$list" ./slurm_wgs.sbatch
    ;;
  rna)
    sample=${1:?sample을 지정하세요}
    ./run_rna_sample.sh "$sample"
    ;;
  rna-list)
    list=${1:?sample_list를 지정하세요}
    count=${2:-${SAMPLE_COUNT:-all}}
    ./runner_local.sh rna "$list" "$count"
    ;;
  slurm-rna)
    list=${1:?sample_list를 지정하세요}
    [[ -s "$list" ]] || die "sample list가 없거나 비어 있음: $list"
    list=$(readlink -f "$list")
    valid_count=$(awk 'NF && $1 !~ /^#/ {n++} END {print n+0}' "$list")
    [[ "$valid_count" -gt 0 ]] || die "유효 sample이 없음: $list"
    array=${2:-1-${valid_count}%4}
    sbatch --array="$array" --export=ALL,SAMPLE_LIST="$list" ./slurm_rna.sbatch
    ;;
  download-plan|legacy-download-plan)
    [[ "$mode" == legacy-* ]] || warn_legacy_download_alias
    run_1000g_download_script 00_make_phase_assignments.sh "$@"
    ;;
  download-list|legacy-download-list)
    [[ "$mode" == legacy-* ]] || warn_legacy_download_alias
    phase=${1:?phase|all을 지정하세요}
    run_1000g_download_script 01_make_phase_url_lists.sh "$phase"
    ;;
  download|legacy-download)
    [[ "$mode" == legacy-* ]] || warn_legacy_download_alias
    phase=${1:?phase를 지정하세요}
    shard=${2:-all}
    run_1000g_download_script 02_download_phase.sh "$phase" "$shard"
    ;;
  verify-download|legacy-verify-download)
    [[ "$mode" == legacy-* ]] || warn_legacy_download_alias
    phase=${1:?phase를 지정하세요}
    shard=${2:-all}
    run_1000g_download_script 03_verify_phase.sh "$phase" "$shard"
    ;;
  download-vcf|legacy-download-vcf)
    [[ "$mode" == legacy-* ]] || warn_legacy_download_alias
    run_1000g_download_script 04_download_common_vcf.sh
    ;;
  preprocess-phase|legacy-preprocess-phase)
    [[ "$mode" == legacy-* ]] || warn_legacy_download_alias
    run_1000g_preprocess_phase "$@"
    ;;
  current-validate-package)
    run_current_share_script 00_validate_package.sh
    ;;
  current-download-list)
    phase=${1:?phase|all을 지정하세요}
    run_current_share_script 01_make_phase_url_lists.sh "$phase"
    ;;
  current-download)
    phase=${1:?phase를 지정하세요}
    run_current_share_script 02_download_phase.sh "$phase"
    ;;
  current-verify-download)
    phase=${1:?phase를 지정하세요}
    run_current_share_script 03_verify_phase.sh "$phase"
    ;;
  current-download-vcf)
    run_current_share_script 04_download_illumina_vcf.sh
    ;;
  current-preprocess-phase)
    phase=${1:?phase|all을 지정하세요}
    run_current_share_script 05_preprocess_phase.sh "$phase"
    ;;
  smoke)
    run_fast_smoke
    ;;
  smoke-wgs)
    run_wgs_smoke "$@"
    ;;
  help|-h|--help)
    usage
    ;;
  *)
    usage >&2
    die "unknown mode: $mode"
    ;;
esac
