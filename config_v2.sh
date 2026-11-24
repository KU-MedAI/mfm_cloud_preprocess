#!/usr/bin/env bash
# integrated WGS/RNA preprocessing v2 공통 설정.
#
# 이 파일은 process_code/v2의 config_v2.sh 역할을 integrated v2에 맞춘 것이다.
# root-level 단계별 shell wrapper는 모두 이 파일을 source한다.
set -euo pipefail
# native/NumPy 보조 thread는 기본1로 두고 단계별 병렬화 예산을 따름.
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNNER="$PIPELINE_DIR/scripts/run_integrated_preprocess.sh"
source "$PIPELINE_DIR/lib/data_paths.sh"
# 모든 단계 wrapper가 config를 직접 source하므로, memory profile도 여기서 한 번
# 초기화해야 한다. run_pipeline만 초기화한다고 가정하면 Stage 10 단독 shell에서
# MFM_MEMORY_PROFILE_RESOLVED가 unset이 되어 set -u로 중단한다.
export MFM_RAM_LIMIT_BYTES=${MFM_RAM_LIMIT_BYTES:-1650000000000}
source "$PIPELINE_DIR/lib/memory_profile.sh"
mfm_memory_profile_init
# Login shell이나 `conda run ... bash -l`이 PATH를 다시 정렬하더라도 core
# Python이 /usr/bin/python으로 떨어지지 않게 환경 경로를 명시적으로 고정한다.
MFM_CORE_ENV_NAME=${MFM_CORE_ENV_NAME:-jj_mfm}
MFM_CORE_ENV=${MFM_CORE_ENV:-}
if [[ -z "$MFM_CORE_ENV" ]]; then
  if [[ -n "${CONDA_PREFIX:-}" && -x "$CONDA_PREFIX/bin/python" ]]; then
    MFM_CORE_ENV=$CONDA_PREFIX
  fi
fi
if [[ -n "$MFM_CORE_ENV" && -d "$MFM_CORE_ENV/bin" ]]; then
  export PATH="$MFM_CORE_ENV/bin:$PATH"
fi
find_tool() {
  if [[ -n "${CONDA_PREFIX:-}" && -x "$CONDA_PREFIX/bin/$1" ]]; then
    printf '%s\n' "$CONDA_PREFIX/bin/$1"
  else
    command -v "$1" 2>/dev/null || true
  fi
}
CONDA_BIN=${CONDA_BIN:-$(find_tool conda)}
PYTHON_BIN=${PYTHON_BIN:-${MFM_CORE_ENV:+$MFM_CORE_ENV/bin/python}}
[[ -x "$PYTHON_BIN" ]] || PYTHON_BIN=$(find_tool python)
[[ -n "$PYTHON_BIN" ]] || PYTHON_BIN=$(find_tool python3)
if [[ -z "${MFM_REPO_ROOT:-}" ]]; then
  for candidate in     "$PIPELINE_DIR/vendor/MFM-DNARNAFM"
  do
    if [[ -d "$candidate/pre_training/preprocess_raw_readset" ]]; then
      MFM_REPO_ROOT=$candidate
      break
    fi
  done
  MFM_REPO_ROOT=${MFM_REPO_ROOT:-$PIPELINE_DIR/vendor/MFM-DNARNAFM}
fi

export MFM_CORE_ENV_NAME MFM_CORE_ENV PYTHON_BIN MFM_REPO_ROOT
if [[ "${MFM_RNA_BACKEND:-cpu}" == parabricks ]]; then
  export STAR_BIN="$PIPELINE_DIR/scripts/star_gpu.sh"
fi
# ── run/data root ────────────────────────────────────────────────────────────
RUN_ID=${RUN_ID:-$(date +%y%m%d)}
OUTPUT_ROOT=${OUTPUT_ROOT:-${DATASET_ROOT:-$MFM_DEFAULT_OUTPUT_ROOT}}
INPUT_ROOT=${INPUT_ROOT:-}
OUTPUT_ROOT_OVERRIDE=${OUTPUT_ROOT_OVERRIDE:-0}
DATASET_ROOT=${DATASET_ROOT:-$OUTPUT_ROOT}
# 선택적 local NVMe scratch. 비어 있으면 기존 output/work 경로를 사용한다.
MFM_SCRATCH_ROOT=${MFM_SCRATCH_ROOT:-}
if [[ -n "$MFM_SCRATCH_ROOT" ]]; then
  MFM_STAGE11_SCRATCH_ROOT=${MFM_STAGE11_SCRATCH_ROOT:-$MFM_SCRATCH_ROOT/stage11}
  MFM_STAGE12_SCRATCH_ROOT=${MFM_STAGE12_SCRATCH_ROOT:-$MFM_SCRATCH_ROOT/stage12}
else
  MFM_STAGE11_SCRATCH_ROOT=${MFM_STAGE11_SCRATCH_ROOT:-}
  MFM_STAGE12_SCRATCH_ROOT=${MFM_STAGE12_SCRATCH_ROOT:-}
fi
export MFM_SCRATCH_ROOT MFM_STAGE11_SCRATCH_ROOT MFM_STAGE12_SCRATCH_ROOT
PIPELINE_MANIFEST=${PIPELINE_MANIFEST:-}
SAMPLE_LIST=${SAMPLE_LIST:-}
case "$OUTPUT_ROOT_OVERRIDE" in 0|1) ;; *) echo "[FATAL] OUTPUT_ROOT_OVERRIDE must be 0|1" >&2; exit 1 ;; esac

# ── cohort/sample defaults. Manifest 값이 있으면 manifest가 우선한다. ────────
DEFAULT_COHORT=${DEFAULT_COHORT:-igsr_ont_wgs_only}
DEFAULT_PROJECT_ID=${DEFAULT_PROJECT_ID:-igsr_ont_wgs_only}
DEFAULT_PROJECT_LABEL=${DEFAULT_PROJECT_LABEL:-IGSR ONT WGS-only}
DEFAULT_DATASET_ROLE=${DEFAULT_DATASET_ROLE:-pre_training}
DEFAULT_READ_TYPE=${DEFAULT_READ_TYPE:-long_read}
DEFAULT_PLATFORM=${DEFAULT_PLATFORM:-ont}
DEFAULT_SAMPLE_STATUS=${DEFAULT_SAMPLE_STATUS:-healthy}
DEFAULT_DISEASE=${DEFAULT_DISEASE:-NaN}
REFERENCE_BUILD=${REFERENCE_BUILD:-hg38}

# ── reference roles ──────────────────────────────────────────────────────────
# REF_DECODE: CRAM decode, Clair3, WhatsHap, samtools view/index 등에 사용.
# REF_SEQUENCE: bcftools consensus로 personal FASTA를 만들 때 사용.
REF_DECODE=${REF_DECODE:-$MFM_REFERENCE_DIR/long_read_decode.fa}
REF_SEQUENCE=${REF_SEQUENCE:-$MFM_REFERENCE_DIR/hg38_autosomes.fa}
ILLUMINA_REF_DECODE=${ILLUMINA_REF_DECODE:-$MFM_REFERENCE_DIR/illumina_decode.fa}

# ── 1000G panel/common VCF ───────────────────────────────────────────────────
# 원격 주소는 안내/다운로드용이고 Stage 02/03의 입력은 로컬 파일이다.
# 과거 로컬 PANEL_BASE 설정은 호환한다. URL을 로컬 경로로 검사하지 않는다.
PANEL_SOURCE_URL=${PANEL_SOURCE_URL:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV}
PANEL_LOCAL_DIR=${PANEL_LOCAL_DIR:-}
if [[ -n "${PANEL_BASE:-}" ]]; then
  if [[ "$PANEL_BASE" == *://* ]]; then
    PANEL_SOURCE_URL=$PANEL_BASE
  else
    [[ -z "$PANEL_LOCAL_DIR" || "$PANEL_LOCAL_DIR" == "$PANEL_BASE" ]] || {
      echo "[FATAL] PANEL_BASE와 PANEL_LOCAL_DIR의 로컬 경로가 다릅니다" >&2; exit 1;
    }
    PANEL_LOCAL_DIR=$PANEL_BASE
  fi
fi
PANEL_LOCAL_DIR=${PANEL_LOCAL_DIR:-$MFM_REFERENCE_DIR/illumina_panel}
[[ "$PANEL_LOCAL_DIR" != *://* ]] || {
  echo "[FATAL] PANEL_LOCAL_DIR에는 URL이 아닌 로컬 디렉토리를 지정하세요" >&2; exit 1;
}
PANEL_BASE=$PANEL_LOCAL_DIR
export PANEL_SOURCE_URL PANEL_LOCAL_DIR PANEL_BASE
if [[ -z "${PANEL_FILE_TEMPLATE:-}" ]]; then
  PANEL_FILE_TEMPLATE='1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz'
fi
COMMON_VCF_DIR=${COMMON_VCF_DIR:-}
PERSAMPLE_VCF_DIR=${PERSAMPLE_VCF_DIR:-}
CREATE_COMMON_VCF=${CREATE_COMMON_VCF:-0}
VCF_JOBS=${VCF_JOBS:-8}
VCF_EXTRACTION_MODE=${VCF_EXTRACTION_MODE:-chromosome}
case "$VCF_EXTRACTION_MODE" in chromosome|concat) ;; *) echo "[FATAL] VCF_EXTRACTION_MODE는 chromosome|concat이어야 함" >&2; exit 1 ;; esac
VCF_CHROM_JOBS=${VCF_CHROM_JOBS:-6}
VCF_THREADS_PER_CHROM=${VCF_THREADS_PER_CHROM:-2}
SEQUENCE_VARIANT_POLICY=${SEQUENCE_VARIANT_POLICY:-snv_indel}
COMMON_VCF_FILTER_EXPR=${COMMON_VCF_FILTER_EXPR:-'GT="het"'}
CONSENSUS_FILTER_EXPR=${CONSENSUS_FILTER_EXPR:-'(FILTER="PASS" || FILTER=".")'}
HAPLOTAG_VCF_TYPES=${HAPLOTAG_VCF_TYPES:-snps,indels}

# ── VCF source policy ────────────────────────────────────────────────────────
# platform_matched: 과거 Clair3 self-call 재현용 opt-in 정책이다.
# manifest: manifest의 vcf_source=public을 사용한다. 현재 production 기본값이다.
# public 분기는 원본 VCF를 직접 수정하지 않으며, 검증 실패 시 Clair3로 fallback하지 않는다.
VCF_SOURCE_POLICY=${VCF_SOURCE_POLICY:-manifest}
VCF_SOURCE_POLICY=$(echo "$VCF_SOURCE_POLICY" | tr '[:upper:]' '[:lower:]' | tr '-' '_' | tr -d ' ')
case "$VCF_SOURCE_POLICY" in
  platform_matched|manifest) ;;
  *) echo "[FATAL] VCF_SOURCE_POLICY must be platform_matched|manifest: $VCF_SOURCE_POLICY" >&2; exit 1 ;;
esac
ENFORCE_VCF_SOURCE_POLICY=1

# ── tools ────────────────────────────────────────────────────────────────────
BCFTOOLS_BIN=${BCFTOOLS_BIN:-$(find_tool bcftools)}
MINIMAP2_BIN=${MINIMAP2_BIN:-$(find_tool minimap2)}
BWA_MEM2_BIN=${BWA_MEM2_BIN:-$(find_tool bwa-mem2)}
SAMTOOLS_BIN=${SAMTOOLS_BIN:-$(find_tool samtools)}
WHATSHAP_BIN=${WHATSHAP_BIN:-$(find_tool whatshap)}

# ── Clair3 env. mamba는 사용하지 않는다. ─────────────────────────────────────
CLAIR3_ENV_NAME=${CLAIR3_ENV_NAME:-clair3_v2}
if [[ -z "${CONDA_ENVS_DIR:-}" ]]; then
  if [[ -n "$CONDA_BIN" ]]; then
    CONDA_BASE=$($CONDA_BIN info --base 2>/dev/null || true)
    CONDA_ENVS_DIR=${CONDA_BASE:+$CONDA_BASE/envs}
  fi
  CONDA_ENVS_DIR=${CONDA_ENVS_DIR:-$PIPELINE_DIR/envs}
fi
CLAIR3_ENV=${CLAIR3_ENV:-$CONDA_ENVS_DIR/$CLAIR3_ENV_NAME}
CLAIR3_PYTHON_BIN=${CLAIR3_PYTHON_BIN:-$CLAIR3_ENV/bin/python}
CLAIR3_PY=${CLAIR3_PY:-$CLAIR3_ENV/bin/run_clair3.py}
CLAIR3_WHATSHAP_BIN=${CLAIR3_WHATSHAP_BIN:-$CLAIR3_ENV/bin/whatshap}
CLAIR3_SAMTOOLS_BIN=${CLAIR3_SAMTOOLS_BIN:-$CLAIR3_ENV/bin/samtools}
CLAIR3_PYPY_BIN=${CLAIR3_PYPY_BIN:-$CLAIR3_ENV/bin/pypy3}
CLAIR3_PARALLEL_BIN=${CLAIR3_PARALLEL_BIN:-$CLAIR3_ENV/bin/parallel}
CLAIR3_TABIX_BIN=${CLAIR3_TABIX_BIN:-$CLAIR3_ENV/bin/tabix}
CLAIR3_MODELS_DIR=${CLAIR3_MODELS_DIR:-$CLAIR3_ENV/bin/models}
# CLAIR3_MODEL_PATH를 직접 지정하면 전체 플랫폼 공통 override로 사용한다.
# 기본 운영은 sample platform에 따라 ONT/HiFi 모델을 자동 선택한다.
CLAIR3_MODEL_PATH_OVERRIDE=${CLAIR3_MODEL_PATH:-}
CLAIR3_MODEL_PATH_ONT=${CLAIR3_MODEL_PATH_ONT:-$CLAIR3_MODELS_DIR/ont}
CLAIR3_MODEL_PATH_HIFI=${CLAIR3_MODEL_PATH_HIFI:-$CLAIR3_MODELS_DIR/hifi}
CLAIR3_MODEL_PATH=${CLAIR3_MODEL_PATH_OVERRIDE:-$CLAIR3_MODEL_PATH_ONT}
CLAIR3_THREADS=${CLAIR3_THREADS:-32}
CLAIR3_OUTPUT_THREADS=${CLAIR3_OUTPUT_THREADS:-8}
CLAIR3_JOBS=${CLAIR3_JOBS:-1}
CLAIR3_CHUNK_SIZE=${CLAIR3_CHUNK_SIZE:-5000000}
CLAIR3_MIN_MQ=${CLAIR3_MIN_MQ:-5}
CLAIR3_MIN_COVERAGE=${CLAIR3_MIN_COVERAGE:-2}
CLAIR3_QUAL=${CLAIR3_QUAL:-2}
CLAIR3_BED_FN=${CLAIR3_BED_FN:-}
CLAIR3_DRY_RUN=${CLAIR3_DRY_RUN:-0}
CLAIR3_ALLOW_EMPTY_HAPLOTAG=${CLAIR3_ALLOW_EMPTY_HAPLOTAG:-0}
CLAIR3_ALLOW_PILEUP_FALLBACK=${CLAIR3_ALLOW_PILEUP_FALLBACK:-0}
CLAIR3_QC_INSPECT_READS=${CLAIR3_QC_INSPECT_READS:-50000}
CLAIR3_QC_MIN_HP_FRACTION=${CLAIR3_QC_MIN_HP_FRACTION:-0.50}


# ── sequence source policy ───────────────────────────────────────────────────
# SEQUENCE_SOURCE controls how HAP1/HAP2 FASTA is obtained.
# - vcf_consensus: bcftools consensus. v2 production default for both short/long reads.
# - auto: resolves to vcf_consensus for all platforms; kept for one-command compatibility.
# - denovo: opt-in research/QC mode for long-read de novo assembly; not production default.
SEQUENCE_SOURCE=${SEQUENCE_SOURCE:-vcf_consensus}
# 현재 선정 cohort는 모두 제공 VCF를 사용한다. Clair3 self-call은 과거 결과
# 재현/연구 실행에서 아래 opt-in과 sample별 vcf_source=clair3를 함께 명시할 때만 허용한다.
ENABLE_CLAIR3_SELF_CALL=${ENABLE_CLAIR3_SELF_CALL:-0}
case "$ENABLE_CLAIR3_SELF_CALL" in 0|1) ;; *) echo "[FATAL] ENABLE_CLAIR3_SELF_CALL must be 0|1" >&2; exit 1 ;; esac
# manifest를 쓰지 않는 단일 sample 실행의 호환용 override다.
# auto는 long-read에서 더 이상 Clair3로 해석하지 않고 오류 처리한다.
LONGREAD_SEQUENCE_VCF_SOURCE=${LONGREAD_SEQUENCE_VCF_SOURCE:-auto}
LONGREAD_SEQUENCE_VCF_SOURCE=$(echo "$LONGREAD_SEQUENCE_VCF_SOURCE" | tr '[:upper:]' '[:lower:]' | tr '-' '_' | tr -d ' ')
case "$LONGREAD_SEQUENCE_VCF_SOURCE" in auto|clair3|public) ;; *) echo "[FATAL] LONGREAD_SEQUENCE_VCF_SOURCE must be auto|clair3|public" >&2; exit 1 ;; esac
LONGREAD_HAPLOTAG_VCF_SOURCE=${LONGREAD_HAPLOTAG_VCF_SOURCE:-auto}
LONGREAD_HAPLOTAG_VCF_SOURCE=$(echo "$LONGREAD_HAPLOTAG_VCF_SOURCE" | tr '[:upper:]' '[:lower:]' | tr '-' '_' | tr -d ' ')
case "$LONGREAD_HAPLOTAG_VCF_SOURCE" in auto|clair3|public) ;; *) echo "[FATAL] LONGREAD_HAPLOTAG_VCF_SOURCE must be auto|clair3|public" >&2; exit 1 ;; esac
# sequence와 haplotag에는 동일한 VCF evidence를 사용한다. 한쪽만 지정한
# 기존 단일-sample 실행은 그 값을 양쪽에 전파하고, 상충하는 지정은 차단한다.
if [[ "$LONGREAD_SEQUENCE_VCF_SOURCE" == auto && "$LONGREAD_HAPLOTAG_VCF_SOURCE" != auto ]]; then
  LONGREAD_SEQUENCE_VCF_SOURCE=$LONGREAD_HAPLOTAG_VCF_SOURCE
elif [[ "$LONGREAD_SEQUENCE_VCF_SOURCE" != auto && "$LONGREAD_HAPLOTAG_VCF_SOURCE" == auto ]]; then
  LONGREAD_HAPLOTAG_VCF_SOURCE=$LONGREAD_SEQUENCE_VCF_SOURCE
elif [[ "$LONGREAD_SEQUENCE_VCF_SOURCE" != auto && "$LONGREAD_HAPLOTAG_VCF_SOURCE" != auto && "$LONGREAD_SEQUENCE_VCF_SOURCE" != "$LONGREAD_HAPLOTAG_VCF_SOURCE" ]]; then
  echo "[FATAL] long-read sequence/haplotag VCF source가 서로 다름: sequence=$LONGREAD_SEQUENCE_VCF_SOURCE haplotag=$LONGREAD_HAPLOTAG_VCF_SOURCE" >&2
  exit 1
fi
PUBLIC_VCF_PATH=${PUBLIC_VCF_PATH:-}
PUBLIC_VCF_SAMPLE_ID=${PUBLIC_VCF_SAMPLE_ID:-}
VCF_REFERENCE_FASTA=${VCF_REFERENCE_FASTA:-}
VCF_CALLSET=${VCF_CALLSET:-}
PHASE_SOURCE=${PHASE_SOURCE:-provided}
LONGREAD_SEQ_VCF_EXPORT_MODE=${LONGREAD_SEQ_VCF_EXPORT_MODE:-copy}  # copy|symlink
LONGREAD_SEQ_VCF_EXPORT_MODE=$(echo "$LONGREAD_SEQ_VCF_EXPORT_MODE" | tr '[:upper:]' '[:lower:]' | tr '-' '_' | tr -d ' ')
LONGREAD_ALLOW_UNPHASED_SEQUENCE_VCF=${LONGREAD_ALLOW_UNPHASED_SEQUENCE_VCF:-0}
if [[ "${LONGREAD_PREPARE_PANEL_VCF_FOR_QC:-0}" == 1 ]]; then
  echo "[FATAL] LONGREAD_PREPARE_PANEL_VCF_FOR_QC=1 is removed. Long-read panel VCF preparation is not part of production v2." >&2
  exit 1
fi
LONGREAD_PREPARE_PANEL_VCF_FOR_QC=0
LONGREAD_RUN_CLAIR3_FOR_QC=0
DENOVO_ASSEMBLER=${DENOVO_ASSEMBLER:-hifiasm}
HIFIASM_BIN=${HIFIASM_BIN:-$(find_tool hifiasm)}
DENOVO_THREADS=${DENOVO_THREADS:-48}
DENOVO_FASTQ_THREADS=${DENOVO_FASTQ_THREADS:-8}
DENOVO_CONTIGS=${DENOVO_CONTIGS:-${CHROMS:-chr1-chr22}}
DENOVO_LABEL=${DENOVO_LABEL:-}
DENOVO_EXCLUDE_FLAGS=${DENOVO_EXCLUDE_FLAGS:-2308}  # unmapped(4)+secondary(256)+supplementary(2048)
DENOVO_MIN_CONTIG_LEN=${DENOVO_MIN_CONTIG_LEN:-100000}
DENOVO_KEEP_FASTQ=${DENOVO_KEEP_FASTQ:-0}
DENOVO_KEEP_GFA=${DENOVO_KEEP_GFA:-1}
DENOVO_KEEP_BIN=${DENOVO_KEEP_BIN:-0}
DENOVO_OVERWRITE=${DENOVO_OVERWRITE:-0}
DENOVO_DRY_RUN=${DENOVO_DRY_RUN:-0}
ASSEMBLY_ANCHOR_REF=${ASSEMBLY_ANCHOR_REF:-$REF_SEQUENCE}
ASSEMBLY_ANCHOR_PRESET=${ASSEMBLY_ANCHOR_PRESET:-asm20}
ASSEMBLY_ANCHOR_THREADS=${ASSEMBLY_ANCHOR_THREADS:-48}
ASSEMBLY_MIN_ANCHOR_MAPQ=${ASSEMBLY_MIN_ANCHOR_MAPQ:-30}
ASSEMBLY_ANCHOR_OVERWRITE=${ASSEMBLY_ANCHOR_OVERWRITE:-0}
ASSEMBLY_ANCHOR_DRY_RUN=${ASSEMBLY_ANCHOR_DRY_RUN:-0}

# ── WGS preprocess ───────────────────────────────────────────────────────────
# 현재 Stage 01은 모든 원본 alignment/read 형식을 canonical analysis BAM으로 통일한다.
# STAGE04_* 변수명은 예전 내부 단계 번호와 기존 실행 설정의 호환을 위해 유지한다.
STAGE04_ALIGN_THREADS=${STAGE04_ALIGN_THREADS:-32}
STAGE04_SORT_THREADS=${STAGE04_SORT_THREADS:-8}
# memory_profile.sh가 선택한 값(균형 2G, high 16G)을 기본값으로 쓴다.
# 명시한 STAGE04_SORT_MEM은 그대로 보존한다.
STAGE04_SORT_MEM=${STAGE04_SORT_MEM:-$MFM_SORT_MEM_DEFAULT}
STAGE04_EXTRACT_THREADS=${STAGE04_EXTRACT_THREADS:-4}
# 전수 record/tag 통계: auto는64MiB 이상에서 native samtools를 사용함.
MFM_RECORD_STATS_BACKEND=${MFM_RECORD_STATS_BACKEND:-auto}
case "$MFM_RECORD_STATS_BACKEND" in auto|python|native) ;; *) die "MFM_RECORD_STATS_BACKEND는 auto/python/native여야 함" ;; esac
export MFM_RECORD_STATS_BACKEND
STAGE04_MIN_QUAL_VALUES=${STAGE04_MIN_QUAL_VALUES:-2}
STAGE04_OVERWRITE=${STAGE04_OVERWRITE:-0}
STAGE04_KEEP_RUN_BAMS=${STAGE04_KEEP_RUN_BAMS:-0}
# Raw Illumina는 일반 human paired-end WGS 관행에 맞춰 BWA-MEM2가 기본이다.
# minimap2는 이전 결과 재현/ablation 때만 명시적으로 선택한다. 현재 Stage 10 개인 HAP 재정렬과는 별개다.
STAGE04_ILLUMINA_ALIGNER=${STAGE04_ILLUMINA_ALIGNER:-bwa-mem2}
STAGE04_ILLUMINA_ALIGNER=$(echo "$STAGE04_ILLUMINA_ALIGNER" | tr '[:upper:]' '[:lower:]')
case "$STAGE04_ILLUMINA_ALIGNER" in
  bwa-mem2|minimap2) ;;
  *) echo "[FATAL] STAGE04_ILLUMINA_ALIGNER must be bwa-mem2|minimap2: $STAGE04_ILLUMINA_ALIGNER" >&2; exit 1 ;;
esac
# Illumina PCR duplicate는 canonical BAM 좌표에서 한 번만 mark한다. auto는
# source Picard provenance가 검증되면 보존하고, raw/unprovenanced 입력만 mark한다.
STAGE04_ILLUMINA_DUPLICATE_POLICY=${STAGE04_ILLUMINA_DUPLICATE_POLICY:-auto}
STAGE04_ILLUMINA_DUPLICATE_POLICY=$(echo "$STAGE04_ILLUMINA_DUPLICATE_POLICY" | tr '[:upper:]' '[:lower:]')
case "$STAGE04_ILLUMINA_DUPLICATE_POLICY" in
  auto|mark|preserve) ;;
  *) echo "[FATAL] STAGE04_ILLUMINA_DUPLICATE_POLICY must be auto|mark|preserve: $STAGE04_ILLUMINA_DUPLICATE_POLICY" >&2; exit 1 ;;
esac
STAGE04_MARKDUP_HEAP_GB=${STAGE04_MARKDUP_HEAP_GB:-$MFM_STAGE04_MARKDUP_HEAP_GB_DEFAULT}
[[ "$STAGE04_MARKDUP_HEAP_GB" =~ ^[1-9][0-9]*$ ]] || {
  echo "[FATAL] STAGE04_MARKDUP_HEAP_GB must be a positive integer: $STAGE04_MARKDUP_HEAP_GB" >&2; exit 1;
}
# BWA-MEM2 index는 sample work dir가 아닌 reference fingerprint별 공유 cache에 한 번만 생성한다.
BWA_MEM2_INDEX_ROOT=${BWA_MEM2_INDEX_ROOT:-$DATASET_ROOT/_reference_cache/bwa_mem2}
CHROMS=${CHROMS:-chr1-chr22}
WINDOW_SIZE=${WINDOW_SIZE:-1048576}
OFFSET_UNIT=${OFFSET_UNIT:-}
FIXED_WINDOW_START=${FIXED_WINDOW_START:-0}
WINDOW_SEED=${WINDOW_SEED:-42}
MIN_MAPQ=${MIN_MAPQ:-30}
WINDOW_JOBS=${WINDOW_JOBS:-2}
WINDOW_DROP_N=${WINDOW_DROP_N:-0}
case "$WINDOW_DROP_N" in 0|1) ;; *) echo "[FATAL] N-window: WINDOW_DROP_N은 0 또는 1이어야 함" >&2; exit 1 ;; esac
# 내부 N은 보존하고 실제 구간이 전부 N인 창만 제외한다. 5%는 학습 선택용이다.
WINDOW_STORAGE_POLICY=${WINDOW_STORAGE_POLICY:-annotated_acgtn_v1}
# N 원값만 보존함. 학습용 N 비율 기준은 소비자에서 결정함.
export WINDOW_STORAGE_POLICY WINDOW_DROP_N WINDOW_MAX_BLACKLIST_FRAC
BLACKLIST_BED=${BLACKLIST_BED:-$MFM_REFERENCE_DIR/hg38_autosomes_blacklist_union.bed}
WINDOW_MAX_BLACKLIST_FRAC=${WINDOW_MAX_BLACKLIST_FRAC:-1.0}
# config는 run_pipeline→phased→runner→sample처럼 여러 process에서 다시 source된다.
# 첫 source의 사용자 요청을 별도로 export해, 부모 fallback(예: align=6)이 child에서
# 명시값으로 오인되어 CPU 예산 재계산을 막지 않게 한다.
if [[ -z "${MFM_WGS_ALIGN_REQUESTED+x}" ]]; then MFM_WGS_ALIGN_REQUESTED=${WGS_REALIGN_THREADS:-auto}; fi
if [[ -z "${MFM_WGS_SORT_REQUESTED+x}" ]]; then MFM_WGS_SORT_REQUESTED=${WGS_SORT_THREADS:-auto}; fi
if [[ -z "${MFM_WGS_SORT_MEM_REQUESTED+x}" ]]; then MFM_WGS_SORT_MEM_REQUESTED=${WGS_SORT_MEM:-auto}; fi
if [[ -z "${MFM_WGS_HAP_REQUESTED+x}" ]]; then MFM_WGS_HAP_REQUESTED=${WGS_HAP_JOBS:-auto}; fi
if [[ -z "${MFM_WGS_ROUTE_REQUESTED+x}" ]]; then MFM_WGS_ROUTE_REQUESTED=${WGS_ROUTE_JOBS:-auto}; fi
case "$MFM_WGS_ALIGN_REQUESTED" in
  auto) MFM_WGS_ALIGN_EXPLICIT=0; WGS_REALIGN_THREADS=6 ;;
  *) MFM_WGS_ALIGN_EXPLICIT=1; WGS_REALIGN_THREADS=$MFM_WGS_ALIGN_REQUESTED ;;
esac
case "$MFM_WGS_SORT_REQUESTED" in
  auto) MFM_WGS_SORT_EXPLICIT=0; WGS_SORT_THREADS=2 ;;
  *) MFM_WGS_SORT_EXPLICIT=1; WGS_SORT_THREADS=$MFM_WGS_SORT_REQUESTED ;;
esac
case "$MFM_WGS_SORT_MEM_REQUESTED" in
  auto) WGS_SORT_MEM=$MFM_SORT_MEM_DEFAULT ;;
  *) WGS_SORT_MEM=$MFM_WGS_SORT_MEM_REQUESTED ;;
esac
WGS_REALIGN_JOBS=${WGS_REALIGN_JOBS:-1}
# sample 내부 HAP1/HAP2 병렬도. auto는 CPU 예산 16 이상이면 기본 2이며,
# runner가 자원 부족을 명시한 경우 또는 CPU 예산 미달이면 1로 fallback한다.
WGS_HAP_JOBS=$MFM_WGS_HAP_REQUESTED
# all_to_both Illumina source routing은 indexed contig별 process로 병렬화한다.
# HAP 정렬과 동시에 실행되지 않으므로 sample CPU 예산 전체를 사용할 수 있다.
WGS_ROUTE_JOBS=$MFM_WGS_ROUTE_REQUESTED
export MFM_WGS_ALIGN_REQUESTED MFM_WGS_SORT_REQUESTED MFM_WGS_SORT_MEM_REQUESTED MFM_WGS_HAP_REQUESTED MFM_WGS_ROUTE_REQUESTED
WGS_AUTO_THREADS=${WGS_AUTO_THREADS:-1}
# 과거 production 운영의 0.90 참고선이다. record_only에서는 비교·기록에만 쓰고,
# strict를 명시한 경우에만 Stage 10 중단 기준으로 쓴다.
WGS_MIN_MAPPED_FRAC=${WGS_MIN_MAPPED_FRAC:-0.90}
awk -v value="$WGS_MIN_MAPPED_FRAC" 'BEGIN { exit !(value ~ /^[0-9]+([.][0-9]+)?$/ && value >= 0 && value <= 1) }' || {
  echo "[FATAL] WGS_MIN_MAPPED_FRAC must be a number in [0,1]: $WGS_MIN_MAPPED_FRAC" >&2
  exit 1
}
# 비율 QC는 기본적으로 표본 탈락/필터에 쓰지 않고 값만 기록한다.
# 과거 0.90 hard-fail 재현이 필요한 경우에만 strict를 명시한다.
QC_RATIO_POLICY=${QC_RATIO_POLICY:-record_only}
case "$QC_RATIO_POLICY" in
  record_only|strict) ;;
  *) echo "[FATAL] QC_RATIO_POLICY must be record_only|strict: $QC_RATIO_POLICY" >&2; exit 1 ;;
esac
# 현재 Stage 10(legacy 내부 Stage22): short-read personal-HAP realignment defaults to BWA-MEM2.
# Long-read platforms continue to use the platform-specific minimap2 preset.
STAGE22_ILLUMINA_ALIGNER=${STAGE22_ILLUMINA_ALIGNER:-bwa-mem2}
STAGE22_ILLUMINA_ALIGNER=$(echo "$STAGE22_ILLUMINA_ALIGNER" | tr '[:upper:]' '[:lower:]')
case "$STAGE22_ILLUMINA_ALIGNER" in
  bwa-mem2|minimap2) ;;
  *) echo "[FATAL] STAGE22_ILLUMINA_ALIGNER must be bwa-mem2|minimap2: $STAGE22_ILLUMINA_ALIGNER" >&2; exit 1 ;;
esac
STAGE22_BWA_INDEX_ROOT=${STAGE22_BWA_INDEX_ROOT:-$BWA_MEM2_INDEX_ROOT/personal_hap}
# 개인 HAP index는 대부분 sample마다 고유하므로 기본은 scratch에서 만들고 성공 후 제거한다.
# shared는 반복 실험에서 같은 HAP index를 재사용할 때만 명시적으로 사용한다.
STAGE11_BWA_INDEX_MODE=${STAGE11_BWA_INDEX_MODE:-ephemeral}
case "$STAGE11_BWA_INDEX_MODE" in
  ephemeral|shared) ;;
  *) echo "[FATAL] STAGE11_BWA_INDEX_MODE must be ephemeral|shared: $STAGE11_BWA_INDEX_MODE" >&2; exit 1 ;;
esac
STAGE22_ILLUMINA_DUPLICATE_MODE=${STAGE22_ILLUMINA_DUPLICATE_MODE:-preserve_source}
case "$STAGE22_ILLUMINA_DUPLICATE_MODE" in
  preserve_source|none) ;;
  *) echo "[FATAL] STAGE22_ILLUMINA_DUPLICATE_MODE must be preserve_source|none: $STAGE22_ILLUMINA_DUPLICATE_MODE" >&2; exit 1 ;;
esac
WGS_REALIGN_FORMAT=${WGS_REALIGN_FORMAT:-cram}
if [[ "$WGS_REALIGN_FORMAT" != cram ]]; then
  echo "[FATAL] WGS_REALIGN_FORMAT must be cram: Stage 12(legacy Stage70)가 personal target.fa와 CRAM @SQ M5를 검증합니다." >&2
  exit 1
fi
# 현재 Stage 10의 high-level routing mode가 legacy Stage22 저수준 변수보다 우선한다.
# all_to_both: coverage/evidence; hp_routed_with_unknown: retain HP0 twice;
# hap_specific: retain only HP-consistent read/fragment on one HAP.
HAP_READ_MODE=${HAP_READ_MODE:-all_to_both}
case "$HAP_READ_MODE" in
  all_to_both|hp_routed_with_unknown|hap_specific) ;;
  *) echo "[FATAL] HAP_READ_MODE must be all_to_both|hp_routed_with_unknown|hap_specific: $HAP_READ_MODE" >&2; exit 1 ;;
esac
REALIGNED_SUBDIR=${REALIGNED_SUBDIR:-"realigned/$HAP_READ_MODE"}
COMPACT_LAYOUT_BASE=${COMPACT_LAYOUT_BASE:-full_shared_hap_4bit_acgtnp_depth}
COMPACT_LAYOUT_NAME=${COMPACT_LAYOUT_NAME:-"${COMPACT_LAYOUT_BASE}__${HAP_READ_MODE}"}
STAGE12_TRACE_POLICY=${STAGE12_TRACE_POLICY:-required}
case "$STAGE12_TRACE_POLICY" in
  required|deferred|off) ;;
  *) echo "[FATAL] STAGE12_TRACE_POLICY must be required|deferred|off: $STAGE12_TRACE_POLICY" >&2; exit 1 ;;
esac
ALLOW_CLEANUP_WITHOUT_TRACE=${ALLOW_CLEANUP_WITHOUT_TRACE:-0}
READ_ROUTING_POLICY=${READ_ROUTING_POLICY:-all_to_both}
UNKNOWN_READ_POLICY=${UNKNOWN_READ_POLICY:-include_both}
if [[ -z "${MFM_COMPACT_JOBS_REQUESTED+x}" ]]; then MFM_COMPACT_JOBS_REQUESTED=${COMPACT_JOBS:-auto}; fi
COMPACT_JOBS=$MFM_COMPACT_JOBS_REQUESTED
export MFM_COMPACT_JOBS_REQUESTED
# trace는 compact writer와 독립적인 읽기 전용 검증이다. contig×HAP의 평면
# task(기본 chr1–22 × HAP1/2 = 최대 44개)만 병렬화하며, 출력 Zarr나 filter
# policy를 바꾸지 않는다. sample CPU 예산이 있으면 trace는 compact와 별도로
# 그 예산 안에서 자원을 배분한다.
if [[ -z "${MFM_TRACE_JOBS_REQUESTED+x}" ]]; then MFM_TRACE_JOBS_REQUESTED=${TRACE_JOBS:-auto}; fi
TRACE_JOBS=$MFM_TRACE_JOBS_REQUESTED
export MFM_TRACE_JOBS_REQUESTED
MFM_TRACE_MAX_JOBS=${MFM_TRACE_MAX_JOBS:-44}
MFM_TRACE_CPU_RESERVE=${MFM_TRACE_CPU_RESERVE:-4}
if [[ -z "${MFM_TRACE_BLOSC_THREADS_REQUESTED+x}" ]]; then
  MFM_TRACE_BLOSC_THREADS_REQUESTED=${MFM_TRACE_BLOSC_THREADS:-auto}
fi
MFM_TRACE_BLOSC_THREADS=$MFM_TRACE_BLOSC_THREADS_REQUESTED
MFM_TRACE_PREFETCH_RECORDS=${MFM_TRACE_PREFETCH_RECORDS:-65536}
MFM_TRACE_SQLITE_BATCH=${MFM_TRACE_SQLITE_BATCH:-8192}
# 임시 SQLite page cache는 worker당64MiB. DB 전체를 메모리에 올리지 않는다.
MFM_TRACE_SQLITE_CACHE_MIB=${MFM_TRACE_SQLITE_CACHE_MIB:-64}
MFM_TRACE_MULTISET_BACKEND=${MFM_TRACE_MULTISET_BACKEND:-auto}
case "$MFM_TRACE_MULTISET_BACKEND" in auto|partitioned|sqlite) ;; *) echo "[FATAL] MFM_TRACE_MULTISET_BACKEND는 auto|partitioned|sqlite이어야 함" >&2; exit 1 ;; esac
MFM_TRACE_SORT_LEAF_MIB=${MFM_TRACE_SORT_LEAF_MIB:-64}
MFM_TRACE_DECODED_CACHE_MIB=${MFM_TRACE_DECODED_CACHE_MIB:-64}
MFM_TRACE_PREFETCH_PACKED_MIB=${MFM_TRACE_PREFETCH_PACKED_MIB:-16}
export MFM_TRACE_MULTISET_BACKEND MFM_TRACE_SORT_LEAF_MIB MFM_TRACE_DECODED_CACHE_MIB MFM_TRACE_PREFETCH_PACKED_MIB
MFM_VALIDATION_RECORD_CHUNK=${MFM_VALIDATION_RECORD_CHUNK:-262144}
MFM_VALIDATION_PACKED_BYTE_CHUNK=${MFM_VALIDATION_PACKED_BYTE_CHUNK:-16777216}
MFM_WGS_CPU_BALANCE=${MFM_WGS_CPU_BALANCE:-auto}
export MFM_WGS_CPU_BALANCE
MFM_WINDOW_CONTIG_JOBS=${MFM_WINDOW_CONTIG_JOBS:-auto}
MFM_BLOSC_THREADS_REQUESTED=${MFM_BLOSC_THREADS_REQUESTED:-${MFM_BLOSC_THREADS:-auto}}
export MFM_WINDOW_CONTIG_JOBS MFM_BLOSC_THREADS_REQUESTED
MFM_COMPACT_RESERVATION_MODEL=${MFM_COMPACT_RESERVATION_MODEL:-input_v2}
export MFM_COMPACT_RESERVATION_MODEL
MFM_VALIDATION_CACHE_CHUNKS=${MFM_VALIDATION_CACHE_CHUNKS:-4}
export MFM_TRACE_MAX_JOBS MFM_TRACE_CPU_RESERVE MFM_TRACE_BLOSC_THREADS_REQUESTED
export MFM_TRACE_PREFETCH_RECORDS MFM_TRACE_SQLITE_BATCH MFM_VALIDATION_RECORD_CHUNK MFM_VALIDATION_PACKED_BYTE_CHUNK MFM_VALIDATION_CACHE_CHUNKS
export MFM_TRACE_SQLITE_CACHE_MIB
WINDOWS_PER_CHUNK=${WINDOWS_PER_CHUNK:-2}
COMPACT_COMPRESSOR=${COMPACT_COMPRESSOR:-zstd}
COMPACT_COMPRESSION_LEVEL=${COMPACT_COMPRESSION_LEVEL:-5}
COMPACT_VALIDATION_LEVEL=${COMPACT_VALIDATION_LEVEL:-full}
COMPACT_MAX_WINDOWS=${COMPACT_MAX_WINDOWS:-}
COMPACT_ALLOW_UNPHASED_LONGREAD=${COMPACT_ALLOW_UNPHASED_LONGREAD:-0}
# original은 중간 산출물을 기본 보존한다. 기존 cleanup은 명시적 opt-in에서 사용한다.
CLEAN_INTERMEDIATES_ON_SUCCESS=${CLEAN_INTERMEDIATES_ON_SUCCESS:-0}
CLEAN_INTERMEDIATES_DRY_RUN=${CLEAN_INTERMEDIATES_DRY_RUN:-0}
CLEAN_AFTER_RNA_ONLY=${CLEAN_AFTER_RNA_ONLY:-1}
CLEAN_SHARED_ASSETS_ON_SUCCESS=${CLEAN_SHARED_ASSETS_ON_SUCCESS:-0}

log() { echo "[$(date '+%F %T')] $*" >&2; }
warn() { echo "[$(date '+%F %T')] [WARN] $*" >&2; }
die() { echo "[$(date '+%F %T')] [FATAL] $*" >&2; exit 1; }

manifest_field() {
  local sample=$1 col=$2
  [[ -n "$PIPELINE_MANIFEST" && -s "$PIPELINE_MANIFEST" ]] || return 3
  awk -F'\t' -v S="$sample" -v C="$col" '
    NR==1 && $1=="sample_id" { next }
    $1==S { print $C; found=1; exit }
    END { if (!found) exit 3 }
  ' "$PIPELINE_MANIFEST"
}

manifest_named_field() {
  local sample=$1 name=$2
  [[ -n "$PIPELINE_MANIFEST" && -s "$PIPELINE_MANIFEST" ]] || return 3
  awk -F'\t' -v S="$sample" -v N="$name" '
    NR==1 { for (i=1; i<=NF; i++) if ($i==N) col=i; next }
    $1==S && col { print $col; found=1; exit }
    END { if (!found) exit 3 }
  ' "$PIPELINE_MANIFEST"
}

manifest_has_sample() {
  local sample=$1
  [[ -n "${PIPELINE_MANIFEST:-}" && -s "$PIPELINE_MANIFEST" ]] || return 1
  awk -F'\t' -v S="$sample" 'NR>1 && $1==S {found=1; exit} END {exit(found ? 0 : 1)}' "$PIPELINE_MANIFEST"
}

field_or_default() {
  local sample=$1 col=$2 default=$3 value
  value=$(manifest_field "$sample" "$col" 2>/dev/null || true)
  if [[ -n "$value" && "$value" != "-" ]]; then
    echo "$value"
  else
    echo "$default"
  fi
}

named_field_or_default() {
  local sample=$1 name=$2 default=$3 value
  value=$(manifest_named_field "$sample" "$name" 2>/dev/null || true)
  if [[ -n "$value" && "$value" != "-" ]]; then
    printf '%s\n' "$value"
  else
    printf '%s\n' "$default"
  fi
}

resolve_input_path() {
  local value=${1:-} manifest_dir
  [[ -n "$value" && "$value" != "-" ]] || return 0
  if [[ "$value" == /* ]]; then
    printf '%s\n' "$value"
  elif [[ -n "${INPUT_ROOT:-}" ]]; then
    readlink -m "$INPUT_ROOT/$value"
  elif [[ -n "${PIPELINE_MANIFEST:-}" ]]; then
    manifest_dir=$(cd "$(dirname "$PIPELINE_MANIFEST")" && pwd)
    readlink -m "$manifest_dir/$value"
  else
    readlink -m "$value"
  fi
}

manifest_positional_path() {
  local value
  value=$(field_or_default "$1" "$2" "$3")
  resolve_input_path "$value"
}

named_path_or_default() {
  local sample=$1 name=$2 default=$3 value
  value=$(named_field_or_default "$sample" "$name" "$default")
  [[ -n "$value" && "$value" != "-" ]] || return 0
  if [[ "$value" == /* ]]; then
    printf '%s\n' "$value"
  else
    resolve_input_path "$value"
  fi
}

sample_platform() { field_or_default "$1" 2 "$DEFAULT_PLATFORM"; }
sample_read_type() { field_or_default "$1" 3 "$DEFAULT_READ_TYPE"; }

# 서열전용은 short-read에만 적용하며 long-read의 기존 경로를 유지함.
sample_sequence_only() {
  local mode=${MFM_SHORT_READ_OUTPUT_MODE:-full}
  [[ "$(sample_read_type "$1")" != long_read ]] || mode=${MFM_LONG_READ_OUTPUT_MODE:-full}
  [[ "$mode" == sequence_only ]]
}
sample_output_mode() {
  local mode=${MFM_SHORT_READ_OUTPUT_MODE:-full}
  [[ "$(sample_read_type "$1")" != long_read ]] || mode=${MFM_LONG_READ_OUTPUT_MODE:-full}
  printf '%s\n' "$mode"
}
sample_pileup_feature() { [[ "$(sample_output_mode "$1")" == pileup_feature ]]; }
sample_hap_realign_required() {
  # 실측에서만 HAP 재정렬 비용을 추가함. reference pileup 입력은 변경하지 않음.
  if [[ "${MFM_PILEUP_BENCHMARK_REALIGN:-0}" == 1 ]] && sample_pileup_feature "$1"; then return 0; fi
  [[ "$(sample_output_mode "$1")" == full ]] || { [[ -n "$(sample_rna "$1")" && "$(sample_rna "$1")" != - ]] && ! sample_rna_deferred "$1" && [[ "${MFM_RNA_OUTPUT_MODE:-full}" != coverage_only ]]; }
}
sample_rna_deferred() {
  [[ "$(sample_output_mode "$1")" != full ]] && [[ "${MFM_SEQUENCE_DEFER_RNA:-0}" == 1 ]] && [[ -n "$(sample_rna "$1")" && "$(sample_rna "$1")" != - ]]
}
sample_wgs_reads_required() {
  ! sample_sequence_only "$1" || { [[ -n "$(sample_rna "$1")" && "$(sample_rna "$1")" != - ]] && ! sample_rna_deferred "$1" && [[ "${MFM_RNA_OUTPUT_MODE:-full}" != coverage_only ]]; }
}
apply_sample_payload() {
  local sample=$1
  case "${MFM_SEQUENCE_DEFER_RNA:-0}" in 0|1) ;; *) die "invalid MFM_SEQUENCE_DEFER_RNA" ;; esac
  case "${MFM_SHORT_READ_OUTPUT_MODE:-full}" in full|sequence_only|pileup_feature) ;; *) die "invalid MFM_SHORT_READ_OUTPUT_MODE" ;; esac
  case "${MFM_LONG_READ_OUTPUT_MODE:-full}" in full|sequence_only|pileup_feature) ;; *) die "invalid MFM_LONG_READ_OUTPUT_MODE" ;; esac
  case "${MFM_RNA_OUTPUT_MODE:-full}" in full|coverage_only) ;; *) die "invalid MFM_RNA_OUTPUT_MODE" ;; esac
  case "${MFM_PILEUP_BENCHMARK_REALIGN:-0}" in 0|1) ;; *) die "invalid MFM_PILEUP_BENCHMARK_REALIGN" ;; esac
  if [[ "${MFM_RNA_OUTPUT_MODE:-full}" == coverage_only ]] && sample_rna_deferred "$sample"; then die "RNA coverage_only와RNA보류는동시설정할수없음"; fi
  if sample_pileup_feature "$sample"; then
    [[ "$(sequence_source_for_sample "$sample")" == vcf_consensus ]] || die "pileup_feature requires vcf_consensus chain"
    export MFM_WGS_PAYLOAD=pileup_feature
    export COMPACT_LAYOUT_NAME=hap_sequence_reference_pileup_v1
  elif sample_sequence_only "$sample"; then
    [[ "$(sequence_source_for_sample "$sample")" == vcf_consensus ]] || die "sequence_only requires vcf_consensus"
    local expected_vcf=illumina_panel
    [[ "$(sample_read_type "$sample")" != long_read ]] || expected_vcf=public
    [[ "$(sample_vcf_source "$sample")" == "$expected_vcf" && "$(sample_phase_source "$sample")" == provided ]] || die "sequence_only requires provided phased VCF: $expected_vcf"
    if [[ -n "$(sample_rna "$sample")" && "$(sample_rna "$sample")" != - ]] && ! sample_rna_deferred "$sample" && [[ "${MFM_RNA_OUTPUT_MODE:-full}" != coverage_only ]]; then
      [[ "${MFM_SEQUENCE_RNA_WGS_EVIDENCE:-}" == temporary ]] || die "RNA ASE/CNV 유지에는 WGS 임시근거 정책 확인이 필요함"
    fi
    export MFM_WGS_PAYLOAD=sequence_only
    export COMPACT_LAYOUT_NAME=sequence_4bit_acgtnp__windowed
  fi
}

sample_dataset_role() { field_or_default "$1" 4 "$DEFAULT_DATASET_ROLE"; }
sample_project_id() { field_or_default "$1" 5 "$DEFAULT_PROJECT_ID"; }
sample_cohort() { field_or_default "$1" 6 "$DEFAULT_COHORT"; }
sample_source_alignment() {
  local sample=$1 value
  value=$(manifest_field "$sample" 7 2>/dev/null || true)
  [[ -n "$value" && "$value" != "-" ]] || value=${ALIGNMENT_PATH:-}
  resolve_input_path "$value"
}
sample_read_files_manifest() {
  local sample=$1 value manifest_dir
  value=$(named_field_or_default "$sample" read_files_manifest "${READ_FILES_MANIFEST:-}")
  [[ -n "$value" && "$value" != "-" ]] || return 0
  if [[ "$value" == /* ]]; then
    printf '%s\n' "$value"
  elif [[ -n "${PIPELINE_MANIFEST:-}" ]]; then
    manifest_dir=$(cd "$(dirname "$PIPELINE_MANIFEST")" && pwd)
    readlink -m "$manifest_dir/$value"
  else
    readlink -m "$value"
  fi
}
sample_reusable_file_alignment() {
  local sample=$1 file_manifest path manifest_dir
  file_manifest=$(sample_read_files_manifest "$sample")
  [[ -n "$file_manifest" && -s "$file_manifest" ]] || return 0
  path=$(awk -F'\t' -v S="$sample" '
    NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
    $h["sample_id"]==S && $h["file_role"]=="alignment" && $h["alignment_state"]=="aligned" {
      n++; value=$h["path"]
    }
    END {if(n==1) print value}
  ' "$file_manifest")
  [[ -n "$path" && "$path" != "-" ]] || return 0
  if [[ "$path" == /* ]]; then
    printf '%s\n' "$path"
  elif [[ -n "${INPUT_ROOT:-}" ]]; then
    readlink -m "$INPUT_ROOT/$path"
  else
    manifest_dir=$(cd "$(dirname "$file_manifest")" && pwd)
    readlink -m "$manifest_dir/$path"
  fi
}
sample_prepared_alignment() {
  # Stage 10 이후는 원본 포맷을 보지 않는다. Stage 04 canonical BAM만 소비한다.
  local sample=$1 explicit work_root
  explicit=$(named_path_or_default "$sample" prepared_alignment_path "${PREPARED_ALIGNMENT_PATH:-}")
  if [[ -n "$explicit" ]]; then
    [[ "$explicit" == *.bam ]] || die "$sample prepared_alignment_path must be BAM: $explicit"
    printf '%s\n' "$explicit"
    return 0
  fi
  work_root=$(sample_work_root "$sample")
  printf '%s\n' "$work_root/prepared_input/$sample.analysis.bam"
}
# Stage 10 이후 소비자는 Stage 04가 검증·표준화한 BAM만 본다.
sample_alignment() { sample_prepared_alignment "$1"; }
sample_raw_root() { manifest_positional_path "$1" 8 "${RAW_ROOT:-}"; }
sample_rna() { manifest_positional_path "$1" 9 "${RNA_PATH:-}"; }
sample_data_group() {
  local sample=$1 value rna
  value=$(named_field_or_default "$sample" data_group "")
  value=$(echo "$value" | tr '[:upper:]' '[:lower:]' | tr '-' '_' | tr -d ' ')
  if [[ -z "$value" ]]; then
    rna=$(sample_rna "$sample" 2>/dev/null || true)
    [[ -n "$rna" && "$rna" != "-" ]] && value=wgs_rna_pair || value=wgs
  fi
  case "$value" in
    wgs|wgs_only) printf 'wgs\n' ;;
    wgs_rna_pair|rna_pair|pair) printf 'wgs_rna_pair\n' ;;
    *) die "invalid data_group: sample=$sample value=$value" ;;
  esac
}
sample_final_root() {
  local sample=$1 project group explicit
  project=$(sample_project_id "$sample")
  group=$(sample_data_group "$sample")
  if [[ "$OUTPUT_ROOT_OVERRIDE" == 1 ]]; then
    printf '%s/%s/%s\n' "$OUTPUT_ROOT" "$group" "$project"
    return 0
  fi
  explicit=$(field_or_default "$sample" 10 "${FINAL_ROOT:-}")
  if [[ -n "$explicit" && "$explicit" != "-" ]]; then
    printf '%s\n' "$explicit"
  else
    printf '%s/%s/%s\n' "$OUTPUT_ROOT" "$group" "$project"
  fi
}
sample_work_root() {
  local sample=$1 project group explicit
  project=$(sample_project_id "$sample")
  group=$(sample_data_group "$sample")
  if [[ "$OUTPUT_ROOT_OVERRIDE" == 1 ]]; then
    printf '%s/_work/%s/%s/%s/%s\n' "$OUTPUT_ROOT" "$RUN_ID" "$group" "$project" "$sample"
    return 0
  fi
  explicit=$(field_or_default "$sample" 11 "${WORK_ROOT:-}")
  if [[ -n "$explicit" && "$explicit" != "-" ]]; then
    printf '%s\n' "$explicit"
  else
    printf '%s/_work/%s/%s/%s/%s\n' "$OUTPUT_ROOT" "$RUN_ID" "$group" "$project" "$sample"
  fi
}
sample_status_root() {
  local sample=$1 project group
  project=$(sample_project_id "$sample")
  group=$(sample_data_group "$sample")
  printf '%s/_run_status/%s/%s/%s/%s\n' "$OUTPUT_ROOT" "$RUN_ID" "$group" "$project" "$sample"
}
runner_status_root() {
  printf '%s/_run_status/%s/runner/%s\n' "$OUTPUT_ROOT" "$RUN_ID" "$1"
}
sample_stage04_qc() {
  local sample=$1
  printf '%s/prepared_input/%s.stage04.json\n' "$(sample_work_root "$sample")" "$sample"
}
sample_common_vcf_dir() {
  local sample=$1 cohort fallback value
  cohort=$(sample_cohort "$sample")
  fallback="$OUTPUT_ROOT/_work/$RUN_ID/_shared/$cohort/common_vcf"
  value=$(field_or_default "$sample" 12 "${COMMON_VCF_DIR:-$fallback}")
  resolve_input_path "$value"
}
sample_persample_vcf_dir() {
  local sample=$1 cohort fallback
  cohort=$(sample_cohort "$sample")
  fallback="$OUTPUT_ROOT/_work/$RUN_ID/_shared/$cohort/persample_vcf"
  field_or_default "$sample" 13 "${PERSAMPLE_VCF_DIR:-$fallback}"
}
sample_rna_bigwigs() { field_or_default "$1" 14 "${RNA_BIGWIGS:-}"; }

normalize_platform() {
  local p
  p=$(echo "$1" | tr '[:upper:]' '[:lower:]' | tr '-' '_' | tr -d ' ')
  case "$p" in
    ont|nanopore|oxfordnanopore|oxford_nanopore) echo ont ;;
    hifi|pacbio_hifi|pacbiohifi|ccs|pacbio_ccs|pacbioccs) echo hifi ;;
    clr|pacbio_clr|pacbioclr|pb_clr) echo clr ;;
    illumina|ilmn|short_read|shortread) echo illumina ;;
    *) echo "$p" ;;
  esac
}

sample_vcf_source() {
  local sample=$1 value is_long=0
  is_long_read_sample "$sample" && is_long=1
  value=$(named_field_or_default "$sample" vcf_source "")
  if [[ -z "$value" ]]; then
    if [[ "$is_long" == 1 ]]; then
      [[ -z "${PIPELINE_MANIFEST:-}" ]] || die "long-read manifest에는 vcf_source=public이 필요함: sample=$sample"
      value=$LONGREAD_SEQUENCE_VCF_SOURCE
    else
      value=illumina_panel
    fi
  fi
  value=$(echo "$value" | tr '[:upper:]' '[:lower:]' | tr '-' '_' | tr -d ' ')
  if [[ -z "$value" || "$value" == auto ]]; then
    [[ "$is_long" == 0 ]] || die "long-read vcf_source=auto는 금지됨. 제공 VCF는 public을 명시: sample=$sample"
    value=illumina_panel
  fi
  case "$value" in
    clair3|public|illumina_panel) ;;
    *) die "지원하지 않는 vcf_source: sample=$sample value=$value" ;;
  esac
  if [[ "$is_long" == 1 && "$value" == public && "$VCF_SOURCE_POLICY" != manifest ]]; then
    die "vcf_source=public은 VCF_SOURCE_POLICY=manifest 명시가 필요함: sample=$sample"
  fi
  if [[ "$is_long" == 1 && "$value" == clair3 && "$ENABLE_CLAIR3_SELF_CALL" != 1 ]]; then
    die "Clair3 self-call은 기본 비활성임: sample=$sample (명시적 재현 실행은 ENABLE_CLAIR3_SELF_CALL=1 필요)"
  fi
  if [[ "$is_long" == 1 && "$value" == illumina_panel ]]; then
    die "long-read vcf_source=illumina_panel은 지원하지 않음: sample=$sample"
  fi
  if [[ "$is_long" == 0 && "$value" != illumina_panel ]]; then
    die "short-read는 vcf_source=illumina_panel만 지원함: sample=$sample value=$value"
  fi
  printf '%s\n' "$value"
}

sample_public_vcfs() {
  local sample=$1 raw item
  raw=$(named_field_or_default "$sample" public_vcf_path "$PUBLIC_VCF_PATH")
  [[ -n "$raw" && "$raw" != "-" ]] || return 0
  IFS=';' read -r -a items <<< "$raw"
  for item in "${items[@]}"; do
    item=$(echo "$item" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    [[ -n "$item" ]] || continue
    resolve_input_path "$item"
  done
}
sample_public_vcf() { sample_public_vcfs "$1" | head -1; }
sample_public_vcf_sample_id() { named_field_or_default "$1" public_vcf_sample_id "${PUBLIC_VCF_SAMPLE_ID:-$1}"; }
# 모델 산출물 ID와 source BAM의 @RG SM이 다른 공개 alias를 명시적으로 보존한다.
# public VCF의 donor sample ID는 sample_public_vcf_sample_id로 별도 관리한다.
# manifest 외 단일 샘플 실행에서도 source BAM의 @RG SM alias를 명시할 수 있다.
sample_source_sample_id() { named_field_or_default "$1" source_sample_id "${SOURCE_SAMPLE_ID:-$1}"; }
sample_vcf_reference() { named_path_or_default "$1" vcf_reference_fasta "$VCF_REFERENCE_FASTA"; }
sample_vcf_callset() { named_field_or_default "$1" vcf_callset "${VCF_CALLSET:-unknown}"; }
sample_phase_source() {
  local value
  value=$(named_field_or_default "$1" phase_source "$PHASE_SOURCE")
  value=$(echo "$value" | tr '[:upper:]' '[:lower:]' | tr '-' '_' | tr -d ' ')
  case "$value" in provided) ;; *) die "현재 phase_source는 provided만 지원함: sample=$1 value=$value" ;; esac
  printf '%s\n' "$value"
}

sample_source_decode_reference() {
  local sample=$1 platform explicit
  explicit=$(named_path_or_default "$sample" source_decode_reference_fasta "${SOURCE_DECODE_REFERENCE_FASTA:-}")
  if [[ -n "$explicit" ]]; then
    printf '%s\n' "$explicit"
    return 0
  fi
  platform=$(normalize_platform "$(sample_platform "$sample")")
  if [[ "$platform" == illumina ]]; then
    printf '%s\n' "$ILLUMINA_REF_DECODE"
  else
    printf '%s\n' "$REF_DECODE"
  fi
}

sample_analysis_reference() {
  local sample=$1 platform source explicit public_ref
  explicit=$(named_path_or_default "$sample" analysis_reference_fasta "${ANALYSIS_REFERENCE_FASTA:-}")
  if [[ -n "$explicit" ]]; then
    printf '%s\n' "$explicit"
    return 0
  fi
  platform=$(normalize_platform "$(sample_platform "$sample")")
  source=$(sample_vcf_source "$sample")
  if [[ "$source" == public ]]; then
    public_ref=$(sample_vcf_reference "$sample")
    [[ -n "$public_ref" ]] || die "public VCF reference가 비어 있음: sample=$sample"
    printf '%s\n' "$public_ref"
    return 0
  fi
  if [[ "$platform" == illumina ]]; then
    printf '%s\n' "$ILLUMINA_REF_DECODE"
  else
    printf '%s\n' "$REF_DECODE"
  fi
}

# 기존 단계의 함수 이름은 분석 reference 의미로 유지한다.
sample_decode_reference() { sample_analysis_reference "$1"; }

# canonical BAM의 실제 정렬 좌표를 만든 reference를 선택한다. uBAM/FASTQ를
# 새로 정렬한 BAM은 analysis reference, 기존 정렬을 단순 변환한 BAM은 source
# decode reference다. 재사용 시에도 Stage 01의 alignment_policy를 보존한다.
sample_alignment_reference() {
  local sample=$1 qc mode
  qc=$(sample_stage04_qc "$sample")
  mode=source
  if [[ -s "$qc" ]]; then
    mode=$("$PYTHON" - "$qc" <<'PY'
import json
import sys
data = json.load(open(sys.argv[1]))
mode = data.get("alignment_policy", {}).get("mode", "")
raw = mode in {"raw_long_read", "raw_illumina", "raw_short_read"} or data.get("status") in {
    "created_raw_longread", "created_raw_illumina",
}
print("analysis" if raw else "source")
PY
    ) || die "Stage 01 reference 정책을 읽지 못했습니다: $qc"
  fi
  if [[ "$mode" == analysis ]]; then
    sample_analysis_reference "$sample"
  else
    sample_source_decode_reference "$sample"
  fi
}

sample_sequence_reference() {
  local sample=$1
  if [[ "$(sample_vcf_source "$sample")" == public ]]; then
    sample_vcf_reference "$sample"
  else
    printf '%s\n' "$REF_SEQUENCE"
  fi
}

# RNA 결과 디렉터리는 전역 PHASE_MODE가 아니라 샘플별 VCF evidence를 표시한다.
# manifest가 없는 구형 단일-sample 실행만 PHASE_MODE 호환값을 유지한다.
sample_rna_evidence_label() {
  local sample=$1 source
  if is_long_read_sample "$sample"; then
    if [[ -n "${PIPELINE_MANIFEST:-}" ]]; then
      source=$(sample_vcf_source "$sample")
      case "$source" in
        public|clair3) printf '%s\n' "$source" ;;
        *) die "long-read RNA evidence label을 결정할 수 없음: sample=$sample source=$source" ;;
      esac
    else
      printf '%s\n' "$PHASE_MODE"
    fi
  else
    printf 'panel\n'
  fi
}
is_long_read_sample() {
  local p
  p=$(normalize_platform "$(sample_platform "$1")")
  [[ "$p" == ont || "$p" == hifi || "$p" == clr ]]
}

require_manifest_sample_contract() {
  local sample=$1 platform read_type role source_alignment file_manifest alignment row_count alignment_count data_group rna
  [[ -n "${PIPELINE_MANIFEST:-}" ]] || return 0
  [[ -s "$PIPELINE_MANIFEST" ]] || die "PIPELINE_MANIFEST가 없거나 비어 있음: $PIPELINE_MANIFEST"
  awk -F'\t' 'NR==1 {
    expected[1]="sample_id"; expected[2]="platform"; expected[3]="read_type";
    expected[4]="dataset_role"; expected[5]="project_id"; expected[6]="cohort";
    expected[7]="alignment_path";
    for (i=1; i<=7; i++) if ($i != expected[i]) exit 1;
    exit 0
  }' "$PIPELINE_MANIFEST" || die "manifest 필수 1~7열 schema 오류: $PIPELINE_MANIFEST"
  manifest_has_sample "$sample" || die "manifest에 sample 행이 없음: sample=$sample manifest=$PIPELINE_MANIFEST"
  row_count=$(awk -F'\t' -v S="$sample" 'NR>1 && $1==S{n++}END{print n+0}' "$PIPELINE_MANIFEST")
  [[ "$row_count" == 1 ]] || die "manifest sample 행 중복: sample=$sample count=$row_count manifest=$PIPELINE_MANIFEST"
  platform=$(normalize_platform "$(sample_platform "$sample")")
  read_type=$(sample_read_type "$sample")
  role=$(sample_dataset_role "$sample")
  source_alignment=$(sample_source_alignment "$sample")
  file_manifest=$(sample_read_files_manifest "$sample")
  alignment=$(sample_alignment "$sample")
  case "$platform:$read_type" in
    illumina:short_read|ont:long_read|hifi:long_read|clr:long_read) ;;
    *) die "manifest platform/read_type 불일치: sample=$sample platform=$platform read_type=$read_type" ;;
  esac
  case "$role" in
    pre_training|post_training) ;;
    *) die "manifest dataset_role 오류: sample=$sample dataset_role=$role" ;;
  esac
  data_group=$(sample_data_group "$sample")
  rna=$(sample_rna "$sample")
  if [[ "$data_group" == wgs_rna_pair ]]; then
    [[ -n "$rna" && "$rna" != "-" ]] || die "manifest RNA pair 입력 누락: sample=$sample data_group=$data_group rna_path=$rna"
  elif [[ -n "$rna" && "$rna" != "-" ]]; then
    die "manifest data_group/rna_path 불일치: sample=$sample data_group=$data_group rna_path=$rna"
  fi
  if ! sample_wgs_reads_required "$sample"; then return 0; fi
  if [[ -n "$source_alignment" && "$source_alignment" != "-" ]]; then
    [[ -z "$file_manifest" ]] || die "alignment_path와 read_files_manifest를 함께 지정할 수 없음: sample=$sample"
    alignment_count=$(awk -F'\t' -v P="$(manifest_field "$sample" 7)" 'NR>1 && $7==P{n++}END{print n+0}' "$PIPELINE_MANIFEST")
    [[ "$alignment_count" == 1 ]] || die "manifest alignment_path 중복: sample=$sample path=$source_alignment count=$alignment_count"
  else
    [[ -n "$file_manifest" ]] || die "alignment_path 또는 read_files_manifest가 필요함: sample=$sample"
    if ! "$RUNNER" validate-wgs-input-manifest --manifest "$PIPELINE_MANIFEST" --sample "$sample" >/dev/null; then
      die "raw/file manifest 입력 계약 검증 실패: sample=$sample"
    fi
  fi
  [[ -n "$alignment" && "$alignment" != "-" ]] || die "prepared alignment 경로 해석 실패: sample=$sample"
}

clair3_model_path_for_platform() {
  local platform
  platform=$(normalize_platform "$1")
  if [[ -n "${CLAIR3_MODEL_PATH_OVERRIDE:-}" ]]; then
    echo "$CLAIR3_MODEL_PATH_OVERRIDE"
    return 0
  fi
  case "$platform" in
    ont) echo "$CLAIR3_MODEL_PATH_ONT" ;;
    hifi) echo "$CLAIR3_MODEL_PATH_HIFI" ;;
    clr) die "PacBio CLR용 Clair3 model은 구성하지 않음: 공개 VCF(vcf_source=public)를 사용해야 함" ;;
    *) die "Clair3 model selection only supports long-read platform ont|hifi: $platform" ;;
  esac
}

clair3_model_path_for_sample() {
  clair3_model_path_for_platform "$(sample_platform "$1")"
}

longread_sequence_vcf_source_for_sample() {
  if is_long_read_sample "$1"; then
    sample_vcf_source "$1"
  else
    echo illumina_panel
  fi
}

longread_haplotag_vcf_source_for_sample() {
  if is_long_read_sample "$1"; then
    sample_vcf_source "$1"
  else
    echo illumina_panel
  fi
}

longread_needs_clair3() {
  is_long_read_sample "$1" && [[ "$(sample_vcf_source "$1")" == clair3 ]]
}

sample_haplotagged_alignment() {
  local sample=$1 root platform source label
  root=$(sample_final_root "$sample")
  platform=$(normalize_platform "$(sample_platform "$sample")")
  if [[ "$platform" == ont || "$platform" == hifi || "$platform" == clr ]]; then
    source=$(sample_vcf_source "$sample")
    if [[ "$source" == public ]]; then label=public_vcf; else label=$source; fi
    echo "$root/$sample/haplotagged/$sample.${label}_longread.haplotagged.bam"
  else
    echo "$root/$sample/haplotagged/$sample.illumina.hg38.haplotagged.bam"
  fi
}

sample_hap_alignment() {
  local sample=$1 hap=$2 root ext
  root=$(sample_final_root "$sample")
  ext=$WGS_REALIGN_FORMAT
  echo "$root/$sample/$REALIGNED_SUBDIR/HAP$hap/$sample.personal.HAP$hap.$ext"
}

sample_hap_target() {
  local sample=$1 hap=$2 root
  root=$(sample_final_root "$sample")
  echo "$root/$sample/$REALIGNED_SUBDIR/HAP$hap/target.fa"
}

have_index() {
  local p=$1
  [[ -s "$p.bai" || -s "$p.csi" || -s "${p%.bam}.bai" || -s "$p.crai" || -s "${p%.cram}.crai" ]]
}

expand_chroms() {
  local raw=${1:-$CHROMS} part left right i
  raw=${raw//,/ }
  for part in $raw; do
    [[ -n "$part" ]] || continue
    if [[ "$part" == *-* ]]; then
      left=${part%-*}; right=${part#*-}
      left=${left#chr}; right=${right#chr}
      [[ "$left" =~ ^[0-9]+$ && "$right" =~ ^[0-9]+$ && "$left" -le "$right" ]] \
        || die "invalid chromosome range: $part"
      for ((i=left; i<=right; i++)); do printf 'chr%s\n' "$i"; done
    else
      part=${part#chr}
      [[ "$part" =~ ^[0-9]+$ ]] || die "invalid chromosome token: $part"
      printf 'chr%s\n' "$part"
    fi
  done
}

chroms_space() {
  expand_chroms "${1:-$CHROMS}" | paste -sd ' ' -
}

first_chrom() {
  expand_chroms "${1:-$CHROMS}" | sed -n '1p'
}

require_hap_contig_assets() {
  local root=$1 chroms=${2:-$CHROMS} require_chain=${3:-0} h c
  [[ -d "$root" ]] || die "personal haplotype root 없음: $root"
  while IFS= read -r c; do
    [[ -n "$c" ]] || continue
    for h in 1 2; do
      require_file "$root/HAP$h/$c.fa"
      require_file "$root/HAP$h/$c.fa.fai"
      if [[ "$require_chain" == 1 ]]; then
        require_file "$root/HAP$h/$c.chain"
      fi
    done
  done < <(expand_chroms "$chroms")
}

require_file() { [[ -s "$1" ]] || die "missing file: $1"; }
require_exe() { [[ -x "$1" ]] || die "missing executable: $1"; }

run_cli() {
  "$RUNNER" "$@"
}


sequence_source_for_sample() {
  local sample=$1 value platform
  value=$(echo "${SEQUENCE_SOURCE:-auto}" | tr '[:upper:]' '[:lower:]' | tr '-' '_' | tr -d ' ')
  if [[ "$value" != "auto" ]]; then
    echo "$value"
    return 0
  fi
  platform=$(normalize_platform "$(sample_platform "$sample")")
  # v2 default is VCF-consensus for every platform. Long-read de novo assembly
  # remains available only when SEQUENCE_SOURCE=denovo is explicitly requested.
  echo vcf_consensus
}

denovo_label() {
  local sample=$1 platform contig_label label
  if [[ -n "${DENOVO_LABEL:-}" ]]; then
    echo "$DENOVO_LABEL"
    return 0
  fi
  platform=$(normalize_platform "$(sample_platform "$sample")")
  contig_label=$(echo "${DENOVO_CONTIGS:-$CHROMS}" | tr ',' '_' | tr ' ' '_' | tr ':' '_' | tr '-' '_')
  label="${DENOVO_ASSEMBLER}_${platform}_${contig_label}"
  echo "$label"
}

maybe_fixed_window_arg() {
  if [[ "$FIXED_WINDOW_START" == 1 ]]; then
    echo --fixed-window-start
  elif [[ -n "$OFFSET_UNIT" ]]; then
    echo --offset-unit "$OFFSET_UNIT"
  fi
}

maybe_max_windows_arg() {
  [[ -n "$COMPACT_MAX_WINDOWS" ]] && echo --max-windows "$COMPACT_MAX_WINDOWS"
}

print_config() {
  cat <<EOF
PIPELINE_DIR=$PIPELINE_DIR
RUN_ID=$RUN_ID
DATASET_ROOT=$DATASET_ROOT
PIPELINE_MANIFEST=$PIPELINE_MANIFEST
DEFAULT_COHORT=$DEFAULT_COHORT
DEFAULT_PROJECT_ID=$DEFAULT_PROJECT_ID
DEFAULT_PLATFORM=$DEFAULT_PLATFORM
REF_DECODE=$REF_DECODE
REF_SEQUENCE=$REF_SEQUENCE
CLAIR3_ENV=$CLAIR3_ENV
CLAIR3_PY=$CLAIR3_PY
CLAIR3_MODEL_PATH=$CLAIR3_MODEL_PATH
CLAIR3_MODEL_PATH_OVERRIDE=${CLAIR3_MODEL_PATH_OVERRIDE:-}
CLAIR3_MODEL_PATH_ONT=$CLAIR3_MODEL_PATH_ONT
CLAIR3_MODEL_PATH_HIFI=$CLAIR3_MODEL_PATH_HIFI
CHROMS=$CHROMS
WINDOW_SIZE=$WINDOW_SIZE
WGS_MIN_MAPPED_FRAC=$WGS_MIN_MAPPED_FRAC
QC_RATIO_POLICY=$QC_RATIO_POLICY
WGS_ROUTE_JOBS=$WGS_ROUTE_JOBS
WGS_HAP_JOBS=$WGS_HAP_JOBS
WGS_REALIGN_THREADS=$WGS_REALIGN_THREADS
WGS_SORT_THREADS=$WGS_SORT_THREADS
STAGE11_BWA_INDEX_MODE=$STAGE11_BWA_INDEX_MODE
MFM_SCRATCH_ROOT=${MFM_SCRATCH_ROOT:-}
MFM_STAGE11_SCRATCH_ROOT=${MFM_STAGE11_SCRATCH_ROOT:-}
MFM_STAGE12_SCRATCH_ROOT=${MFM_STAGE12_SCRATCH_ROOT:-}
HAP_READ_MODE=$HAP_READ_MODE
REALIGNED_SUBDIR=$REALIGNED_SUBDIR
COMPACT_LAYOUT_NAME=$COMPACT_LAYOUT_NAME
COMPACT_VALIDATION_LEVEL=$COMPACT_VALIDATION_LEVEL
STAGE12_TRACE_POLICY=$STAGE12_TRACE_POLICY
ALLOW_CLEANUP_WITHOUT_TRACE=$ALLOW_CLEANUP_WITHOUT_TRACE
UNKNOWN_READ_POLICY=$UNKNOWN_READ_POLICY
READ_ROUTING_POLICY=$READ_ROUTING_POLICY
SEQUENCE_SOURCE=$SEQUENCE_SOURCE
VCF_SOURCE_POLICY=$VCF_SOURCE_POLICY
LONGREAD_SEQUENCE_VCF_SOURCE=$LONGREAD_SEQUENCE_VCF_SOURCE
LONGREAD_HAPLOTAG_VCF_SOURCE=$LONGREAD_HAPLOTAG_VCF_SOURCE
ENFORCE_VCF_SOURCE_POLICY=$ENFORCE_VCF_SOURCE_POLICY
DENOVO_ASSEMBLER=$DENOVO_ASSEMBLER
HIFIASM_BIN=$HIFIASM_BIN
CLEAN_INTERMEDIATES_ON_SUCCESS=$CLEAN_INTERMEDIATES_ON_SUCCESS
CLEAN_INTERMEDIATES_DRY_RUN=$CLEAN_INTERMEDIATES_DRY_RUN
CLEAN_AFTER_RNA_ONLY=$CLEAN_AFTER_RNA_ONLY
CLEAN_SHARED_ASSETS_ON_SUCCESS=$CLEAN_SHARED_ASSETS_ON_SUCCESS
RNA_COLLATE_THREADS=$RNA_COLLATE_THREADS
RNA_FASTQ_THREADS=$RNA_FASTQ_THREADS
RNA_BAM_INDEX_THREADS=$RNA_BAM_INDEX_THREADS
RNA_MIN_PAIRED_FRACTION=$RNA_MIN_PAIRED_FRACTION
RNA_REQUIRE_SM_MATCH=$RNA_REQUIRE_SM_MATCH
HAP_VARIANT_MIN_MAPQ=$HAP_VARIANT_MIN_MAPQ
HAP_VARIANT_MIN_DEPTH=$HAP_VARIANT_MIN_DEPTH
HAP_VARIANT_MIN_ALT_FRAC=$HAP_VARIANT_MIN_ALT_FRAC
HAP_VARIANT_THREADS=$HAP_VARIANT_THREADS
RNA_CNV_MIN_RD=$RNA_CNV_MIN_RD
RNA_CNV_MAX_RD=$RNA_CNV_MAX_RD
RNA_CNV_FILTER=$RNA_CNV_FILTER
RNA_SAMPLE_MEMORY_GIB=${RNA_SAMPLE_MEMORY_GIB:-48}
RNA_SAMPLE_TMP_GIB=${RNA_SAMPLE_TMP_GIB:-80}
EOF
}


# ── process_code/v2 RNA branch 호환 alias ───────────────────────────────
# 13~19 RNA 스크립트는 v2_final/v2_work/v2_rna/PYTHON/SAMTOOLS 등
# 구형 변수명을 사용한다. WGS는 integrated naming을 쓰되 RNA 구현을 보존하기 위해
# 여기에서만 호환 shim을 제공한다.
MFM_ROOT=${MFM_ROOT:-${INPUT_ROOT:-$PIPELINE_DIR}}
V2_ROOT=${V2_ROOT:-$DATASET_ROOT}
V2_RUN=${V2_RUN:-$RUN_ID}
V2_MANIFEST=${V2_MANIFEST:-${PIPELINE_MANIFEST:-}}
if [[ -z "${PIPELINE_MANIFEST:-}" && -n "${V2_MANIFEST:-}" ]]; then
  PIPELINE_MANIFEST="$V2_MANIFEST"
fi
CONDA=${CONDA:-$CONDA_ENVS_DIR}
PYTHON=${PYTHON:-$PYTHON_BIN}
BCFTOOLS=${BCFTOOLS:-$BCFTOOLS_BIN}
MINIMAP2=${MINIMAP2:-$MINIMAP2_BIN}
SAMTOOLS=${SAMTOOLS:-$SAMTOOLS_BIN}
WHATSHAP=${WHATSHAP:-$WHATSHAP_BIN}
AUTOSOMES=${AUTOSOMES:-"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"}
PHASE_MODE=${PHASE_MODE:-clair3}
STAR_BIN=${STAR_BIN:-$(find_tool STAR)}
PICARD_BIN=${PICARD_BIN:-$(find_tool picard)}
STAR_INDEX_THREADS=${STAR_INDEX_THREADS:-4}
STAR_ALIGN_THREADS=${STAR_ALIGN_THREADS:-6}
RNA_COLLATE_THREADS=${RNA_COLLATE_THREADS:-8}
RNA_FASTQ_THREADS=${RNA_FASTQ_THREADS:-4}
RNA_BAM_INDEX_THREADS=${RNA_BAM_INDEX_THREADS:-4}
# record_only에서는 비교·기록용 참고선이며 strict에서만 중단 기준이다.
RNA_MIN_PAIRED_FRACTION=${RNA_MIN_PAIRED_FRACTION:-0.90}
RNA_REQUIRE_SM_MATCH=${RNA_REQUIRE_SM_MATCH:-1}
STAR_SJDB_OVERHANG=${STAR_SJDB_OVERHANG:-74}
STAR_KEEP_INDEX=${STAR_KEEP_INDEX:-0}
# RNA 재정렬 뒤 Picard MarkDuplicates로 duplicate flag만 표기한다. 원본 BAM과
# MarkDuplicates BAM을 모두 보존하여 downstream이 duplicate 정책을 명시적으로 선택한다.
RNA_MARKDUP_HEAP_GB=${RNA_MARKDUP_HEAP_GB:-$MFM_RNA_MARKDUP_HEAP_GB_DEFAULT}
RNA_STAR_INDEX_HAP_JOBS=${RNA_STAR_INDEX_HAP_JOBS:-2}
RNA_STAR_INDEX_BENCHMARK=${RNA_STAR_INDEX_BENCHMARK:-0}
RNA_ALLELIC_JOBS=${RNA_ALLELIC_JOBS:-auto}
RNA_ALLELIC_BENCHMARK=${RNA_ALLELIC_BENCHMARK:-0}
RNA_ALLELIC_BATCH_SIZE=${RNA_ALLELIC_BATCH_SIZE:-50000}
export RNA_ALLELIC_BATCH_SIZE
[[ "$RNA_ALLELIC_BATCH_SIZE" =~ ^[1-9][0-9]*$ ]] || { echo "[FATAL] RNA_ALLELIC_BATCH_SIZE는 양수이어야 함" >&2; exit 1; }
case "$RNA_ALLELIC_BENCHMARK" in 0|1) ;; *) echo "[FATAL] RNA_ALLELIC_BENCHMARK는 0|1이어야 함" >&2; exit 1 ;; esac
[[ "$RNA_ALLELIC_JOBS" == auto || "$RNA_ALLELIC_JOBS" =~ ^[1-9][0-9]*$ ]] || { echo "[FATAL] RNA_ALLELIC_JOBS는 auto 또는 양수이어야 함" >&2; exit 1; }
export RNA_ALLELIC_BENCHMARK
export RNA_STAR_INDEX_HAP_JOBS RNA_STAR_INDEX_BENCHMARK RNA_ALLELIC_JOBS
case "$RNA_STAR_INDEX_HAP_JOBS" in 1|2) ;; *) echo '[FATAL] RNA_STAR_INDEX_HAP_JOBS는 1|2이어야 함' >&2; exit 1 ;; esac
case "$RNA_STAR_INDEX_BENCHMARK" in 0|1) ;; *) echo '[FATAL] RNA_STAR_INDEX_BENCHMARK는 0|1이어야 함' >&2; exit 1 ;; esac
RNA_PIPELINE_OVERLAP=${RNA_PIPELINE_OVERLAP:-auto}
case "$RNA_PIPELINE_OVERLAP" in auto|0|1) ;; *) echo '[FATAL] RNA_PIPELINE_OVERLAP는 auto/0/1임' >&2; exit 1 ;; esac
export RNA_PIPELINE_OVERLAP
RNA_ALIGN_HAP_JOBS=${RNA_ALIGN_HAP_JOBS:-2}
RNA_COVERAGE_HAP_JOBS=${RNA_COVERAGE_HAP_JOBS:-2}
export RNA_ALIGN_HAP_JOBS RNA_COVERAGE_HAP_JOBS
for rna_hap_jobs in "$RNA_ALIGN_HAP_JOBS" "$RNA_COVERAGE_HAP_JOBS"; do
  [[ "$rna_hap_jobs" == 1 || "$rna_hap_jobs" == 2 ]] || { echo '[FATAL] RNA HAP jobs는1/2임' >&2; exit 1; }
done
RNA_MARKDUP_HAP_JOBS=${RNA_MARKDUP_HAP_JOBS:-2}
export RNA_MARKDUP_HEAP_GB RNA_MARKDUP_HAP_JOBS
RNA_GPU_MARKDUP_MODE=${RNA_GPU_MARKDUP_MODE:-integrated}
export RNA_GPU_MARKDUP_MODE
case "$RNA_GPU_MARKDUP_MODE" in integrated|separate) ;; *) echo '[FATAL] RNA_GPU_MARKDUP_MODE는 integrated|separate이어야 함' >&2; exit 1 ;; esac
case "$RNA_MARKDUP_HAP_JOBS" in 1|2) ;; *) echo "[FATAL] RNA_MARKDUP_HAP_JOBS는 1|2이어야 함" >&2; exit 1 ;; esac
# 2026-09-28: 기본 annotation 을 GENCODE v47(Ensembl 113, GRCh38)로 올렸다. v29 는 references/ 에 남겨 둔다.
GENCODE_GTF=${GENCODE_GTF:-$MFM_REFERENCE_DIR/gencode.v47.annotation.gtf.gz}
# GTF 가 바뀌면 lift GTF·STAR index·정렬이 모두 달라지므로 캐시 디렉터리 이름에 이 tag 를 넣는다.
# 파일명이 gencode.vNN.* 이면 gencodevNN, 아니면 경로 md5 앞 8자리.
gencode_tag() {
  local b; b=$(basename "$1")
  if [[ "$b" =~ ^gencode\.v([0-9]+[A-Za-z0-9]*)\. ]]; then printf 'gencodev%s\n' "${BASH_REMATCH[1]}"
  else printf 'gtf%s\n' "$(printf '%s' "$1" | md5sum | cut -c1-8)"; fi
}
GENCODE_TAG=${GENCODE_TAG:-$(gencode_tag "$GENCODE_GTF")}
RNA_ZARR_CHUNK=${RNA_ZARR_CHUNK:-4194304}
AG_NORM_POWER=${AG_NORM_POWER:-0.75}
AG_NORM_CLIP=${AG_NORM_CLIP:-10.0}
# ── RNA index 추가 contig · CNV 필터 ──────────────────────────────────────
# STAR index 에만 넣는 비-autosome contig. GTEx index(noALT·noHLA·noDecoy) 중 개인 HAP(chr1–22)에 없는 것.
# 갈 곳 없는 read(chrM·chrX·rRNA 반복 등)를 흡수시킨다. coverage·count 에는 쓰지 않는다(chr1–22 만).
RNA_INDEX_ADDNL_SEQNS=${RNA_INDEX_ADDNL_SEQNS:-1}
RNA_ADDNL_SEQNS_SOURCE=${RNA_ADDNL_SEQNS_SOURCE:-$MFM_REFERENCE_DIR/Homo_sapiens_assembly38.fasta}
# GTEx index 와 같은 비-autosome contig (chrX·chrY·chrM·chrUn·*_random·chrEBV; decoy·alt·HLA 제외).
RNA_ADDNL_SEQNS_FASTA=${RNA_ADDNL_SEQNS_FASTA:-$MFM_REFERENCE_DIR/rna_index_nonautosome_gtex.fa}
# Stage 20은 HAP FASTA를 만든 이형접합 variant가 실제 WGS read에서도 지지되는지 확인한다.
# HG00096 EDA에서 alt fraction 0 부근 오류군과 0.4–0.5 지지군이 0.1에서 분리되어 이 값을 기본으로 둔다.
HAP_VARIANT_MIN_MAPQ=${HAP_VARIANT_MIN_MAPQ:-20}
HAP_VARIANT_MIN_DEPTH=${HAP_VARIANT_MIN_DEPTH:-10}
HAP_VARIANT_MIN_ALT_FRAC=${HAP_VARIANT_MIN_ALT_FRAC:-0.1}
# autosome별 독립 작업이므로 22를 요청하되 sample budget과 동적 stage lease가 실제 상한을 정한다.
HAP_VARIANT_THREADS=${HAP_VARIANT_THREADS:-22}
# CNV 필터 임계. AlleleSeq2 filter_cnv_sites.py와 같은 구간이다.
RNA_CNV_MIN_RD=${RNA_CNV_MIN_RD:-0.5}
RNA_CNV_MAX_RD=${RNA_CNV_MAX_RD:-1.5}
RNA_CNV_FILTER=${RNA_CNV_FILTER:-1}
case "$RNA_CNV_FILTER" in 0|1) ;; *) echo "[FATAL] RNA_CNV_FILTER는 0|1이어야 합니다" >&2; exit 1 ;; esac

# local runner가 정한 sample별 CPU 상한. 사용자가 LOCAL_JOBS를 늘려도 각
# 도구 thread를 함께 줄여 72코어 서버 oversubscription을 방지한다.
if [[ -n "${MFM_SAMPLE_CPU_BUDGET:-}" ]]; then
  [[ "$MFM_SAMPLE_CPU_BUDGET" =~ ^[0-9]+$ && "$MFM_SAMPLE_CPU_BUDGET" -ge 1 ]] || {
    echo "[FATAL] invalid MFM_SAMPLE_CPU_BUDGET: $MFM_SAMPLE_CPU_BUDGET" >&2; exit 1;
  }
  mfm_cap_threads() {
    local name=$1 value=${!1}
    if (( value > MFM_SAMPLE_CPU_BUDGET )); then printf -v "$name" '%d' "$MFM_SAMPLE_CPU_BUDGET"; fi
  }
  case "$WGS_HAP_JOBS" in
    auto)
      # runner의 명시적인 자원 거부는 존중한다. 직접 실행처럼 flag가 없으면
      # sample CPU 예산으로 판단해 16 CPU 이상에서 HAP 두 개를 병렬화한다.
      if [[ "${MFM_HAP_PARALLEL_ALLOWED:-1}" == 1 ]] && (( MFM_SAMPLE_CPU_BUDGET >= 16 )); then
        WGS_HAP_JOBS=2
      else
        WGS_HAP_JOBS=1
      fi
      ;;
    1|2) ;;
    *) echo "[FATAL] WGS_HAP_JOBS must be auto|1|2: $WGS_HAP_JOBS" >&2; exit 1 ;;
  esac
  case "$WGS_ROUTE_JOBS" in
    auto)
      WGS_ROUTE_JOBS=$MFM_SAMPLE_CPU_BUDGET
      (( WGS_ROUTE_JOBS <= 16 )) || WGS_ROUTE_JOBS=16
      ;;
    *)
      [[ "$WGS_ROUTE_JOBS" =~ ^[0-9]+$ && "$WGS_ROUTE_JOBS" -ge 1 ]] || {
        echo "[FATAL] WGS_ROUTE_JOBS must be auto or a positive integer: $WGS_ROUTE_JOBS" >&2; exit 1;
      }
      (( WGS_ROUTE_JOBS > MFM_SAMPLE_CPU_BUDGET )) && WGS_ROUTE_JOBS=$MFM_SAMPLE_CPU_BUDGET
      ;;
  esac
  # aligner와 samtools sort가 동시에 실행되고, HAP 병렬 시 그 파이프가
  # WGS_HAP_JOBS개 생긴다. 따라서 hap별 합산 예산으로 나눠 oversubscription을 막는다.
  wgs_per_hap_budget=$((MFM_SAMPLE_CPU_BUDGET / WGS_HAP_JOBS))
  (( wgs_per_hap_budget >= 2 )) || { WGS_HAP_JOBS=1; wgs_per_hap_budget=$MFM_SAMPLE_CPU_BUDGET; }
  if [[ "$WGS_AUTO_THREADS" == 1 ]]; then
    if (( MFM_WGS_SORT_EXPLICIT == 0 )); then
      WGS_SORT_THREADS=$((wgs_per_hap_budget / 4))
      (( WGS_SORT_THREADS >= 1 )) || WGS_SORT_THREADS=1
      (( WGS_SORT_THREADS <= 4 )) || WGS_SORT_THREADS=4
    fi
    if (( MFM_WGS_ALIGN_EXPLICIT == 0 )); then
      WGS_REALIGN_THREADS=$((wgs_per_hap_budget - WGS_SORT_THREADS))
      (( WGS_REALIGN_THREADS >= 1 )) || WGS_REALIGN_THREADS=1
    fi
  fi
  (( WGS_SORT_THREADS >= wgs_per_hap_budget )) && WGS_SORT_THREADS=$((wgs_per_hap_budget - 1))
  (( WGS_SORT_THREADS >= 1 )) || WGS_SORT_THREADS=1
  if (( WGS_REALIGN_THREADS + WGS_SORT_THREADS > wgs_per_hap_budget )); then
    WGS_REALIGN_THREADS=$((wgs_per_hap_budget - WGS_SORT_THREADS))
    (( WGS_REALIGN_THREADS >= 1 )) || WGS_REALIGN_THREADS=1
  fi
  for name in CLAIR3_THREADS DENOVO_THREADS DENOVO_FASTQ_THREADS ASSEMBLY_ANCHOR_THREADS \
    STAGE04_ALIGN_THREADS STAGE04_SORT_THREADS STAGE04_EXTRACT_THREADS WINDOW_JOBS \
    STAR_INDEX_THREADS STAR_ALIGN_THREADS HAP_VARIANT_THREADS \
    RNA_COLLATE_THREADS RNA_FASTQ_THREADS RNA_BAM_INDEX_THREADS; do
    mfm_cap_threads "$name"
  done
  if (( RNA_COLLATE_THREADS + RNA_FASTQ_THREADS > MFM_SAMPLE_CPU_BUDGET )); then
    RNA_FASTQ_THREADS=$((MFM_SAMPLE_CPU_BUDGET / 3))
    (( RNA_FASTQ_THREADS >= 1 )) || RNA_FASTQ_THREADS=1
    RNA_COLLATE_THREADS=$((MFM_SAMPLE_CPU_BUDGET - RNA_FASTQ_THREADS))
    (( RNA_COLLATE_THREADS >= 1 )) || RNA_COLLATE_THREADS=1
  fi
  case "$COMPACT_JOBS" in
    auto)
      if [[ "${MFM_CPU_UTILIZATION:-maximum}" == maximum ]]; then
        # 염색체별 read 변환은 주로 단일 CPU 계산이다. 예산의 약수 조건은
        # 소수 배정에서 jobs=1로 축소되므로 실제 염색체 수까지만 제한한다.
        COMPACT_JOBS=$(expand_chroms "$CHROMS" | wc -l)
        (( COMPACT_JOBS > MFM_SAMPLE_CPU_BUDGET )) && COMPACT_JOBS=$MFM_SAMPLE_CPU_BUDGET
      else
      # contig process와 각 process의 Blosc thread 곱이 sample 예산 안에
      # 들도록 하며, writer RSS 급증을 피하기 위해 process는 최대 8개다.
      COMPACT_JOBS=$((MFM_SAMPLE_CPU_BUDGET / 4))
      (( COMPACT_JOBS >= 1 )) || COMPACT_JOBS=1
      (( COMPACT_JOBS <= 8 )) || COMPACT_JOBS=8
      # 동일한 Blosc thread 수를 쓰는 process들의 곱이 CPU 예산과 정확히
      # 맞도록 가능한 가장 큰 약수를 선택한다(36코어면 6×6).
      while (( COMPACT_JOBS > 1 && MFM_SAMPLE_CPU_BUDGET % COMPACT_JOBS != 0 )); do
        COMPACT_JOBS=$((COMPACT_JOBS - 1))
      done
      fi
      ;;
    *)
      [[ "$COMPACT_JOBS" =~ ^[0-9]+$ && "$COMPACT_JOBS" -ge 1 ]] || {
        echo "[FATAL] COMPACT_JOBS must be auto or a positive integer: $COMPACT_JOBS" >&2; exit 1;
      }
      (( COMPACT_JOBS > MFM_SAMPLE_CPU_BUDGET )) && COMPACT_JOBS=$MFM_SAMPLE_CPU_BUDGET
      ;;
  esac
  MFM_BLOSC_THREADS=$((MFM_SAMPLE_CPU_BUDGET / COMPACT_JOBS))
  if [[ "${MFM_CPU_UTILIZATION:-maximum}" == maximum && "$MFM_COMPACT_JOBS_REQUESTED" == auto ]]; then
    MFM_BLOSC_THREADS=1
  fi
  (( MFM_BLOSC_THREADS >= 1 )) || MFM_BLOSC_THREADS=1
  export MFM_BLOSC_THREADS
else
  [[ "$COMPACT_JOBS" == auto ]] && COMPACT_JOBS=4
  # sample별 CPU budget 없이 단독 실행할 때의 보수적 기본값. trace가 이
  # 값을 상속하되, 명시적 MFM_TRACE_BLOSC_THREADS가 있으면 그것을 우선한다.
  MFM_BLOSC_THREADS=${MFM_BLOSC_THREADS:-8}
  [[ "$MFM_BLOSC_THREADS" =~ ^[0-9]+$ && "$MFM_BLOSC_THREADS" -ge 1 ]] || {
    echo "[FATAL] MFM_BLOSC_THREADS must be a positive integer: $MFM_BLOSC_THREADS" >&2; exit 1;
  }
  [[ "$WGS_ROUTE_JOBS" == auto ]] && WGS_ROUTE_JOBS=8
  case "$WGS_HAP_JOBS" in
    auto) WGS_HAP_JOBS=2 ;;
    1|2) ;;
    *) echo "[FATAL] WGS_HAP_JOBS must be auto|1|2: $WGS_HAP_JOBS" >&2; exit 1 ;;
  esac
fi
[[ "$MFM_TRACE_MAX_JOBS" =~ ^[0-9]+$ && "$MFM_TRACE_MAX_JOBS" -ge 1 ]] || {
  echo "[FATAL] MFM_TRACE_MAX_JOBS must be a positive integer: $MFM_TRACE_MAX_JOBS" >&2; exit 1;
}
[[ "$MFM_TRACE_CPU_RESERVE" =~ ^[0-9]+$ ]] || {
  echo "[FATAL] MFM_TRACE_CPU_RESERVE must be a non-negative integer: $MFM_TRACE_CPU_RESERVE" >&2; exit 1;
}
trace_contig_count=$(expand_chroms "$CHROMS" | wc -l)
trace_task_cap=$((trace_contig_count * 2))
(( trace_task_cap >= 1 )) || { echo "[FATAL] no trace task for CHROMS=$CHROMS" >&2; exit 1; }
(( MFM_TRACE_MAX_JOBS < trace_task_cap )) && trace_task_cap=$MFM_TRACE_MAX_JOBS
if [[ -n "${MFM_SAMPLE_CPU_BUDGET:-}" ]]; then
  # 부모 shell·I/O가 쓸 최소 여유를 남긴다. 작은 예산에서는 적어도 한
  # worker는 보장한다. 각 worker의 Blosc thread도 이 나눗셈 안에 포함한다.
  trace_available=$((MFM_SAMPLE_CPU_BUDGET - MFM_TRACE_CPU_RESERVE))
  (( trace_available >= 1 )) || trace_available=1
  (( trace_available < trace_task_cap )) && trace_task_cap=$trace_available
fi
case "$TRACE_JOBS" in
  auto) TRACE_JOBS=$trace_task_cap ;;
  *)
    [[ "$TRACE_JOBS" =~ ^[0-9]+$ && "$TRACE_JOBS" -ge 1 ]] || {
      echo "[FATAL] TRACE_JOBS must be auto or a positive integer: $TRACE_JOBS" >&2; exit 1;
    }
    (( TRACE_JOBS > trace_task_cap )) && TRACE_JOBS=$trace_task_cap
    ;;
esac
if [[ "$MFM_BLOSC_THREADS_REQUESTED" != auto ]]; then
  [[ "$MFM_BLOSC_THREADS_REQUESTED" =~ ^[1-9][0-9]*$ ]] || { echo "[FATAL] MFM_BLOSC_THREADS_REQUESTED는 auto 또는 양수여야 함" >&2; exit 1; }
  MFM_BLOSC_THREADS=$MFM_BLOSC_THREADS_REQUESTED
  export MFM_BLOSC_THREADS
fi
case "$MFM_TRACE_BLOSC_THREADS" in
  auto)
    if [[ -n "${MFM_SAMPLE_CPU_BUDGET:-}" ]]; then
      MFM_TRACE_BLOSC_THREADS=$((MFM_SAMPLE_CPU_BUDGET / TRACE_JOBS))
      if [[ "${MFM_CPU_UTILIZATION:-maximum}" == maximum && "$MFM_TRACE_JOBS_REQUESTED" == auto ]]; then
        MFM_TRACE_BLOSC_THREADS=1
      fi
    else
      MFM_TRACE_BLOSC_THREADS=$MFM_BLOSC_THREADS
    fi
    (( MFM_TRACE_BLOSC_THREADS >= 1 )) || MFM_TRACE_BLOSC_THREADS=1
    ;;
  *)
    [[ "$MFM_TRACE_BLOSC_THREADS" =~ ^[0-9]+$ && "$MFM_TRACE_BLOSC_THREADS" -ge 1 ]] || {
      echo "[FATAL] MFM_TRACE_BLOSC_THREADS must be auto or a positive integer: $MFM_TRACE_BLOSC_THREADS" >&2; exit 1;
    }
    if [[ -n "${MFM_SAMPLE_CPU_BUDGET:-}" ]]; then
      trace_blosc_cap=$((MFM_SAMPLE_CPU_BUDGET / TRACE_JOBS))
      (( trace_blosc_cap >= 1 )) || trace_blosc_cap=1
      (( MFM_TRACE_BLOSC_THREADS > trace_blosc_cap )) && MFM_TRACE_BLOSC_THREADS=$trace_blosc_cap
    fi
    ;;
esac
export WGS_HAP_JOBS WGS_ROUTE_JOBS WGS_REALIGN_THREADS WGS_SORT_THREADS TRACE_JOBS MFM_TRACE_BLOSC_THREADS

v2_final() { echo "$(sample_final_root "$1")/$1"; }
v2_work()  { echo "$(sample_work_root "$1")"; }
v2_platform() { sample_platform "$1"; }
v2_wgs() { sample_alignment "$1"; }
v2_rna() {
  local value
  value=$(sample_rna "$1" 2>/dev/null || true)
  [[ -n "$value" ]] && echo "$value" || echo "-"
}
v2_is_long() { is_long_read_sample "$1"; }

v2_require() {
  local miss=0
  for t in "$PYTHON" "$SAMTOOLS" "$BCFTOOLS" "$MINIMAP2"; do
    [[ -x "$t" ]] || { warn "도구 없음: $t"; miss=1; }
  done
  if [[ -n "${PIPELINE_MANIFEST:-}" ]]; then
    [[ -s "$PIPELINE_MANIFEST" ]] || { warn "매니페스트 없음: $PIPELINE_MANIFEST"; miss=1; }
  fi
  (( miss == 0 )) || die "필수 도구/입력 누락"
}
