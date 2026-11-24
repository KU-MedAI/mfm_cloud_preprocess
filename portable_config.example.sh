#!/usr/bin/env bash
# 이 파일을 portable_config.sh로 복사한 뒤 새 서버 경로만 수정한다.
# 실행 전: source ./portable_config.sh

export INPUT_ROOT="/data/mfm_input"
export OUTPUT_ROOT="/data/mfm_output"
export RUN_ID="$(date +%y%m%d)"

# 새 서버의 jj_mfm conda 환경. pipeline은 이 Python을 우선 사용한다.
export MFM_CORE_ENV="/opt/conda/envs/jj_mfm"
export PYTHON_BIN="$MFM_CORE_ENV/bin/python"
# Stage 11/12 임시 파일용 노드 로컬 scratch(선택). 충분한 공간이 없으면 비워 둔다.
export MFM_SCRATCH_ROOT=""  # 예: /local_scratch/mfm

# manifest의 상대 입력은 INPUT_ROOT 아래에서 찾는다.
export PIPELINE_MANIFEST="$PWD/manifests/preprocess_manifest.tsv"

# manifest에 명시되지 않은 공통 reference/RNA resource의 기본 위치다.
export REF_DECODE="$INPUT_ROOT/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
export REF_SEQUENCE="$INPUT_ROOT/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
export ILLUMINA_REF_DECODE="$INPUT_ROOT/references/illumina_decode.fa"
export BLACKLIST_BED="$INPUT_ROOT/references/hg38_autosomes_blacklist_union.bed"
export GENCODE_GTF="$INPUT_ROOT/references/gencode.v47.annotation.gtf.gz"
export RNA_ADDNL_SEQNS_SOURCE="$INPUT_ROOT/references/Homo_sapiens_assembly38.fasta"
export RNA_ADDNL_SEQNS_FASTA="$INPUT_ROOT/references/addnl_seqns_unplaced_unlocalized_EBV_decoy.fa"

# 현재 선정 cohort는 제공 VCF를 사용하며 Clair3 self-call은 기본 비활성이다.
export VCF_SOURCE_POLICY=manifest
export ENABLE_CLAIR3_SELF_CALL=0

# 서버 자원에 맞게 조정한다. 내부 도구 thread는 sample별 CPU 예산으로 자동 제한된다.
export MFM_CPU_BUDGET="$(nproc)"
export LOCAL_JOBS=2
# phased WGS 기본: 준비 4, HAP 재정렬 4, compact 2 samples
export WGS_PREP_SAMPLE_JOBS=4
export WGS_ALIGN_SAMPLE_JOBS=4
export WGS_COMPACT_SAMPLE_JOBS=2
export WGS_ROUTE_JOBS=auto
export WGS_HAP_JOBS=auto
export STAGE11_BWA_INDEX_MODE=ephemeral
export DOWNLOAD_JOBS=2
export RNA_SAMPLE_MEMORY_GIB=48
export RNA_SAMPLE_TMP_GIB=80
