#!/usr/bin/env bash
# 02 — 실제 WGS stage 실행 전 sample별 입력 계약 검사.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh
export PYTHONPATH="$PIPELINE_DIR:${PYTHONPATH:-}"

S=${1:?사용법: $0 <sample>}
require_manifest_sample_contract "$S"
ALIGN=$(sample_alignment "$S")
SOURCE_SAMPLE_ID=$(sample_source_sample_id "$S" 2>/dev/null || true)
[[ -n "$SOURCE_SAMPLE_ID" ]] || SOURCE_SAMPLE_ID="$S"
PLATFORM=$(normalize_platform "$(sample_platform "$S")")
# raw 입력을 새로 정렬한 canonical BAM은 analysis reference, 기존 정렬을
# 변환한 BAM은 source decode reference와 대조한다. VCF REF 계약은 별도다.
DECODE_REF=$(sample_alignment_reference "$S")
READ_TYPE=$(sample_read_type "$S")
ROLE=$(sample_dataset_role "$S")
COMMON=$(sample_common_vcf_dir "$S")
PERSAMPLE=$(sample_persample_vcf_dir "$S")
VCF_SOURCE=$(sample_vcf_source "$S")
SEQUENCE_REF=$(sample_sequence_reference "$S")
PANEL_ROOT=$PANEL_LOCAL_DIR
STATUS_ROOT=${RUN_STATUS_DIR:-$(sample_status_root "$S")}
OUT_DIR="$STATUS_ROOT/preflight"
OUT="$OUT_DIR/$S.preflight.tsv"
OUT_JSON="$OUT_DIR/$S.preflight.json"
mkdir -p "$OUT_DIR"
TMP="$OUT.part.$$"
FAIL=0
printf 'check\tstatus\tdetail\n' > "$TMP"

record() {
  local check=$1 status=$2 detail=${3:-}
  printf '%s\t%s\t%s\n' "$check" "$status" "$detail" >> "$TMP"
  [[ "$status" == PASS ]] || FAIL=1
}
STAGE04_QC=$(sample_stage04_qc "$S")
STAGE04_PROVENANCE=0
if [[ "$ALIGN" == *.bam && -s "$STAGE04_QC" ]]; then
  stage04_check=$("$PYTHON" - "$STAGE04_QC" "$ALIGN" "$PLATFORM" <<'PY'
import json
import pathlib
import sys
try:
    data = json.load(open(sys.argv[1]))
    if data.get("status") not in {
        "created_raw_longread", "created_raw_illumina", "converted_aligned", "reused_canonical_bam",
    } or data.get("canonical_format") != "bam":
        print("not_generated")
    else:
        alignment = pathlib.Path(sys.argv[2]).resolve()
        platform = sys.argv[3]
        fingerprint = data.get("normalization_reference", {})
        source_reference = pathlib.Path(fingerprint.get("path", ""))
        duplicate = data.get("duplicate_marking", {})
        duplicate_valid = True
        if platform == "illumina":
            action = duplicate.get("action")
            duplicate_valid = (
                duplicate.get("status") == "pass"
                and duplicate.get("remove_duplicates") is False
                and action in {"picard_marked", "preserved_source_picard"}
            )
            if action == "picard_marked":
                metrics_path = duplicate.get("metrics_path")
                duplicate_valid = duplicate_valid and bool(metrics_path) and pathlib.Path(metrics_path).is_file()
            elif action == "preserved_source_picard":
                provenance = duplicate.get("source_provenance", [])
                duplicate_valid = duplicate_valid and bool(provenance) and all(
                    item.get("trusted_picard_marking") is True for item in provenance
                )
        valid = (
            pathlib.Path(data.get("analysis_alignment", "")).resolve() == alignment
            and source_reference.is_file()
            and pathlib.Path(str(source_reference) + ".fai").is_file()
            and duplicate_valid
        )
        print("valid" if valid else "invalid")
except Exception:
    print("invalid")
PY
  )
  case "$stage04_check" in
    valid) STAGE04_PROVENANCE=1; record stage04_reference_provenance PASS "$STAGE04_QC" ;;
    invalid) record stage04_reference_provenance FAIL "$STAGE04_QC" ;;
  esac
fi
if [[ "$PLATFORM" == illumina ]]; then
  if [[ "$STAGE04_PROVENANCE" == 1 ]]; then
    record illumina_duplicate_provenance PASS "$STAGE04_QC"
  else
    record illumina_duplicate_provenance FAIL "Stage 01 QC v3의 Picard metrics 또는 trusted source provenance 필요"
  fi
fi

# 실행 단계 02는 source alignment decode reference와 personal-FASTA/VCF reference를 별도
# 계약으로 검사한다. 두 FASTA는 의도적으로 sequence가 다를 수 있으므로 서로의 M5 일치를
# 요구하지 않는다. BAM은 @SQ M5가 없을 수 있어 Stage 01이 기록한 source reference의 22개
# autosome M5를 decode FASTA와 직접 대조한다. CRAM/M5 보유 BAM은 header M5를 우선 사용한다.
AUTOSOME_REGIONS=${AUTOSOMES// /,}
CRAM_HEADER=$({ "$SAMTOOLS_BIN" view -H "$ALIGN" 2>/dev/null || true; })

fai_length() {
  local fasta=$1 contig=$2
  awk -F'\t' -v c="$contig" '$1==c {print $2; exit}' "$fasta.fai" 2>/dev/null || true
}

header_sq_value() {
  local contig=$1 key=$2
  awk -F'\t' -v c="$contig" -v k="$key" '
    $1=="@SQ" {
      found=0
      for (i=2; i<=NF; i++) if ($i=="SN:" c) found=1
      if (!found) next
      for (i=2; i<=NF; i++) if ($i ~ ("^" k ":")) {sub("^" k ":", "", $i); print $i; exit}
    }
  ' <<< "$CRAM_HEADER"
}

# 공유 성공 cache는 실제 입력·index·reference·검사기 변경 시 무효화됨.
PREFLIGHT_CACHE_ROOT=${PREFLIGHT_CACHE_ROOT:-$DATASET_ROOT/_validation_cache/input_preflight}
PREFLIGHT_JOBS=${PREFLIGHT_JOBS:-${MFM_CPU_BUDGET:-1}}
[[ "$PREFLIGHT_JOBS" =~ ^[1-9][0-9]*$ ]] || die "PREFLIGHT_JOBS는 양수여야 함"
VCF_BATCH_RESULTS="$OUT_DIR/$S.vcf_ref_results.tsv"

build_reference_m5_cache() {
  local fasta=$1 cache=$2
  PYTHONPATH="$PIPELINE_DIR/src${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON" -m integrate_wgs_rna.input_preflight m5 \
    --reference "$fasta" --contigs "$AUTOSOME_REGIONS" --output "$cache" \
    --cache-root "$PREFLIGHT_CACHE_ROOT" --jobs "$PREFLIGHT_JOBS"
}

prepare_vcf_batch() {
  local items="$OUT_DIR/$S.vcf_ref_inputs.tsv"
  : > "$items"
  while (( $# )); do
    printf '%s\t%s\n' "$1" "$2" >> "$items"
    shift 2
  done
  PYTHONPATH="$PIPELINE_DIR/src${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON" -m integrate_wgs_rna.input_preflight vcf \
    --items "$items" --reference "$SEQUENCE_REF" --regions "$AUTOSOME_REGIONS" \
    --bcftools "$BCFTOOLS_BIN" --cache-root "$PREFLIGHT_CACHE_ROOT" \
    --error-dir "$OUT_DIR" --sample "$S" --output "$VCF_BATCH_RESULTS" --jobs "$PREFLIGHT_JOBS"
}

validate_vcf_reference() {
  local vcf=$1 label=$2 line ok cached seconds
  VCF_VALIDATE_DETAIL=
  line=$(awk -F'\t' -v label="$label" '$1==label {print; exit}' "$VCF_BATCH_RESULTS")
  if [[ -z "$line" ]]; then VCF_VALIDATE_DETAIL="VCF 검사 결과 없음: $vcf"; return 1; fi
  IFS=$'\t' read -r _ ok VCF_VALIDATE_DETAIL cached seconds <<< "$line"
  [[ "$ok" == 1 ]]
}

write_preflight_json() {
  local status=$1 json_tmp="$OUT_JSON.part.$$"
  "$PYTHON" - "$OUT" "$json_tmp" "$S" "$DECODE_REF" "$SEQUENCE_REF" "$status" <<'PY'
import csv
import json
import sys

tsv, out, sample, decode_ref, sequence_ref, status = sys.argv[1:]
with open(tsv, newline="") as fh:
    checks = list(csv.DictReader(fh, delimiter="\t"))
with open(out, "w") as fh:
    json.dump({
        "schema": "mfm_preflight_qc_v2",
        "sample": sample,
        "status": status,
        "decode_reference": decode_ref,
        "sequence_reference": sequence_ref,
        "checks": checks,
    }, fh, ensure_ascii=False, indent=2)
    fh.write("\n")
PY
  mv "$json_tmp" "$OUT_JSON"
}

for tool in "$SAMTOOLS_BIN" "$BCFTOOLS_BIN" "$WHATSHAP_BIN"; do
  [[ -x "$tool" ]] && record tool PASS "$tool" || record tool FAIL "$tool"
done
if [[ "$PLATFORM" == illumina ]]; then
  [[ -x "$BWA_MEM2_BIN" ]] && record aligner PASS "$BWA_MEM2_BIN" || record aligner FAIL "$BWA_MEM2_BIN"
else
  [[ -x "$MINIMAP2_BIN" ]] && record aligner PASS "$MINIMAP2_BIN" || record aligner FAIL "$MINIMAP2_BIN"
fi
case "$PLATFORM:$READ_TYPE" in
  illumina:short_read|ont:long_read|hifi:long_read|clr:long_read) platform_contract=PASS ;;
  *) platform_contract=FAIL ;;
esac
case "$ROLE" in pre_training|post_training) role_contract=PASS ;; *) role_contract=FAIL ;; esac
[[ "$platform_contract" == PASS && "$role_contract" == PASS ]] \
  && record manifest_contract PASS "$PLATFORM/$READ_TYPE/$ROLE" \
  || record manifest_contract FAIL "$PLATFORM/$READ_TYPE/$ROLE"
record vcf_source PASS "$VCF_SOURCE"
[[ -s "$ALIGN" ]] && record alignment PASS "$ALIGN" || record alignment FAIL "$ALIGN"
have_index "$ALIGN" && record alignment_index PASS "$ALIGN" || record alignment_index FAIL "$ALIGN"
if [[ -s "$ALIGN" ]] && "$SAMTOOLS_BIN" quickcheck "$ALIGN" 2>/dev/null; then
  record cram_quickcheck PASS "$ALIGN"
else
  record cram_quickcheck FAIL "$ALIGN"
fi

sms=$({ "$SAMTOOLS_BIN" view -H "$ALIGN" 2>/dev/null || true; } | awk -F'\t' '$1=="@RG"{for(i=2;i<=NF;i++)if($i~/^SM:/){sub(/^SM:/,"",$i);print $i}}' | sort -u | paste -sd, -)
allowed_sms="$S"
[[ "$SOURCE_SAMPLE_ID" != "$S" ]] && allowed_sms="$allowed_sms,$SOURCE_SAMPLE_ID"
sms_valid=1
IFS=, read -r -a observed_sms <<< "$sms"
for observed_sm in "${observed_sms[@]}"; do
  [[ -n "$observed_sm" && ",$allowed_sms," == *",$observed_sm,"* ]] || { sms_valid=0; break; }
done
[[ -n "$sms" && "$sms_valid" == 1 ]] \
  && record cram_sample_id PASS "observed=$sms allowed=$allowed_sms" \
  || record cram_sample_id FAIL "observed=${sms:-none} allowed=$allowed_sms"

for ref in "$DECODE_REF" "$SEQUENCE_REF"; do
  [[ -s "$ref" && -s "$ref.fai" ]] && record reference PASS "$ref" || record reference FAIL "$ref"
done
ref_cache_key=$( { printf '%s\n' "$DECODE_REF"; stat -c '%s:%y:%z:%i' "$DECODE_REF" 2>/dev/null || true; sha256sum "$DECODE_REF.fai" 2>/dev/null || true; } \
  | sha256sum | awk '{print $1}')
REF_MD5_CACHE="$STATUS_ROOT/reference_autosomes.${ref_cache_key}.tsv"
if [[ ! -s "$REF_MD5_CACHE" ]]; then
  if build_reference_m5_cache "$DECODE_REF" "$REF_MD5_CACHE"; then
    record reference_autosome_cache PASS "$REF_MD5_CACHE"
  else
    record reference_autosome_cache FAIL "$DECODE_REF"
  fi
else
  record reference_autosome_cache PASS "$REF_MD5_CACHE"
fi

# canonical BAM은 원 BAM/CRAM header의 M5를 보존하지 않을 수 있다. 그 경우 source
# reference의 염기서열을 decode reference와 22개 autosome에서 비교한다. 경로·mtime 비교만으로
# 통과시키지 않아 서로 다른 reference를 사용하는 BAM은 계속 차단된다.
ALIGNMENT_HEADER_HAS_M5=$(awk -F'\t' '$1=="@SQ" {for(i=2;i<=NF;i++) if($i~/^M5:/) {print "1"; exit}}' <<< "$CRAM_HEADER")
SOURCE_REFERENCE_EQUIVALENT=0
SOURCE_REFERENCE_PATH=
if [[ "$ALIGN" == *.bam && "$STAGE04_PROVENANCE" == 1 && -z "$ALIGNMENT_HEADER_HAS_M5" ]]; then
  SOURCE_REFERENCE_PATH=$("$PYTHON" - "$STAGE04_QC" <<'PY'
import json
import sys
try:
    print(json.load(open(sys.argv[1])).get("normalization_reference", {}).get("path", ""))
except Exception:
    pass
PY
)
  if [[ -n "$SOURCE_REFERENCE_PATH" && -s "$SOURCE_REFERENCE_PATH" && -s "$SOURCE_REFERENCE_PATH.fai" ]]; then
    source_ref_key=$( { printf '%s\n' "$SOURCE_REFERENCE_PATH"; stat -c '%s:%y:%z:%i' "$SOURCE_REFERENCE_PATH" 2>/dev/null || true; sha256sum "$SOURCE_REFERENCE_PATH.fai" 2>/dev/null || true; } \
      | sha256sum | awk '{print $1}')
    SOURCE_REF_MD5_CACHE="$STATUS_ROOT/source_reference_autosomes.${source_ref_key}.tsv"
    if [[ ! -s "$SOURCE_REF_MD5_CACHE" ]] && ! build_reference_m5_cache "$SOURCE_REFERENCE_PATH" "$SOURCE_REF_MD5_CACHE"; then
      record bam_source_reference_autosomes FAIL "$SOURCE_REFERENCE_PATH"
    elif [[ -s "$SOURCE_REF_MD5_CACHE" ]]; then
      SOURCE_REFERENCE_EQUIVALENT=1
      for contig in $AUTOSOMES; do
        decode_m5=$(awk -F'\t' -v c="$contig" '$1==c {print $3; exit}' "$REF_MD5_CACHE")
        source_m5=$(awk -F'\t' -v c="$contig" '$1==c {print $3; exit}' "$SOURCE_REF_MD5_CACHE")
        [[ -n "$decode_m5" && "$decode_m5" == "$source_m5" ]] || { SOURCE_REFERENCE_EQUIVALENT=0; break; }
      done
      [[ "$SOURCE_REFERENCE_EQUIVALENT" == 1 ]] \
        && record bam_source_reference_autosomes PASS "$SOURCE_REFERENCE_PATH" \
        || record bam_source_reference_autosomes FAIL "$SOURCE_REFERENCE_PATH"
    fi
  else
    record bam_source_reference_autosomes FAIL "Stage 01 source reference 없음: ${SOURCE_REFERENCE_PATH:-unset}"
  fi
fi
for contig in $AUTOSOMES; do
  decode_len=$(fai_length "$DECODE_REF" "$contig")
  sequence_len=$(fai_length "$SEQUENCE_REF" "$contig")
  cram_len=$(header_sq_value "$contig" LN)
  [[ -n "$cram_len" && "$cram_len" == "$decode_len" && "$decode_len" == "$sequence_len" ]] \
    && record "autosome_length_$contig" PASS "$cram_len" \
    || record "autosome_length_$contig" FAIL "cram=$cram_len decode=$decode_len sequence=$sequence_len"
  ref_m5=$(awk -F'\t' -v c="$contig" '$1==c {print $3; exit}' "$REF_MD5_CACHE" 2>/dev/null || true)
  cram_m5=$(header_sq_value "$contig" M5)
  if [[ -n "$cram_m5" && "$cram_m5" == "$ref_m5" ]]; then
    record "cram_decode_m5_$contig" PASS "$cram_m5"
  elif [[ "$ALIGN" == *.bam && "$STAGE04_PROVENANCE" == 1 && -z "$cram_m5" && "$SOURCE_REFERENCE_EQUIVALENT" == 1 ]]; then
    record "cram_decode_m5_$contig" PASS "BAM M5 없음; source reference autosome M5 검증 통과"
  else
    record "cram_decode_m5_$contig" FAIL "alignment=$cram_m5 decode=$ref_m5"
  fi
done

qual_metrics=$({ "$SAMTOOLS_BIN" view -T "$DECODE_REF" -F 4 "$ALIGN" 2>/dev/null || true; } | awk 'BEGIN{n=miss=bad=0} NR<=1000{n++; if($10=="*"||$11=="*")miss++; else if(length($10)!=length($11))bad++} NR==1000{exit} END{print n,miss,bad}')
read -r qual_reads qual_missing qual_bad <<< "$qual_metrics"
[[ "$qual_reads" -gt 0 && "$qual_missing" == 0 && "$qual_bad" == 0 ]] \
  && record base_qual PASS "reads=$qual_reads missing=$qual_missing length_mismatch=$qual_bad" \
  || record base_qual FAIL "reads=$qual_reads missing=$qual_missing length_mismatch=$qual_bad"

if [[ "$PLATFORM" == illumina ]]; then
  for region in chr1:10000-20000 chr2:38143382-38182146 chr21:5010012-5072204; do
    label=${region/:/_}
    if "$SAMTOOLS_BIN" view -T "$DECODE_REF" "$ALIGN" "$region" >/dev/null 2>"$OUT_DIR/$S.$label.decode.err"; then
      record "cram_decode_$label" PASS "$region"
    else
      detail=$(head -1 "$OUT_DIR/$S.$label.decode.err" 2>/dev/null || true)
      record "cram_decode_$label" FAIL "$region $detail"
    fi
  done
fi

if [[ "$VCF_SOURCE" == illumina_panel ]]; then
  # Stage 03과 같은 template resolver를 사용한다. {chrom}/%s/custom 이름 지원.
  # Python을 염색체마다 다시 시작하지 않고 22개 이름을 한 번에 해석한다.
  PANEL_NAMES=$(PYTHONPATH="$PIPELINE_DIR/src${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON" - "$PANEL_FILE_TEMPLATE" <<'PY'
import sys
from integrate_wgs_rna.cohort_vcf import panel_source
for chrom in range(1, 23):
    print(panel_source("", sys.argv[1], f"chr{chrom}"))
PY
  ) || die "panel filename template 해석 실패: $PANEL_FILE_TEMPLATE"
  mapfile -t PANEL_NAMES_ARRAY <<< "$PANEL_NAMES"
  vcf_fail=0
  for chrom in $(seq 1 22); do
    name=${PANEL_NAMES_ARRAY[chrom-1]}
    vcf="$PANEL_ROOT/$name"
    [[ -s "$vcf" && -s "$vcf.tbi" ]] || { vcf_fail=1; break; }
  done
  [[ "$vcf_fail" == 0 ]] && record panel_files PASS "$PANEL_ROOT" || record panel_files FAIL "local=$PANEL_ROOT download_source=$PANEL_SOURCE_URL (자동 다운로드 없음)"
  if [[ "$vcf_fail" == 0 ]]; then
    VCF_BATCH_ARGS=()
    for chrom in $(seq 1 22); do
      VCF_BATCH_ARGS+=("$PANEL_ROOT/${PANEL_NAMES_ARRAY[chrom-1]}" "panel_chr$chrom")
    done
    prepare_vcf_batch "${VCF_BATCH_ARGS[@]}"
    for chrom in $(seq 1 22); do
      name=${PANEL_NAMES_ARRAY[chrom-1]}
      vcf="$PANEL_ROOT/$name"
      if validate_vcf_reference "$vcf" "panel_chr$chrom"; then
        record "vcf_reference_chr$chrom" PASS "$VCF_VALIDATE_DETAIL"
      else
        record "vcf_reference_chr$chrom" FAIL "$VCF_VALIDATE_DETAIL"
      fi
    done
  else
    record vcf_reference_autosomes FAIL "panel VCF 파일 또는 index 누락"
  fi
  mkdir -p "$COMMON"
  [[ -d "$COMMON" && -w "$COMMON" ]] && record common_vcf_output PASS "$COMMON" || record common_vcf_output FAIL "$COMMON"
  panel_sample_fail=0
  for chrom in $(seq 1 22); do
    name=${PANEL_NAMES_ARRAY[chrom-1]}
    if "$BCFTOOLS_BIN" query -l "$PANEL_ROOT/$name" 2>/dev/null \
        | awk -v s="$S" '$1==s{found=1}END{exit(found?0:1)}'; then
      record "panel_sample_chr$chrom" PASS "$S"
    else
      record "panel_sample_chr$chrom" FAIL "$S"
      panel_sample_fail=1
    fi
  done
  [[ "$panel_sample_fail" == 0 ]] && record panel_sample_autosomes PASS "$S" \
    || record panel_sample_autosomes FAIL "$S"
elif [[ "$VCF_SOURCE" == public ]]; then
  mapfile -t PUBLIC_VCFS < <(sample_public_vcfs "$S")
  PUBLIC_SAMPLE=$(sample_public_vcf_sample_id "$S")
  PUBLIC_REF=$(sample_vcf_reference "$S")
  PUBLIC_CALLSET=$(sample_vcf_callset "$S")
  PHASE_SOURCE_RESOLVED=$(sample_phase_source "$S")
  [[ "${#PUBLIC_VCFS[@]}" -gt 0 ]] && record public_vcf_count PASS "${#PUBLIC_VCFS[@]}" || record public_vcf_count FAIL 0
  public_sample_fail=0
  VCF_BATCH_ARGS=()
  for PUBLIC_VCF in "${PUBLIC_VCFS[@]}"; do
    VCF_BATCH_ARGS+=("$PUBLIC_VCF" "$(basename "$PUBLIC_VCF")")
  done
  if [[ "${#PUBLIC_VCFS[@]}" -gt 0 ]]; then
    "$PYTHON" -m handoff.vcf_indexes --bcftools "$BCFTOOLS_BIN" "${PUBLIC_VCFS[@]}"
  fi

  prepare_vcf_batch "${VCF_BATCH_ARGS[@]}"
  for PUBLIC_VCF in "${PUBLIC_VCFS[@]}"; do
    [[ -s "$PUBLIC_VCF" ]] && record public_vcf PASS "$PUBLIC_VCF" || record public_vcf FAIL "$PUBLIC_VCF"
    [[ -s "$PUBLIC_VCF.tbi" || -s "$PUBLIC_VCF.csi" ]] && record public_vcf_index PASS "$PUBLIC_VCF" || record public_vcf_index FAIL "$PUBLIC_VCF"
    if validate_vcf_reference "$PUBLIC_VCF" "$(basename "$PUBLIC_VCF")"; then
      record public_vcf_ref_alleles PASS "$VCF_VALIDATE_DETAIL"
    else
      record public_vcf_ref_alleles FAIL "$VCF_VALIDATE_DETAIL"
    fi
    if [[ ! -s "$PUBLIC_VCF" ]] || ! "$BCFTOOLS_BIN" query -l "$PUBLIC_VCF" 2>/dev/null | awk -v s="$PUBLIC_SAMPLE" '$1==s{found=1}END{exit(found?0:1)}'; then
      public_sample_fail=1
    fi
  done
  [[ -n "$PUBLIC_CALLSET" && "$PUBLIC_CALLSET" != unknown ]] && record public_vcf_callset PASS "$PUBLIC_CALLSET" || record public_vcf_callset FAIL "$PUBLIC_CALLSET"
  [[ "$PHASE_SOURCE_RESOLVED" == provided ]] && record public_phase_source PASS "$PHASE_SOURCE_RESOLVED" || record public_phase_source FAIL "$PHASE_SOURCE_RESOLVED"
  [[ -s "$PUBLIC_REF" && -s "$PUBLIC_REF.fai" ]] && record public_vcf_reference_file PASS "$PUBLIC_REF" || record public_vcf_reference_file FAIL "$PUBLIC_REF"
  if [[ "$public_sample_fail" == 0 && "${#PUBLIC_VCFS[@]}" -gt 0 ]]; then
    record public_vcf_sample PASS "$PUBLIC_SAMPLE"
  else
    record public_vcf_sample FAIL "$PUBLIC_SAMPLE"
  fi
  mkdir -p "$PERSAMPLE"
  [[ -d "$PERSAMPLE" && -w "$PERSAMPLE" ]] && record public_vcf_output PASS "$PERSAMPLE" || record public_vcf_output FAIL "$PERSAMPLE"
elif [[ "$VCF_SOURCE" == clair3 ]]; then
  MODEL_PATH=$(clair3_model_path_for_sample "$S")
  [[ -d "$MODEL_PATH" ]] && record clair3_model PASS "$MODEL_PATH" || record clair3_model FAIL "$MODEL_PATH"
  record vcf_reference_deferred PASS "Clair3 VCF는 stage 04 생성 직후 stage 09 전에 동일 REF allele 검사를 수행"
else
  record vcf_source_supported FAIL "$VCF_SOURCE"
fi

# RNA pair는 WGS가 끝난 뒤가 아니라 Stage 02에서 원본 자체를 먼저 검사한다.
# 개인 HAP FASTA/chain은 아직 생성 전이므로 여기서는 source RNA·index·GTF·header·
# pair 계약을 동결하고, Stage 13 진입 시 같은 contract와 새 HAP 자산을 재검증한다.
RNA=$(sample_rna "$S" 2>/dev/null || true)
RNA_PREFLIGHT_QC="$(sample_final_root "$S")/$S/qc/rna_preflight.json"
if sample_rna_deferred "$S"; then
  record rna_pair_status PASS "WGS_RNA_PAIR_RNA_DEFERRED"
  record rna_input_contract PASS "DEFERRED"
elif [[ -n "$RNA" && "$RNA" != "-" ]]; then
  mkdir -p "$(dirname "$RNA_PREFLIGHT_QC")"
  PUBLIC_ID=$(sample_public_vcf_sample_id "$S" 2>/dev/null || true)
  SOURCE_ID=$(sample_source_sample_id "$S" 2>/dev/null || true)
  EXPECTED_SM="$S"
  [[ -n "$PUBLIC_ID" && "$PUBLIC_ID" != "$S" ]] && EXPECTED_SM="$S,$PUBLIC_ID"
  [[ -n "$SOURCE_ID" && "$SOURCE_ID" != "$S" && "$SOURCE_ID" != "$PUBLIC_ID" ]] && EXPECTED_SM="$EXPECTED_SM,$SOURCE_ID"
  RNA_DECODE_REF=$(named_path_or_default "$S" rna_reference_fasta "${RNA_DECODE_REFERENCE:-}")
  RNA_REF_ARGS=()
  [[ -n "$RNA_DECODE_REF" ]] && RNA_REF_ARGS=(--decode-reference "$RNA_DECODE_REF")
  RNA_SM_ARGS=(--require-sm-match)
  [[ "$RNA_REQUIRE_SM_MATCH" == 1 ]] || RNA_SM_ARGS=(--no-require-sm-match)
  RNA_ADDNL_READY=1
  if [[ "$RNA_INDEX_ADDNL_SEQNS" == 1 ]]; then
    if [[ ! -s "$RNA_ADDNL_SEQNS_SOURCE" || ! -s "$RNA_ADDNL_SEQNS_SOURCE.fai" || ! -s "$RNA_ADDNL_SEQNS_FASTA" ]]; then
      record rna_additional_reference FAIL "FASTA/FAI/additional FASTA 필요: source=$RNA_ADDNL_SEQNS_SOURCE fai=$RNA_ADDNL_SEQNS_SOURCE.fai addnl=$RNA_ADDNL_SEQNS_FASTA"
      RNA_ADDNL_READY=0
    elif awk -F'\t' 'NF >= 5 && $1 != "" && $2 ~ /^[0-9]+$/ && $2 > 0 {ok=1} END {exit !ok}' "$RNA_ADDNL_SEQNS_SOURCE.fai"; then
      record rna_additional_reference PASS "$RNA_ADDNL_SEQNS_SOURCE"
    else
      record rna_additional_reference FAIL "invalid FASTA index: $RNA_ADDNL_SEQNS_SOURCE.fai"
      RNA_ADDNL_READY=0
    fi
  else
    record rna_additional_reference PASS "RNA_INDEX_ADDNL_SEQNS=0"
  fi
  if [[ "$RNA_ADDNL_READY" == 1 ]] && "$PYTHON" ./rna_preflight.py \
      --sample "$S" --expected-samples "$EXPECTED_SM" --rna "$RNA" \
      --quickcheck-json "$(sample_final_root "$S")/$S/qc/stage01_rna_quickcheck.json" \
      --gtf "$GENCODE_GTF" --contigs "$(chroms_space "$CHROMS")" \
      --samtools "$SAMTOOLS_BIN" --hap-read-mode "$HAP_READ_MODE" \
      --threads "$RNA_BAM_INDEX_THREADS" \
      --min-paired-fraction "$RNA_MIN_PAIRED_FRACTION" \
      --ratio-qc-policy "$QC_RATIO_POLICY" \
      --sm-placeholders "${RNA_SM_PLACEHOLDERS:-sm1}" --skip-hap-assets \
      "${RNA_REF_ARGS[@]}" "${RNA_SM_ARGS[@]}" --output "$RNA_PREFLIGHT_QC" \
      > "$OUT_DIR/$S.rna_preflight.log" 2>&1; then
    record rna_pair_status PASS "WGS_RNA_PAIR"
    record rna_input_contract PASS "$RNA_PREFLIGHT_QC"
    record rna_quickcheck PASS "$RNA"
  else
    if [[ "$RNA_ADDNL_READY" != 1 ]]; then
      detail="RNA 추가 reference 검사 실패; rna_additional_reference 항목 참조"
    elif [[ -f "$OUT_DIR/$S.rna_preflight.log" ]]; then
      detail=$(tr '\n' ' ' < "$OUT_DIR/$S.rna_preflight.log" | cut -c1-500)
    else
      detail="RNA preflight 실패; 로그가 생성되지 않음"
    fi
    record rna_pair_status FAIL "manifest RNA input failed"
    record rna_input_contract FAIL "$detail"
  fi
else
  record rna_pair_status PASS "WGS_ONLY"
  record rna_input_contract PASS "NOT_APPLICABLE"
fi

mkdir -p "$OUTPUT_ROOT"
free_kib=$(df -Pk "$OUTPUT_ROOT" | awk 'NR==2{print $4}')
# 고정 운영 예산은 기본 미적용이다. 가용량 관찰과 실제 I/O 오류는 유지한다.
# 사용자가 필요할 때만 PREFLIGHT_MIN_FREE_GIB로 명시적인 최소값을 설정할 수 있다.
min_free_gib=${PREFLIGHT_MIN_FREE_GIB:-0}
[[ "$min_free_gib" =~ ^[0-9]+$ ]] || die "PREFLIGHT_MIN_FREE_GIB는 0 이상의 정수여야 함"
min_kib=$(( min_free_gib * 1024 * 1024 ))
if [[ "$min_kib" == 0 ]]; then
  record free_space PASS "policy=observe_only free_kib=$free_kib min_kib=0"
else
  [[ "$free_kib" -ge "$min_kib" ]] && record free_space PASS "policy=explicit_user_minimum free_kib=$free_kib min_kib=$min_kib" || record free_space FAIL "policy=explicit_user_minimum free_kib=$free_kib min_kib=$min_kib"
fi

mv "$TMP" "$OUT"
if [[ "$FAIL" == 0 ]]; then
  write_preflight_json pass
else
  write_preflight_json fail
fi
# 상태 root가 cleanup으로 사라져도 입력 계약을 최종 package에 보존한다.
FINAL_QC_DIR="$(sample_final_root "$S")/$S/qc"
mkdir -p "$FINAL_QC_DIR"
cp -a "$OUT_JSON" "$FINAL_QC_DIR/wgs_preflight.json"
cp -a "$OUT" "$FINAL_QC_DIR/wgs_preflight.tsv"
if [[ "$FAIL" != 0 ]]; then
  echo "[FATAL] preflight 실패: sample=$S report=$OUT qc=$OUT_JSON" >&2
  exit 1
fi
echo "[OK] preflight 통과: sample=$S report=$OUT qc=$OUT_JSON" >&2
