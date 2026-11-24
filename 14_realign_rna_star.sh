#!/usr/bin/env bash
# 구번호 RNA 통합 STAR 실행기. 현재 v2 runner는 Stage 16 STAR 정렬과 Stage 17 MarkDuplicates를 사용한다.
# 이 파일은 보존된 이전 실행 경로이며 새 전처리 entrypoint가 아니다.
#
# 14 — RNA 리드를 개인 유전체 H1/H2 에 STAR 로 재정렬하고 duplicate flag를 기록한다.
#
#   ./14_realign_rna_star.sh <sample>
#
#  새 규칙: *"wgs 에서 만든 개인 유전체 h1, h2 에다가 다시 rna-seq bam 의 각 read 를 star 로
#  alignment 해서 h1, h2 마다의 각 1bp 별 rna-seq signal 을 만들면 됨. alphagenome 의 전처리
#  방법을 따르면 돼."*
#
#  ★ 산출물의 정체를 오해하지 말 것 (파일럿에서 확인)
#     H1 에 정렬한 커버리지는 **H1 좌표에서 측정한 총 발현**이다 — haplotype 특이 발현이 아니다.
#     H2 분자에서 나온 리드도 H1 에 잘 붙는다(두 haplotype 은 위치의 0.1% 에서만 다르다).
#     haplotype 특이 신호는 두 정렬의 **리드별 AS 비교**로만 나온다(→ 16_allelic_counts.py).
#
#  ⚠ index 는 haplotype 마다 새로 만든다(STAR 는 증분 index 가 없다).
#    파일럿 실측(HG00149 HAP1, 2.8 Gb): 4스레드 25분 43초 / index 26 GB → 정렬 직후 삭제한다.
#  ⚠ splice 주석은 **개인 좌표로 lift 한 GTF** 를 쓴다(13_lift_annotation.py).
#    참조 좌표 GTF 를 그대로 주면 indel 만큼 어긋난 junction 이 삽입된다.
#
#  ★ 2026-09-17 AlleleSeq2 대조 반영
#    --scoreGenomicLengthLog2scale 0 을 명시한다(STAR 기본 -0.25). AlleleSeq2 는 4개 모드
#    전부에 0.0 을 준다. HAP 간 indel 로 genomic span 이 달라지면 염기 일치와 무관한 AS 차이가
#    생기는데, 16 은 AS 1 점 차이로 배정하므로 이 항이 곧바로 배정을 바꾼다.
#    ⚠ --sjdbScore 0 은 넣지 않는다. AlleleSeq2 가 그걸 주는 것은 splicing 이 없는
#      ASB/ASCA(DNA) 모드뿐이고, RNA 용 ASE-params 는 기본값을 쓴다.
#    index 에는 ADDNL_SEQNS(unplaced/unlocalized/EBV/decoy)를 함께 넣어 갈 곳 없는 read 를
#    흡수시킨다. 카운트 대상 contig 는 CHROMS 그대로다.
set -uo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
RNA=$(v2_rna "$S") || die "$S 매니페스트에 없다"
[[ "$RNA" != "-" && -s "$RNA" ]] || die "$S RNA BAM 없음: $RNA"
F=$(v2_final "$S"); W=$(v2_work "$S"); mkdir -p "$W/logs" "$W/rna"
MTAG=$(sample_rna_evidence_label "$S")
RNA_CONTIGS=$(chroms_space "$CHROMS")
[[ -n "$RNA_CONTIGS" ]] || die "$S RNA contig 목록이 비어 있음: CHROMS=$CHROMS"
require_hap_contig_assets "$F/individual_fasta" "$CHROMS" 1
RNA_PREFLIGHT_QC="$F/qc/rna_preflight.json"
require_file "$RNA_PREFLIGHT_QC"
RNA_CONTRACT_SHA256=$("$PYTHON" -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d.get("status")=="pass"; print(d["contract_sha256"])' "$RNA_PREFLIGHT_QC") || die "$S RNA preflight contract 읽기 실패"

count_fastq_records() {
  local fastq=$1
  [[ -e "$fastq" && -s "$fastq" ]] || { printf '0\n'; return 0; }
  awk 'END { if (NR % 4) exit 2; print NR / 4 }' "$fastq"
}

contract_matches() {
  local path=$1
  [[ -s "$path" ]] && "$PYTHON" -c 'import json,sys; d=json.load(open(sys.argv[1])); raise SystemExit(0 if d.get("rna_preflight_contract_sha256")==sys.argv[2] else 1)' "$path" "$RNA_CONTRACT_SHA256"
}

metric_value() {
  local metrics=$1 key=$2
  awk -F'\t' -v key="$key" '
    /^#/ || NF == 0 { next }
    !header {
      for (i=1; i<=NF; i++) col[$i]=i
      header=1
      next
    }
    key in col { print $(col[key]); exit }
  ' "$metrics"
}

metric_number() {
  local value=${1:-0}
  [[ "$value" =~ ^[0-9]+([.][0-9]+)?$ ]] && printf '%s\n' "$value" || printf '0\n'
}

# ── ① BAM → FASTQ (paired 보존) ──
#  ⚠ 원본은 coordinate-sorted 다. collate 로 이름 순서를 맞추지 않으면 페어가 깨진다.
#  ⚠ 파일럿 실측: 이 단계가 정렬보다 6배 오래 걸린다(14분 vs 2.3분). 병목이다.
R1=$W/rna/R1.fq; R2=$W/rna/R2.fq
FASTQ_QC=$W/rna/fastq_qc.json
UNPAIRED=$W/rna/discard_unpaired.fq
SINGLETON=$W/rna/discard_singleton.fq
if [[ "${RNA_FORCE:-0}" == 1 || ! -s "$R1" || ! -s "$R2" ]] || ! contract_matches "$FASTQ_QC"; then
  log "$S: RNA BAM → FASTQ"
  rm -f "$R1.part" "$R2.part" "$UNPAIRED.part" "$SINGLETON.part"
  "$SAMTOOLS" collate -u -O -@ "$RNA_COLLATE_THREADS" "$RNA" "$W/rna/collate_tmp" 2>/dev/null \
    | "$SAMTOOLS" fastq -@ "$RNA_FASTQ_THREADS" -1 "$R1.part" -2 "$R2.part" \
        -0 "$UNPAIRED.part" -s "$SINGLETON.part" -n - \
        2>"$W/logs/14.fastq.err"
  st=("${PIPESTATUS[@]}"); [[ "${st[0]}" == 0 && "${st[1]}" == 0 ]] || die "$S FASTQ 추출 실패"
  r1_pairs=$(count_fastq_records "$R1.part") || die "$S R1 FASTQ 레코드 형식 오류"
  r2_pairs=$(count_fastq_records "$R2.part") || die "$S R2 FASTQ 레코드 형식 오류"
  singleton_discard=$(count_fastq_records "$SINGLETON.part") || die "$S singleton FASTQ 레코드 형식 오류"
  unpaired_discard=$(count_fastq_records "$UNPAIRED.part") || die "$S unpaired FASTQ 레코드 형식 오류"
  [[ "$r1_pairs" == "$r2_pairs" && "$r1_pairs" -gt 0 ]] \
    || die "$S FASTQ pair 계약 위반: R1=$r1_pairs R2=$r2_pairs"
  mv "$R1.part" "$R1"; mv "$R2.part" "$R2"
  "$PYTHON" - "$FASTQ_QC.part" "$RNA" "$R1" "$R2" "$r1_pairs" "$r2_pairs" \
    "$singleton_discard" "$unpaired_discard" "$RNA_CONTRACT_SHA256" <<'PY'
import json
import sys
from pathlib import Path

out, alignment, r1, r2, r1_pairs, r2_pairs, singleton, unpaired, contract_sha256 = sys.argv[1:]
Path(out).write_text(json.dumps({
    "rna_preflight_contract_sha256": contract_sha256,
    "schema": "mfm_rna_fastq_qc_v1",
    "input_alignment": alignment,
    "r1_fastq": r1,
    "r2_fastq": r2,
    "r1_pairs": int(r1_pairs),
    "r2_pairs": int(r2_pairs),
    "singleton_discard_records": int(singleton),
    "unpaired_discard_records": int(unpaired),
    "pair_contract_pass": r1_pairs == r2_pairs,
}, ensure_ascii=False, indent=2) + "\n")
PY
  mv "$FASTQ_QC.part" "$FASTQ_QC"
  rm -f "$UNPAIRED.part" "$SINGLETON.part"
fi

for h in 1 2; do
  TD=$F/rna_personal/$MTAG/HAP$h; mkdir -p "$TD"
  COV=$TD/coverage.zarr
  FINAL_BAM=$TD/$S.rna.personal.HAP$h.bam
  HD=$F/individual_fasta/HAP$h
  GTF=$TD/annotation.gtf
  LIFT_QC=$TD/gtf_lift_qc.json
  EXONS=$TD/exons.bygene.tsv
  # Stage 13과 14의 경계를 강제한다. 독립 재시작 시 이전 단계가 없으면 묵시적으로
  # 다시 계산하지 않고 즉시 실패해야 잘못된 좌표계의 GTF가 STAR에 들어가지 않는다.
  require_file "$GTF"
  require_file "$EXONS"
  require_file "$LIFT_QC"
  if [[ "${RNA_FORCE:-0}" != 1 && -s "$FINAL_BAM" && -s "$FINAL_BAM.bai" && -s "$TD/rna_realign.json" ]] \
      && grep -q '"schema": "mfm_v2_rna_realign_v3"' "$TD/rna_realign.json" \
      && contract_matches "$TD/rna_realign.json"; then
    log "  [skip] HAP$h"
    continue
  fi

  # ── target concat + STAR index ──
  TFA=$W/rna/HAP$h.fa
  if [[ ! -s "$TFA.fai" ]]; then
    rm -f "$TFA"
    for c in $RNA_CONTIGS; do cat "$HD/$c.fa" >> "$TFA.part" || die "$S HAP$h $c FASTA concat 실패"; done
    mv "$TFA.part" "$TFA" || die "$S HAP$h target FASTA 이동 실패"
    "$SAMTOOLS" faidx "$TFA" || die "$S HAP$h target FASTA index 실패"
  fi
  # ── ③-b ADDNL_SEQNS — index 에만 넣는 흡수용 contig (AlleleSeq2 makePersonalGenome.mk) ──
  #   unplaced/unlocalized/EBV/decoy 를 index 에 함께 넣어 autosome 에 갈 곳 없는 read 가
  #   억지로 붙는 것을 막는다. alt/HLA 는 제외한다(AlleleSeq2 도 제외한다).
  #   exon 주석이 없으므로 16 의 gene 카운트에는 들어가지 않고, 15 는 --contigs 로 제한된다.
  ADDNL_ARGS=()
  if [[ "$RNA_INDEX_ADDNL_SEQNS" == 1 ]]; then
    mkdir -p "$(dirname "$RNA_ADDNL_SEQNS_FASTA")"
    exec {ADDNL_LOCK_FD}>"$RNA_ADDNL_SEQNS_FASTA.lock"
    flock "$ADDNL_LOCK_FD" || die "ADDNL_SEQNS lock 획득 실패"
    if [[ ! -s "$RNA_ADDNL_SEQNS_FASTA" ]]; then
      [[ -s "$RNA_ADDNL_SEQNS_SOURCE" && -s "$RNA_ADDNL_SEQNS_SOURCE.fai" ]] \
        || die "ADDNL_SEQNS 원본 또는 .fai 없음: $RNA_ADDNL_SEQNS_SOURCE (RNA_INDEX_ADDNL_SEQNS=0으로 끌 수 있음)"
      log "  ADDNL_SEQNS 생성(최초 1회): $RNA_ADDNL_SEQNS_FASTA"
      mapfile -t ADDNL_CONTIGS < <(
        awk '$1 ~ /^chrUn_/ || $1 ~ /_random$/ || $1 == "chrEBV" {print $1}' "$RNA_ADDNL_SEQNS_SOURCE.fai")
      [[ "${#ADDNL_CONTIGS[@]}" -gt 0 ]] || die "ADDNL_SEQNS 대상 contig를 찾지 못함: $RNA_ADDNL_SEQNS_SOURCE.fai"
      ADDNL_TMP="$RNA_ADDNL_SEQNS_FASTA.part.$$"
      "$SAMTOOLS" faidx "$RNA_ADDNL_SEQNS_SOURCE" "${ADDNL_CONTIGS[@]}" \
        > "$ADDNL_TMP" || die "ADDNL_SEQNS 추출 실패"
      mv "$ADDNL_TMP" "$RNA_ADDNL_SEQNS_FASTA" || die "ADDNL_SEQNS 이동 실패"
      log "  ADDNL_SEQNS contig ${#ADDNL_CONTIGS[@]}개"
    fi
    flock -u "$ADDNL_LOCK_FD"
    exec {ADDNL_LOCK_FD}>&-
    ADDNL_ARGS=("$RNA_ADDNL_SEQNS_FASTA")
  fi

  IDX=$W/rna/idx_hap$h
  if [[ ! -s "$IDX/SAindex" ]]; then
    log "  HAP$h: STAR index (파일럿 실측 4스레드 26분 / 26 GB)"
    mkdir -p "$IDX"
    "$STAR_BIN" --runMode genomeGenerate --runThreadN "$STAR_INDEX_THREADS" \
      --genomeDir "$IDX" --genomeFastaFiles "$TFA" ${ADDNL_ARGS[@]+"${ADDNL_ARGS[@]}"} \
      --sjdbGTFfile "$GTF" --sjdbOverhang "$STAR_SJDB_OVERHANG" \
      --outFileNamePrefix "$IDX/" > "$W/logs/14.index.HAP$h.log" 2>&1 \
      || die "$S HAP$h STAR index 실패: $(tail -3 "$W/logs/14.index.HAP$h.log" | tr '\n' ' ')"
  fi

  # ── ④ 정렬 (TOPMed/GTEx와 동일한 filter·junction 기준, 아래 JSON에 차이 명시) ──
  AD=$W/rna/aln_hap$h; mkdir -p "$AD"
  STAR_BAM=$AD/Aligned.sortedByCoord.out.bam
  if [[ ! -s "$STAR_BAM" ]]; then
    log "  HAP$h: STAR 정렬"
    "$STAR_BIN" --runMode alignReads --runThreadN "$STAR_ALIGN_THREADS" \
      --genomeDir "$IDX" --readFilesIn "$R1" "$R2" --readFilesCommand cat \
      --outFileNamePrefix "$AD/" \
      --twopassMode Basic --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI AS nM NM --outSAMattrRGline "ID:${S}.HAP${h}" "SM:${S}" \
      --outSAMunmapped Within --outSAMstrandField intronMotif \
      --outFilterType BySJout --outFilterMultimapNmax 20 \
      --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 \
      --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 \
      --scoreGenomicLengthLog2scale 0 \
      --limitSjdbInsertNsj 1200000 --outFilterIntronMotifs None \
      --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
      --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
      --alignSoftClipAtReferenceEnds Yes --genomeLoad NoSharedMemory \
      --quantMode GeneCounts > "$W/logs/14.align.HAP$h.log" 2>&1 \
      || die "$S HAP$h STAR 정렬 실패"
  fi
  cp -f "$AD/Log.final.out" "$TD/Log.final.out" 2>/dev/null
  cp -f "$AD/ReadsPerGene.out.tab" "$TD/ReadsPerGene.out.tab" 2>/dev/null

  # ── ⑤ MarkDuplicates — TOPMed와 같이 duplicate를 제거하지 않고 flag/metrics를 남긴다. ──
  MD_BAM=$AD/Aligned.sortedByCoord.out.md.bam
  MD_METRICS=$AD/$S.HAP$h.markduplicates.metrics.txt
  [[ -x "$PICARD_BIN" ]] || die "$S HAP$h Picard 없음: $PICARD_BIN (./00_environment.sh install-rna 또는 PICARD_BIN 지정 필요)"
  [[ "$RNA_MARKDUP_HEAP_GB" =~ ^[1-9][0-9]*$ ]] || die "RNA_MARKDUP_HEAP_GB는 양의 정수여야 함: $RNA_MARKDUP_HEAP_GB"
  [[ "$RNA_MARKDUP_REMOVE_DUPLICATES" == 0 ]] || die "RNA_MARKDUP_REMOVE_DUPLICATES는 0만 허용: flag만 표기해야 함"
  if [[ "${RNA_FORCE:-0}" == 1 || ! -s "$MD_BAM" || ! -s "$MD_BAM.bai" || ! -s "$MD_METRICS" ]]; then
    log "  HAP$h: Picard MarkDuplicates"
    PICARD_TMP=$W/rna/picard_tmp_hap$h; mkdir -p "$PICARD_TMP"
    # -Djava.io.tmpdir 는 TMP_DIR 와 다르다. TMP_DIR 은 Picard 가 정렬 spill 을 쓰는 곳이고,
    # java.io.tmpdir 은 JVM 이 snappy 같은 native .so 를 풀어 mmap 하는 곳이다.
    # 이 호스트의 /tmp 는 noexec 라 후자를 옮기지 않으면 실행 즉시
    # UnsatisfiedLinkError(failed to map segment from shared object)로 죽는다.
    "$PICARD_BIN" "-Xmx${RNA_MARKDUP_HEAP_GB}g" "-Djava.io.tmpdir=$PICARD_TMP" MarkDuplicates \
      I="$STAR_BAM" O="$MD_BAM.part.bam" M="$MD_METRICS.part.txt" \
      ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
      REMOVE_DUPLICATES=false TMP_DIR="$PICARD_TMP" \
      > "$W/logs/14.markduplicates.HAP$h.log" 2>&1 \
      || die "$S HAP$h Picard MarkDuplicates 실패"
    "$SAMTOOLS" index -@ "$RNA_BAM_INDEX_THREADS" "$MD_BAM.part.bam" || die "$S HAP$h MarkDuplicates BAM index 실패"
    mv "$MD_BAM.part.bam" "$MD_BAM" || die "$S HAP$h MarkDuplicates BAM 이동 실패"
    mv "$MD_BAM.part.bam.bai" "$MD_BAM.bai" || die "$S HAP$h MarkDuplicates BAI 이동 실패"
    mv "$MD_METRICS.part.txt" "$MD_METRICS" || die "$S HAP$h MarkDuplicates metrics 이동 실패"
  fi
  cp -f "$MD_BAM" "$FINAL_BAM" || die "$S HAP$h final BAM 복사 실패"
  cp -f "$MD_BAM.bai" "$FINAL_BAM.bai" || die "$S HAP$h final BAM index 복사 실패"

  uniq=$(awk -F'|' '/Uniquely mapped reads %/{gsub(/[ \t%]/,"",$2); print $2}' "$TD/Log.final.out")
  pairs_examined=$(metric_number "$(metric_value "$MD_METRICS" READ_PAIRS_EXAMINED)")
  pair_duplicates=$(metric_number "$(metric_value "$MD_METRICS" READ_PAIR_DUPLICATES)")
  unpaired_examined=$(metric_number "$(metric_value "$MD_METRICS" UNPAIRED_READS_EXAMINED)")
  unpaired_duplicates=$(metric_number "$(metric_value "$MD_METRICS" UNPAIRED_READ_DUPLICATES)")
  percent_duplication=$(metric_number "$(metric_value "$MD_METRICS" PERCENT_DUPLICATION)")
  duplicate_primary=$($SAMTOOLS view -c -F 2308 -f 1024 "$MD_BAM") || die "$S HAP$h duplicate flag count 실패"
  primary_records=$($SAMTOOLS view -c -F 2308 "$MD_BAM") || die "$S HAP$h primary record count 실패"
  "$PYTHON" - "$TD/rna_realign.json.part" "$FASTQ_QC" "$LIFT_QC" "$S" "HAP$h" "$MTAG" \
    "$STAR_BIN" "$STAR_SJDB_OVERHANG" "$uniq" "$COV" "$MD_BAM" "$MD_METRICS" \
    "$pairs_examined" "$pair_duplicates" "$unpaired_examined" "$unpaired_duplicates" \
    "$percent_duplication" "$duplicate_primary" "$primary_records" "$RNA_CONTIGS" \
    "$RNA_INDEX_ADDNL_SEQNS" "$RNA_CONTRACT_SHA256" <<'PY'
import json
import subprocess
import sys
from pathlib import Path

(out, fastq_qc, lift_qc, sample, hap, phase_mode, star_bin, sjdb_overhang, unique_pct,
 coverage_zarr, marked_bam, metrics_path, pairs_examined, pair_duplicates, unpaired_examined,
 unpaired_duplicates, percent_duplication, duplicate_primary, primary_records, contigs,
 addnl_seqns_flag, contract_sha256) = sys.argv[1:]
addnl_seqns = addnl_seqns_flag == "1"
fastq = json.loads(Path(fastq_qc).read_text())
lift = json.loads(Path(lift_qc).read_text())
payload = {
    "rna_preflight_contract_sha256": contract_sha256,
    "schema": "mfm_v2_rna_realign_v3",
    "sample": sample,
    "hap": hap,
    "phase_mode": phase_mode,
    "aligner": f"STAR {subprocess.getoutput(f'{star_bin} --version')}",
    "sjdb_overhang": int(sjdb_overhang),
    "fastq_reconstruction": fastq,
    "gtf_lift_qc": lift,
    "markduplicates": {
        "tool": "Picard MarkDuplicates",
        "remove_duplicates": False,
        "marked_bam": marked_bam,
        "metrics_path": metrics_path,
        "read_pairs_examined": int(float(pairs_examined)),
        "read_pair_duplicates": int(float(pair_duplicates)),
        "unpaired_reads_examined": int(float(unpaired_examined)),
        "unpaired_read_duplicates": int(float(unpaired_duplicates)),
        "duplicate_read_records_from_metrics": int(2 * float(pair_duplicates) + float(unpaired_duplicates)),
        "primary_duplicate_records": int(float(duplicate_primary)),
        "primary_records": int(float(primary_records)),
        "percent_duplication": float(percent_duplication),
    },
    "alleleseq2_parameter_policy": {
        "reference": "gersteinlab/AlleleSeq2 1f3d424 (= trgaleev/AlleleSeq2 cfe8acf, 내용 동일) PIPELINE.mk ASE-params",
        "adopted": {
            "scoreGenomicLengthLog2scale": (
                "0 — AlleleSeq2는 custom/ASE/ASB/ASCA 4개 모드 전부에 0.0을 명시한다. "
                "STAR 기본값 -0.25는 정렬 genomic length의 log2 비례 항을 AS에 섞는데, "
                "HAP1/HAP2는 indel로 같은 fragment의 span이 달라질 수 있어 "
                "염기 일치와 무관한 AS 차이가 된다. 33이 AS 1점 차이로 배정하므로 critical하다."),
            "addnl_seqns": (
                "unplaced/unlocalized/EBV/decoy contig를 index에만 포함 — "
                "AlleleSeq2 makePersonalGenome.mk의 ADDNL_SEQNS. alt/HLA는 양쪽 모두 제외."
                if addnl_seqns else "비활성(RNA_INDEX_ADDNL_SEQNS=0)"),
        },
        "deliberately_not_adopted": {
            "sjdbScore": (
                "AlleleSeq2가 --sjdbScore 0을 주는 것은 ASB/ASCA 즉 splicing이 없는 DNA 모드뿐이고 "
                "(--scoreGap -100이 함께 붙는다), RNA-seq용 ASE-params는 기본값(+2)을 쓴다. "
                "따라서 RNA 경로에 0을 넣는 것은 AlleleSeq2를 따르는 것이 아니다."),
            "single_diploid_index": (
                "AlleleSeq2는 hap1+hap2를 합친 단일 index에 1회 정렬해 MAPQ 255가 "
                "diploid 전체에서 unique를 뜻한다. v2는 독립 2회 정렬로 AS를 비교하는 "
                "phASER/WASP 계열 설계라 구조가 다르다. twopassMode Basic이 두 index에서 "
                "따로 도는 데서 오는 novel junction 집합 차이는 미해결로 남아 있다."),
            "hetsnv_site_counting": (
                "AlleleSeq2는 hetSNV site 단위 mpileup으로 세고 filter_phase_warnings.py로 "
                "hap allele이 최다가 아닌 site를 버린다. v2는 fragment의 AS 비교 후 gene 단위로 "
                "세므로 site 단위 판정 지점이 없다."),
        },
    },
    "topmed_parameter_policy": {
        "same_as_topmed": {
            "twopass_mode": "Basic",
            "out_filter_multimap_nmax": 20,
            "out_filter_mismatch_nmax": 999,
            "out_filter_mismatch_nover_lmax": 0.1,
            "out_filter_score_min_over_lread": 0.33,
            "out_filter_match_nmin_over_lread": 0.33,
            "intron_limits": [20, 1000000],
            "mate_gap_max": 1000000,
            "junction_overhangs": [8, 1],
            "out_filter_intron_motifs": "None",
            "align_soft_clip_at_reference_ends": "Yes",
            "markduplicates": "Picard; duplicate flag only, no removal",
        },
        "personal_hap_purpose_deviations": {
            "reference": "HAP별 SNV+indel personal FASTA (TOPMed의 공통 GRCh38+ERCC 대신)",
            "annotation": "chain으로 lift한 HAP별 GTF",
            "sjdb_overhang": "75 bp RNA read의 read_length-1 = 74 (TOPMed 2x101 bp index의 100 대신)",
            "contigs": (f"카운트 대상은 요청 CHROMS만: {contigs}; ERCC/ALT/HLA 미포함. "
                        + ("index에는 AlleleSeq2 ADDNL_SEQNS(unplaced/unlocalized/EBV/decoy)를 함께 넣는다"
                           if addnl_seqns else "decoy도 미포함(RNA_INDEX_ADDNL_SEQNS=0)")),
            "score_genomic_length_log2scale": "0 (STAR 기본 -0.25 대신; AlleleSeq2 전 모드 공통. TOPMed는 미지정)",
            "fastq_command": "임시 FASTQ가 비압축이므로 zcat 대신 cat",
            "bam_sorting": "coverage·ASE 직전 MarkDuplicates를 위해 STAR coordinate-sorted BAM 사용",
            "quantification": "RSEM 입력용 TranscriptomeSAM은 만들지 않고 GeneCounts만 생성",
            "chimeric_output": "fusion/구조변이 산출물이 목적이 아니므로 TOPMed의 chimera 출력 옵션은 사용하지 않음",
        },
    },
    "annotation": "GENCODE v29 lifted to personal coords (13_lift_annotation.py)",
    "uniquely_mapped_pct": float(unique_pct) if unique_pct else None,
    "coverage_zarr": coverage_zarr,
    "coverage_duplicate_policy": "samtools depth 기본 필터로 DUP flag 제외; marked BAM은 보존",
    "note": "커버리지는 이 haplotype 좌표에서 측정한 총 발현이다. haplotype 특이 신호는 16_allelic_counts.py의 read별 AS 비교를 사용한다.",
}
Path(out).write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
PY
  mv "$TD/rna_realign.json.part" "$TD/rna_realign.json" || die "$S HAP$h RNA QC JSON 이동 실패"
  log "  [ok] HAP$h uniquely mapped ${uniq:-?}%"
  (( STAR_KEEP_INDEX )) || rm -rf "$IDX"
done
rm -f "$W"/rna/collate_tmp*.bam 2>/dev/null
log "14 완료: $S"
