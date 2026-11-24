#!/usr/bin/env bash
# 구번호 Stage 15 직접 실행기. 현재 RNA runner의 Stage 18은 `18_coverage_zarr.py`를 사용한다.
# 이 wrapper의 filtering/output 계약은 현재 단계 계약과 구분한다.
# 15 — Stage 14 marked BAM에서 HAP별 1 bp coverage Zarr와 통합 RNA QC를 만든다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
v2_require
require_manifest_sample_contract "$S"
F=$(v2_final "$S")
W=$(v2_work "$S")
MTAG=$(sample_rna_evidence_label "$S")
RNA_CONTIGS=$(chroms_space "$CHROMS")
FASTQ_QC=$W/rna/fastq_qc.json
RNA_PREFLIGHT_QC=$F/qc/rna_preflight.json
mkdir -p "$W/logs" "$W/rna"
require_file "$FASTQ_QC"
require_file "$RNA_PREFLIGHT_QC"
require_hap_contig_assets "$F/individual_fasta" "$CHROMS" 1

RNA_CONTRACT_SHA256=$("$PYTHON" -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d.get("status")=="pass"; print(d["contract_sha256"])' "$RNA_PREFLIGHT_QC") \
  || die "$S RNA preflight contract 읽기 실패"

contract_matches() {
  local path=$1
  [[ -s "$path" ]] && "$PYTHON" -c 'import json,sys; d=json.load(open(sys.argv[1])); raise SystemExit(0 if d.get("rna_preflight_contract_sha256")==sys.argv[2] else 1)' "$path" "$RNA_CONTRACT_SHA256"
}

for h in 1 2; do
  HD=$F/individual_fasta/HAP$h
  TD=$F/rna_personal/$MTAG/HAP$h
  FINAL_BAM=$TD/$S.rna.personal.HAP$h.bam
  REALIGN_QC=$TD/rna_realign.json
  COV=$TD/coverage.zarr
  TFA=$W/rna/HAP$h.fa
  require_file "$FINAL_BAM"
  require_file "$FINAL_BAM.bai"
  require_file "$REALIGN_QC"
  contract_matches "$REALIGN_QC" || die "$S HAP$h RNA realign contract 불일치"

  if [[ "${RNA_FORCE:-0}" != 1 && -d "$COV" ]]; then
    log "  [skip] HAP$h coverage"
    continue
  fi

  # Stage 15만 재시작해도 동작하도록 Stage 14 작업 FASTA가 정리된 경우 동일한
  # contig 순서로 재구성한다. 원본 HAP FASTA에는 손대지 않는다.
  if [[ ! -s "$TFA" || ! -s "$TFA.fai" ]]; then
    rm -f "$TFA.part" "$TFA.fai"
    for c in $RNA_CONTIGS; do
      cat "$HD/$c.fa" >> "$TFA.part" || die "$S HAP$h $c FASTA concat 실패"
    done
    mv "$TFA.part" "$TFA" || die "$S HAP$h target FASTA 이동 실패"
    "$SAMTOOLS" faidx "$TFA" || die "$S HAP$h target FASTA index 실패"
  fi

  "$PYTHON" ./15_coverage_zarr.py --bam "$FINAL_BAM" --fai "$TFA.fai" --out "$COV" \
    --contigs "$RNA_CONTIGS" --chunk "$RNA_ZARR_CHUNK" \
    --power "$AG_NORM_POWER" --clip "$AG_NORM_CLIP" \
    --sample "$S" --hap "HAP$h" --samtools "$SAMTOOLS" \
    2>&1 | tee "$W/logs/15.cov.HAP$h.log" || die "$S HAP$h coverage 생성 실패"
  [[ -d "$COV" ]] || die "$S HAP$h coverage Zarr 없음: $COV"
done

SUMMARY=$F/rna_personal/$MTAG/rna_qc.json
"$PYTHON" - "$SUMMARY.part" "$S" "$MTAG" "$FASTQ_QC" \
  "$F/rna_personal/$MTAG/HAP1/rna_realign.json" \
  "$F/rna_personal/$MTAG/HAP2/rna_realign.json" <<'PY'
import json
import sys
from pathlib import Path

out, sample, phase_mode, fastq_qc, hap1, hap2 = sys.argv[1:]
payload = {
    "schema": "mfm_rna_qc_v2",
    "sample": sample,
    "phase_mode": phase_mode,
    "fastq_reconstruction": json.loads(Path(fastq_qc).read_text()),
    "haplotypes": [json.loads(Path(path).read_text()) for path in (hap1, hap2)],
}
Path(out).write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
PY
mv "$SUMMARY.part" "$SUMMARY"
log "15 완료: $S"
