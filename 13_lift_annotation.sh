#!/usr/bin/env bash
# 구번호 Stage 13 직접 실행기. 현재 RNA runner에서는 Stage 14 `14_lift_annotation.py`가 실행된다.
# 이전 실행을 재현하기 위해 보존한다.
# 13 — GENCODE 주석을 HAP1/HAP2 개인 좌표로 lift한다.
set -euo pipefail
cd "$(dirname "$0")"
source ./config_v2.sh

S=${1:?사용법: $0 <sample>}
v2_require
require_manifest_sample_contract "$S"
if [[ "$(sequence_source_for_sample "$S")" == assembly_direct ]]; then
  die "RNA assembly_direct는 공식 multi-chain 기반 GTF lift 검증 전까지 지원하지 않음: sample=$S"
fi

F=$(v2_final "$S")
W=$(v2_work "$S")
MTAG=$(sample_rna_evidence_label "$S")
RNA_CONTIGS=$(chroms_space "$CHROMS")
[[ -n "$RNA_CONTIGS" ]] || die "$S RNA contig 목록이 비어 있음: CHROMS=$CHROMS"
mkdir -p "$W/logs"
require_hap_contig_assets "$F/individual_fasta" "$CHROMS" 1

normalize_lift_qc_paths() {
  local qc=$1 gtf=$2 exons=$3 realign=${4:-}
  "$PYTHON" - "$qc" "$gtf" "$exons" "$realign" <<'PY'
import json
import os
import sys
from pathlib import Path

qc_path, gtf, exons, realign_path = map(Path, sys.argv[1:])
qc = json.loads(qc_path.read_text())
qc["output_gtf"] = str(gtf)
qc["output_exons"] = str(exons)
tmp = Path(str(qc_path) + ".pathfix.part")
tmp.write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n")
os.replace(tmp, qc_path)
if str(realign_path) and realign_path.is_file():
    payload = json.loads(realign_path.read_text())
    payload["gtf_lift_qc"] = qc
    tmp = Path(str(realign_path) + ".pathfix.part")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    os.replace(tmp, realign_path)
PY
}

for h in 1 2; do
  HD=$F/individual_fasta/HAP$h
  TD=$F/rna_personal/$MTAG/HAP$h
  GTF=$TD/annotation.gtf
  EXONS=$TD/exons.bygene.tsv
  LIFT_QC=$TD/gtf_lift_qc.json
  mkdir -p "$TD"

  if [[ "${RNA_FORCE:-0}" != 1 && -s "$GTF" && -s "$EXONS" && -s "$LIFT_QC" ]]; then
    normalize_lift_qc_paths "$LIFT_QC" "$GTF" "$EXONS" "$TD/rna_realign.json"
    log "  [skip] HAP$h GTF lift"
    continue
  fi

  log "  HAP$h: GTF lift"
  "$PYTHON" ./13_lift_annotation.py --chain-dir "$HD" --gtf "$GENCODE_GTF" \
    --out-gtf "$GTF.part" --out-exons "$EXONS.part" --qc-json "$LIFT_QC.part" \
    --on-unlifted "$RNA_LIFT_ON_UNLIFTED" --contigs "$RNA_CONTIGS" \
    2>"$W/logs/13.lift.HAP$h.log" || die "$S HAP$h GTF lift 실패"
  mv "$GTF.part" "$GTF" || die "$S HAP$h GTF 출력 이동 실패"
  mv "$EXONS.part" "$EXONS" || die "$S HAP$h exon 출력 이동 실패"
  mv "$LIFT_QC.part" "$LIFT_QC" || die "$S HAP$h GTF QC 출력 이동 실패"
  normalize_lift_qc_paths "$LIFT_QC" "$GTF" "$EXONS" "$TD/rna_realign.json"
done

log "13 완료: $S"
