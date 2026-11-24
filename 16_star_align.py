#!/usr/bin/env python3
"""16 — RNA read 를 HAP 하나의 STAR index 에 정렬 (GTEx/TOPMed 파라미터).

    python 16_star_align.py --sample S --hap 1 --index-json <W>/rna/hap1.index.json \
        --fastq-qc <W>/rna/fastq_qc.json --out-dir <R>/HAP1 --work-dir <W>/rna --log-dir <W>/logs

입력   13 의 R1/R2 FASTQ(fastq_qc.json 에 경로), 15 의 index(hap{h}.index.json)
출력   <work-dir>/aln_hap{h}_oh{N}_{cache-tag}/Aligned.sortedByCoord.out.bam (있으면 재사용)
       <out-dir>/Log.final.out, ReadsPerGene.out.tab
       <work-dir>/hap{h}.align.json      다음 단계(17)가 읽는 BAM 경로 기록

파라미터 (TOPMed/GTEx run_STAR.py 와 같음; 차이는 19 의 rna_realign.json 에 기록)
  --twopassMode Basic --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999
  --outFilterMismatchNoverLmax 0.1 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33
  --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8
  --alignSJDBoverhangMin 1 --outFilterIntronMotifs None --alignSoftClipAtReferenceEnds Yes
  --limitSjdbInsertNsj 1200000 --outSAMunmapped Within --outSAMstrandField intronMotif --quantMode GeneCounts
★ --scoreGenomicLengthLog2scale -0.25 : STAR 기본값을 명시만 한다(GTEx 는 지정하지 않음 = 같은 값). 이 항이 HAP 간 AS 비교에
  섞이는 것은 21이 AS에서 계산해 뺀다.
★ --outSAMattributes NH HI AS nM NM : 21이 AS를 읽는다.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import uuid
from pathlib import Path
from rna_cache import cache_fields, cache_valid, file_fingerprint, fingerprint, read_cache, tool_fingerprint

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rna_common import die, log, nonempty, read_json, run, write_json  # noqa: E402
import rna_profile  # noqa: E402

STAR_TOPMED_ARGS = rna_profile.STAR_ALIGN_ARGS    # GTEx run_STAR.py 기본값 (rna_profile.py)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--hap", required=True, choices=("1", "2"))
    ap.add_argument("--index-json", required=True, type=Path)
    ap.add_argument("--fastq-qc", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--work-dir", required=True, type=Path)
    ap.add_argument("--log-dir", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=int(os.environ.get("STAR_ALIGN_THREADS", "8")))
    ap.add_argument("--star", default=os.environ.get("STAR_BIN", "STAR"))
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    S, h = a.sample, a.hap
    ix = read_json(a.index_json)
    fq = read_json(a.fastq_qc)
    r1, r2 = fq["r1_fastq"], fq["r2_fastq"]
    if not nonempty(Path(r1), Path(r2)):
        die(f"{S} FASTQ 없음: {r1} {r2} (13 먼저)")
    idx = Path(ix["index_dir"])
    if not cache_valid(ix):
        die(f"{S} HAP{h} STAR index 입력/완료 계약이 없거나 변경됨 (15부터 재생성 필요)")
    if not nonempty(idx / "SAindex"):
        die(f"{S} HAP{h} STAR index 없음: {idx} (15 먼저)")
    inputs = {"stage": 16, "sample": S, "hap": h,
              "index_json": file_fingerprint(a.index_json, content=True),
              "r1": file_fingerprint(r1), "r2": file_fingerprint(r2),
              "star_args": list(STAR_TOPMED_ARGS), "star": tool_fingerprint(a.star),
              "script": file_fingerprint(__file__, content=True)}
    sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))
    from integrate_wgs_rna import gpu_rna
    integrated=gpu_rna.integrated_markdup_enabled()
    if os.environ.get('MFM_RNA_BACKEND')=='parabricks':
        inputs['gpu_adapter']=file_fingerprint(gpu_rna.__file__,content=True)
        inputs['integrated_markdup']=integrated
    record_path = a.work_dir / f"hap{h}.align.json"
    old = read_cache(record_path)
    reuse = not a.force and cache_valid(old, inputs)
    if reuse:
        ad = Path(old["align_dir"])
        bam = Path(old["star_bam"])
        log(f"  HAP{h}: STAR 정렬 재사용 ({bam})")
    else:
        ad = a.work_dir / f"aln_hap{h}_oh{ix['sjdb_overhang']}_{ix['cache_tag']}_{fingerprint(inputs)[:16]}_{uuid.uuid4().hex[:12]}"
        ad.mkdir(parents=True, exist_ok=True)
        bam = ad / "Aligned.sortedByCoord.out.bam"
        log(f"  HAP{h}: STAR 정렬")
        run([a.star, "--runMode", "alignReads", "--runThreadN", str(a.threads), "--genomeDir", str(idx),
             "--readFilesIn", r1, r2, "--readFilesCommand", "cat", "--outFileNamePrefix", f"{ad}/",
             *STAR_TOPMED_ARGS, "--outSAMattrRGline", f"ID:{S}.HAP{h}", f"SM:{S}"],
            log_path=a.log_dir / f"16.align.HAP{h}.log", what=f"{S} HAP{h} STAR 정렬")
    if not nonempty(bam):
        die(f"{S} HAP{h} STAR BAM 없음: {bam}")
    a.out_dir.mkdir(parents=True, exist_ok=True)
    for name in ("Log.final.out", "ReadsPerGene.out.tab"):
        if not nonempty(ad / name):
            die(f"{S} HAP{h} STAR {name} 없음: {ad}")
        shutil.copyfile(ad / name, a.out_dir / name)
    uniq = None
    for line in (ad / "Log.final.out").read_text().splitlines():
        if "Uniquely mapped reads %" in line:
            uniq = float(line.split("|")[1].strip().rstrip("%"))
    if reuse:
        return 0
    integrated_record=None
    outputs=[bam, ad / "Log.final.out", ad / "ReadsPerGene.out.tab"]
    if os.environ.get('MFM_RNA_BACKEND')=='parabricks':
        backend=ad/'gpu_backend.json'
        execution=json.loads(backend.read_text())
        integrated_record=execution.get('integrated_markdup')
        if integrated != (integrated_record is not None):die('GPU 중복표시 실행 모드와 완료 계약이 다름')
        outputs.append(backend)
        if integrated:
            outputs.append(gpu_rna.validate_integrated_markdup(integrated_record,bam))
    write_json(record_path, {
        "schema": "mfm_rna_star_align_v1", "sample": S, "hap": f"HAP{h}", "align_dir": str(ad),
        "star_bam": str(bam), "integrated_markdup": integrated_record, "score_genomic_length_log2scale": rna_profile.STAR_GENOMIC_LENGTH_LOG2SCALE, "uniquely_mapped_pct": uniq,
        "index_json": str(a.index_json), "fastq_qc": str(a.fastq_qc),
        **cache_fields(inputs, outputs)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
