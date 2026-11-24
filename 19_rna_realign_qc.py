#!/usr/bin/env python3
"""19 — HAP 별 RNA 재정렬 QC와 sample summary.

    python 19_rna_realign_qc.py hap --sample S --hap 1 --label clair3 --work-dir <W>/rna --out-dir <R>/HAP1 \
        --fastq-qc <W>/rna/fastq_qc.json --lift-qc <R>/HAP1/gtf_lift_qc.json --coverage <R>/HAP1/coverage.zarr \
        --gencode-gtf <gtf> --gencode-tag gencodev47
    python 19_rna_realign_qc.py summary --sample S --label clair3 --fastq-qc … --rna-dir <R>
    python 19_rna_realign_qc.py check --realign-json <R>/HAP1/rna_realign.json --gencode-tag …

입력   13 fastq_qc.json, 14 gtf_lift_qc.json, 15 hap{h}.index.json, 16 hap{h}.align.json, 17 hap{h}.markdup.json,
       18 coverage.zarr
출력   <out-dir>/rna_realign.json (schema mfm_v2_rna_realign_v6) — HAP별 14–18 산출물과 profile fingerprint를
       한곳에 기록한다. `check` 명령은 수동/호환성 검증용이며, 현재 runner의 단계별 재사용은 각 단계가 직접 판정한다.
       <R>/rna_realign_summary.json (schema mfm_rna_qc_v2) — 두 HAP summary (canonical)
       <R>/rna_qc.json — 기존 consumer를 위한 동일 payload 호환 사본
"""
from __future__ import annotations

import os
import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rna_common import die, log, nonempty, read_json, write_json  # noqa: E402
import rna_profile  # noqa: E402

REALIGN_SCHEMA = "mfm_v2_rna_realign_v6"
PROFILE = f"{rna_profile.PROFILE_NAME}: GTEx/TOPMed V11 (broadinstitute/gtex-pipeline run_STAR.py, STAR 2.7.11b) on HAP reference"
if os.environ.get('MFM_RNA_BACKEND','cpu')=='parabricks':
    PROFILE = 'gtex_v11_hap_gpu_star272a: Parabricks4.8; outFilterType=Normal; CPU GTEx profile과 비동등함'



def realign_payload(a, fastq: dict, lift: dict, ix: dict, al: dict, mdj: dict) -> dict:
    md = {k: mdj[k] for k in ("tool", "remove_duplicates", "marked_bam", "metrics_path", "read_pairs_examined",
                              "read_pair_duplicates", "unpaired_reads_examined", "unpaired_read_duplicates",
                              "duplicate_read_records_from_metrics", "primary_duplicate_records", "primary_records",
                              "percent_duplication")}
    gpu={}
    if os.environ.get('MFM_RNA_BACKEND','cpu')=='parabricks':
        gpu={'gpu_execution':read_json(Path(al['align_dir'])/'gpu_backend.json')}
    return {
        **gpu,
        "schema": REALIGN_SCHEMA,
        "sample": a.sample,
        "hap": f"HAP{a.hap}",
        "phase_mode": a.label,
        "aligner": f"STAR {subprocess.getoutput(f'{a.star} --version')}",
        "profile": PROFILE,
        "intentional_modifications_vs_gtex": [
            "M1 reference = HAP 개인 chr1–22 (+ GTEx index 의 나머지 contig chrX·chrY·chrM·chrUn·random·EBV 를 reference 서열로)",
            "M2 annotation = HAP 좌표로 lift 한 GENCODE (14, trim)",
            "M3 coverage = GTEx bam2coverage와 같이 duplicate 포함(-g 1024)",
        ],
        "output_only_differences_vs_gtex": ["TranscriptomeSAM·chimeric 출력 없음", "좌표 정렬 BAM 을 STAR 가 바로 씀",
                                            "비압축 FASTQ(cat)", "RG = <S>.HAP<h>"],
        "sjdb_overhang": int(ix["sjdb_overhang"]),
        "score_genomic_length_log2scale": al.get("score_genomic_length_log2scale", -0.25),
        "index_contigs": ix["contigs"] + (["GTEx 비-autosome contig (15)"] if ix.get("addnl_seqns") else []),
        "fastq_reconstruction": fastq,
        "gtf_lift_qc": lift,
        "markduplicates": md,
        "annotation": f"{a.gencode_tag} lifted to HAP coordinates (14_lift_annotation.py, trim)",
        "gencode_gtf": a.gencode_gtf,
        "gencode_tag": a.gencode_tag,
        "uniquely_mapped_pct": al.get("uniquely_mapped_pct"),
        "coverage_zarr": str(a.coverage),
        "coverage_definition": "GTEx-compatible samtools depth 정책: " + rna_profile.COVERAGE_POLICY,
        "reuse_fingerprint": rna_profile.fingerprint(a.gencode_gtf, a.extra_contig_fai, bool(ix.get("addnl_seqns"))),
        "note": "커버리지는 이 HAP 좌표에서 측정한 총 발현이다. HAP 특이 신호는 21의 fragment별 AS 비교를 쓴다.",
        "step_records": {"index": ix, "align": al},
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    h = sub.add_parser("hap")
    for k in ("--sample", "--hap", "--label", "--work-dir", "--out-dir", "--fastq-qc", "--lift-qc", "--coverage",
              "--gencode-gtf", "--gencode-tag"):
        h.add_argument(k, required=True)
    h.add_argument("--star", default="STAR")
    h.add_argument("--extra-contig-fai", default="", help="15 의 추가 contig 원본 FASTA .fai (fingerprint 용)")
    s = sub.add_parser("summary")
    for k in ("--sample", "--label", "--fastq-qc", "--rna-dir"):
        s.add_argument(k, required=True)
    c = sub.add_parser("check")
    for k in ("--realign-json", "--gencode-tag", "--gencode-gtf"):
        c.add_argument(k, required=True)
    c.add_argument("--extra-contig-fai", default="")
    c.add_argument("--addnl", type=int, default=1)
    a = ap.parse_args()

    if a.cmd == "check":        # 0 = 이 HAP의 QC·fingerprint·coverage 계약이 현재 설정과 일치함
        p = Path(a.realign_json)
        if not nonempty(p):
            return 1
        d = read_json(p)
        want = rna_profile.fingerprint(a.gencode_gtf, a.extra_contig_fai, bool(a.addnl))["sha256"]
        ok = (d.get("schema") == REALIGN_SCHEMA and d.get("gencode_tag") == a.gencode_tag
              and (d.get("reuse_fingerprint") or {}).get("sha256") == want
              and Path(d.get("coverage_zarr", "")).is_dir())
        return 0 if ok else 1

    if a.cmd == "summary":
        rd = Path(a.rna_dir)
        # Stage 25의 final package 계약은 재정렬 summary라는 의미가 드러나는
        # rna_realign_summary.json을 canonical로 사용한다. 기존 분석 노트와
        # 재개 도구가 rna_qc.json을 참조하므로 같은 immutable payload를 함께
        # 남긴다. 둘 중 하나가 다른 내용을 갖는 migration 상태를 만들지 않는다.
        payload = {
            "schema": "mfm_rna_qc_v2", "sample": a.sample, "phase_mode": a.label,
            "fastq_reconstruction": read_json(Path(a.fastq_qc)),
            "haplotypes": [read_json(rd / hp / "rna_realign.json") for hp in ("HAP1", "HAP2")]}
        canonical = rd / "rna_realign_summary.json"
        legacy = rd / "rna_qc.json"
        write_json(canonical, payload)
        write_json(legacy, payload)
        log(f"{a.sample}: RNA QC 요약 → {canonical} (compatibility: {legacy.name})")
        return 0

    wd, od = Path(a.work_dir), Path(a.out_dir)
    need = [Path(a.fastq_qc), Path(a.lift_qc), wd / f"hap{a.hap}.index.json", wd / f"hap{a.hap}.align.json",
            wd / f"hap{a.hap}.markdup.json"]
    miss = [str(p) for p in need if not nonempty(p)]
    if miss or not Path(a.coverage).is_dir():
        die(f"{a.sample} HAP{a.hap} QC 입력 없음: {miss or a.coverage}")
    payload = realign_payload(a, read_json(need[0]), read_json(need[1]), read_json(need[2]),
                              read_json(need[3]), read_json(need[4]))
    write_json(od / "rna_realign.json", payload)
    log(f"  [ok] HAP{a.hap} uniquely mapped {payload['uniquely_mapped_pct']}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
