#!/usr/bin/env python3
"""15 — HAP 하나의 STAR index 생성 (HAP FASTA 연결 + 흡수용 contig + lift GTF).

    python 15_star_index.py --sample S --hap 1 --hap-fasta-dir <F>/individual_fasta/HAP1 --contigs "chr1 … chr22" \
        --gtf <R>/HAP1/annotation.gtf --fastq-qc <W>/rna/fastq_qc.json --work-dir <W>/rna --log-dir <W>/logs \
        --cache-tag lift4trim_gencodev47 [--addnl-fasta <…addnl_seqns….fa> --addnl-source <hs38dh.fa>]

입력   HAP FASTA(<chr>.fa, 09 영구 자산), 14 가 만든 HAP 좌표 GTF, 13 의 fastq_qc.json(최대 read 길이)
출력   <work-dir>/HAP{h}.fa(.fai)        요청 contig 를 이어 붙인 index 입력
       <work-dir>/idx_hap{h}_oh{N}_{cache-tag}/   STAR index (SAindex 가 있으면 재사용)
       <work-dir>/hap{h}.index.json      다음 단계(16)가 읽는 index 경로·sjdbOverhang 기록

sjdbOverhang (--sjdb-overhang auto)
  STAR 매뉴얼: "Ideally, this length should be equal to the ReadLength-1 … In case of reads of varying length,
  the ideal value is max(ReadLength)-1." 이전에는 75 bp 를 가정해 74 로 고정했는데 76 bp RNA 도 있다(HG00373).
  기본은 FASTQ 최대 read 길이 - 1, 정수를 주면 그 값으로 고정.
흡수용 contig (2026-09-28: GTEx 와 같게)
  GTEx index reference = Broad Homo_sapiens_assembly38 에서 ALT·HLA·decoy 를 뺀 것(TOPMed_RNAseq_pipeline.md 49·65–70줄).
  개인 HAP 은 chr1–22 만 있으므로, 그 나머지(chrX, chrY, chrM, chrUn_* 중 decoy 아닌 것, *_random, chrEBV)를 reference
  서열 그대로 index 에 넣는다. 이 read 들이 autosome 의 NUMT·paralog 에 억지로 붙지 않게 하는 용도이며, coverage 와
  count 는 chr1–22 만 쓴다. decoy(*_decoy), alt(*_alt), HLA 는 GTEx 처럼 넣지 않는다.
  --addnl-fasta 가 없으면 --addnl-source(.fai 포함)에서 위 규칙으로 뽑아 한 번만 만든다.
캐시 키에 sjdbOverhang·lift 정책·GENCODE tag 가 들어가므로, 이 중 하나라도 바뀌면 index 를 새로 만든다.
파일럿 실측(HG00149 HAP1, 2.8 Gb): 4 스레드 25분 43초 / index 26 GB.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import time
import sys
import uuid
from pathlib import Path
from rna_cache import cache_fields, cache_valid, file_fingerprint, fingerprint, read_cache, tool_fingerprint

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rna_common import die, log, nonempty, read_json, run, write_json  # noqa: E402
from rna_profile import select_nonautosome  # noqa: E402,F401


def sjdb_overhang(spec: str, fastq_qc: Path) -> tuple[int, str]:
    if spec == "auto":
        qc = read_json(fastq_qc)
        mx = qc.get("max_read_length")
        if not mx:
            out = subprocess.run(["awk", "NR%4==2 && length($0)>m{m=length($0)} END{print m+0}",
                                  qc["r1_fastq"], qc["r2_fastq"]], capture_output=True, text=True, check=True)
            mx = int(out.stdout.strip())
        if int(mx) < 2:
            die(f"RNA 최대 read 길이 산출 실패: {mx}")
        return int(mx) - 1, f"max(read length) - 1 = {mx} - 1"
    if not spec.isdigit() or int(spec) < 1:
        die(f"RNA_SJDB_OVERHANG은 auto 또는 양의 정수: {spec}")
    return int(spec), "RNA_SJDB_OVERHANG 명시값"


def ensure_addnl(fasta: Path, source: str, samtools: str) -> None:
    if nonempty(fasta):
        return
    if not (source and nonempty(Path(source)) and nonempty(Path(source + ".fai"))):
        die(f"ADDNL_SEQNS 원본 또는 .fai 없음: {source} (RNA_INDEX_ADDNL_SEQNS=0으로 끌 수 있음)")
    contigs = select_nonautosome([l.split("\t")[0] for l in open(source + ".fai")])
    if not contigs:
        die(f"ADDNL_SEQNS 대상 contig를 찾지 못함: {source}.fai")
    log(f"  ADDNL_SEQNS 생성(최초 1회): {fasta} ({len(contigs)} contig)")
    tmp = fasta.with_name(f"{fasta.name}.{os.getpid()}.part")   # 여러 sample 이 동시에 만들어도 반쪽 파일이 안 남게
    with tmp.open("w") as fh:
        if subprocess.run([samtools, "faidx", source, *contigs], stdout=fh).returncode != 0:
            tmp.unlink(missing_ok=True); die("ADDNL_SEQNS 추출 실패")
    os.replace(tmp, fasta)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--hap", required=True, choices=("1", "2"))
    ap.add_argument("--hap-fasta-dir", required=True, type=Path)
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--gtf", required=True, type=Path)
    ap.add_argument("--fastq-qc", required=True, type=Path)
    ap.add_argument("--work-dir", required=True, type=Path)
    ap.add_argument("--log-dir", required=True, type=Path)
    ap.add_argument("--cache-tag", required=True, help="예: lift4trim_gencodev47")
    ap.add_argument("--sjdb-overhang", default=os.environ.get("RNA_SJDB_OVERHANG", "auto"))
    ap.add_argument("--addnl-fasta", default="", help="비우면 ADDNL_SEQNS 를 index 에 넣지 않는다")
    ap.add_argument("--addnl-source", default="")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("STAR_INDEX_THREADS", "8")))
    ap.add_argument("--star", default=os.environ.get("STAR_BIN", "STAR"))
    ap.add_argument("--samtools", default=os.environ.get("SAMTOOLS", "samtools"))
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    S, h, wd = a.sample, a.hap, a.work_dir
    contigs = a.contigs.split()
    if not contigs:
        die(f"{S} RNA contig 목록이 비어 있음")
    if not nonempty(a.gtf):
        die(f"{S} HAP{h} lift GTF 없음: {a.gtf} (14 먼저)")
    oh, oh_src = sjdb_overhang(a.sjdb_overhang, a.fastq_qc)
    log(f"{S}: STAR sjdbOverhang={oh} ({oh_src})")

    # ── target FASTA: 요청 contig 를 이어 붙인다 ──
    wd.mkdir(parents=True, exist_ok=True)
    target_inputs = {"sample": S, "hap": h, "contigs": contigs,
                     "sources": [file_fingerprint(a.hap_fasta_dir / f"{c}.fa") for c in contigs],
                     "samtools": tool_fingerprint(a.samtools),
                     "script": file_fingerprint(__file__, content=True)}
    target_record_path = wd / f"hap{h}.target.json"
    target_record = read_cache(target_record_path)
    if cache_valid(target_record, target_inputs):
        tfa = Path(target_record["target_fasta"])
    else:
        tfa = wd / f"HAP{h}.{fingerprint(target_inputs)[:16]}.{uuid.uuid4().hex[:12]}.fa"
        part = tfa.with_name(tfa.name + ".part")
        with part.open("wb") as out:
            for c in contigs:
                src = a.hap_fasta_dir / f"{c}.fa"
                if not nonempty(src):
                    die(f"{S} HAP{h} {c} FASTA 없음: {src}")
                with src.open("rb") as fh:
                    shutil.copyfileobj(fh, out, 16 << 20)
        os.replace(part, tfa)
        run([a.samtools, "faidx", str(tfa)], what=f"{S} HAP{h} target FASTA index")
        write_json(target_record_path, {"target_fasta": str(tfa),
                   **cache_fields(target_inputs, [tfa, Path(str(tfa) + ".fai")])})

    addnl: list[str] = []
    if a.addnl_fasta:
        ensure_addnl(Path(a.addnl_fasta), a.addnl_source, a.samtools)
        addnl = [a.addnl_fasta]

    inputs = {"stage": 15, "sample": S, "hap": h, "contigs": contigs,
              "target": file_fingerprint(tfa), "target_fai": file_fingerprint(str(tfa) + ".fai", content=True),
              "gtf": file_fingerprint(a.gtf, content=True),
              "addnl": [file_fingerprint(path) for path in addnl],
              "sjdb_overhang": oh, "cache_tag": a.cache_tag,
              "star": tool_fingerprint(a.star), "script": file_fingerprint(__file__, content=True)}
    record_path = wd / f"hap{h}.index.json"
    old = read_cache(record_path)
    required = ("Genome", "SA", "SAindex", "chrName.txt", "chrLength.txt", "chrStart.txt", "genomeParameters.txt")
    old_dir = Path(old.get("index_dir", wd / "missing_index"))
    if not a.force and cache_valid(old, inputs) and nonempty(*(old_dir / name for name in required)):
        log(f"  HAP{h}: STAR index 재사용 ({old_dir})")
        return 0
    else:
        idx = wd / f"idx_hap{h}_oh{oh}_{a.cache_tag}_{fingerprint(inputs)[:16]}_{uuid.uuid4().hex[:12]}"
        log(f"  HAP{h}: STAR index (파일럿 실측 4스레드 26분 / 26 GB)")
        idx.mkdir(parents=True, exist_ok=True)
        index_started = time.time()
        run([a.star, "--runMode", "genomeGenerate", "--runThreadN", str(a.threads), "--genomeDir", str(idx),
             "--genomeFastaFiles", str(tfa), *addnl, "--sjdbGTFfile", str(a.gtf), "--sjdbOverhang", str(oh),
             "--outFileNamePrefix", f"{idx}/"],
            log_path=a.log_dir / f"15.index.HAP{h}.log", what=f"{S} HAP{h} STAR index")
        index_finished = time.time()
    if not nonempty(*(idx / name for name in required)):
        die(f"{S} HAP{h} STAR index 필수 출력 누락: {idx}")
    write_json(record_path, {
        "schema": "mfm_rna_star_index_v1", "sample": S, "hap": f"HAP{h}", "target_fasta": str(tfa),
        "target_fai": str(tfa) + ".fai", "index_dir": str(idx), "sjdb_overhang": oh,
        "index_execution": {"started_at": index_started, "finished_at": index_finished, "threads": a.threads},
        "sjdb_overhang_source": oh_src, "gtf": str(a.gtf), "addnl_fasta": a.addnl_fasta or None,
        "addnl_seqns": bool(a.addnl_fasta), "contigs": contigs, "cache_tag": a.cache_tag,
        **cache_fields(inputs, sorted(p for p in idx.iterdir() if p.is_file()))})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
