#!/usr/bin/env python3
"""20 — ASE 근거로 쓸 HAP1–HAP2 차이(het variant)를 WGS read 로 검증한다.

    python 20_hap_variants.py --vcf <HAP FASTA 를 만든 VCF> --wgs-bam <06 haplotagged BAM> \
        --hap-dir1 <individual_fasta/HAP1> --hap-dir2 <individual_fasta/HAP2> --gtf <GENCODE> \
        --contigs "chr1 ..." --out <hap_variants.tsv> --qc-json <hap_variants_qc.json>

★ 왜
  21 은 같은 fragment 를 HAP1·HAP2 에 정렬한 결과로 allele 을 정한다. 이 방법은 "두 정렬이 다르면 그 fragment 가
  HAP 사이 차이를 지났다"는 전제에 기대는데, HG00096 EDA 에서 이 전제가 깨지는 두 경우가 나왔다.
    ① 한쪽 HAP 에서만 unique: 다른 염색체 pseudogene 의 변이 때문에 유일성이 바뀐다. gene 위엔 차이가 없다.
       (RPS15 3:707 — 같은 read 의 AS 가 두 HAP 에서 같았다. 이진 ASE 적격 211 중 82 개가 이 경로에 절반 이상 의존)
    ② VCF 오류: 실제로 없는 variant 가 HAP 서열에 들어가 그 위 read 의 점수를 깎는다.
       (PKN1 chr19:14,469,962 TC>T — HiFi 0/35, RNA 0/400 read 에 결실 없음. 발현 gene exon het indel 의 7.4%)
  그래서 "실제 WGS read 로 확인된 HAP variant 를 덮는 fragment 만 ASE 근거로 인정"하도록, 그 variant 목록을 만든다.

★ 판정 (variant 마다)
  not_applied  HAP FASTA 의 그 위치 서열이 VCF allele 과 다르다(consensus 가 겹침 등으로 적용하지 않음) → 근거 아님
  low_coverage 덮는 WGS read < --min-depth
  unsupported  alt 를 가진 read 비율 < --min-alt-frac
  validated    그 외
  read 가 alt 를 가짐: SNV/MNV 는 그 위치 염기가 alt, indel 은 같은 종류·같은 길이 indel 이 VCF 위치 ±10 bp 안에 있음
  (반복 서열에서 정렬기마다 indel 위치가 달라지는 것을 허용). indel 은 앞뒤 10 bp 까지 덮는 read 만 센다.
  기준값 근거: HG00096 HiFi 에서 exon het indel alt 비율이 0 근처 봉우리와 0.4–0.5 중심 분포로 갈라졌다(0.1 에서 분리).
★ 범위: GENCODE exon ± --flank 안의 het variant (fragment 의 두 mate 가 닿는 범위). 범위 밖 variant 는 표에 없고,
  21 은 표에 없는 variant 를 근거로 쓰지 않는다.
★ 출력 좌표: 각 variant 가 HAP1·HAP2 개인 좌표에서 차지하는 구간(0-based half-open, chain 으로 변환).
  21 이 fragment 정렬 블록과 바로 겹침을 본다.
"""
from __future__ import annotations

import argparse, bisect, gzip, hashlib, json, os, re, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,str(Path(__file__).resolve().parent/"src"))
from rna_common import log, die, write_json  # noqa: E402

HET = {"0|1": (0, 1), "1|0": (1, 0), "1|.": (1, 0), ".|1": (0, 1)}
INDEL_SLOP = 10
COLS = ["contig", "pos", "ref", "alt", "gt", "kind", "status", "wgs_depth", "wgs_alt", "wgs_alt_frac",
        "hp1_alt", "hp1_depth", "hp2_alt", "hp2_depth", "h1_start", "h1_end", "h1_allele", "h2_start", "h2_end", "h2_allele", "genes"]


def sha256_file(path: str, chunk: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def load_chain(path: Path) -> np.ndarray:
    """bcftools consensus chain → (ref_start, hap_start, size) block 배열."""
    out = []
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if not p:
                continue
            if p[0] == "chain":
                t, q = int(p[5]), int(p[10]); continue
            s = int(p[0]); out.append((t, q, s))
            if len(p) == 3:
                t += s + int(p[1]); q += s + int(p[2])
    return np.array(out, dtype=np.int64)


def ref2hap(blocks: np.ndarray, x: int) -> int:
    i = max(0, int(np.searchsorted(blocks[:, 0], x, side="right")) - 1)
    t, q, s = (int(v) for v in blocks[i])
    return q + (x - t) if x < t + s else q + s


def exon_regions(gtf: str, contigs: set[str], flank: int):
    """contig → (merged regions [(s,e)], exon list [(s,e,gene)])."""
    ex: dict[str, list] = {}
    with gzip.open(gtf, "rt") if gtf.endswith(".gz") else open(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 9)
            if f[2] != "exon" or f[0] not in contigs:
                continue
            gid = re.search(r'gene_id "([^"]+)"', f[8])[1]
            ex.setdefault(f[0], []).append((int(f[3]) - 1, int(f[4]), gid))
    out = {}
    for c, lst in ex.items():
        lst.sort()
        merged = []
        for s, e, _ in lst:
            s, e = max(0, s - flank), e + flank
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[c] = (merged, lst)
    return out


def genes_near(exons, pos0: int, flank: int, starts) -> str:
    i = bisect.bisect_right(starts, pos0 + flank)
    hit = {g for s, e, g in exons[max(0, i - 200):i] if s - flank <= pos0 < e + flank}
    return ",".join(sorted(hit))


def count_support(bam, ref, c, region, variants, min_mapq):
    """region 안 variant 들의 WGS 지지 read 수. read 를 한 번만 훑는다."""
    s0, e0 = region
    vpos = [v["p0"] for v in variants]
    for r in bam.fetch(c, s0, e0):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail \
                or r.mapping_quality < min_mapq:
            continue
        rs, re_ = r.reference_start, r.reference_end
        i = bisect.bisect_left(vpos, rs - INDEL_SLOP); j = bisect.bisect_right(vpos, re_)
        if i >= j:
            continue
        hp = r.get_tag("HP") if r.has_tag("HP") else 0
        events, base_at = [], None
        rp = rs
        for op, ln in r.cigartuples:
            if op in (0, 7, 8):
                rp += ln
            elif op == 2:
                events.append(("D", rp, ln)); rp += ln
            elif op == 3:
                events.append(("N", rp, ln)); rp += ln
            elif op == 1:
                events.append(("I", rp, ln))
        for v in variants[i:j]:
            p0 = v["p0"]
            if v["kind"] == "INDEL":
                if rs > p0 - INDEL_SLOP or re_ < p0 + len(v["ref"]) + INDEL_SLOP:
                    continue
                if any(k == "N" and a <= p0 + 1 < a + ln for k, a, ln in events):
                    continue
                kind = "D" if len(v["ref"]) > len(v["alt"]) else "I"
                ln_v = abs(len(v["ref"]) - len(v["alt"]))
                has = any(k == kind and ln == ln_v and abs(a - (p0 + 1)) <= INDEL_SLOP for k, a, ln in events)
            else:
                if base_at is None:
                    base_at = {rpos: qpos for qpos, rpos in r.get_aligned_pairs(matches_only=True)}
                qs = [base_at.get(p0 + k) for k in range(len(v["ref"]))]
                if any(q is None for q in qs):
                    continue
                seq = r.query_sequence
                has = "".join(seq[q] for q in qs).upper() == v["alt"].upper()
            v["wgs_depth"] += 1; v["wgs_alt"] += has
            if hp in (1, 2):
                v[f"hp{hp}_depth"] += 1; v[f"hp{hp}_alt"] += has


def process_contig(job):
    (c, vcf, bam_path, ref_fa, hap_dirs, regions, exons, flank, min_mapq, min_depth, min_alt) = job[:11]
    tool,decode=job[11:] if len(job)>11 else (None,0)
    groups=[]
    chains = {h: load_chain(Path(hap_dirs[h]) / f"{c}.chain") for h in (1, 2)}
    fastas = {h: pysam.FastaFile(str(Path(hap_dirs[h]) / f"{c}.fa")) for h in (1, 2)}
    starts = [s for s, _, _ in exons]
    rows = []
    with pysam.VariantFile(vcf) as vf, pysam.AlignmentFile(bam_path, reference_filename=ref_fa or None) as bam:
        sample = list(vf.header.samples)[0]
        for reg in regions:
            variants = []
            for rec in vf.fetch(c, reg[0], reg[1]):
                gt = rec.samples[sample]["GT"]
                key = "|".join("." if x is None else str(x) for x in gt)
                if key not in HET or not rec.alts or len(rec.alts) != 1:
                    continue
                ref, alt = rec.ref, rec.alts[0]
                if alt.startswith("<") or not re.fullmatch(r"[ACGTNacgtn]+", alt):
                    continue
                kind = "SNV" if len(ref) == len(alt) == 1 else ("MNV" if len(ref) == len(alt) else "INDEL")
                v = dict(contig=c, pos=rec.pos, p0=rec.pos - 1, ref=ref, alt=alt, gt=key, kind=kind,
                         wgs_depth=0, wgs_alt=0, hp1_alt=0, hp1_depth=0, hp2_alt=0, hp2_depth=0)
                applied = True
                for h in (1, 2):
                    allele = alt if HET[key][h - 1] == 1 else ref
                    hs = ref2hap(chains[h], rec.pos - 1)
                    he = hs + len(allele)
                    v[f"h{h}_start"], v[f"h{h}_end"], v[f"h{h}_allele"] = hs, he, allele.upper()
                    if fastas[h].fetch(c, hs, he).upper() != allele.upper():
                        applied = False
                v["applied"] = applied
                variants.append(v)
            variants.sort(key=lambda d: d["p0"])
            if tool:groups.append((reg,variants))
            elif variants:count_support(bam, ref_fa, c, reg, variants, min_mapq)
            rows.extend(variants)
        if tool:
            from integrate_wgs_rna.rna_wgs_native import count_support as native_support
            native_support(tool,bam_path,ref_fa,c,groups,min_mapq,decode)
        for v in rows:
            d, a = v["wgs_depth"], v["wgs_alt"]
            v["wgs_alt_frac"] = round(a / d, 4) if d else ""
            v["status"] = ("not_applied" if not v["applied"] else "low_coverage" if d < min_depth
                           else "unsupported" if a / d < min_alt else "validated")
            v["genes"] = genes_near(exons, v["p0"], flank, starts)
    return c, rows


_VARIANT_JOBS={}


def process_chunk(contig,indices,tool,decode):
    import time
    from integrate_wgs_rna.rna_task_chunks import cpu_seconds,timing
    started=time.monotonic();cpu=cpu_seconds();base=_VARIANT_JOBS[contig]
    job=base[:5]+([base[5][i] for i in indices],)+base[6:]+(tool,decode)
    c,rows=process_contig(job)
    stat=timing(started,cpu,contig=c,regions=len(indices),variants=len(rows),begin=min(base[5][i][0] for i in indices),end=max(base[5][i][1] for i in indices))
    return c,rows,stat


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="09 가 HAP FASTA 를 만든 phased VCF (personal_fasta_qc.json 의 refok_vcf)")
    ap.add_argument("--wgs-bam", required=True, help="06 haplotagged WGS BAM/CRAM (GRCh38)")
    ap.add_argument("--reference", default="", help="WGS 가 CRAM 이면 decode reference")
    ap.add_argument("--hap-dir1", required=True); ap.add_argument("--hap-dir2", required=True)
    ap.add_argument("--gtf", required=True); ap.add_argument("--contigs", required=True)
    ap.add_argument("--flank", type=int, default=1000)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-depth", type=int, default=10)
    ap.add_argument("--min-alt-frac", type=float, default=0.1)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", required=True); ap.add_argument("--qc-json", required=True)
    a = ap.parse_args()
    for p in (a.vcf, a.wgs_bam, a.gtf):
        if not Path(p).is_file():
            die(f"입력 없음: {p}")
    contigs = a.contigs.split()
    log(f"20 HAP variant 검증: exon ±{a.flank} bp, depth ≥ {a.min_depth}, alt ≥ {a.min_alt_frac}")
    regs = exon_regions(a.gtf, set(contigs), a.flank)
    hap_dirs = {1: a.hap_dir1, 2: a.hap_dir2}
    jobs = [(c, a.vcf, a.wgs_bam, a.reference, hap_dirs, regs[c][0], regs[c][1], a.flank, a.min_mapq,
             a.min_depth, a.min_alt_frac) for c in contigs if c in regs]
    if not jobs:
        die("요청 contig와 겹치는 GTF exon이 없어 HAP variant를 검증할 수 없음")
    from integrate_wgs_rna.rna_wgs_native import select,plan
    from integrate_wgs_rna.cpu_executor import cpu_results
    tool=select('RNA_HAP_VARIANT_BACKEND','rna_variant_support')
    out={}
    execution={}
    if tool:
        from multiprocessing import get_context
        from integrate_wgs_rna.rna_task_chunks import partition_rows,task_layout
        global _VARIANT_JOBS
        _VARIANT_JOBS={job[0]:job for job in jobs}
        groups=[(job[0],indices) for job in jobs for indices in partition_rows(job[5],lambda row:row)]
        layout=task_layout();budget=plan(len(groups),cpu_request=a.threads,decode_threads=0 if layout=='chunks' else None)
        tasks=[(c,indices,tool,budget['decode_threads']) for c,indices in groups]
        completed={};traces=[]
        try:
            for index,result,error in cpu_results(process_chunk,tasks,budget['workers_initial'],elastic=True,
                    mp_context=get_context('fork'),max_workers=budget['workers_upper'],cpus_per_worker=budget['cpu_per_worker'],
                    memory_per_worker_bytes=budget['worker_address_space_bytes'],retain_worker_memory=True,parent_cpus=budget['parent_cpus']):
                if error:raise error
                c,rows,trace=result;completed[index]=(c,rows);traces.append(trace)
            for index in range(len(tasks)):
                c,rows=completed[index];out.setdefault(c,[]).extend(rows)
        finally:_VARIANT_JOBS={}
        execution=dict(layout=layout,budget=budget,tasks=traces)
        log('RNA-VARIANT-TASKS '+json.dumps(execution))
    else:
        out = {}
        # autosome 수보다 많은 process를 만들지 않는다. 결과 순서는 ex.map과 아래 contig 순회로 고정된다.
        with ProcessPoolExecutor(max_workers=min(len(jobs), max(1, a.threads))) as ex:
            for c, rows in ex.map(process_contig, jobs):
                out[c] = rows
                log(f"  {c}: het variant {len(rows):,}")
    tmp = Path(a.out + ".part")
    stat: dict = {}
    with open(tmp, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for c in contigs:
            for v in out.get(c, []):
                fh.write("\t".join(str(v[k]) for k in COLS) + "\n")
                stat.setdefault(v["kind"], {}).setdefault(v["status"], 0)
                stat[v["kind"]][v["status"]] += 1
    os.replace(tmp, a.out)
    n = sum(sum(d.values()) for d in stat.values())
    write_json(Path(a.qc_json), {
        "schema": "mfm_hap_variants_qc_v1",
        "inputs": {"vcf": a.vcf, "vcf_sha256": sha256_file(a.vcf), "wgs_bam": a.wgs_bam,
                   "hap_dir1": a.hap_dir1, "hap_dir2": a.hap_dir2, "gtf": a.gtf},
        "thresholds": {"flank": a.flank, "min_mapq": a.min_mapq, "min_depth": a.min_depth,
                       "min_alt_frac": a.min_alt_frac, "indel_slop_bp": INDEL_SLOP},
        "resource_execution": execution,
        "variants": n, "by_kind_status": stat,
        "validated_fraction": sum(d.get("validated", 0) for d in stat.values()) / max(1, n),
        "output": a.out, "output_sha256": sha256_file(a.out),
    })
    log(f"20 완료: {n:,} variant, " + ", ".join(f"{k} {json.dumps(v)}" for k, v in stat.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
