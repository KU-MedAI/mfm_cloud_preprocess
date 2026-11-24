#!/usr/bin/env python3
"""보정별 기여도 분해 — BAM 을 한 번만 읽고 모든 arm 의 카운트를 동시에 만든다.

    python tools/ablation.py --bam1 <HAP1.md.bam> --bam2 <HAP2.md.bam> \
        --exons1 <e1.tsv> --exons2 <e2.tsv> --sample <ID> [--cnv-rd <gene_rd.tsv>] \
        --out <ablation.json>

정책마다 16_allelic_counts.py 를 다시 돌리면 BAM 을 arm 수만큼 읽게 된다(샘플당 10 GB × N).
여기서는 mate 별 (AS, blocks, dup, uniq) 를 한 번만 수집해 두고, 배정 규칙만 바꿔 가며
여러 카운트를 만든다. 결과는 16 과 같은 규칙을 쓴다.

arm 정의 — 누적이 아니라 **baseline 대비 단일 변경**이다. 무엇이 무엇을 움직였는지 보려면
한 번에 하나만 바꿔야 한다.
  baseline        기존 동작 (duplicate 포함, 한쪽 hap 포함, uniq only, binomial)
  dup             + duplicate 제외
  both_haps       + 한쪽 hap 만 정렬된 fragment 제외
  mmap            + multimapper 보정으로 방향이 뒤집히는 유전자 제외
  cnv             + CNV 의심 유전자 제외 (rd 파일이 있을 때만)
  betabinom       + overdispersion 보정 (통계만)
  all             위 전부 동시 적용
"""
from __future__ import annotations
import argparse, collections, hashlib, json, math, statistics, sys
from pathlib import Path
import pysam

CODE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CODE))
from rna_ase_stats import (  # noqa: E402
    bh_fdr, betabinom_two_sided, binom_two_sided, fit_overdispersion,
)
from rna_allelic_common import build_exon_index, genes_for_indexed, load_exons  # noqa: E402


def legacy_fit_overdispersion(rows):
    """Ablation baseline은 구번호 Stage 18의 저깊이 null 동작을 유지한다."""
    return fit_overdispersion(rows, include_low_depth_null=True)

MATES = (0, 1)


def collect(bam_path: str, mapq_min: int):
    """query_name → [mate0, mate1]; 각 원소 (AS, contig, blocks, is_dup, is_uniq) 또는 None."""
    out: dict = {}
    n_dup = 0
    with pysam.AlignmentFile(bam_path, "rb") as bf:
        for r in bf.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            mq = r.mapping_quality
            if mq == 0:
                continue
            try:
                a = int(r.get_tag("AS"))
            except KeyError:
                continue
            if r.is_duplicate:
                n_dup += 1
            i = 0 if r.is_read1 else 1
            slot = out.get(r.query_name)
            if slot is None:
                slot = out[r.query_name] = [None, None]
            prev = slot[i]
            if prev is None or a > prev[0]:
                slot[i] = (a, r.reference_name, r.get_blocks(), r.is_duplicate, mq >= mapq_min)
    return out, n_dup


def assign(A, B, ex1, ex2, *, drop_dup: bool, require_both: bool, uniq_only: bool):
    aC, bC = collections.Counter(), collections.Counter()
    st = dict(assigned=0, tie=0, single=0, nocommon=0)

    def keep(m):
        if m is None:
            return None
        if drop_dup and m[3]:
            return None
        if uniq_only and not m[4]:
            return None
        return m

    for key in set(A) | set(B):
        ra = [keep(m) for m in A.get(key, [None, None])]
        rb = [keep(m) for m in B.get(key, [None, None])]
        ha, hb = any(m is not None for m in ra), any(m is not None for m in rb)
        if not ha and not hb:
            continue
        if ha and hb:
            common = [i for i in MATES if ra[i] is not None and rb[i] is not None]
            if not common:
                st["nocommon"] += 1; continue
            s1 = sum(ra[i][0] for i in common); s2 = sum(rb[i][0] for i in common)
            if s1 == s2:
                st["tie"] += 1; continue
            tgt, idx, cnt = (ra, ex1, aC) if s1 > s2 else (rb, ex2, bC)
        else:
            if require_both:
                continue
            st["single"] += 1
            tgt, idx, cnt = (ra, ex1, aC) if ha else (rb, ex2, bC)
        hits = set()
        for i in MATES:
            if tgt[i] is not None:
                hits |= genes_for_indexed(idx, tgt[i][1], tgt[i][2])
        for g in hits:
            cnt[g] += 1
        st["assigned"] += 1
    return aC, bC, st


def evaluate(aC, bC, *, betabinom: bool, exclude: set[str] | None, min_total: int = 20):
    genes = sorted((set(aC) | set(bC)) - (exclude or set()))
    deep = [(g, aC[g], bC[g]) for g in genes if aC[g] + bC[g] >= min_total]
    if not deep:
        return dict(genes=len(genes), genes_deep=0)
    A = sum(x for _, x, _ in deep); B = sum(y for _, _, y in deep)
    rho, diag = (legacy_fit_overdispersion([(aC[g], aC[g] + bC[g]) for g in genes])
                 if betabinom else (0.0, {"status": "not_used"}))
    p = [(betabinom_two_sided(x, y, rho) if betabinom else binom_two_sided(x, y))
         for _, x, y in deep]
    q = bh_fdr(p)
    sig = [g for (g, _, _), v in zip(deep, q) if v <= 0.05]
    return dict(genes=len(genes), genes_deep=len(deep),
                global_ratio=A / (A + B) if A + B else None,
                global_bias_abs=abs(A / (A + B) - 0.5) if A + B else None,
                rho=rho, rho_status=diag.get("status"),
                n_significant=len(sig), is_ase=len(sig) / len(deep),
                excess_over_truth=(len(sig) / len(deep)) / 0.067,
                sig_genes=sig)


def main() -> int:
    ap = argparse.ArgumentParser()
    for k in ("bam1", "bam2", "exons1", "exons2", "sample", "out"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--cnv-rd", default="")
    ap.add_argument("--min-mapq", type=int, default=255)
    a = ap.parse_args()

    print(f"[{a.sample}] BAM 읽는 중...", flush=True)
    A, dupA = collect(a.bam1, a.min_mapq)
    B, dupB = collect(a.bam2, a.min_mapq)
    ex1 = build_exon_index(load_exons(a.exons1))
    ex2 = build_exon_index(load_exons(a.exons2))
    print(f"[{a.sample}] fragment H1 {len(A):,} / H2 {len(B):,} · DUP 레코드 {dupA+dupB:,}", flush=True)

    cnv_bad: set[str] = set()
    if a.cnv_rd and Path(a.cnv_rd).is_file():
        import csv
        with open(a.cnv_rd) as fh:
            cnv_bad = {r["gene_id"] for r in csv.DictReader(fh, delimiter="\t")
                       if r.get("pass_cnv") == "0"}

    base = dict(drop_dup=False, require_both=False, uniq_only=True)
    arms = {
        "baseline":  (dict(base),                          dict(betabinom=False, exclude=None)),
        "dup":       (dict(base, drop_dup=True),           dict(betabinom=False, exclude=None)),
        "both_haps": (dict(base, require_both=True),       dict(betabinom=False, exclude=None)),
        "cnv":       (dict(base),                          dict(betabinom=False, exclude=cnv_bad)),
        "betabinom": (dict(base),                          dict(betabinom=True,  exclude=None)),
        "all":       (dict(drop_dup=True, require_both=True, uniq_only=True),
                                                           dict(betabinom=True,  exclude=cnv_bad)),
    }
    # mmap arm: uniq(MAPQ>=min) 과 mmap(0<MAPQ<min) 을 각각 배정해 방향이 뒤집히는 유전자를 뺀다
    res: dict = {"sample": a.sample, "duplicate_records": dupA + dupB,
                 "cnv_excluded_genes": len(cnv_bad), "arms": {}}
    cache: dict = {}
    for name, (apol, epol) in arms.items():
        key = tuple(sorted(apol.items()))
        if key not in cache:
            cache[key] = assign(A, B, ex1, ex2, **apol)
        aC, bC, st = cache[key]
        out = evaluate(aC, bC, **epol)
        out["assignment"] = st
        res["arms"][name] = out
        print(f"  {name:<10} 유전자 {out.get('genes_deep',0):>6,}  "
              f"ratio {out.get('global_ratio',0):.4f}  "
              f"is_ase {out.get('is_ase',0):.3f}  "
              f"({out.get('excess_over_truth',0):.1f}x)", flush=True)

    # mmap 보정
    uq = cache[tuple(sorted(base.items()))]
    mmA, mmB, _ = assign(A, B, ex1, ex2, drop_dup=False, require_both=False, uniq_only=False)
    flip = set()
    for g in set(uq[0]) | set(uq[1]):
        x, y = uq[0][g], uq[1][g]
        mx, my = max(mmA[g] - x, 0), max(mmB[g] - y, 0)
        if x == y:
            continue
        ax, ay = (x, y + mx + my) if x > y else (x + mx + my, y)
        if not ((x > y and ax > ay) or (x < y and ax < ay)):
            flip.add(g)
    out = evaluate(uq[0], uq[1], betabinom=False, exclude=flip)
    out["assignment"] = uq[2]; out["mmap_flipped_genes"] = len(flip)
    res["arms"]["mmap"] = out
    print(f"  {'mmap':<10} 유전자 {out['genes_deep']:>6,}  ratio {out['global_ratio']:.4f}  "
          f"is_ase {out['is_ase']:.3f}  (뒤집힌 유전자 {len(flip):,})", flush=True)

    Path(a.out).write_text(json.dumps(res, ensure_ascii=False, indent=2) + "\n")
    print(f"[{a.sample}] → {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
