#!/usr/bin/env python3
"""H1/H2 배정 품질을 세 축으로 재는 공통 지표. 어느 실행이든 같은 코드로 통과시킨다.

    python tools/hap_assignment_report.py --label ONT_old --counts <allelic_counts.tsv> ... \
        [--hap1-dir <HAP1>] [--hap2-dir <HAP2>] [--out report.json]

★ 세 질문에 각각 대응한다

  ① H1/H2 read 를 잘 붙였는가 — 배정의 양과 근거
     informative_fraction   AS 비교가 가능했던 fragment 비율. het 밀도 상한이 있어 5–12%가 정상.
     tie_fraction           두 hap 에서 AS 가 같아 버린 비율.
     single_hap_fraction    배정분 중 "한쪽에만 정렬"에서 온 몫. AS 비교가 아니므로 낮을수록 좋다.
     global_ratio           전체 aCount/(aCount+bCount). HAP1=GT[0] 이고 phase block 모자이크이므로
                            **0.5 에서 벗어나면 그 자체가 mapping bias 신호다.** 이 지표가 핵심이다.
     duplicates_excluded    카운트가 duplicate 를 뺐는지.

  ② haplotype 특이 신호를 잡았는가
     coverage_gene_r        유전자 단위 HAP1/HAP2 커버리지 상관. 설계상 ~0.99 여야 정상이고,
                            **커버리지에는 대립 신호가 없다는 사실의 확인**이다(낮아지면 좌표 오정렬 의심).
     ratio_sd / ratio_iqr   대립 비율의 흩어짐. 신호가 있으면 0.5 주변에서 넓어진다.
     overdispersion_rho     잡음 크기. 낮을수록 카운트가 깨끗하다.
     split_half_r           같은 유전자를 fragment 절반씩으로 나눠 잰 비율의 상관. 재현성.

  ③ ASE 가 제대로 나오는가
     is_ase_binom / _betabin   depth>=20 유전자 중 유의 비율. 참조축 진리값 0.067 과 비교한다.
     ase_inflation_ratio       binom / betabin. overdispersion 무시로 부풀린 배수.
"""
from __future__ import annotations
import argparse, csv, json, math, statistics, collections, sys
from pathlib import Path

CODE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CODE))
from rna_ase_stats import (  # noqa: E402
    bh_fdr, betabinom_two_sided, binom_two_sided, fit_overdispersion,
)


def read_counts(p: Path):
    with p.open() as fh:
        return [(r["gene_id"], int(r["aCount"]), int(r["bCount"]))
                for r in csv.DictReader(fh, delimiter="\t")]


def gene_coverage(hapdir: Path, contig: str):
    """exon 합집합의 평균 depth. 좌표계가 hap 마다 다르므로 유전자 단위로만 비교 가능하다."""
    import numpy as np, zarr
    ex = collections.defaultdict(list)
    f = hapdir / "exons.bygene.tsv"
    if not f.is_file():
        return None
    for line in f.open():
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 4 and parts[1] == contig:
            ex[parts[0]].append((int(parts[2]), int(parts[3])))
    z = hapdir / "coverage.zarr"
    if not z.is_dir():
        return None
    arr = np.asarray(zarr.open(str(z), mode="r")[contig][:]).astype(np.float64)
    out = {}
    for g, ivs in ex.items():
        merged = []
        for s, e in sorted(ivs):
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(e, merged[-1][1]))
            else:
                merged.append((s, e))
        tot = 0.0; L = 0
        for s, e in merged:
            s = max(0, s); e = min(len(arr), e)
            if e > s:
                tot += arr[s:e].sum(); L += e - s
        if L >= 200:
            out[g] = tot / L
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--meta", default="")
    ap.add_argument("--split-half", default="")
    ap.add_argument("--hap1-dir", default="")
    ap.add_argument("--hap2-dir", default="")
    ap.add_argument("--contig", default="chr22")
    ap.add_argument("--min-total", type=int, default=20)
    ap.add_argument("--restrict-contig", default="",
                    help="이 contig 의 유전자만 본다. --exons 로 gene→contig 를 읽는다")
    ap.add_argument("--exons", default="", help="gene→contig 매핑용 exons.bygene.tsv")
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    counts = read_counts(Path(a.counts))
    if a.restrict_contig:
        if not a.exons:
            raise SystemExit("--restrict-contig 에는 --exons 가 필요합니다")
        want = {c.strip() for c in a.restrict_contig.replace(",", " ").split() if c.strip()}
        keep = set()
        for line in Path(a.exons).open():
            f = line.rstrip("\n").split("\t")
            if len(f) == 4 and f[1] in want:
                keep.add(f[0])
        before = len(counts)
        counts = [(g, x, y) for g, x, y in counts if g in keep]
        print(f"[restrict] {a.restrict_contig}: 유전자 {before:,} → {len(counts):,}", file=sys.stderr)
    deep = [(g, x, y) for g, x, y in counts if x + y >= a.min_total]
    rep: dict = {"label": a.label, "counts_path": a.counts,
                 "genes_total": len(counts), "genes_deep": len(deep)}

    # ── ① 배정 ──
    meta = {}
    if a.meta and Path(a.meta).is_file():
        meta = json.loads(Path(a.meta).read_text())
    tot_f = int(meta.get("fragments_total", 0))
    asg = int(meta.get("fragments_assigned", 0))
    single = int(meta.get("fragments_only_hap1", 0)) + int(meta.get("fragments_only_hap2", 0))
    A = sum(x for _, x, _ in counts); B = sum(y for _, _, y in counts)
    # 전역비는 소수의 고발현 유전자가 지배할 수 있다(chr22 의 IGL 좌위가 대표적이다.
    # LCL 은 준클론성이라 면역글로불린이 대립 배제로 단일대립 발현한다 — 진짜 생물학이지만
    # 전역비로는 "배정 편향"처럼 보인다). 그래서 견고한 중심값을 함께 낸다.
    deep_sorted = sorted(deep, key=lambda r: -(r[1] + r[2]))
    ratios_all = [x / (x + y) for _, x, y in deep] if deep else []
    top5 = deep_sorted[:5]
    t5 = sum(x + y for _, x, y in top5)
    tot_deep = sum(x + y for _, x, y in deep)
    rest = deep_sorted[5:]
    rest_a = sum(x for _, x, _ in rest); rest_b = sum(y for _, _, y in rest)
    rep["robust"] = {
        "median_gene_ratio": statistics.median(ratios_all) if ratios_all else None,
        "median_gene_bias_abs": abs(statistics.median(ratios_all) - 0.5) if ratios_all else None,
        "top5_count_share": t5 / tot_deep if tot_deep else None,
        "top5_genes": [g for g, _, _ in top5],
        "ratio_excluding_top5": rest_a / (rest_a + rest_b) if (rest_a + rest_b) else None,
        "note": ("median_gene_ratio 가 0.5 에서 벗어나야 배정 편향이다. "
                 "global_ratio 만 벗어나면 소수 유전자의 실제 단일대립 발현일 수 있다."),
    }
    rep["assignment"] = {
        "fragments_total": tot_f, "fragments_assigned": asg,
        "informative_fraction": asg / tot_f if tot_f else None,
        "tie_fraction": int(meta.get("fragments_tie", 0)) / tot_f if tot_f else None,
        "single_hap_fraction_of_assigned": single / asg if asg else None,
        "global_aCount": A, "global_bCount": B,
        "global_ratio": A / (A + B) if A + B else None,
        "global_bias_abs": abs(A / (A + B) - 0.5) if A + B else None,
        "duplicates_excluded": (not meta.get("include_duplicates", True)
                                if "include_duplicates" in meta else None),
        "script_revision": meta.get("script_revision"),
    }

    # ── ② haplotype 특이 신호 ──
    ratios = [x / (x + y) for _, x, y in deep]
    # 이 보고서는 기존 ablation과 같이 구번호 Stage 18의 baseline 적합을 유지한다.
    rho, rho_diag = fit_overdispersion(
        [(x, x + y) for _, x, y in counts], include_low_depth_null=True
    )
    sig = {
        "ratio_sd": statistics.pstdev(ratios) if len(ratios) > 1 else None,
        "ratio_iqr": (statistics.quantiles(ratios, n=4)[2] - statistics.quantiles(ratios, n=4)[0])
                     if len(ratios) > 3 else None,
        "median_abs_effect": statistics.median(abs(r - 0.5) for r in ratios) if ratios else None,
        "overdispersion_rho": rho, "overdispersion_status": rho_diag.get("status"),
    }
    if a.split_half and Path(a.split_half).is_file():
        with Path(a.split_half).open() as fh:
            xs, ys = [], []
            for r in csv.DictReader(fh, delimiter="\t"):
                t0, t1 = int(r["split0_totalCount"]), int(r["split1_totalCount"])
                if t0 >= 5 and t1 >= 5:
                    xs.append(float(r["split0_ratio"])); ys.append(float(r["split1_ratio"]))
        if len(xs) > 2:
            mx, my = statistics.fmean(xs), statistics.fmean(ys)
            dx = [v - mx for v in xs]; dy = [v - my for v in ys]
            den = math.sqrt(sum(v * v for v in dx) * sum(v * v for v in dy))
            sig["split_half_r"] = (sum(p * q for p, q in zip(dx, dy)) / den) if den else None
            sig["split_half_genes"] = len(xs)
    if a.hap1_dir and a.hap2_dir:
        try:
            c1 = gene_coverage(Path(a.hap1_dir), a.contig)
            c2 = gene_coverage(Path(a.hap2_dir), a.contig)
            if c1 and c2:
                import numpy as np
                common = sorted(set(c1) & set(c2))
                u = np.array([c1[g] for g in common]); v = np.array([c2[g] for g in common])
                m = (u > 1) & (v > 1)
                if m.sum() > 2:
                    sig["coverage_gene_r"] = float(np.corrcoef(u[m], v[m])[0, 1])
                    sig["coverage_median_abs_log2"] = float(np.median(np.abs(np.log2(u[m] / v[m]))))
                    sig["coverage_genes"] = int(m.sum())
                    sig["coverage_contig"] = a.contig
        except Exception as e:                       # zarr/numpy 없거나 산출물이 없을 수 있다
            sig["coverage_error"] = str(e)
    rep["haplotype_signal"] = sig

    # ── ③ ASE ──
    pb = [binom_two_sided(x, y) for _, x, y in deep]
    pbb = [betabinom_two_sided(x, y, rho) for _, x, y in deep]
    qb, qbb = bh_fdr(pb), bh_fdr(pbb)
    nb = sum(v <= 0.05 for v in qb); nbb = sum(v <= 0.05 for v in qbb)
    rep["ase"] = {
        "significant_binom": nb, "significant_betabinom": nbb,
        "is_ase_binom": nb / len(deep) if deep else None,
        "is_ase_betabinom": nbb / len(deep) if deep else None,
        "ase_inflation_ratio": nb / nbb if nbb else None,
        "reference_truth": 0.067,
        "excess_over_truth": (nbb / len(deep) / 0.067) if deep and nbb else None,
    }

    txt = json.dumps(rep, ensure_ascii=False, indent=2)
    if a.out:
        Path(a.out).write_text(txt + "\n")
    print(txt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
