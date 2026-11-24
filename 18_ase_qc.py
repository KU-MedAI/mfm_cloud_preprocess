#!/usr/bin/env python3
"""[구번호] fragment 기반 ASE label의 깊이·재현성·민감도 QC.

현재 v2 Stage 23은 `23_ase_qc.py`이며 Stage 20 variant validation 및 현행 Stage 21
allelic schema를 요구한다. 이 파일은 이전 ASE QC 입력·통계를 재현하도록 보존한다.

연속형 ratio 학습과 이진 ASE 판정을 분리한다. 연속형 label에는 FDR을 gate로
쓰지 않는다. 그렇게 하면 균형 발현 유전자가 제거되어 학습 표본이 편향되기 때문이다.

★ 2026-09-17 AlleleSeq2 대조 반영 — beta-binomial overdispersion
  plain binomial은 overdispersion 때문에 ASE를 과대 호출한다. AlleleSeq2/AlleleDB가
  `alleledb_calcOverdispersion.R` + `alleledb_alleleseqBetabinomial.R`를 두는 이유가 이것이고,
  두 산출물(binom / betabinom)을 **둘 다** 기본으로 만든다(PIPELINE.mk의 `all:` 타깃).

  여기서도 같은 구조를 만든다. 추정 절차는 AlleleSeq2와 같다 — MLE가 아니라
  **비율 히스토그램에 대한 가중 SSE 적합**이다:
    1. 관측 ratio를 40구간 히스토그램으로 만든다(minN=6 이상, maxN=min(max,2500) 이하).
    2. 후보 rho마다, 관측된 depth 분포 w[n]으로 가중한 beta-binomial 혼합분포를
       같은 구간에 넣어 이론 히스토그램을 만든다.
    3. 가운데(ratio 0.5 근처)에 큰 가중을 주는 삼각 가중으로 SSE를 재고 최소가 되는 rho를 고른다.
       양 끝을 낮게 누르는 것은 그쪽이 ASE 신호이고 지금 맞추려는 것은 null이기 때문이다.
    4. 굵은 격자(0.1)로 훑고 그 주변을 조밀하게 다시 훑는다.
  p-value도 AlleleSeq2와 같은 형태다: `min(1, 2 * pbetabinom(minor, total, 0.5, rho))`.
  parameterization은 VGAM과 같다 — alpha = p(1-rho)/rho, beta = (1-p)(1-rho)/rho.

  ⚠ AlleleSeq2와 다른 점: multiple testing 보정을 AlleleSeq2의 시뮬레이션 FDR(FalsePos.py,
    500회) 대신 BH로 한다. 난수 없이 재현되고 기존 산출물과 컬럼이 이어지기 때문이다.
    과대 호출의 주된 원인인 overdispersion 자체는 위 모델로 교정된다.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


from rna_ase_stats import (
    BB_BINS, BB_MAX_N, BB_MIN_GENES_FOR_FIT, BB_MIN_N,
    _binned_empirical, _binned_null, _graded_weights, bh_fdr,
    betabinom_pmf_series, betabinom_two_sided, binom_two_sided,
    flag, pearson, read_tsv, wilson_interval,
    fit_overdispersion as _fit_overdispersion,
)

def fit_overdispersion(rows: list[tuple[int, int]]) -> tuple[float, dict]:
    """구번호 Stage 18의 저깊이 null 포함 방식을 보존한다."""
    return _fit_overdispersion(rows, include_low_depth_null=True)

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--split-half")
    ap.add_argument("--both-haps")
    ap.add_argument("--meta")
    ap.add_argument("--mmap-adjusted", help="33의 .mmap_adjusted.tsv. 있으면 multimapper로 "
                                            "설명되는 불균형을 reason code로 표시한다")
    ap.add_argument("--cnv-rd", help="현재 v2 Stage 22 CNV의 gene_rd.tsv. 있으면 의심 유전자를 제외한다")
    ap.add_argument("--min-cnv-rd", type=float, default=0.5)
    ap.add_argument("--max-cnv-rd", type=float, default=1.5)
    ap.add_argument("--out-qc")
    ap.add_argument("--out-json")
    ap.add_argument("--out-eligibility")
    ap.add_argument("--min-total-count", type=int, default=20)
    ap.add_argument("--min-split-count", type=int, default=5)
    ap.add_argument("--max-split-ratio-delta", type=float, default=0.25)
    ap.add_argument("--max-fdr", type=float, default=0.05)
    ap.add_argument("--min-informative-fraction", type=float, default=0.01)
    ap.add_argument("--min-eligible-genes", type=int, default=100)
    ap.add_argument("--min-split-correlation", type=float, default=0.30)
    args = ap.parse_args()

    counts_path = Path(args.counts)
    prefix = counts_path.with_suffix("")
    split_path = Path(args.split_half) if args.split_half else prefix.with_name(prefix.name + ".split_half.tsv")
    both_path = Path(args.both_haps) if args.both_haps else prefix.with_name(prefix.name + ".both_haps.tsv")
    meta_path = Path(args.meta) if args.meta else prefix.with_name(prefix.name + ".meta.json")
    out_qc = Path(args.out_qc) if args.out_qc else counts_path.parent / "ase_qc.tsv"
    out_json = Path(args.out_json) if args.out_json else counts_path.parent / "ase_qc.json"
    out_elig = Path(args.out_eligibility) if args.out_eligibility else counts_path.parent / "ase_training_eligibility.tsv"

    for required in (counts_path, split_path, both_path, meta_path):
        if not required.is_file():
            raise FileNotFoundError(f"ASE QC 필수 입력 없음: {required}")

    counts = read_tsv(counts_path)
    split = read_tsv(split_path)
    both = read_tsv(both_path)
    meta = json.loads(meta_path.read_text())
    if meta.get("counting_unit") != "fragment":
        raise ValueError("read 단위 또는 단위 미상 산출물은 학습 QC에 사용할 수 없습니다")

    mmap_rows: dict[str, dict[str, str]] = {}
    if args.mmap_adjusted:
        mmap_path = Path(args.mmap_adjusted)
        if not mmap_path.is_file():
            raise FileNotFoundError(f"mmap 보정 파일 없음: {mmap_path}")
        mmap_rows = read_tsv(mmap_path)
    cnv_rows: dict[str, dict[str, str]] = {}
    if args.cnv_rd:
        cnv_path = Path(args.cnv_rd)
        if not cnv_path.is_file():
            raise FileNotFoundError(f"CNV rd 파일 없음: {cnv_path}")
        cnv_rows = read_tsv(cnv_path)

    gene_ids = sorted(counts)
    p_values = [binom_two_sided(int(counts[g]["aCount"]), int(counts[g]["bCount"]))
                for g in gene_ids]
    q_values = bh_fdr(p_values)
    # beta-binomial — overdispersion을 sample 단위로 추정한 뒤 같은 형태의 양측 p를 만든다.
    rho, rho_diag = fit_overdispersion(
        [(int(counts[g]["aCount"]), int(counts[g]["aCount"]) + int(counts[g]["bCount"]))
         for g in gene_ids])
    bb_p_values = [betabinom_two_sided(int(counts[g]["aCount"]), int(counts[g]["bCount"]), rho)
                   for g in gene_ids]
    bb_q_values = bh_fdr(bb_p_values)
    records: list[dict[str, object]] = []
    corr_x: list[float] = []
    corr_y: list[float] = []

    for gene, p_value, q_value, bb_p, bb_q in zip(
            gene_ids, p_values, q_values, bb_p_values, bb_q_values):
        row = counts[gene]
        s = split.get(gene, {})
        b = both.get(gene, {})
        a_count, b_count = int(row["aCount"]), int(row["bCount"])
        total = a_count + b_count
        ratio = a_count / total if total else 0.5
        s0_total = int(s.get("split0_totalCount", 0))
        s1_total = int(s.get("split1_totalCount", 0))
        s0_ratio = float(s.get("split0_ratio", 0.5))
        s1_ratio = float(s.get("split1_ratio", 0.5))
        split_delta = abs(s0_ratio - s1_ratio)
        split_direction = ((s0_ratio >= 0.5 and s1_ratio >= 0.5) or
                           (s0_ratio <= 0.5 and s1_ratio <= 0.5))
        pass_depth = total >= args.min_total_count
        pass_split_depth = s0_total >= args.min_split_count and s1_total >= args.min_split_count
        pass_stability = pass_split_depth and split_delta <= args.max_split_ratio_delta
        # CNV 의심 구간 — AlleleSeq2 filter_cnv_sites.py 와 같은 [0.5, 1.5] 구간.
        cnv_row = cnv_rows.get(gene, {})
        cnv_rd = float(cnv_row["rd"]) if cnv_row.get("rd") not in (None, "", "NA") else None
        pass_cnv = cnv_rd is None or args.min_cnv_rd <= cnv_rd <= args.max_cnv_rd
        # multimapper 로 설명되는 불균형 — AlleleSeq2 filter_sites_w_mmaps.py adjust 대응.
        mmap_row = mmap_rows.get(gene, {})
        mmap_survives = mmap_row.get("mmap_survives")
        pass_mmap = mmap_survives is None or mmap_survives == "1"
        continuous_eligible = pass_depth and pass_stability and pass_cnv
        # AlleleSeq2는 binom과 betabinom을 둘 다 내지만 overdispersion을 반영한 쪽이 본선이다.
        significant = bb_q <= args.max_fdr
        significant_binom = q_value <= args.max_fdr
        binary_ase_eligible = (continuous_eligible and significant
                               and split_direction and pass_mmap)
        if pass_split_depth:
            corr_x.append(s0_ratio)
            corr_y.append(s1_ratio)
        both_total = int(b.get("totalCount", 0))
        both_ratio = float(b.get("ratio", 0.5)) if both_total else 0.5
        low, high = wilson_interval(a_count, total)
        reasons = []
        if not pass_depth: reasons.append("low_total_count")
        if not pass_split_depth: reasons.append("low_split_count")
        elif not pass_stability: reasons.append("unstable_split_ratio")
        if not pass_cnv: reasons.append("cnv_region")
        if not pass_mmap: reasons.append("mmap_explains_imbalance")
        records.append({
            "gene_id": gene, "aCount": a_count, "bCount": b_count,
            "totalCount": total, "ratio": ratio, "ratio_ci95_low": low,
            "ratio_ci95_high": high, "effect_abs": abs(ratio - 0.5),
            "binom_p": p_value, "fdr_q": q_value,
            "betabinom_p": bb_p, "betabinom_fdr_q": bb_q,
            "ase_significant": significant,
            "ase_significant_binom": significant_binom,
            "cnv_rd": cnv_rd if cnv_rd is not None else "",
            "pass_cnv": pass_cnv, "pass_mmap_adjustment": pass_mmap,
            "split0_totalCount": s0_total, "split0_ratio": s0_ratio,
            "split1_totalCount": s1_total, "split1_ratio": s1_ratio,
            "split_ratio_delta": split_delta,
            "split_direction_consistent": split_direction,
            "both_haps_totalCount": both_total, "both_haps_ratio": both_ratio,
            "single_hap_sensitivity_delta": abs(ratio - both_ratio),
            "pass_depth": pass_depth, "pass_split_depth": pass_split_depth,
            "pass_split_stability": pass_stability,
            "continuous_eligible": continuous_eligible,
            "binary_ase_eligible": binary_ase_eligible,
            "reason_codes": ";".join(reasons) if reasons else "PASS",
        })

    split_corr = pearson(corr_x, corr_y)
    fragments_total = int(meta.get("fragments_total", 0))
    fragments_assigned = int(meta.get("fragments_assigned", 0))
    informative_fraction = fragments_assigned / fragments_total if fragments_total else 0.0
    continuous_count = sum(bool(r["continuous_eligible"]) for r in records)
    binary_count = sum(bool(r["binary_ase_eligible"]) for r in records)
    provenance_ok = (meta.get("schema") == "mfm_allelic_counts_meta_v3" and
                     meta.get("script_revision") == "2026-09-17_alleleseq2_dup_mmap")
    # AlleleSeq2는 카운트 전에 duplicate를 제거한다. flag만 남은 산출물을 그대로 쓰면
    # totalCount가 부풀어 검정이 과도하게 유의해지므로 sample 단위로 막는다.
    duplicates_excluded = not bool(meta.get("include_duplicates", False))
    sample_reasons = []
    if not provenance_ok: sample_reasons.append("invalid_provenance")
    if not duplicates_excluded: sample_reasons.append("duplicates_included")
    if rho_diag.get("status") != "ok": sample_reasons.append("overdispersion_not_estimable")
    if informative_fraction < args.min_informative_fraction: sample_reasons.append("low_informative_fraction")
    if continuous_count < args.min_eligible_genes: sample_reasons.append("too_few_eligible_genes")
    if split_corr is None or split_corr < args.min_split_correlation: sample_reasons.append("low_split_correlation")
    sample_ready = not sample_reasons

    qc_fields = list(records[0]) if records else []
    with out_qc.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=qc_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow({k: flag(v) if isinstance(v, bool) else v for k, v in record.items()})

    eligibility_fields = ["gene_id", "sample_training_ready", "continuous_eligible",
                          "binary_ase_eligible", "ase_significant", "reason_codes"]
    with out_elig.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=eligibility_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow({
                "gene_id": record["gene_id"], "sample_training_ready": flag(sample_ready),
                "continuous_eligible": flag(sample_ready and bool(record["continuous_eligible"])),
                "binary_ase_eligible": flag(sample_ready and bool(record["binary_ase_eligible"])),
                "ase_significant": flag(bool(record["ase_significant"])),
                "reason_codes": record["reason_codes"],
            })

    single_hap = int(meta.get("fragments_only_hap1", 0)) + int(meta.get("fragments_only_hap2", 0))
    summary = {
        "schema": "mfm_ase_qc_v2",
        "inputs": {"counts": str(counts_path), "split_half": str(split_path),
                   "both_haps": str(both_path), "meta": str(meta_path),
                   "mmap_adjusted": args.mmap_adjusted or None,
                   "cnv_rd": args.cnv_rd or None},
        "thresholds": {k: getattr(args, k) for k in (
            "min_total_count", "min_split_count", "max_split_ratio_delta", "max_fdr",
            "min_informative_fraction", "min_eligible_genes", "min_split_correlation",
            "min_cnv_rd", "max_cnv_rd")},
        "statistics": {
            "primary_test": "beta-binomial (overdispersion 보정) + BH FDR",
            "secondary_test": "exact binomial + BH FDR (비교용으로 함께 보고)",
            "overdispersion_rho": rho,
            "overdispersion_fit": rho_diag,
            "reference": ("AlleleSeq2 alleledb_calcOverdispersion.R / "
                          "alleledb_alleleseqBetabinomial.R"),
            "deviation": ("AlleleSeq2는 FalsePos.py 시뮬레이션 FDR(500회)을 쓰지만 "
                          "여기서는 난수 없이 재현되는 BH를 쓴다"),
            "genes_significant_betabinom": sum(bool(r["ase_significant"]) for r in records),
            "genes_significant_binom": sum(bool(r["ase_significant_binom"]) for r in records),
        },
        "sample_qc": {
            "sample_training_ready": sample_ready,
            "reason_codes": sample_reasons or ["PASS"],
            "provenance_ok": provenance_ok,
            "duplicates_excluded": duplicates_excluded,
            "genes_excluded_cnv": sum(not bool(r["pass_cnv"]) for r in records),
            "genes_mmap_explained": sum(not bool(r["pass_mmap_adjustment"]) for r in records),
            "fragments_total": fragments_total, "fragments_assigned": fragments_assigned,
            "informative_fraction": informative_fraction,
            "tie_fraction": int(meta.get("fragments_tie", 0)) / max(fragments_total, 1),
            "no_common_mate_fraction": int(meta.get("fragments_no_common_mate", 0)) / max(fragments_total, 1),
            "single_hap_fraction": single_hap / max(fragments_total, 1),
            "genes_total": len(records), "genes_continuous_eligible": continuous_count,
            "genes_binary_ase_eligible": binary_count,
            "split_half_pearson": split_corr,
            "genes_in_split_correlation": len(corr_x),
            "median_single_hap_sensitivity_delta": (statistics.median(
                float(r["single_hap_sensitivity_delta"]) for r in records) if records else None),
            "genes_single_hap_sensitivity_delta_gt_0_10": sum(
                float(r["single_hap_sensitivity_delta"]) > 0.10 for r in records),
        },
        "interpretation": {
            "continuous_eligible": "최소 깊이·split-half 안정성·비CNV를 통과한 연속형 ASE ratio label",
            "binary_ase_eligible": ("연속형 QC + beta-binomial FDR + split 방향 일관성 + "
                                    "multimapper 보정 생존을 통과한 이진 ASE-positive label"),
            "ase_significant": "beta-binomial q <= max_fdr (본선)",
            "ase_significant_binom": "exact binomial q <= max_fdr (비교용; 보통 더 많이 잡힌다)",
            "warning": "FDR로 연속형 label을 거르면 균형 발현 음성군이 제거되므로 사용하지 않는다",
        },
    }
    out_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(f"[ase-qc] sample_ready={sample_ready} informative={informative_fraction:.4f} "
          f"continuous={continuous_count}/{len(records)} binary_ASE={binary_count} "
          f"split_r={split_corr if split_corr is not None else 'NA'}")
    print(f"[ase-qc] overdispersion rho={rho:.4f} ({rho_diag.get('status')}, "
          f"적합에 쓴 유전자 {rho_diag.get('genes_used_for_fit', 0):,}) / "
          f"유의 betabinom={sum(bool(r['ase_significant']) for r in records):,} "
          f"vs binom={sum(bool(r['ase_significant_binom']) for r in records):,}")
    print(f"[ase-qc] {out_qc} / {out_elig} / {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
