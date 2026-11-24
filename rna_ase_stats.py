#!/usr/bin/env python3
"""RNA ASE 단계들이 공유하는 통계 primitive.

단계별 입력 계약과 QC 판정은 각 stage 파일에 둔다. 여기에는 수식이 동일한
저수준 계산만 모아 두고, 구버전 Stage 18의 overdispersion depth 집계 차이는
`include_low_depth_null=True`로 명시해 결과를 보존한다.
"""
from __future__ import annotations

import csv
import math
import statistics
from pathlib import Path

BB_MIN_N = 6
BB_MAX_N = 2500
BB_BINS = 40
BB_MIN_GENES_FOR_FIT = 200


def read_tsv(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "gene_id" not in reader.fieldnames:
            raise ValueError(f"유효한 gene TSV가 아닙니다: {path}")
        return {row["gene_id"]: row for row in reader}


def binom_two_sided(a_count: int, b_count: int) -> float:
    """p=0.5인 exact two-sided binomial p-value(대칭 꼬리 2배)."""
    n = a_count + b_count
    if n == 0:
        return 1.0
    k = min(a_count, b_count)
    logs = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
            - n * math.log(2.0) for i in range(k + 1)]
    pivot = max(logs)
    cdf = math.exp(pivot) * sum(math.exp(x - pivot) for x in logs)
    return min(1.0, 2.0 * cdf)


def _lbeta(x: float, y: float) -> float:
    return math.lgamma(x) + math.lgamma(y) - math.lgamma(x + y)


def betabinom_pmf_series(n: int, p: float, rho: float, upto: int) -> list[float]:
    """k=0..upto의 beta-binomial PMF를 점화식으로 계산한다."""
    upto = min(upto, n)
    if rho <= 0.0:
        out = [math.exp(n * math.log1p(-p))] if p < 1.0 else [0.0]
        for k in range(upto):
            out.append(out[k] * (n - k) / (k + 1) * (p / (1.0 - p)))
        return out
    s = (1.0 - rho) / rho
    al, be = p * s, (1.0 - p) * s
    out = [math.exp(_lbeta(al, n + be) - _lbeta(al, be))]
    for k in range(upto):
        out.append(out[k] * (n - k) / (k + 1) * (al + k) / (be + n - k - 1))
    return out


def betabinom_two_sided(a_count: int, b_count: int, rho: float) -> float:
    """AlleleSeq2 형태의 대칭 beta-binomial p-value."""
    n = a_count + b_count
    if n == 0:
        return 1.0
    return min(1.0, 2.0 * math.fsum(
        betabinom_pmf_series(n, 0.5, rho, min(a_count, b_count))
    ))


def _binned_empirical(totals: list[int], ratios: list[float], nbins: int) -> list[float] | None:
    counts = [0.0] * nbins
    n_used = 0
    for ratio in ratios:
        idx = math.ceil(ratio * nbins) - 1
        idx = min(max(idx, 0), nbins - 1)
        counts[idx] += 1.0
        n_used += 1
    if n_used == 0:
        return None
    return [count / n_used for count in counts]


def _binned_null(depth_weights: dict[int, int], rho: float, nbins: int) -> list[float]:
    counts = [0.0] * nbins
    total_w = 0.0
    for n, weight in depth_weights.items():
        pmf = betabinom_pmf_series(n, 0.5, rho, n)
        for k, probability in enumerate(pmf):
            idx = math.ceil((k / n) * nbins) - 1
            idx = min(max(idx, 0), nbins - 1)
            counts[idx] += probability * weight
        total_w += weight
    if total_w <= 0:
        return counts
    total = math.fsum(counts)
    return [count / total for count in counts] if total > 0 else counts


def _graded_weights(nbins: int) -> list[float]:
    half = nbins // 2
    up = [(i + 1) / half for i in range(half)]
    return up + up[::-1]


def fit_overdispersion(
    rows: list[tuple[int, int]], *, include_low_depth_null: bool = False
) -> tuple[float, dict]:
    """AlleleSeq2식 가중 SSE 적합.

    현재 Stage 22/23은 경험·null 모두 total>=6만 사용한다. 구번호 Stage 18은
    과거 계산을 재현하도록 null depth 가중치에 total 1..5를 포함한다.
    """
    depth_weights: dict[int, int] = {}
    fit_ratios: list[float] = []
    fit_totals: list[int] = []
    max_total = 0
    for a_count, total in rows:
        if total <= 0:
            continue
        max_total = max(max_total, total)
        if total <= BB_MAX_N and (include_low_depth_null or total >= BB_MIN_N):
            depth_weights[total] = depth_weights.get(total, 0) + 1
        if BB_MIN_N <= total <= BB_MAX_N:
            fit_ratios.append(a_count / total)
            fit_totals.append(total)
    diag = {
        "method": "AlleleSeq2 calcOverdispersion.R 가중 SSE 적합",
        "min_total": BB_MIN_N, "max_total_used": min(max_total, BB_MAX_N),
        "bins": BB_BINS, "genes_used_for_fit": len(fit_ratios),
    }
    empirical = _binned_empirical(fit_totals, fit_ratios, BB_BINS)
    if empirical is None or len(fit_ratios) < BB_MIN_GENES_FOR_FIT:
        diag["status"] = "insufficient_data"
        diag["min_genes_required"] = BB_MIN_GENES_FOR_FIT
        return 0.0, diag
    weights = _graded_weights(BB_BINS)

    def sse(rho: float) -> float:
        null = _binned_null(depth_weights, rho, BB_BINS)
        return math.fsum(w * (e - t) ** 2 for w, e, t in zip(weights, empirical, null))

    best_rho, best_sse = 0.0, sse(0.0)
    for i in range(1, 10):
        rho = i * 0.1
        score = sse(rho)
        if score < best_sse:
            best_rho, best_sse = rho, score
        else:
            break
    lo, hi = max(0.0, best_rho - 0.1), min(0.99, best_rho + 0.1)
    rho = lo
    while rho <= hi + 1e-12:
        score = sse(rho)
        if score < best_sse:
            best_rho, best_sse = rho, score
        rho += 0.01
    diag["status"] = "ok"
    diag["sse"] = best_sse
    return best_rho, diag


def bh_fdr(p_values: list[float]) -> list[float]:
    order = sorted(range(len(p_values)), key=p_values.__getitem__)
    out = [1.0] * len(p_values)
    running = 1.0
    m = len(p_values)
    for rank0 in range(m - 1, -1, -1):
        idx = order[rank0]
        running = min(running, p_values[idx] * m / (rank0 + 1))
        out[idx] = min(1.0, running)
    return out


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total == 0:
        return 0.0, 1.0
    p = successes / total
    den = 1.0 + z * z / total
    center = (p + z * z / (2.0 * total)) / den
    half = z * math.sqrt(p * (1.0 - p) / total + z * z / (4.0 * total * total)) / den
    return max(0.0, center - half), min(1.0, center + half)


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2:
        return None
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    dx, dy = [x - mx for x in xs], [y - my for y in ys]
    den = math.sqrt(sum(x * x for x in dx) * sum(y * y for y in dy))
    return sum(x * y for x, y in zip(dx, dy)) / den if den else None


def flag(value: bool) -> str:
    return "1" if value else "0"
