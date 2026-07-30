"""Small statistical helpers for validation reports."""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

from statsmodels.stats.contingency_tables import StratifiedTable


@dataclass(frozen=True)
class EnrichmentResult:
    name: str
    benign_observed: int
    pathogenic_observed: int
    benign_not_observed: int
    pathogenic_not_observed: int
    odds_ratio: float
    ci_low: float
    ci_high: float
    fisher_p: float


@dataclass(frozen=True)
class AdjustedResult:
    odds_ratio_mh: float
    ci_low: float
    ci_high: float
    cmh_chi2: float
    cmh_p: float


def odds_ratio_and_ci(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    """Return raw OR and approximate 95% CI for [[a, b], [c, d]].

    Cells are:
    a = benign observed, b = pathogenic observed,
    c = benign not observed, d = pathogenic not observed.
    """

    denominator = b * c
    numerator = a * d
    if denominator == 0:
        raw_or = float("inf") if numerator > 0 else float("nan")
    else:
        raw_or = numerator / denominator

    aa, bb, cc, dd = float(a), float(b), float(c), float(d)
    if min(aa, bb, cc, dd) == 0.0:
        aa += 0.5
        bb += 0.5
        cc += 0.5
        dd += 0.5

    corrected_or = (aa * dd) / (bb * cc)
    se = math.sqrt((1.0 / aa) + (1.0 / bb) + (1.0 / cc) + (1.0 / dd))
    log_or = math.log(corrected_or)
    return raw_or, math.exp(log_or - 1.96 * se), math.exp(log_or + 1.96 * se)


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    try:
        from scipy.stats import fisher_exact

        return float(fisher_exact([[a, b], [c, d]], alternative="two-sided").pvalue)
    except Exception:
        return _fisher_exact_two_sided_fallback(a, b, c, d)


def enrichment_result(name: str, a: int, b: int, c: int, d: int) -> EnrichmentResult:
    odds_ratio, ci_low, ci_high = odds_ratio_and_ci(a, b, c, d)
    return EnrichmentResult(
        name=name,
        benign_observed=a,
        pathogenic_observed=b,
        benign_not_observed=c,
        pathogenic_not_observed=d,
        odds_ratio=odds_ratio,
        ci_low=ci_low,
        ci_high=ci_high,
        fisher_p=fisher_exact_two_sided(a, b, c, d),
    )


def mantel_haenszel_adjusted(results: list[EnrichmentResult]) -> AdjustedResult | None:
    """Return a pooled stratified OR and CMH p-value across 2x2 strata."""
    tables = [
        [
            [result.benign_observed, result.pathogenic_observed],
            [result.benign_not_observed, result.pathogenic_not_observed],
        ]
        for result in results
        if (
            result.benign_observed
            + result.pathogenic_observed
            + result.benign_not_observed
            + result.pathogenic_not_observed
        )
        > 1
    ]
    if not tables:
        return None

    table = StratifiedTable(tables, shift_zeros=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        odds_ratio_mh = float(table.oddsratio_pooled)
        if math.isfinite(odds_ratio_mh) and odds_ratio_mh > 0:
            ci_low, ci_high = map(float, table.oddsratio_pooled_confint(alpha=0.05))
        else:
            ci_low, ci_high = float("nan"), float("nan")
        cmh = table.test_null_odds(correction=False)
    return AdjustedResult(
        odds_ratio_mh,
        ci_low,
        ci_high,
        float(cmh.statistic),
        float(cmh.pvalue),
    )


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    """Return Benjamini-Hochberg adjusted p-values, preserving NaNs."""
    adjusted = [float("nan")] * len(pvalues)
    valid = [
        (index, float(value))
        for index, value in enumerate(pvalues)
        if value is not None and math.isfinite(float(value)) and 0 <= float(value) <= 1
    ]
    if not valid:
        return adjusted

    ordered = sorted(valid, key=lambda item: item[1])
    count = len(ordered)
    running_minimum = 1.0
    for rank_index in range(count - 1, -1, -1):
        original_index, pvalue = ordered[rank_index]
        rank = rank_index + 1
        running_minimum = min(running_minimum, pvalue * count / rank)
        adjusted[original_index] = min(1.0, running_minimum)
    return adjusted


def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _hypergeom_log_prob(k: int, row1: int, row2: int, col1: int, total: int) -> float:
    return _log_comb(row1, k) + _log_comb(row2, col1 - k) - _log_comb(total, col1)


def _fisher_exact_two_sided_fallback(a: int, b: int, c: int, d: int) -> float:
    row1 = a + b
    row2 = c + d
    col1 = a + c
    total = row1 + row2
    if total == 0:
        return float("nan")

    low = max(0, col1 - row2)
    high = min(row1, col1)
    observed_log = _hypergeom_log_prob(a, row1, row2, col1, total)
    included_logs = []
    for k in range(low, high + 1):
        log_prob = _hypergeom_log_prob(k, row1, row2, col1, total)
        if log_prob <= observed_log + 1e-12:
            included_logs.append(log_prob)
    if not included_logs:
        return 0.0
    max_log = max(included_logs)
    log_p = max_log + math.log(sum(math.exp(value - max_log) for value in included_logs))
    return min(1.0, math.exp(log_p))
