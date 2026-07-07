"""Small statistical helpers for validation reports."""

from __future__ import annotations

import math
from dataclasses import dataclass


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

    numerator = 0.0
    denominator = 0.0
    observed_minus_expected = 0.0
    variance_sum = 0.0
    weighted_log_or = 0.0
    weight_sum = 0.0

    for result in results:
        a = result.benign_observed
        b = result.pathogenic_observed
        c = result.benign_not_observed
        d = result.pathogenic_not_observed
        n = a + b + c + d
        if n <= 1:
            continue

        numerator += (a * d) / n
        denominator += (b * c) / n

        alt_total = a + b
        no_alt_total = c + d
        benign_total = a + c
        pathogenic_total = b + d
        expected_a = alt_total * benign_total / n
        variance_a = alt_total * no_alt_total * benign_total * pathogenic_total / (n * n * (n - 1))
        observed_minus_expected += a - expected_a
        variance_sum += variance_a

        aa, bb, cc, dd = float(a), float(b), float(c), float(d)
        if min(aa, bb, cc, dd) == 0.0:
            aa += 0.5
            bb += 0.5
            cc += 0.5
            dd += 0.5
        var_log_or = (1.0 / aa) + (1.0 / bb) + (1.0 / cc) + (1.0 / dd)
        if var_log_or > 0:
            weight = 1.0 / var_log_or
            weighted_log_or += weight * math.log((aa * dd) / (bb * cc))
            weight_sum += weight

    if numerator == 0.0 and denominator == 0.0:
        odds_ratio_mh = float("nan")
    elif denominator == 0.0:
        odds_ratio_mh = float("inf")
    else:
        odds_ratio_mh = numerator / denominator

    if weight_sum > 0:
        pooled_log_or = weighted_log_or / weight_sum
        se = math.sqrt(1.0 / weight_sum)
        ci_low = math.exp(pooled_log_or - 1.96 * se)
        ci_high = math.exp(pooled_log_or + 1.96 * se)
    else:
        ci_low = float("nan")
        ci_high = float("nan")

    if variance_sum > 0:
        cmh_chi2 = (observed_minus_expected * observed_minus_expected) / variance_sum
        cmh_p = math.erfc(math.sqrt(cmh_chi2 / 2.0))
    else:
        cmh_chi2 = float("nan")
        cmh_p = float("nan")

    return AdjustedResult(odds_ratio_mh, ci_low, ci_high, cmh_chi2, cmh_p)


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
