from __future__ import annotations

import os
import warnings
from functools import lru_cache
from math import ceil
from typing import Sequence

import numpy as np
from scipy import stats


METHOD_ALIASES = {
    "fisher": "fisher",
    "mannwhitney": "mannwhitney_greater",
    "mannwhitney_greater": "mannwhitney_greater",
    "brunnermunzel": "brunnermunzel_greater",
    "brunnermunzel_greater": "brunnermunzel_greater",
    "permutation": "permutation_one_sided",
    "permutation_one_sided": "permutation_one_sided",
}

METHOD_HEADERS = {
    "fisher": "fisher.exact.pVal",
    "mannwhitney_greater": "mannwhitney.greater.pVal",
    "brunnermunzel_greater": "brunnermunzel.greater.pVal",
    "permutation_one_sided": "permutation.oneSided.pVal",
}


class StatisticUnavailable(ValueError):
    """The test is undefined for these observations; callers report NA with this reason, never p=1."""


def supported_stat_methods() -> tuple[str, ...]:
    return (
        "fisher",
        "mannwhitney_greater",
        "brunnermunzel_greater",
        "permutation_one_sided",
    )


def normalize_stat_method(stat_method: str | None) -> str:
    method = (stat_method or "fisher").strip().lower()
    return METHOD_ALIASES.get(method, "fisher")


def pvalue_header_label(stat_method: str | None) -> str:
    method = normalize_stat_method(stat_method)
    return METHOD_HEADERS.get(method, METHOD_HEADERS["fisher"])


def _fisher_greater_pvalue(values_one: Sequence[float], values_two: Sequence[float], fisher_scale: float) -> float:
    scale = fisher_scale if fisher_scale and fisher_scale > 0 else 1.0

    succ1 = int(ceil(sum(values_one) / scale))
    fail1 = int(ceil(max(0.0, len(values_one) - sum(values_one)) / scale))
    succ2 = int(ceil(sum(values_two) / scale))
    fail2 = int(ceil(max(0.0, len(values_two) - sum(values_two)) / scale))

    return float(stats.fisher_exact([[succ1, fail1], [succ2, fail2]], "greater")[1])


def _mannwhitney_greater_pvalue(values_one: Sequence[float], values_two: Sequence[float]) -> float:
    return float(stats.mannwhitneyu(values_one, values_two, alternative="greater")[1])


def _brunnermunzel_greater_pvalue(values_one: Sequence[float], values_two: Sequence[float]) -> float:
    # AUDIT F14: undefined rank-test variance is an error, never evidence of p=1.
    if len(values_one) < 2 or len(values_two) < 2:
        raise StatisticUnavailable("Brunner-Munzel requires at least two observations per group")

    if len(set(values_one)) < 2 or len(set(values_two)) < 2:
        raise StatisticUnavailable("Brunner-Munzel requires nonconstant observations in each group")

    pooled = list(values_one) + list(values_two)
    if len(set(pooled)) < 2:
        raise StatisticUnavailable("Brunner-Munzel pooled observations are constant")

    # Reject effectively undefined variance instead of manufacturing a p-value.
    arr1 = np.asarray(values_one, dtype=float)
    arr2 = np.asarray(values_two, dtype=float)
    if np.std(arr1) < 1e-10 or np.std(arr2) < 1e-10:
        raise StatisticUnavailable("Brunner-Munzel group variance is too small")

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        # Use distribution='normal' for robustness; 't' can fail with 0/0 degrees of freedom
        pvalue = stats.brunnermunzel(values_one, values_two, alternative="greater", distribution="normal")[1]
    return float(pvalue)


def _sanitize_values(values: Sequence[float]) -> list[float]:
    # AUDIT F14: callers must remove ineligible observations explicitly.
    if any(not np.isfinite(v) for v in values):
        raise ValueError("Statistical input contains nonfinite values")
    return [float(v) for v in values]


def _int_from_env(name: str, default: int, minimum: int) -> int:
    raw = os.environ.get(name)
    if not raw:
        return default
    try:
        parsed = int(raw)
    except ValueError:
        return default
    return parsed if parsed >= minimum else default


@lru_cache(maxsize=1)
def _permutation_config() -> tuple[int, int]:
    n_permutations = _int_from_env("RMAPS_STAT_PERMUTATIONS", default=500, minimum=50)
    seed = _int_from_env("RMAPS_STAT_SEED", default=1337, minimum=0)
    return n_permutations, seed


def _permutation_one_sided_pvalue(values_one: Sequence[float], values_two: Sequence[float]) -> float:
    n_permutations, seed = _permutation_config()

    first = np.asarray(values_one, dtype=float)
    second = np.asarray(values_two, dtype=float)
    n_first = first.size

    observed = float(np.mean(first) - np.mean(second))
    pooled = np.concatenate([first, second])
    rng = np.random.default_rng(seed)

    at_least_observed = 0
    for _ in range(n_permutations):
        permuted = rng.permutation(pooled)
        perm_first = permuted[:n_first]
        perm_second = permuted[n_first:]
        stat = float(np.mean(perm_first) - np.mean(perm_second))
        if stat >= observed:
            at_least_observed += 1

    return (at_least_observed + 1.0) / (n_permutations + 1.0)


def _clamp_pvalue(pvalue: float) -> float:
    # AUDIT F14: invalid statistical results propagate with a diagnostic.
    if not np.isfinite(pvalue) or not 0.0 <= pvalue <= 1.0:
        raise ValueError(f"Invalid statistical p-value: {pvalue!r}")
    return float(pvalue)


def compute_locus_pvalue(
    values_one: Sequence[float],
    values_two: Sequence[float],
    stat_method: str | None,
    fisher_scale: float = 1.0,
) -> float:
    method = normalize_stat_method(stat_method)

    clean_one = _sanitize_values(values_one)
    clean_two = _sanitize_values(values_two)

    if not clean_one or not clean_two:
        raise StatisticUnavailable("Statistical comparison has an empty eligible group")

    try:
        if method == "mannwhitney_greater":
            return _clamp_pvalue(_mannwhitney_greater_pvalue(clean_one, clean_two))
        if method == "brunnermunzel_greater":
            return _clamp_pvalue(_brunnermunzel_greater_pvalue(clean_one, clean_two))
        if method == "permutation_one_sided":
            return _clamp_pvalue(_permutation_one_sided_pvalue(clean_one, clean_two))
        return _clamp_pvalue(_fisher_greater_pvalue(clean_one, clean_two, fisher_scale))
    except StatisticUnavailable as exc:
        raise StatisticUnavailable(f"{method}: {exc}") from exc
    except Exception as exc:
        raise ValueError(f"{method} statistical comparison failed: {exc}") from exc


def locus_pvalue_or_na(
    values_one: Sequence[float],
    values_two: Sequence[float],
    stat_method: str | None,
    fisher_scale: float = 1.0,
) -> tuple[float, str]:
    """compute_locus_pvalue with an undefined test returned as (NaN, reason); every other failure still raises."""
    try:
        return compute_locus_pvalue(values_one, values_two, stat_method, fisher_scale=fisher_scale), ""
    except StatisticUnavailable as exc:
        return float("nan"), str(exc)


def format_pvalue(pvalue: float) -> str:
    return str(pvalue) if np.isfinite(pvalue) else "NA"
