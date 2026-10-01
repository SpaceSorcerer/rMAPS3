"""Tests for tools/count_aware_stats.py: span->count reconstruction and the two tests.

Synthetic inputs only; no lab run directory is read.
"""

import importlib.util
import tempfile
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu, norm

MODULE_PATH = Path(__file__).resolve().parents[1] / "tools" / "count_aware_stats.py"
SPEC = importlib.util.spec_from_file_location("count_aware_stats", MODULE_PATH)
cas = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cas)

SEED = 20260920


def write_npz(path, n_events, window, step, length, spans, elig_lo, elig_hi, labels):
    event_index = np.array([s[0] for s in spans], dtype=np.int32)
    hit_start = np.array([s[1] for s in spans], dtype=np.int16)
    hit_end = np.array([s[2] for s in spans], dtype=np.int16)
    np.savez(path, schema_version=np.int16(2), window=np.int32(window),
             step=np.int32(step), region_length=np.int32(length),
             event_index=event_index, hit_start=hit_start, hit_end=hit_end,
             elig_lo=np.asarray(elig_lo, dtype=np.int32),
             elig_hi=np.asarray(elig_hi, dtype=np.int32),
             set_label=np.asarray(labels, dtype=np.uint8))


def brute_force_counts(n_events, window, step, length, spans, elig_lo, elig_hi):
    positions = np.arange(0, length - window + 1, step)
    counts = np.zeros((n_events, positions.size), dtype=np.int64)
    eligible = np.zeros_like(counts, dtype=bool)
    for event in range(n_events):
        for j, pos in enumerate(positions):
            eligible[event, j] = (elig_lo[event] >= 0 and pos >= elig_lo[event]
                                  and pos < elig_hi[event])
            if not eligible[event, j]:
                continue
            for ev, start, end in spans:
                if ev == event and start < pos + window and end > pos:
                    counts[event, j] += 1
    return counts, eligible


def histograms(fg_values, bg_values):
    kmax = int(max(fg_values.max(), bg_values.max()))
    counts = np.concatenate([fg_values, bg_values])[:, None].astype(np.int16)
    eligible = np.ones_like(counts, dtype=bool)
    fg_rows = np.arange(fg_values.size)
    bg_rows = np.arange(fg_values.size, counts.shape[0])
    return (cas.group_histograms(counts, eligible, fg_rows, kmax),
            cas.group_histograms(counts, eligible, bg_rows, kmax))


def test_counts_match_brute_force():
    rng = np.random.default_rng(SEED)
    n_events, window, step, length = 40, 10, 2, 60
    spans = []
    for event in range(n_events):
        for _ in range(int(rng.integers(0, 5))):
            start = int(rng.integers(0, length - 6))
            spans.append((event, start, start + int(rng.integers(3, 7))))
    elig_lo = rng.integers(0, 10, n_events).astype(np.int32)
    elig_hi = (length - window + 1) - rng.integers(0, 5, n_events).astype(np.int32)
    elig_lo[0] = -1
    labels = rng.integers(0, 3, n_events).astype(np.uint8)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "t.npz"
        write_npz(path, n_events, window, step, length, spans, elig_lo, elig_hi, labels)
        counts, eligible, got_labels = cas.load_counts(path)
    expected_counts, expected_elig = brute_force_counts(
        n_events, window, step, length, spans, elig_lo, elig_hi)
    np.testing.assert_array_equal(eligible, expected_elig)
    np.testing.assert_array_equal(counts, expected_counts * expected_elig)
    np.testing.assert_array_equal(got_labels, labels)
    assert expected_counts.max() > 1, "synthetic case must contain multi-hit exons"


def test_mw_matches_scipy_asymptotic():
    rng = np.random.default_rng(SEED)
    fg = rng.poisson(1.5, 120)
    bg = rng.poisson(1.0, 400)
    hist1, hist0 = histograms(fg, bg)
    got = float(cas.mw_counts_pvalues(hist1, hist0)[0])
    expected = float(mannwhitneyu(fg, bg, alternative="greater",
                                  method="asymptotic", use_continuity=True).pvalue)
    assert abs(got - expected) < 1e-12


def test_planted_density_enrichment_is_significant():
    rng = np.random.default_rng(SEED + 1)
    fg = rng.poisson(3.0, 200)
    bg = rng.poisson(1.0, 2000)
    hist1, hist0 = histograms(fg, bg)
    p_mw = float(cas.mw_counts_pvalues(hist1, hist0)[0])
    p_pois, ratio = cas.poisson_rate_pvalues(hist1, hist0)
    assert p_mw < 1e-10
    assert float(p_pois[0]) < 1e-10
    assert float(ratio[0]) > 2.0


def test_null_data_is_not_significant():
    rng = np.random.default_rng(SEED + 2)
    fg = rng.poisson(1.0, 200)
    bg = rng.poisson(1.0, 2000)
    hist1, hist0 = histograms(fg, bg)
    assert float(cas.mw_counts_pvalues(hist1, hist0)[0]) > 0.05
    assert float(cas.poisson_rate_pvalues(hist1, hist0)[0][0]) > 0.05


def test_binary_data_matches_proportion_test():
    """With 0/1 observations mw_counts reduces to a two-proportion test."""
    rng = np.random.default_rng(SEED + 3)
    fg = (rng.random(300) < 0.30).astype(np.int64)
    bg = (rng.random(1500) < 0.20).astype(np.int64)
    hist1, hist0 = histograms(fg, bg)
    p_mw = float(cas.mw_counts_pvalues(hist1, hist0)[0])
    n1, n0 = fg.size, bg.size
    p1, p0 = fg.mean(), bg.mean()
    pooled = (fg.sum() + bg.sum()) / (n1 + n0)
    z = (p1 - p0) / np.sqrt(pooled * (1 - pooled) * (1 / n1 + 1 / n0))
    p_prop = float(norm.sf(z))
    p_scipy_mw = float(mannwhitneyu(fg, bg, alternative="greater",
                                    method="asymptotic").pvalue)
    p_fisher = float(fisher_exact([[int(fg.sum()), int((1 - fg).sum())],
                                   [int(bg.sum()), int((1 - bg).sum())]],
                                  alternative="greater").pvalue)
    assert abs(p_mw - p_scipy_mw) < 1e-12
    assert abs(p_mw / p_prop - 1.0) < 0.02
    assert abs(np.log10(p_mw) - np.log10(p_fisher)) < 1.0


def test_poisson_ignores_ineligible_and_reports_rate_ratio():
    counts = np.array([[2], [0], [3], [5]], dtype=np.int16)
    eligible = np.array([[True], [True], [True], [False]])
    counts[~eligible] = 0
    hist_fg = cas.group_histograms(counts, eligible, np.array([0, 1]), 3)
    hist_bg = cas.group_histograms(counts, eligible, np.array([2, 3]), 3)
    p, ratio = cas.poisson_rate_pvalues(hist_fg, hist_bg)
    assert int(hist_bg.sum()) == 1
    assert abs(float(ratio[0]) - (2 / 2) / (3 / 1)) < 1e-12
    assert np.isfinite(p[0])


def test_parse_probe_rejects_malformed_specifications():
    assert cas.parse_probe("ARM:QKI.ACTAAC[ACG]:DownstreamIntron:dn") == (
        "ARM", "QKI.ACTAAC[ACG]", "DownstreamIntron", "dn")
    for bad in ("ARM:MOTIF:REGION", "ARM:MOTIF:REGION:sideways"):
        try:
            cas.parse_probe(bad)
        except ValueError:
            continue
        raise AssertionError("parse_probe accepted " + bad)
