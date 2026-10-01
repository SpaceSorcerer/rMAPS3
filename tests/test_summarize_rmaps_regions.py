import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control, fisher_exact


MODULE_PATH = Path(__file__).resolve().parents[1] / "tools" / "summarize_rmaps_regions.py"
SPEC = importlib.util.spec_from_file_location("summarize_rmaps_regions", MODULE_PATH)
summ = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(summ)


def test_planted_enrichment_and_direct_fisher_minimum():
    labels = np.array([1] * 20 + [0] * 20, dtype=bool)
    hits = np.zeros((40, 5), dtype=float)
    hits[:16, 0] = 1
    hits[20:22, 0] = 1
    rng = np.random.default_rng(149)
    hits[:, 1:] = rng.binomial(1, 0.2, size=(40, 4))
    perms = summ.generate_permutations(labels, 999, np.random.default_rng(149))

    result = summ.calibrate_min_p(hits, labels, perms, "greater")
    direct = min(
        fisher_exact(
            [[int(hits[:20, j].sum()), 20 - int(hits[:20, j].sum())],
             [int(hits[20:, j].sum()), 20 - int(hits[20:, j].sum())]],
            alternative="greater",
        ).pvalue
        for j in range(hits.shape[1])
    )
    np.testing.assert_allclose(result.observed_min_p, direct, rtol=1e-12, atol=0)
    assert result.calib_p < 0.02


def test_exchangeable_random_hits_are_uniformish():
    labels = np.array([1] * 30 + [0] * 30, dtype=bool)
    perm_rng = np.random.default_rng(149)
    perms = summ.generate_permutations(labels, 299, perm_rng)
    data_rng = np.random.default_rng(811)
    calibrated = []
    for _ in range(60):
        hits = data_rng.binomial(1, 0.18, size=(60, 12)).astype(float)
        calibrated.append(summ.calibrate_min_p(hits, labels, perms, "greater").calib_p)
    assert 0.4 <= float(np.mean(calibrated)) <= 0.6


def test_na_cells_never_enter_counts():
    hits = np.array([[1, np.nan], [np.nan, 1], [1, 0], [0, np.nan]], dtype=float)
    labels = np.array([1, 1, 0, 0], dtype=bool)
    counts = summ.window_counts(hits, labels)
    assert counts["fg_hit"].tolist() == [1, 1]
    assert counts["fg_elig"].tolist() == [1, 1]
    assert counts["bg_hit"].tolist() == [1, 0]
    assert counts["bg_elig"].tolist() == [2, 1]


def test_all_na_permutations_still_count_in_requested_B_denominator():
    hits = np.array([[1, np.nan], [0, np.nan], [np.nan, 1], [np.nan, 0]], dtype=float)
    labels = np.array([1, 0, 1, 0], dtype=bool)
    permutations = np.array(
        [[1, 1, 0, 0], [1, 0, 1, 0], [1, 0, 0, 1],
         [0, 1, 1, 0], [0, 1, 0, 1], [0, 0, 1, 1]],
        dtype=np.uint8,
    )
    result = summ.calibrate_min_p(hits, labels, permutations, "greater")
    assert result.observed_min_p == 0.5
    assert np.isnan(result.perm_min_p[[0, 5]]).all()
    assert result.calib_p == 4 / 7


def test_bh_matches_scipy():
    p = np.array([0.001, 0.02, 0.2, 0.9, 0.041, 0.06])
    np.testing.assert_allclose(summ.bh_adjust(p), false_discovery_control(p, method="bh"))


def test_region_pooling_covers_all_eight_exactly_once():
    members = [region for regions in summ.POOL_TO_REGIONS.values() for region in regions]
    assert len(members) == 8
    assert len(set(members)) == 8
    assert set(members) == set(summ.REGIONS)


def test_excel_preserves_infinite_enrichment_ratio_as_text():
    assert summ._excel_value(math.inf) == "Inf"
    assert summ._excel_value(-math.inf) == "-Inf"
    assert summ._excel_value(math.nan) is None


def _slow_reference(hits, eligible, labels, permutations):
    """Scalar Fisher/min-P oracle, deliberately independent of the BLAS path."""
    assignments = np.vstack([labels.astype(np.uint8), permutations.astype(np.uint8)])
    minima = []
    for assignment in assignments:
        per_window = []
        for column in range(hits.shape[1]):
            fg = assignment.astype(bool)
            bg = ~fg
            fg_elig = int(eligible[fg, column].sum())
            bg_elig = int(eligible[bg, column].sum())
            if not fg_elig or not bg_elig:
                continue
            fg_hit = int(hits[fg, column].sum())
            bg_hit = int(hits[bg, column].sum())
            per_window.append(fisher_exact(
                [[fg_hit, fg_elig - fg_hit], [bg_hit, bg_elig - bg_hit]],
                alternative="greater",
            ).pvalue)
        minima.append(min(per_window) if per_window else np.nan)
    return np.asarray(minima)


def test_vectorized_matches_slow_reference_with_ineligible_cells():
    rng = np.random.default_rng(8021)
    labels = np.array([1] * 5 + [0] * 6, dtype=bool)
    eligible = rng.random((11, 7)) >= 0.25
    hits = (rng.random((11, 7)) < 0.3) & eligible
    permutations = summ.generate_permutations(labels, 31, np.random.default_rng(149))

    result = summ.calibrate_min_p(
        hits, labels, permutations, "greater", chunk_size=8, eligible=eligible
    )
    reference = _slow_reference(hits, eligible, labels, permutations)
    observed_per_window = []
    for column in range(hits.shape[1]):
        fg_elig = int(eligible[labels, column].sum())
        bg_elig = int(eligible[~labels, column].sum())
        if not fg_elig or not bg_elig:
            observed_per_window.append(np.nan)
            continue
        fg_hit = int(hits[labels, column].sum())
        bg_hit = int(hits[~labels, column].sum())
        observed_per_window.append(fisher_exact(
            [[fg_hit, fg_elig - fg_hit], [bg_hit, bg_elig - bg_hit]],
            alternative="greater",
        ).pvalue)

    np.testing.assert_allclose(result.observed_min_p, reference[0], rtol=0, atol=1e-12)
    np.testing.assert_allclose(result.observed_p, observed_per_window, rtol=0, atol=1e-12,
                               equal_nan=True)
    np.testing.assert_allclose(result.perm_min_p, reference[1:], rtol=0, atol=1e-12,
                               equal_nan=True)
    expected_calib = (1 + int(np.sum(reference[1:] <= reference[0]))) / (1 + len(permutations))
    assert result.calib_p == expected_calib


def test_chunked_and_unchunked_permutations_are_identical():
    rng = np.random.default_rng(990)
    labels = np.array([1] * 8 + [0] * 9, dtype=bool)
    eligible = rng.random((17, 13)) >= 0.1
    hits = (rng.random((17, 13)) < 0.2) & eligible
    permutations = summ.generate_permutations(labels, 47, np.random.default_rng(149))
    np.testing.assert_array_equal(
        permutations,
        summ.generate_permutations(labels, 47, np.random.default_rng(149)),
    )

    chunked = summ.calibrate_min_p(
        hits, labels, permutations, "greater", chunk_size=7, eligible=eligible
    )
    unchunked = summ.calibrate_min_p(
        hits, labels, permutations, "greater", chunk_size=47, eligible=eligible
    )

    assert chunked.observed_min_p == unchunked.observed_min_p
    assert chunked.calib_p == unchunked.calib_p
    np.testing.assert_array_equal(chunked.perm_min_p, unchunked.perm_min_p)


def test_manifest_schema_2_is_required(tmp_path):
    manifest = tmp_path / "run_manifest.json"
    manifest.write_text(json.dumps({"schema_version": 1}), encoding="utf-8")
    with np.testing.assert_raises_regex(ValueError, "schema_version 2"):
        summ.read_manifest(manifest)


def random_cells():
    rng = np.random.default_rng(1492026)
    size = 20000
    N = rng.integers(2, 90001, size)
    n = rng.integers(1, np.minimum(N - 1, 1000) + 1)
    K = rng.integers(0, np.minimum(N, 20000) + 1)
    lo = np.maximum(0, n + K - N)
    hi = np.minimum(K, n)
    k = rng.integers(lo, hi + 1)
    # Include near-null, tiny tails, degenerate support, and exact endpoints.
    k[::5] = rng.hypergeometric(K[::5], N[::5] - K[::5], n[::5])
    k[1::11] = lo[1::11]
    k[2::11] = hi[2::11]
    n0 = N - n
    c = K - k
    n[3::101] = 0
    n0[4::101] = 0
    k[5::101] = -1
    c[6::101] = -1
    k[7::101] = n[7::101] + 1
    c[8::101] = n0[8::101] + 1
    return k, n, c, n0


def test_20000_random_logspace_equivalence():
    cells = random_cells()
    actual = summ.greater_pvalues(*cells)
    expected = summ.fisher_pvalues(*cells, "greater")
    np.testing.assert_array_equal(np.isnan(actual), np.isnan(expected))
    normal = (actual >= np.finfo(float).tiny) & (expected >= np.finfo(float).tiny)
    errors = np.abs(np.log(actual[normal]) - np.log(expected[normal]))
    assert errors.max() <= 1e-9
    # SciPy can round near-one tails and subnormals differently. Resolve those
    # disagreements against integer combinatorics, not approximate SciPy sf.
    endpoints = (actual == 0) | (actual == 1) | (expected == 0) | (expected == 1)
    subnormal = ((actual > 0) & (actual < np.finfo(float).tiny)) | (
        (expected > 0) & (expected < np.finfo(float).tiny))
    disputed = (endpoints | subnormal) & (actual != expected) & np.isfinite(actual)
    for i in np.flatnonzero(disputed):
        a, n1, c, n0 = [int(x[i]) for x in cells]
        total, successes, draws = n1 + n0, a + c, n1
        numerator = sum(math.comb(successes, j) * math.comb(total - successes, draws - j)
                        for j in range(a, min(successes, draws) + 1))
        exact_rounded = numerator / math.comb(total, draws)
        if endpoints[i]:
            assert actual[i] == exact_rounded
        else:
            # Relative error is ill-conditioned at underflow: one representable
            # float step is the smallest meaningful absolute error allowance.
            tolerance = max(exact_rounded * 1e-9, np.nextafter(0.0, 1.0))
            assert abs(actual[i] - exact_rounded) <= tolerance


def test_hot_path_does_not_call_scipy_stats(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("scipy.stats reference called in hot path")
    monkeypatch.setattr(summ, "fisher_pvalues", forbidden)
    monkeypatch.setattr(summ.hypergeom, "sf", forbidden)
    labels = np.array([1, 1, 0, 0], dtype=bool)
    hits = np.array([[1, 0], [1, 1], [0, 0], [0, 1]], dtype=bool)
    perms = summ.generate_permutations(labels, 20, np.random.default_rng(149))
    result = summ.calibrate_min_p(hits, labels, perms, "greater")
    assert np.isfinite(result.calib_p)


def test_single_hit_tail_preserves_exact_bernoulli_identity():
    N = np.array([374, 398, 1000, 90000])
    n = np.array([187, 199, 127, 917])
    actual = summ.greater_pvalues(np.ones(4, dtype=int), n, np.zeros(4, dtype=int), N-n)
    np.testing.assert_array_equal(actual, n / N)




def test_non_fisher_manifest_is_rejected_before_summary(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    (run / "positional").mkdir()
    (run / "temp").mkdir()
    (run / "pVal.up.vs.bg.RNAmap.txt").touch()
    (run / "pVal.dn.vs.bg.RNAmap.txt").touch()
    (run / "run_manifest.json").write_text(json.dumps({
        "schema_version": 2,
        "statistical_method": "zscore",
        "parameters": {"fisher_alternative": "greater", "genome": "synthetic"},
    }), encoding="utf-8")
    out = tmp_path / "summary"
    with np.testing.assert_raises_regex(ValueError, "requires statistical_method='fisher'"):
        summ.summarize(run, out, "synthetic", 20, 149, 1)
    assert not out.exists()


def test_unsafe_arm_is_rejected_without_writes(tmp_path):
    run = tmp_path / "missing_run"
    out = tmp_path / "summary"
    for arm in ("../outside", "nested/output", "nested\\output", "", ".", "..", "/absolute", "arm:stream"):
        with np.testing.assert_raises_regex(ValueError, "arm must be a filename component"):
            summ.summarize(run, out, arm, 20, 149, 1)
        assert not out.exists()
        assert list(tmp_path.iterdir()) == []
