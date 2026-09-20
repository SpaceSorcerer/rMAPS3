"""Tests for tools/calibrate_ranksum.py and tools/rmaps_countdist_io.py.

Synthetic inputs only. The three checks of the original lab suite that read a
released engine run directory (released root table, countDist text round trip
against the real temp/ files, per-position agreement with a released pVal file)
are not portable and are not reproduced here; the text round trip is covered
against a synthetic countDist table written in the engine's own format.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
from scipy.stats import mannwhitneyu

TOOLS = Path(__file__).resolve().parents[1] / "tools"


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, TOOLS / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


calib = _load("calibrate_ranksum", "calibrate_ranksum.py")
io = _load("rmaps_countdist_io", "rmaps_countdist_io.py")


def synthetic_npz(tmp_path, fg_counts, bg_counts, n_windows, name="TEST.AAA"):
    """An archive with the real position-axis shape but synthetic counts."""
    per_region = n_windows // len(io.REGIONS)
    region_index = np.repeat(np.arange(len(io.REGIONS), dtype=np.int8), per_region)
    position = np.tile(np.arange(per_region, dtype=np.int32), len(io.REGIONS))
    store = {
        "schema_version": np.asarray(io.SCHEMA_VERSION),
        "arm": np.asarray("SYNTH"),
        "motif": np.asarray(name),
        "region_name": np.asarray(io.REGIONS, dtype="U"),
        "region_of_position": region_index,
        "position": position,
    }
    for group, matrix in (("up", fg_counts), ("dn", fg_counts), ("bg", bg_counts)):
        store[group + "_exon_id"] = np.asarray(
            ["exon%d" % i for i in range(matrix.shape[0])], dtype="U")
        io.pack_group(store, group, matrix.astype(np.int16))
    path = tmp_path / (name + ".counts.npz")
    np.savez_compressed(path, **store)
    return path


def write_countdist(path, matrix, positions_per_region):
    """A countDist text file in the exact format printCountDist() emits."""
    with open(path, "w", newline="\n", encoding="utf-8") as handle:
        handle.write("Region\tposition\tsum\tvalues\n")
        column = 0
        for region in io.REGIONS:
            for locus in range(positions_per_region):
                values = matrix[:, column]
                handle.write("{}\t{}\t{}\t[{}]\n".format(
                    region, locus, int(values.sum()),
                    ",".join(str(int(v)) for v in values)))
                column += 1


def test_countdist_text_round_trip_is_self_contained(tmp_path):
    """parse -> pack -> unpack reproduces the text table, with no released input."""
    rng = np.random.default_rng(149)
    positions_per_region = 3
    n_windows = positions_per_region * len(io.REGIONS)
    matrix = rng.poisson(0.3, size=(17, n_windows)).astype(np.int16)
    source = tmp_path / "T.AAA.countDist.bg.txt"
    write_countdist(source, matrix, positions_per_region)

    region_index, position, parsed = io.parse_countdist(source)
    assert np.array_equal(parsed, matrix)
    assert np.array_equal(region_index,
                          np.repeat(np.arange(len(io.REGIONS), dtype=np.int8),
                                    positions_per_region))
    assert np.array_equal(position, np.tile(np.arange(positions_per_region, dtype=np.int32),
                                            len(io.REGIONS)))
    store = {}
    io.pack_group(store, "bg", parsed)
    archive = tmp_path / "T.AAA.counts.npz"
    np.savez_compressed(archive, **store)
    with np.load(archive, allow_pickle=False) as data:
        assert np.array_equal(io.group_csr(data, "bg").toarray(), matrix)


def test_rank_sum_matches_scipy_with_continuity_correction():
    """The kernel equals scipy's asymptotic MWU WITH the continuity correction.

    The released engine calls stats.mannwhitneyu(first, second,
    alternative='greater') with no further keywords, so scipy defaults apply:
    use_continuity=True, and method='auto' resolves to 'asymptotic' whenever
    ties are present, which they always are here.
    """
    rng = np.random.default_rng(149)
    fg = rng.poisson(0.6, size=(11, 5)).astype(np.int16)
    bg = rng.poisson(0.2, size=(37, 5)).astype(np.int16)
    kmax = int(max(fg.max(), bg.max()))
    _, _, _, z, _ = io.rank_statistics(io.count_histograms(fg, kmax),
                                       io.count_histograms(bg, kmax))
    mine = io.p_from_z(z)
    with_continuity = np.array([
        mannwhitneyu(fg[:, w], bg[:, w], alternative="greater", method="asymptotic",
                     use_continuity=True).pvalue for w in range(fg.shape[1])])
    without_continuity = np.array([
        mannwhitneyu(fg[:, w], bg[:, w], alternative="greater", method="asymptotic",
                     use_continuity=False).pvalue for w in range(fg.shape[1])])
    assert np.allclose(mine, with_continuity, rtol=0, atol=0)
    assert not np.allclose(mine, without_continuity, rtol=0, atol=0)


def calibrated_p(tmp_path, fg, bg, n_windows, name, permutations):
    path = synthetic_npz(tmp_path, fg, bg, n_windows, name=name)
    model = calib.MotifModel(path, "up")
    selection = calib.draw_selection(
        np.random.default_rng(np.random.SeedSequence([149, 1, 0])),
        model.n_total, model.n1, permutations)
    return calib.calibrate(model, model.window_masks(), selection, 100)


def test_planted_enrichment_hits_the_floor_and_a_null_motif_does_not(tmp_path):
    """Calibration separates a planted signal from exchangeable null motifs.

    The planted signal sits in windows 20-25, which fall in UpstreamIntron and
    therefore in the pooled Upstream Intron region. Null replicates are built by
    drawing ONE pooled count matrix and splitting it, so the changed/background
    labels are exchangeable by construction and the calibrated p must be roughly
    uniform rather than pinned to the permutation floor.
    """
    rng = np.random.default_rng(149)
    n_windows, n_fg, n_bg, permutations = 80, 40, 800, 500
    floor = 1.0 / (1.0 + permutations)
    pool = "Upstream Intron"

    fg = np.zeros((n_fg, n_windows), dtype=np.int16)
    bg = np.zeros((n_bg, n_windows), dtype=np.int16)
    bg[rng.random((n_bg, n_windows)) < 0.02] = 1
    fg[rng.random((n_fg, n_windows)) < 0.02] = 1
    fg[:, 20:26] = rng.integers(1, 4, size=(n_fg, 6)).astype(np.int16)
    planted = calibrated_p(tmp_path, fg, bg, n_windows, "PLANT.AAA", permutations)
    assert planted[pool]["calibrated_p"] == pytest.approx(floor), (
        "planted calibrated p %g" % planted[pool]["calibrated_p"])

    null_values = []
    for replicate in range(5):
        pooled = np.zeros((n_fg + n_bg, n_windows), dtype=np.int16)
        pooled[rng.random(pooled.shape) < 0.02] = 1
        null_values.append(calibrated_p(tmp_path, pooled[:n_fg], pooled[n_fg:], n_windows,
                                        "NULL%d.AAA" % replicate,
                                        permutations)[pool]["calibrated_p"])
    assert sum(v == pytest.approx(floor) for v in null_values) <= 1, str(null_values)
    assert float(np.median(null_values)) > 10 * floor, str(null_values)


def test_permutation_preserves_the_changed_set_size():
    """Every drawn assignment has exactly n1 distinct changed exons."""
    rng = np.random.default_rng(np.random.SeedSequence([149, 1, 0]))
    n_total, n1, batch = 1000, 37, 200
    selection = calib.draw_selection(rng, n_total, n1, batch)
    assert selection.shape == (batch, n1)
    assert selection.min() >= 0 and selection.max() < n_total
    for row in selection:
        assert np.unique(row).size == n1


def test_permutation_z_matches_a_direct_recomputation(tmp_path):
    """The sparse selection shortcut equals a full rank test on relabelled data."""
    rng = np.random.default_rng(149)
    n_windows, n_fg, n_bg = 24, 15, 120
    fg = rng.poisson(0.4, size=(n_fg, n_windows)).astype(np.int16)
    bg = rng.poisson(0.1, size=(n_bg, n_windows)).astype(np.int16)
    path = synthetic_npz(tmp_path, fg, bg, n_windows, name="CHK.AAA")
    model = calib.MotifModel(path, "up")
    pooled = np.vstack([fg, bg])
    selection = calib.draw_selection(np.random.default_rng(7), model.n_total, model.n1, 5)
    fast = model.permutation_z(selection)
    for i, indices in enumerate(selection):
        mask = np.zeros(model.n_total, dtype=bool)
        mask[indices] = True
        kmax = int(pooled.max())
        _, _, _, direct, _ = io.rank_statistics(
            io.count_histograms(pooled[mask], kmax), io.count_histograms(pooled[~mask], kmax))
        assert np.allclose(fast[i], direct, rtol=1e-10, atol=1e-10, equal_nan=True)


def test_sparse_and_dense_packing_agree(tmp_path):
    """Both storage branches of pack_group reconstruct the same matrix."""
    rng = np.random.default_rng(3)
    dense_matrix = rng.integers(0, 4, size=(40, 30)).astype(np.int16)
    sparse_matrix = np.zeros((40, 30), dtype=np.int16)
    sparse_matrix[rng.random((40, 30)) < 0.02] = 2
    for matrix, expected in ((dense_matrix, "dense"), (sparse_matrix, "csr")):
        store = {}
        io.pack_group(store, "up", matrix)
        assert str(store["up_format"]) == expected
        path = tmp_path / ("pack_" + expected + ".npz")
        np.savez_compressed(path, **store)
        with np.load(path, allow_pickle=False) as data:
            assert np.array_equal(io.group_csr(data, "up").toarray(), matrix)
            assert np.array_equal(io.unpack_group(data, "up"), matrix)


def test_csr_histograms_equal_dense_histograms():
    rng = np.random.default_rng(11)
    matrix = rng.poisson(0.3, size=(200, 17)).astype(np.int16)
    kmax = int(matrix.max())
    csr = sparse.csr_matrix(matrix)
    assert np.array_equal(io.csr_histograms(csr, kmax), io.count_histograms(matrix, kmax))
