"""Synthetic tests for tools/calibrate_ranksum_v2.py and tools/rmaps_calib_v2_lib.py.

Ported from the 2026-09-22 lab calibration tests; the real-output integration tests stay with the
lab run because they read lab result trees.
"""
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import kstest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import calibrate_ranksum as calib  # noqa: E402
import calibrate_ranksum_v2 as v2  # noqa: E402
import rmaps_calib_v2_lib as lib  # noqa: E402
import rmaps_countdist_io as io  # noqa: E402

WINDOWS_PER_REGION = 5
SEED = 149


def make_dataset(tmp_path, rng, fg_sizes, bg_sizes, fg_rate, bg_rate, name="m"):
    """Clustered synthetic archive: rows of one target exon share identical counts."""
    n_windows = WINDOWS_PER_REGION * len(io.REGIONS)
    region = np.repeat(np.arange(len(io.REGIONS), dtype=np.int8), WINDOWS_PER_REGION)
    position = np.tile(np.arange(WINDOWS_PER_REGION, dtype=np.int32), len(io.REGIONS))

    def block(sizes, rate, offset):
        rows, ids = [], []
        for c, size in enumerate(sizes):
            counts = rng.poisson(rate, n_windows).astype(np.int16)
            for r in range(size):
                rows.append(counts)
                ids.append("chr1:+:{}:{}:{}:1:2:3".format(offset + c, offset + c + 50, r))
        return np.vstack(rows), np.asarray(ids)

    fg, fg_ids = block(fg_sizes, fg_rate, 0)
    bg, bg_ids = block(bg_sizes, bg_rate, 10 ** 6)
    store = {"region_of_position": region, "position": position}
    for group, matrix, ids in (("up", fg, fg_ids), ("dn", fg, fg_ids), ("bg", bg, bg_ids)):
        io.pack_group(store, group, matrix)
        store[group + "_exon_id"] = ids
    path = tmp_path / (name + ".counts.npz")
    np.savez(path, **store)
    return path, lib.target_keys(fg_ids), lib.target_keys(bg_ids)


def planted_rate(n_windows):
    rate = np.full(n_windows, 0.2)
    rate[2 * WINDOWS_PER_REGION: 3 * WINDOWS_PER_REGION] = 2.0
    return rate


def test_cluster_integrity(tmp_path):
    rng = np.random.default_rng(1)
    _, fg_keys, bg_keys = make_dataset(tmp_path, rng, [1] * 20 + [2] * 5 + [3] * 2,
                                       [1] * 300 + [2] * 60 + [3] * 20 + [4] * 5, 0.3, 0.3)
    drawer = lib.ClusterDrawer(fg_keys, bg_keys)
    pooled = np.concatenate([fg_keys, bg_keys])
    members = lib.cluster_rows(pooled)
    draws = drawer.draw(np.random.default_rng(7), 300)
    need = dict(drawer.need)
    for draw in draws:
        assert draw.shape[0] == fg_keys.shape[0] == drawer.n1
        assert np.unique(draw).shape[0] == draw.shape[0]
        chosen = set(draw.tolist())
        sizes = {}
        for key in set(pooled[draw].tolist()):
            rows = set(members[key].tolist())
            assert rows <= chosen, "a target exon was split"
            sizes[len(rows)] = sizes.get(len(rows), 0) + 1
        assert sizes == need


def test_stage2_precomputed_draws_equal_chunked_regeneration(tmp_path):
    rng = np.random.default_rng(2)
    _, fg_keys, bg_keys = make_dataset(tmp_path, rng, [1] * 10 + [2] * 3, [1] * 200 + [2] * 40, 0.3, 0.3)
    drawer = lib.ClusterDrawer(fg_keys, bg_keys)
    whole = drawer.draw(lib.stage_rng(SEED, 2, "up"), 1500)
    stream = lib.stage_rng(SEED, 2, "up")
    chunked = np.vstack([drawer.draw(stream, 500) for _ in range(3)])
    assert np.array_equal(whole, chunked)


def test_stage_streams_are_independent_per_stage_and_direction():
    a = lib.stage_rng(SEED, 1, "up").integers(0, 10 ** 9, 5)
    assert not np.array_equal(a, lib.stage_rng(SEED, 2, "up").integers(0, 10 ** 9, 5))
    assert not np.array_equal(a, lib.stage_rng(SEED, 1, "dn").integers(0, 10 ** 9, 5))
    assert np.array_equal(a, lib.stage_rng(SEED, 1, "up").integers(0, 10 ** 9, 5))


def test_drawer_refuses_a_target_exon_in_foreground_and_background():
    keys = np.asarray(["chr1:+:1:2", "chr1:+:3:4"])
    try:
        lib.ClusterDrawer(keys, np.asarray(["chr1:+:3:4", "chr1:+:5:6", "chr1:+:7:8"]))
    except ValueError as exc:
        assert "both foreground and background" in str(exc)
    else:
        raise AssertionError("an exon shared by foreground and background must be refused")


def test_calibrate_keep_matches_v1_calibrate(tmp_path):
    rng = np.random.default_rng(3)
    path, fg_keys, bg_keys = make_dataset(tmp_path, rng, [1] * 15 + [2] * 4, [1] * 200 + [2] * 40, 0.6, 0.3)
    model = calib.MotifModel(path, "up")
    masks = model.window_masks()
    selection = lib.ClusterDrawer(fg_keys, bg_keys).draw(np.random.default_rng(5), 400)
    kept, maxima = lib.calibrate_keep(model, masks, selection, 128)
    reference = calib.calibrate(model, masks, selection, 100)
    for name in masks:
        assert maxima[name].shape == (400,)
        a, b = kept[name]["calibrated_p"], reference[name]["calibrated_p"]
        assert (math.isnan(a) and math.isnan(b)) or a == b


def test_rbp_single_motif_maxz_equals_motif_statistic():
    rng = np.random.default_rng(4)
    null = rng.normal(size=(1, 999))
    o_max, n_max, o_mean, n_mean, arg = lib.rbp_combine(np.asarray([1.7]), null)
    assert o_max == 1.7 and arg == 0
    assert np.array_equal(n_max, null[0])
    assert lib.permutation_p(o_max, n_max) == lib.permutation_p(1.7, null[0])


def test_rbp_single_motif_meanz_equals_maxz():
    rng = np.random.default_rng(5)
    null = rng.normal(size=(1, 999))
    o_max, n_max, o_mean, n_mean, _ = lib.rbp_combine(np.asarray([0.4]), null)
    assert o_mean == o_max
    assert np.array_equal(n_mean, n_max)


def test_rbp_multi_motif_is_elementwise_max_and_mean():
    rng = np.random.default_rng(6)
    null = rng.normal(size=(3, 500))
    observed = np.asarray([0.5, 2.5, 1.0])
    o_max, n_max, o_mean, n_mean, arg = lib.rbp_combine(observed, null)
    assert arg == 1 and o_max == 2.5
    assert np.allclose(n_max, null.max(axis=0))
    assert np.allclose(n_mean, null.mean(axis=0))
    assert math.isclose(o_mean, observed.mean())


def test_rbp_combine_skips_motif_without_usable_windows():
    null = np.asarray([[0.1, 0.2, 0.3], [np.nan, np.nan, np.nan]])
    o_max, n_max, o_mean, n_mean, arg = lib.rbp_combine(np.asarray([1.0, np.nan]), null)
    assert o_max == 1.0 and o_mean == 1.0 and arg == 0
    assert np.array_equal(n_max, null[0])


def test_motif_dedupe_maps_to_every_carrier():
    motifs = ["ESRP1.TGGTGG", "motif_1.TGGTGG", "HNRNPA1.[AGT]TAGGG[AT]",
              "HNRNPA1L2.[AGT]TAGGG[AT]", "HNRNPA2B1.[AGT]TAGGG[AT]", "QKI.ACTAAC[ACG]",
              "HuR.TTTTTT[GT]", "ELAVL1.TTTTTT[GT]", "ELAVL1.ATTTA"]
    groups = lib.unique_motif_groups(motifs)
    assert groups["TGGTGG"] == ["ESRP1.TGGTGG", "motif_1.TGGTGG"]
    assert len(groups["[AGT]TAGGG[AT]"]) == 3
    assert sum(len(v) for v in groups.values()) == len(motifs)
    assert len(groups) == 5
    rbp = lib.rbp_unique_kmers(motifs, {"HuR": "ELAVL1"})
    assert rbp["ELAVL1"] == ["ATTTA", "TTTTTT[GT]"]
    assert rbp["ESRP1"] == ["TGGTGG"] and rbp["motif_1"] == ["TGGTGG"]


def test_duplicate_archives_must_be_bit_identical(tmp_path):
    rng = np.random.default_rng(21)
    first, _, _ = make_dataset(tmp_path, rng, [1] * 5, [1] * 40, 0.3, 0.3, name="ESRP1.TGGTGG")
    (tmp_path / "motif_1.TGGTGG.counts.npz").write_bytes(first.read_bytes())
    audit = v2.verify_duplicate_archives(tmp_path, {"TGGTGG": ["ESRP1.TGGTGG", "motif_1.TGGTGG"]})
    assert audit[0]["archives_identical"] and audit[0]["n_carriers"] == 2
    make_dataset(tmp_path, np.random.default_rng(22), [1] * 5, [1] * 40, 0.3, 0.3, name="motif_1.TGGTGG")
    try:
        v2.verify_duplicate_archives(tmp_path, {"TGGTGG": ["ESRP1.TGGTGG", "motif_1.TGGTGG"]})
    except ValueError as exc:
        assert "differs" in str(exc)
    else:
        raise AssertionError("differing carrier archives must be refused")


def test_planted_enrichment_is_detected_at_motif_and_rbp_level(tmp_path):
    rng = np.random.default_rng(8)
    n_windows = WINDOWS_PER_REGION * len(io.REGIONS)
    path, fg_keys, bg_keys = make_dataset(tmp_path, rng, [1] * 30 + [2] * 5, [1] * 400 + [2] * 80,
                                          planted_rate(n_windows), 0.2)
    null_path, _, _ = make_dataset(tmp_path, np.random.default_rng(8), [1] * 30 + [2] * 5,
                                   [1] * 400 + [2] * 80, 0.2, 0.2, name="null")
    selection = lib.ClusterDrawer(fg_keys, bg_keys).draw(np.random.default_rng(9), 999)
    obs, nulls = [], []
    for p in (path, null_path):
        model = calib.MotifModel(p, "up")
        result, maxima = lib.calibrate_keep(model, model.window_masks(), selection, 250)
        obs.append(result["Upstream Intron"]["observed_max_z"])
        nulls.append(maxima["Upstream Intron"])
        if p == path:
            assert result["Upstream Intron"]["calibrated_p"] == 1 / 1000
            assert result["Downstream Intron"]["calibrated_p"] > 0.01
    o_max, n_max, o_mean, n_mean, arg = lib.rbp_combine(np.asarray(obs), np.vstack(nulls))
    assert arg == 0
    assert lib.permutation_p(o_max, n_max) == 1 / 1000
    assert lib.permutation_p(o_mean, n_mean) <= 0.01
    p_minp, row = lib.rbp_minp(np.asarray(obs), np.vstack(nulls))
    assert row == 0 and p_minp <= 2 / 1000


def test_null_calibrated_p_is_uniform_under_cluster_permutation(tmp_path):
    ps = []
    for seed in range(60):
        rng = np.random.default_rng(100 + seed)
        path, fg_keys, bg_keys = make_dataset(tmp_path, rng, [1] * 12 + [2] * 4 + [3] * 2,
                                              [1] * 150 + [2] * 40 + [3] * 15, 0.5, 0.5,
                                              name="n{}".format(seed))
        model = calib.MotifModel(path, "up")
        selection = lib.ClusterDrawer(fg_keys, bg_keys).draw(np.random.default_rng(seed), 199)
        result, _ = lib.calibrate_keep(model, model.window_masks(), selection, 199)
        ps.append(result["Upstream Intron"]["calibrated_p"])
    ps = np.asarray(ps)
    assert kstest(ps, "uniform").pvalue > 0.01
    assert np.mean(ps <= 0.05) <= 0.2


def test_exceed_counts_with_ties():
    assert lib.exceed_counts(np.asarray([1.0, 3.0, 3.0, 2.0])).tolist() == [4, 2, 2, 3]


def test_minp_single_motif_equals_motif_p():
    rng = np.random.default_rng(11)
    null = np.round(rng.normal(size=(1, 1999)), 1)  # rounded: ties exercised
    for observed in (0.0, 1.3, 2.0, 5.0):
        p, row = lib.rbp_minp(np.asarray([observed]), null)
        assert row == 0
        assert p == lib.permutation_p(observed, null[0])


def test_minp_is_not_penalised_by_heavier_tailed_second_motif():
    rng = np.random.default_rng(12)
    light = rng.normal(size=4999)                 # motif A: narrow null
    heavy = 3.0 * rng.standard_t(2, size=4999)    # motif B: heavy-tailed null
    observed = np.asarray([4.5, 0.0])             # A extreme, B unremarkable
    null = np.vstack([light, heavy])
    p_a = lib.permutation_p(4.5, light)
    p_minp, row = lib.rbp_minp(observed, null)
    o_max, n_max, _, _, _ = lib.rbp_combine(observed, null)
    p_maxz = lib.permutation_p(o_max, n_max)
    assert row == 0
    assert p_minp <= 2 * p_a + 1 / 5000           # Bonferroni-like bound over 2 motifs
    assert p_maxz > 10 * p_minp                   # max-z is dominated by B's tail


def test_minp_null_is_uniform_for_two_exchangeable_motifs():
    rng = np.random.default_rng(13)
    ps = []
    for _ in range(300):
        null = rng.normal(size=(2, 199))
        observed = rng.normal(size=2)
        ps.append(lib.rbp_minp(observed, null)[0])
    assert kstest(ps, "uniform").pvalue > 0.01


def test_permutation_p_convention():
    null = np.asarray([0.0, 1.0, 2.0, np.nan])
    assert lib.permutation_p(1.0, null) == (1 + 2) / 5
    assert lib.permutation_p(5.0, null) == 1 / 5
    assert math.isnan(lib.permutation_p(math.nan, null))


def test_power_label_is_a_named_arm_decision_with_target_exon_counts():
    assert v2.power_label("X_A", [], 30, 20, 25, 18) == v2.ADEQUATE_LABEL
    label = v2.power_label("X_B", ["X_B"], 30, 20, 25, 18)
    assert label.startswith("underpowered") and "30/20 events over 25/18 target exons" in label


def test_cli_defaults_are_the_locked_two_stage_scheme():
    args = v2.build_parser().parse_args(["--arm", "A", "--counts-root", "c", "--released-root", "r",
                                         "--out-root", "o", "--alias-table", "t"])
    assert (args.permutations, args.refine_perms, args.refine_threshold, args.seed) == (2000, 100000, 0.005, 149)
    assert args.rowunit_root is None and not args.stage1_only


def test_readme_carries_the_reportable_null_and_no_review_wording():
    report = {"motif_family_size": 726, "n_unique_motifs": 121, "n_motif_keys": 126, "rbp_family_size": 600,
              "duplicate_kmers": [{"kmer": "TGGTGG", "carrier_keys": "ESRP1.TGGTGG,motif_1.TGGTGG"}],
              "stage1_permutations": 2000, "stage2_permutations": 100000, "refine_threshold": 0.005,
              "seed": 149}
    text = " ".join(v for _, v in v2.readme_rows(report, "alias.tsv", None))
    assert "Reportable p and q" in text and "TARGET EXONS" in text and "never deduplicated" in text
    assert "under review" not in text and "not yet reportable" not in text
    assert "no --rowunit-root" in text


def test_skill_full_mode_reports_events_over_target_exons_and_has_no_review_flag(tmp_path):
    import json
    import pytest
    import rmaps3_skill_run as skill
    report = tmp_path / "refinement_report.json"
    report.write_text(json.dumps({"events": {"up": 102, "dn": 107, "bg": 900},
                                  "target_exons": {"up": 91, "dn": 99, "bg": 800},
                                  "target_exons_in_both_changed_foregrounds": 0}), encoding="utf-8")
    line = skill.target_exon_line(report)
    assert line.startswith("n events (rMATS SE rows) over N target exons")
    assert "included 102 over 91; skipped 107 over 99; background 900 over 800" in line
    base = ["--mode", "full", "--arm", "A_B", "--out", "o", "--up", "u", "--dn", "d", "--bg", "b"]
    assert skill.build_parser().parse_args(base).calib_unit == "cluster"
    for retired in (["--calib-unit-decided"], ["--calib-unit", "row"], ["--unit-status", "x"]):
        with pytest.raises(SystemExit):
            skill.build_parser().parse_args(base + retired)
    source = (ROOT / "tools" / "rmaps3_skill_run.py").read_text(encoding="utf-8")
    assert "under review" not in source and "not yet reportable" not in source
