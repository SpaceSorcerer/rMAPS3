"""Synthetic tests for tools/length_matched.py (length-matched background sensitivity).

Ported from the 2026-09-22 lab tests, which read lab count archives; every fixture here is
synthetic. The port was also checked against the lab outputs on 2026-09-24 (see LESSONS.md).
"""
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import ks_2samp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tests"))

import calibrate_ranksum as calib  # noqa: E402
import length_matched as L  # noqa: E402
import rmaps_countdist_io as io  # noqa: E402
import unit_sensitivity as unit  # noqa: E402
from test_unit_sensitivity import ARM, build_arm  # noqa: E402

FG = (1,) * 30 + (2,) * 6
BG = (1,) * 400 + (2,) * 80 + (3,) * 20


def arm_data(tmp_path):
    cfg, _, _ = build_arm(tmp_path, fg_sizes=FG, dn_sizes=FG, bg_sizes=BG)
    ids, motifs, _ = unit.exon_ids_for_arm(cfg, ARM)
    keys = {g: unit.target_keys(ids[g]) for g in io.GROUPS}
    lens = {g: L.exon_lengths(ids[g]) for g in io.GROUPS}
    return cfg, ids, motifs, keys, lens


def test_exon_length_is_end_minus_zero_based_start():
    ids = np.asarray(["chr1:+:100:160:1:2:3:4", "chr2:-:10:11:1:2:3:4"])
    np.testing.assert_array_equal(L.exon_lengths(ids), [60, 1])


def test_matched_bins_follow_the_changed_set(tmp_path):
    _, _, _, keys, lens = arm_data(tmp_path)
    for d in ("up", "dn"):
        rows, info = L.stratified_background(keys[d], lens[d], keys["bg"], lens["bg"],
                                             L.RESAMPLE_SEEDS[0], d)
        edges = np.asarray(info["edges"])
        fg_u = L.unique_lengths(keys[d], lens[d])
        bg_u = L.unique_lengths(keys["bg"][rows], lens["bg"][rows])
        target = np.bincount(L.assign_bin(fg_u, edges), minlength=10)
        got = np.bincount(L.assign_bin(bg_u, edges), minlength=10)
        assert got.tolist() == info["bg_exons_per_bin_after"]
        assert np.all(np.abs(got - info["ratio_bg_to_changed"] * target) < 1.0 + 1e-9)
        assert ks_2samp(fg_u, bg_u).statistic <= ks_2samp(
            fg_u, L.unique_lengths(keys["bg"], lens["bg"])).statistic


def test_whole_clusters_are_drawn_and_never_split(tmp_path):
    _, _, _, keys, lens = arm_data(tmp_path)
    clusters = unit.cluster_rows(keys["bg"])
    for d in ("up", "dn"):
        rows, _ = L.stratified_background(keys[d], lens[d], keys["bg"], lens["bg"],
                                          L.RESAMPLE_SEEDS[0], d)
        assert np.unique(rows).shape[0] == rows.shape[0]
        chosen = set(keys["bg"][rows].tolist())
        np.testing.assert_array_equal(np.sort(rows),
                                      np.sort(np.concatenate([clusters[k] for k in chosen])))
        drawer = unit.ClusterDrawer(keys[d], keys["bg"][rows])
        pooled = np.concatenate([keys[d], keys["bg"][rows]])
        full = unit.cluster_rows(pooled)
        for draw in drawer.draw(np.random.default_rng(1), 30):
            picked = set(draw.tolist())
            for k in set(pooled[draw].tolist()):
                assert set(full[k].tolist()) <= picked


def test_resample_seed_determinism_and_dependence(tmp_path):
    _, _, _, keys, lens = arm_data(tmp_path)
    a, _ = L.stratified_background(keys["up"], lens["up"], keys["bg"], lens["bg"], 20260922, "up")
    b, _ = L.stratified_background(keys["up"], lens["up"], keys["bg"], lens["bg"], 20260922, "up")
    c, _ = L.stratified_background(keys["up"], lens["up"], keys["bg"], lens["bg"], 20260923, "up")
    np.testing.assert_array_equal(a, b)
    assert not np.array_equal(a, c)


def test_size_control_draws_the_same_number_of_clusters(tmp_path):
    _, _, _, keys, lens = arm_data(tmp_path)
    _, matched = L.stratified_background(keys["up"], lens["up"], keys["bg"], lens["bg"], 7, "up")
    rows, control = L.stratified_background(keys["up"], lens["up"], keys["bg"], lens["bg"], 7, "up",
                                            size_control=True)
    assert control["n_bg_exons_after"] == matched["n_bg_exons_after"]
    assert len(set(keys["bg"][rows].tolist())) == matched["n_bg_exons_after"]


def test_full_background_subset_equals_the_treatment_c_calibration(tmp_path):
    cfg, _, motifs, keys, _ = arm_data(tmp_path)
    npz = Path(cfg.counts_root) / ARM / (motifs[0] + ".counts.npz")
    all_bg = np.arange(keys["bg"].shape[0], dtype=np.int64)
    subset = unit.SubsetModel(npz, "up", None, all_bg)
    full = calib.MotifModel(npz, "up")
    selection = unit.ClusterDrawer(keys["up"], keys["bg"]).draw(unit.stage_rng(cfg, 1, "up"), 199)
    a = calib.calibrate(subset, subset.window_masks(), selection, 50)
    b = calib.calibrate(full, full.window_masks(), selection, 50)
    for name in a:
        pa, pb = a[name]["calibrated_p"], b[name]["calibrated_p"]
        assert (math.isnan(pa) and math.isnan(pb)) or pa == pb


def test_seed_labels_keep_the_names_the_figure_builder_reads():
    assert L.seed_labels((20260922, 20260923)) == {20260922: "lengthmatched",
                                                   20260923: "lengthmatched_seed2"}


def test_cli_run_then_compare(tmp_path):
    cfg, _, _, _, _ = arm_data(tmp_path)
    common = ["--permutations", "199", "--refine-perms", "499", "--refine-threshold", "0.01"]
    assert unit.main(["--mode", "C", "--arms", ARM, "--out-root", str(tmp_path / "unit"),
                      "--counts-root", str(cfg.counts_root), "--released-root",
                      str(cfg.released_root), "--alias-table", str(cfg.alias_table)] + common) == 0
    out = tmp_path / "length"
    assert L.main(["--mode", "run", "--arms", ARM, "--out-root", str(out), "--counts-root",
                   str(cfg.counts_root), "--alias-table", str(cfg.alias_table)] + common) == 0
    for name in ("lengthmatched", "lengthmatched_seed2"):
        report = (out / ARM / (name + "_report.json")).read_text(encoding="utf-8")
        assert '"length_match"' in report
        assert len(unit.read_tsv(out / ARM / (name + "_condensed_per_rbp.tsv"))) == 24
    assert L.main(["--mode", "compare", "--arms", ARM, "--out-root", str(out), "--c-root",
                   str(tmp_path / "unit"), "--positive-control", "QKI", "--control-arms", ARM]) == 0
    panels = unit.read_tsv(out / "panel_summary_C_vs_matched.tsv")
    assert len(panels) == 2 * len(io.PLOT_POOLS)
    control = unit.read_tsv(out / "qki_control_C_vs_matched.tsv")
    assert [(r["pooled_region"], r["direction"]) for r in control] == [
        ("Upstream Intron", "up"), ("Exon Body", "up"), ("Downstream Intron", "dn")]
    assert (out / "comparison_C_vs_matched.xlsx").stat().st_size > 0
