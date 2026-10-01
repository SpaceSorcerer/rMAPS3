"""Synthetic tests for tools/unit_sensitivity.py (row vs target-exon unit sensitivity).

Ported from the 2026-09-21 lab tests, which read lab count archives; every fixture here is
synthetic. The port was also checked against the lab outputs on 2026-09-24 (see LESSONS.md).
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import calibrate_ranksum as calib  # noqa: E402
import rmaps_countdist_io as io  # noqa: E402
import unit_sensitivity as unit  # noqa: E402

WINDOWS_PER_REGION = 4
MOTIFS = ("QKI.ACTAAC", "QKI.TAACTA", "RBFOX2.TGCATG", "SRSF1.GGAGGA")
ARM = "SYN_KO_B"


def exon_id(chrom, start, end, flank):
    return "{}:+:{}:{}:{}:{}:{}:{}".format(chrom, start, end, start - 300 - flank, start - 200,
                                           end + 200, end + 300 + flank)


def group_ids(sizes, offset):
    """Eight-field ids; rows of one target exon differ only in their flanking exons."""
    ids = []
    for c, size in enumerate(sizes):
        start = offset + 1000 * c
        for r in range(size):
            ids.append(exon_id("chr1", start, start + 40 + (c % 7) * 10, r))
    return np.asarray(ids)


def build_arm(root: Path, fg_sizes=(1,) * 12 + (2,) * 4 + (3,) * 2,
              dn_sizes=(1,) * 10 + (2,) * 3, bg_sizes=(1,) * 150 + (2,) * 40 + (3,) * 10,
              seed=11):
    """Counts archives, released root tables, events_*.tsv and an alias table for one arm."""
    rng = np.random.default_rng(seed)
    n_windows = WINDOWS_PER_REGION * len(io.REGIONS)
    region = np.repeat(np.arange(len(io.REGIONS), dtype=np.int8), WINDOWS_PER_REGION)
    position = np.tile(np.arange(WINDOWS_PER_REGION, dtype=np.int32), len(io.REGIONS))
    ids = {"up": group_ids(fg_sizes, 10 ** 6), "dn": group_ids(dn_sizes, 2 * 10 ** 6),
           "bg": group_ids(bg_sizes, 3 * 10 ** 6)}
    sizes = {"up": fg_sizes, "dn": dn_sizes, "bg": bg_sizes}
    counts_dir = root / "counts" / ARM
    counts_dir.mkdir(parents=True)
    roots = {"up": {}, "dn": {}}
    for m, motif in enumerate(MOTIFS):
        store = {"region_of_position": region, "position": position}
        matrices = {}
        for group in io.GROUPS:
            rate = 0.25 + (0.6 if (group == "up" and m == 0) else 0.0)
            rows = []
            for size in sizes[group]:
                base = rng.poisson(rate, n_windows).astype(np.int16)
                for _ in range(size):
                    rows.append(base)
            matrices[group] = np.vstack(rows)
            io.pack_group(store, group, matrices[group])
            store[group + "_exon_id"] = ids[group]
        np.savez(counts_dir / (motif + ".counts.npz"), **store)
        kmax = int(max(matrices[g].max() for g in io.GROUPS))
        for direction in ("up", "dn"):
            z = io.rank_statistics(io.count_histograms(matrices[direction], kmax),
                                   io.count_histograms(matrices["bg"], kmax))[3]
            p = io.p_from_z(z)
            roots[direction][motif] = [float(p[region == r].min()) for r in range(len(io.REGIONS))]
    released = root / "released" / ARM
    released.mkdir(parents=True)
    for direction in ("up", "dn"):
        with open(released / ("pVal." + direction + ".vs.bg.RNAmap.txt"), "w", newline="\n") as h:
            h.write("\t".join(["RBP"] + [io.ROOT_COLUMNS[r] for r in io.REGIONS]) + "\n")
            for motif, values in roots[direction].items():
                h.write("\t".join([motif] + [repr(v) for v in values]) + "\n")
    events = root / "events" / ARM
    events.mkdir(parents=True)
    dpsi = {}
    for group in io.GROUPS:
        dpsi[group] = np.round(rng.uniform(0.1, 0.6, ids[group].shape[0]), 3)
        with open(events / ("events_" + group + ".tsv"), "w", newline="\n") as h:
            h.write("\t".join(("ID",) + unit.EVENT_KEY_COLUMNS + ("IncLevelDifference",)) + "\n")
            for i, (key, value) in enumerate(zip(ids[group], dpsi[group])):
                h.write("\t".join([str(i)] + key.split(":") + [repr(float(value))]) + "\n")
    alias = root / "alias.tsv"
    alias.write_text("table_name\thgnc_symbol\nQKI\tQKI\nRBFOX2\tRBFOX2\nSRSF1\tSRSF1\n",
                     encoding="utf-8")
    cfg = unit.settings(counts_root=root / "counts", released_root=root / "released",
                        alias_table=alias, out_root=root / "out", event_dirs={ARM: str(events)},
                        stage1_perms=199, stage2_perms=499, refine_threshold=0.01)
    return cfg, ids, dpsi


@pytest.fixture()
def arm(tmp_path):
    return build_arm(tmp_path)


def test_full_row_subset_model_equals_motif_model(arm):
    cfg, _, _ = arm
    npz = Path(cfg.counts_root) / ARM / (MOTIFS[0] + ".counts.npz")
    for direction in ("up", "dn"):
        reference = calib.MotifModel(npz, direction)
        mine = unit.SubsetModel(npz, direction)
        assert (mine.n1, mine.n0) == (reference.n1, reference.n0)
        for name in ("observed_z", "baseline", "sigma", "fg_carrying", "bg_carrying"):
            np.testing.assert_array_equal(getattr(mine, name), getattr(reference, name))
        np.testing.assert_array_equal(mine.delta.data, reference.delta.data)


def test_subset_model_matches_a_dense_recomputation(arm):
    cfg, _, _ = arm
    npz = Path(cfg.counts_root) / ARM / (MOTIFS[0] + ".counts.npz")
    rows = np.arange(0, 15, dtype=np.int64)
    model = unit.SubsetModel(npz, "up", rows, rows)
    with np.load(npz, allow_pickle=False) as data:
        fg = io.unpack_group(data, "up")[rows]
        bg = io.unpack_group(data, "bg")[rows]
    kmax = int(max(fg.max(initial=0), bg.max(initial=0)))
    z = io.rank_statistics(io.count_histograms(fg, kmax), io.count_histograms(bg, kmax))[3]
    np.testing.assert_allclose(model.observed_z, z, rtol=0, atol=1e-9, equal_nan=True)


def test_dedup_keeps_one_row_per_target_exon_with_the_largest_dpsi(arm):
    cfg, ids, _ = arm
    for group in io.GROUPS:
        keys = unit.target_keys(ids[group])
        dpsi, source = unit.dpsi_for_group(cfg, ARM, group, ids[group])
        assert source.endswith(".tsv"), source
        keep = unit.dedup_index(keys, dpsi)
        assert keep.shape[0] == len(set(keys.tolist()))
        assert np.all(np.diff(keep) > 0)
        for rows in unit.cluster_rows(keys).values():
            chosen = [r for r in rows.tolist() if r in set(keep.tolist())]
            assert len(chosen) == 1
            assert dpsi[chosen[0]] == float(np.nanmax(dpsi[rows]))


def test_dedup_ties_and_nan_fall_back_to_first_in_file():
    keys = np.asarray(["x", "x", "y"], dtype="U")
    np.testing.assert_array_equal(unit.dedup_index(keys, np.asarray([0.2, 0.2, 0.9])), [0, 2])
    np.testing.assert_array_equal(unit.dedup_index(keys, np.asarray([math.nan, math.nan, 0.9])),
                                  [0, 2])


def test_dpsi_alignment_accepts_only_a_chr_prefix_rename(arm, tmp_path):
    cfg, ids, dpsi = arm
    renamed = np.asarray([s[3:] for s in ids["up"]])
    values, source = unit.dpsi_for_group(cfg, ARM, "up", renamed)
    np.testing.assert_array_equal(values, dpsi["up"])
    assert "chromosome-prefix normalised match" in source
    values, source = unit.dpsi_for_group(cfg, ARM, "up", ids["up"][::-1])
    assert source.startswith("key_sequence_mismatch") and np.isnan(values).all()
    values, source = unit.dpsi_for_group(cfg, "OTHER_ARM", "up", ids["up"])
    assert source.startswith("missing:") and np.isnan(values).all()


def test_zero_duplicate_foreground_keeps_every_row():
    keys = unit.target_keys(group_ids((1,) * 9, 0))
    np.testing.assert_array_equal(unit.dedup_index(keys, np.linspace(0.1, 0.9, 9)), np.arange(9))


def test_power_label_is_a_named_arm_decision():
    assert unit.power_label("A", ("B",), 5, 6) == "adequate"
    assert unit.power_label("B", ("B",), 5, 6).startswith(
        "underpowered (n changed included/skipped = 5/6)")


def test_cli_end_to_end_duplication_b_c_compare(tmp_path):
    cfg, ids, _ = build_arm(tmp_path)
    base = ["--arms", ARM, "--out-root", str(cfg.out_root), "--counts-root", str(cfg.counts_root),
            "--event-dir", "{}={}".format(ARM, cfg.event_dirs[ARM])]
    stats = ["--released-root", str(cfg.released_root), "--alias-table", str(cfg.alias_table),
             "--permutations", "199", "--refine-perms", "499", "--refine-threshold", "0.01"]
    assert unit.main(["--mode", "duplication"] + base) == 0
    duplication = unit.read_tsv(Path(cfg.out_root) / "duplication_table.tsv")
    by_group = {r["group"]: r for r in duplication}
    assert int(by_group["up"]["n_rows"]) == ids["up"].shape[0]
    assert int(by_group["up"]["n_target_exons"]) == len(set(unit.target_keys(ids["up"]).tolist()))
    assert int(by_group["bg"]["max_cluster_size"]) == 3
    assert unit.main(["--mode", "B"] + base + stats) == 0
    assert unit.main(["--mode", "C"] + base + stats) == 0
    out = Path(cfg.out_root) / ARM
    b_rows = unit.read_tsv(out / "treatment_B_condensed_per_rbp.tsv")
    c_rows = unit.read_tsv(out / "treatment_C_condensed_per_rbp.tsv")
    assert len(b_rows) == len(c_rows) == 3 * 2 * len(io.POOL_TO_REGIONS)
    assert {r["n_fg_exons"] for r in b_rows if r["direction"] == "up"} == {
        str(len(set(unit.target_keys(ids["up"]).tolist())))}
    assert {r["n_fg_exons"] for r in c_rows if r["direction"] == "up"} == {str(ids["up"].shape[0])}
    summary_a = tmp_path / "summary_a" / ARM
    summary_a.mkdir(parents=True)
    (summary_a / "condensed_per_rbp.tsv").write_text(
        (out / "treatment_C_condensed_per_rbp.tsv").read_text(encoding="utf-8"), encoding="utf-8")
    assert unit.main(["--mode", "compare", "--arms", ARM, "--out-root", str(cfg.out_root),
                      "--summary-a-root", str(tmp_path / "summary_a"),
                      "--positive-control", "QKI", "--control-arms", ARM]) == 0
    panels = unit.read_tsv(Path(cfg.out_root) / "panel_summary_ABC.tsv")
    assert len(panels) == 2 * len(io.PLOT_POOLS)
    for row in panels:
        assert float(row["spearman_calibrated_p_A_vs_C"]) == pytest.approx(1.0)
    control = unit.read_tsv(Path(cfg.out_root) / "qki_control_ABC.tsv")
    assert [r["pooled_region"] for r in control] == ["Upstream Intron", "Downstream Intron"]
    assert (Path(cfg.out_root) / "comparison_ABC.xlsx").stat().st_size > 0
    assert "exit=0" in (Path(cfg.out_root) / "command.log").read_text(encoding="utf-8")


def test_treatment_b_refuses_without_an_event_dir(tmp_path):
    cfg, _, _ = build_arm(tmp_path)
    with pytest.raises(SystemExit, match="needs --event-dir"):
        unit.main(["--mode", "B", "--arms", ARM, "--out-root", str(cfg.out_root),
                   "--counts-root", str(cfg.counts_root), "--released-root", str(cfg.released_root),
                   "--alias-table", str(cfg.alias_table)])
