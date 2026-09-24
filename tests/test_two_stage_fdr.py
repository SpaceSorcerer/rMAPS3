"""Two-stage calibration v2.1: self-triggered promotion, its validity, and the one contract in every document.

A test reports its stage-2 p only when its OWN stage-1 p <= t (tools/rmaps_calib_v2_lib.two_stage_p). This file checks:
(a) the switched p is super-uniform at alpha below AND above t, on >= 2,000 independent null permutation tests and
    on the pooled tests of real arm-runner calls; the one-sided KS test is shown to have power;
(b) under the global null, with correlated motifs of one RBP, P(any BH rejection) = FDR stays <= nominal + margin at
    the RBP min-P level and the motif level, through the real arm runner;
(c) the old group-triggered rule (a test switched to stage 2 because ANOTHER test of its pair or RBP was promoted) is
    invalid: a constructed group whose partners promote the test while its own stage-1 p > t breaks super-uniformity,
    and the arm runner never switches such a test.
Seeds are fixed; every result is deterministic.
"""
import argparse
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

CONTRACT = "every calibrated p is a valid permutation p"
SELF_TRIGGER = "only when its own stage-1 p"
PER_REGION = 4
MOTIFS = ["RA.AAAC", "RB.CCCA", "RB.GGGA", "RC.TTTG", "RD.ACGT", "RD.CAGT", "RE.GTCA", "RF.TGCA"]
B1, B2, THRESHOLD = 99, 999, 0.05
ALPHAS = (0.01, 0.05, 0.1, 0.2)          # below, at and above the refine threshold t = 0.05
N_ARMS = 200                              # global-null arms through the real arm runner
N_INDEPENDENT = 4000                      # independent null permutation tests for (a)
FDR_ALPHA = 0.05


def margin(alpha, n):
    return 3 * math.sqrt(alpha * (1 - alpha) / n)


# ------------------------------------------------------------------ synthetic arms
def cluster_block(rng, sizes, offset, n_windows, rate, base=None):
    """Rows of one target exon share counts; `base` (one vector per cluster) adds a component shared across motifs."""
    rows, ids = [], []
    for c, size in enumerate(sizes):
        counts = rng.poisson(rate, n_windows)
        if base is not None:
            counts = counts + base[c]
        for r in range(size):
            rows.append(counts.astype(np.int16))
            ids.append("chr1:+:{}:{}:{}:1:2:3".format(offset + 10 * c, offset + 10 * c + 5, r))
    return np.vstack(rows), np.asarray(ids, dtype="U")


def null_arm(root: Path, arm: str, rng, correlated: bool = False):
    """All-null clustered archives plus root tables that reproduce the released regional minima.

    correlated=True gives the motifs of one RBP (RB, RD) a shared per-target-exon count component, so their
    statistics, and the RBP min-P over them, are positively dependent."""
    n_windows = PER_REGION * len(io.REGIONS)
    region = np.repeat(np.arange(len(io.REGIONS), dtype=np.int8), PER_REGION)
    position = np.tile(np.arange(PER_REGION, dtype=np.int32), len(io.REGIONS))
    counts_dir = root / "counts" / arm
    released_dir = root / "released" / arm
    counts_dir.mkdir(parents=True)
    released_dir.mkdir(parents=True)
    sizes = {"up": [1] * 14 + [2] * 3, "dn": [1] * 12 + [2] * 4, "bg": [1] * 150 + [2] * 30 + [3] * 6}
    offsets = {"up": 0, "dn": 10 ** 5, "bg": 10 ** 6}
    shared = {}
    if correlated:
        for rbp in ("RB", "RD"):
            rate = float(rng.uniform(0.3, 0.6))
            shared[rbp] = {g: [rng.poisson(rate, n_windows) for _ in sizes[g]] for g in io.GROUPS}
    roots = {"up": {}, "dn": {}}
    for motif in MOTIFS:
        rbp = lib.table_name_of(motif)
        rate = float(rng.uniform(0.05, 0.2) if rbp in shared else rng.uniform(0.15, 0.6))
        store = {"region_of_position": region, "position": position}
        for group in io.GROUPS:
            base = shared[rbp][group] if rbp in shared else None
            matrix, ids = cluster_block(rng, sizes[group], offsets[group], n_windows, rate, base)
            io.pack_group(store, group, matrix)
            store[group + "_exon_id"] = ids
        path = counts_dir / (motif + ".counts.npz")
        np.savez(path, **store)
        for d in ("up", "dn"):
            model = calib.MotifModel(path, d)
            p = io.p_from_z(model.observed_z)
            roots[d][motif] = {r: float(p[model.region_index == i].min()) for i, r in enumerate(io.REGIONS)}
    for d in ("up", "dn"):
        with open(released_dir / ("pVal." + d + ".vs.bg.RNAmap.txt"), "w", newline="\n", encoding="utf-8") as h:
            h.write("\t".join(["RBP"] + [io.ROOT_COLUMNS[r] for r in io.REGIONS]) + "\n")
            for motif in MOTIFS:
                h.write("\t".join([motif] + [repr(roots[d][motif][r]) for r in io.REGIONS]) + "\n")
    alias = root / "alias.tsv"
    if not alias.exists():
        alias.write_text("table_name\thgnc_symbol\n", encoding="utf-8")
    return argparse.Namespace(
        arm=arm, counts_root=str(root / "counts"), released_root=str(root / "released"),
        out_root=str(root / "out"), alias_table=str(alias), rowunit_root=None, underpowered_arms=[],
        permutations=B1, refine_perms=B2, refine_threshold=THRESHOLD, chunk_size=500,
        seed=int(rng.integers(1, 10 ** 6)), stage1_only=False)


# ------------------------------------------------------------------ independent permutation tests
def permutation_stage_p(data, n1, B, rng):
    """(1 + #{S* >= s}) / (B + 1) per row, S = sum of the first n1 values, B label permutations per row."""
    observed = data[:, :n1].sum(axis=1)
    exceed = np.zeros(data.shape[0], dtype=np.int64)
    for _ in range(B):
        picked = np.argsort(rng.random(data.shape), axis=1)[:, :n1]
        exceed += np.take_along_axis(data, picked, axis=1).sum(axis=1) >= observed - 1e-12
    return (1 + exceed) / (B + 1)


def independent_null_pairs():
    """N_INDEPENDENT null datasets; stage-1 and stage-2 p of each from independent streams on the same data."""
    data = np.random.default_rng(1).standard_normal((N_INDEPENDENT, 32))
    p1 = permutation_stage_p(data, 8, B1, np.random.default_rng(2))
    p2 = permutation_stage_p(data, 8, B2, np.random.default_rng(3))
    return p1, p2


def self_triggered(p1, p2):
    return np.asarray([lib.two_stage_p(float(a), float(b), THRESHOLD)[0] for a, b in zip(p1, p2)])


def assert_super_uniform(p, n_units, label):
    for alpha in ALPHAS:
        rate = float(np.mean(p <= alpha))
        assert rate <= alpha + margin(alpha, n_units), (label, alpha, rate)


def test_a_self_triggered_p_is_super_uniform_below_and_above_t_and_ks_has_power():
    p1, p2 = independent_null_pairs()
    reported = self_triggered(p1, p2)
    switched = int(np.sum(p1 <= THRESHOLD))
    assert switched > 100, switched                      # the mixture is exercised
    assert_super_uniform(reported, N_INDEPENDENT, "independent tests")
    assert kstest(reported, "uniform", alternative="greater").pvalue > 0.01
    # Power: the same one-sided KS on the same n rejects a p that is only 10 % anti-conservative.
    assert kstest(0.9 * reported, "uniform", alternative="greater").pvalue < 1e-3


def group_triggered_construction():
    """(c): p1, p2 independent valid p (the independent-stream framing) and 15 valid partner p-values of the same
    pair/RBP that each fall <= t on a slice (0.2 + j t, 0.2 + (j + 1) t] of the test's OWN p1. Every partner is
    super-uniform: P(partner <= x) = x for x <= t and = t above. The group is promoted when any member's p1 <= t."""
    rng = np.random.default_rng(5)
    lattice = np.arange(1, B1 + 2) / (B1 + 1)
    p1 = rng.choice(lattice, N_INDEPENDENT)
    p2 = rng.choice(np.arange(1, B2 + 2) / (B2 + 1), N_INDEPENDENT)
    partners = []
    for j in range(15):
        low = 0.2 + j * THRESHOLD
        inside = (p1 > low) & (p1 <= low + THRESHOLD)
        partners.append(np.where(inside, p1 - low, 1.0))
    partners = np.vstack(partners)
    assert np.all(np.mean(partners <= THRESHOLD, axis=1) <= THRESHOLD + margin(THRESHOLD, N_INDEPENDENT))
    group_promoted = (p1 <= THRESHOLD) | np.any(partners <= THRESHOLD, axis=0)
    return p1, p2, group_promoted


def test_c_group_triggered_promotion_breaks_super_uniformity_self_triggered_does_not():
    p1, p2, group_promoted = group_triggered_construction()
    promoted_by_group = group_promoted & (p1 > THRESHOLD)
    assert promoted_by_group.sum() > 1000, int(promoted_by_group.sum())
    old_rule = np.where(group_promoted, p2, p1)
    old_rate = float(np.mean(old_rule <= 0.2))
    assert old_rate > 0.2 + margin(0.2, N_INDEPENDENT), old_rate       # the old rule fails here
    assert kstest(old_rule, "uniform", alternative="greater").pvalue < 1e-6
    # Production rule on the SAME draws: stage-2 p exists for every promoted test, yet only own p1 <= t switches.
    reported = np.asarray([lib.two_stage_p(float(a), float(b) if g else None, THRESHOLD)[0]
                           for a, b, g in zip(p1, p2, group_promoted)])
    assert_super_uniform(reported, N_INDEPENDENT, "group construction, self-triggered")
    assert np.array_equal(reported[promoted_by_group], p1[promoted_by_group])


# ------------------------------------------------------------------ the real arm runner
def run_arms(tmp_path, n, correlated, seed):
    rng = np.random.default_rng(seed)
    results = []
    for rep in range(n):
        args = null_arm(tmp_path, "N{}".format(rep), rng, correlated=correlated)
        results.append((args,) + tuple(v2.run_arm(args, lambda message: None)[1:5]))
    return results


def test_c_arm_runner_switches_only_tests_whose_own_stage1_p_is_at_most_t(tmp_path):
    """Stage 1 is the same stream with and without --stage1-only, so a stage-1-only rerun gives every test's own p1."""
    checked_motif = checked_rbp = neighbour_promoted = 0
    for args, per_motif, condensed, _, report in run_arms(tmp_path / "full", 12, True, 11):
        only = argparse.Namespace(**vars(args))
        only.stage1_only, only.out_root = True, str(tmp_path / "stage1" / args.arm)
        _, per1, cond1, _, _, _ = v2.run_arm(only, lambda message: None)
        first = {(r["motif_key"], r["direction"], r["region"]): r for r in per1}
        for r in per_motif:
            own = first[(r["motif_key"], r["direction"], r["region"])]
            for suffix in ("", "_pooled"):
                p1 = own["calibrated_p" + suffix]
                stage = r["calib_stage" + suffix]
                assert stage == (2 if math.isfinite(p1) and p1 <= THRESHOLD else 1), (r["motif_key"], suffix, p1)
                assert r["calib_perms_used" + suffix] == (B2 if stage == 2 else B1)
                if stage == 1:
                    assert r["calibrated_p" + suffix] == p1 or (math.isnan(p1) and math.isnan(r["calibrated_p" + suffix]))
                checked_motif += 1
        rbp1 = {(r["RBP"], r["direction"], r["pooled_region"]): r for r in cond1}
        for r in condensed:
            own = rbp1[(r["RBP"], r["direction"], r["pooled_region"])]
            for p_col, stage_col, perms_col in (
                    ("rbp_calibrated_p_minp", "rbp_calib_stage", "rbp_calib_perms_used"),
                    ("rbp_calibrated_p_maxz", "rbp_calib_stage_maxz", "rbp_calib_perms_used_maxz"),
                    ("rbp_calibrated_p_meanz", "rbp_calib_stage_meanz", "rbp_calib_perms_used_meanz")):
                p1 = own[p_col]
                expected = 2 if math.isfinite(p1) and p1 <= THRESHOLD else 1
                assert r[stage_col] == expected, (r["RBP"], p_col, p1, r[stage_col])
                assert r[perms_col] == (B2 if expected == 2 else B1)
                if expected == 1:
                    assert r[p_col] == p1 or (math.isnan(p1) and math.isnan(r[p_col]))
                    # stage-2 draws exist for this RBP-direction (another statistic or pool promoted it)
                    neighbour_promoted += r["rbp_calib_stage"] == 2 or r["rbp_calib_stage_maxz"] == 2 \
                        or r["rbp_calib_stage_meanz"] == 2
                checked_rbp += 1
        assert report["promotion_rule"].startswith("self-triggered")
    assert checked_motif > 0 and checked_rbp > 0
    assert neighbour_promoted > 0, "no test sat next to a promoted neighbour; the check did not exercise (c)"


def test_b_global_null_fdr_with_correlated_motifs_and_pooled_super_uniformity(tmp_path):
    results = run_arms(tmp_path, N_ARMS, True, 20260924)
    motif_any = rbp_any = 0
    pooled_motif, pooled_rbp, fixed_p, perms = [], [], [], set()
    for _, per_motif, condensed, _, report in results:
        plotted = [r for r in per_motif if r["plot"]]
        motif_any += any(r["calibrated_q"] < FDR_ALPHA for r in plotted)
        rbp_plotted = [r for r in condensed if r["plot"]]
        rbp_any += any(r["rbp_calibrated_q_minp"] < FDR_ALPHA for r in rbp_plotted)
        # one p per (unique motif, direction, pooled region) test and per (RBP, direction, pooled region) test
        seen = {}
        for r in plotted:
            seen[(r["unique_motif_id"], r["direction"], r["pooled_region"])] = r["calibrated_p_pooled"]
            perms.add(r["calib_perms_used_pooled"])
        pooled_motif += [p for p in seen.values() if math.isfinite(p)]
        pooled_rbp += [r["rbp_calibrated_p_minp"] for r in rbp_plotted if math.isfinite(r["rbp_calibrated_p_minp"])]
        fixed_p.append(next(r["rbp_calibrated_p_minp"] for r in condensed if r["RBP"] == "RB"
                            and r["direction"] == "up" and r["pooled_region"] == "Upstream Intron"))
    assert perms == {B1, B2}                                    # both resolutions are reported
    assert len(pooled_motif) >= 2000 and len(pooled_rbp) >= 2000
    # Global null: FDR = P(any BH rejection) <= nominal + Monte Carlo margin at the RBP min-P and motif levels.
    assert rbp_any / N_ARMS <= FDR_ALPHA + margin(FDR_ALPHA, N_ARMS), (rbp_any, N_ARMS)
    assert motif_any / N_ARMS <= FDR_ALPHA + margin(FDR_ALPHA, N_ARMS), (motif_any, N_ARMS)
    # Pooled switched p, below and above t. Tests within an arm are dependent, so the margin uses the arm count.
    assert_super_uniform(np.asarray(pooled_motif), N_ARMS, "arm runner, motif level")
    assert_super_uniform(np.asarray(pooled_rbp), N_ARMS, "arm runner, RBP min-P")
    # One correlated-motif RBP min-P test per independent arm: not anti-conservative (one-sided KS).
    assert kstest(np.asarray(fixed_p), "uniform", alternative="greater").pvalue > 0.01


def test_correlated_arm_really_correlates_the_rbp_motifs(tmp_path):
    args = null_arm(tmp_path, "C0", np.random.default_rng(3), correlated=True)
    with np.load(Path(args.counts_root) / "C0" / "RB.CCCA.counts.npz") as a, \
            np.load(Path(args.counts_root) / "C0" / "RB.GGGA.counts.npz") as b:
        x, y = io.unpack_group(a, "bg").ravel(), io.unpack_group(b, "bg").ravel()
    assert np.corrcoef(x, y)[0, 1] > 0.5


# ------------------------------------------------------------------ one contract in every document
def test_one_inferential_contract_in_lessons_readme_validity_doc_legend_workbook_and_readout():
    lessons = (ROOT / "LESSONS.md").read_text(encoding="utf-8")
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    validity = (ROOT / "docs" / "two_stage_validity.md").read_text(encoding="utf-8")
    builder = (ROOT / "tools" / "build_region_lollipops_v4.py").read_text(encoding="utf-8")
    assert len(validity.splitlines()) <= 25
    for name, text in (("LESSONS.md", lessons), ("README.md", readme), ("figure legend", builder),
                       ("two_stage_validity.md", validity)):
        assert CONTRACT in text.lower(), name
        assert "no bh guarantee" not in text.lower(), name
        assert "PRDS" in text, name
    for name, text in (("LESSONS.md", lessons), ("README.md", readme), ("two_stage_validity.md", validity)):
        assert SELF_TRIGGER in text.lower(), name
        assert "mixed resolutions do not break that" not in text, name
    assert "(b+1)/(B+1)" in validity and "PRDS is assumed" in validity
    report = {"motif_family_size": 726, "n_unique_motifs": 121, "n_motif_keys": 126, "rbp_family_size": 600,
              "motif_bh_divisor": 726, "rbp_bh_divisor": {"minp": 600, "maxz": 600, "meanz": 600},
              "duplicate_kmers": [], "stage1_permutations": 2000, "stage2_permutations": 100000,
              "refine_threshold": 0.005, "seed": 149}
    workbook = " ".join(v for _, v in v2.readme_rows(report, "alias.tsv", None)).lower()
    assert CONTRACT in workbook and SELF_TRIGGER in workbook and "calib_perms_used" in workbook
    assert "1/2,001" in workbook and "1/100,001" in workbook
    assert "no bh guarantee" not in workbook and "prds is assumed" in workbook
