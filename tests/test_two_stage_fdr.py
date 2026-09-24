"""Two-stage calibration v2.1 under the global null, and the one inferential contract in every document.

The simulation runs the real arm runner (calibrate_ranksum_v2.run_arm) on synthetic all-null arms
at small B, with a refine threshold that promotes tests often, so the reported p mixes stage-1 and
stage-2 resolutions exactly as in production. Under the global null the FDR equals the chance of
any BH rejection. Seeds are fixed; the result is deterministic.
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
import rmaps_countdist_io as io  # noqa: E402

CONTRACT = "every calibrated p is a valid permutation p"
PER_REGION = 4
MOTIFS = ["RA.AAAC", "RB.CCCA", "RB.GGGA", "RC.TTTG", "RD.ACGT", "RD.CAGT", "RE.GTCA", "RF.TGCA"]
B1, B2, THRESHOLD = 99, 999, 0.05
N_REPLICATES = 120
ALPHA = 0.05


def cluster_block(rng, sizes, offset, n_windows, rate):
    rows, ids = [], []
    for c, size in enumerate(sizes):
        counts = rng.poisson(rate, n_windows).astype(np.int16)
        for r in range(size):
            rows.append(counts)
            ids.append("chr1:+:{}:{}:{}:1:2:3".format(offset + 10 * c, offset + 10 * c + 5, r))
    return np.vstack(rows), np.asarray(ids, dtype="U")


def null_arm(root: Path, arm: str, rng):
    """All-null clustered archives plus root tables that reproduce the released regional minima."""
    n_windows = PER_REGION * len(io.REGIONS)
    region = np.repeat(np.arange(len(io.REGIONS), dtype=np.int8), PER_REGION)
    position = np.tile(np.arange(PER_REGION, dtype=np.int32), len(io.REGIONS))
    counts_dir = root / "counts" / arm
    released_dir = root / "released" / arm
    counts_dir.mkdir(parents=True)
    released_dir.mkdir(parents=True)
    sizes = {"up": [1] * 14 + [2] * 3, "dn": [1] * 12 + [2] * 4, "bg": [1] * 150 + [2] * 30 + [3] * 6}
    offsets = {"up": 0, "dn": 10 ** 5, "bg": 10 ** 6}
    roots = {"up": {}, "dn": {}}
    for motif in MOTIFS:
        rate = float(rng.uniform(0.15, 0.6))
        store = {"region_of_position": region, "position": position}
        for group in io.GROUPS:
            matrix, ids = cluster_block(rng, sizes[group], offsets[group], n_windows, rate)
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


def simulate(tmp_path):
    rng = np.random.default_rng(20260924)
    results = []
    for rep in range(N_REPLICATES):
        args = null_arm(tmp_path, "N{}".format(rep), rng)
        _, per_motif, condensed, _, report, _ = v2.run_arm(args, lambda message: None)
        results.append((per_motif, condensed, report))
    return results


def test_two_stage_null_controls_fdr_and_p_is_uniform(tmp_path):
    results = simulate(tmp_path)
    margin = 3 * math.sqrt(ALPHA * (1 - ALPHA) / N_REPLICATES)
    motif_any, rbp_any, promoted, fixed_p, stage2_rejections = 0, 0, 0, [], 0
    for per_motif, condensed, report in results:
        plotted = [r for r in per_motif if r["plot"]]
        motif_any += any(r["calibrated_q"] < ALPHA for r in plotted)
        stage2_rejections += sum(1 for r in plotted if r["calibrated_q"] < ALPHA and r["calib_stage_pooled"] == 2)
        rbp_any += any(r["rbp_calibrated_q_minp"] < ALPHA for r in condensed if r["plot"])
        promoted += report["unique_pairs_promoted_motif_level"]
        first = next(r for r in per_motif if r["motif_key"] == MOTIFS[0] and r["direction"] == "up"
                     and r["region"] == "UpstreamIntron")
        fixed_p.append(first["calibrated_p_pooled"])
    # The mixture is exercised: stage 2 ran and every reported p is at one of the two resolutions.
    assert promoted > 0
    perms = {r["calib_perms_used_pooled"] for per_motif, _, _ in results for r in per_motif}
    assert perms == {B1, B2}
    # Global null: FDR = P(any BH rejection) <= nominal + Monte Carlo margin, both families.
    assert motif_any / N_REPLICATES <= ALPHA + margin, (motif_any, N_REPLICATES)
    assert rbp_any / N_REPLICATES <= ALPHA + margin, (rbp_any, N_REPLICATES)
    # One fixed test per independent replicate: super-uniform and not distinguishable from uniform.
    fixed_p = np.asarray(fixed_p)
    assert kstest(fixed_p, "uniform").pvalue > 0.01
    for t in (0.01, 0.05, 0.1, 0.25, 0.5):
        assert np.mean(fixed_p <= t) <= t + 3 * math.sqrt(t * (1 - t) / N_REPLICATES), t


def test_one_inferential_contract_in_lessons_readme_legend_workbook_and_readout():
    lessons = (ROOT / "LESSONS.md").read_text(encoding="utf-8")
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    builder = (ROOT / "tools" / "build_region_lollipops_v4.py").read_text(encoding="utf-8")
    for name, text in (("LESSONS.md", lessons), ("README.md", readme), ("figure legend", builder)):
        assert CONTRACT in text.lower(), name
        assert "no bh guarantee" not in text.lower(), name
        assert "PRDS" in text, name
    report = {"motif_family_size": 726, "n_unique_motifs": 121, "n_motif_keys": 126, "rbp_family_size": 600,
              "motif_bh_divisor": 726, "rbp_bh_divisor": {"minp": 600, "maxz": 600, "meanz": 600},
              "duplicate_kmers": [], "stage1_permutations": 2000, "stage2_permutations": 100000,
              "refine_threshold": 0.005, "seed": 149}
    workbook = " ".join(v for _, v in v2.readme_rows(report, "alias.tsv", None)).lower()
    assert CONTRACT in workbook and "calib_perms_used" in workbook
    assert "1/2,001" in workbook and "1/100,001" in workbook
    assert "no bh guarantee" not in workbook
