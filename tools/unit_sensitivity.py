"""Row-vs-exon unit sensitivity for the released rMAPS3 rank-sum layer (SENSITIVITY, 2026-09-21).

Three treatments of duplicate target exons, all on the same released count archives and the same
Westfall-Young two-stage calibration (seed 149 by default):

A  one row per rMATS SE event                  (the v1 row-unit summary; read, not recomputed)
B  one row per target exon, representative = max |IncLevelDifference|
C  rows kept, permutation labels assigned per target-exon cluster, size-matched

Treatment C was adopted as the reportable null on 2026-09-22 and lives on in
tools/calibrate_ranksum_v2.py; this tool reproduces the A/B/C comparison that motivated it.
Ported from the 2026-09-21 lab script without changing any arithmetic or draw stream; every
lab path became a required argument. Human / GRCh38 (hg38) / GENCODE v49 unless the inputs say
otherwise; SE events only.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
import os
import platform
import random
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import scipy
from scipy import sparse
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import calibrate_ranksum as calib        # noqa: E402
import rmaps_calib_v2_lib as lib         # noqa: E402
import rmaps_countdist_io as io          # noqa: E402

SEED = 149
STAGE1_PERMS = 2000
STAGE2_PERMS = 100000
REFINE_THRESHOLD = 0.005
CHUNK = 500
CONTROL_PANELS = (("Upstream Intron", "up"), ("Downstream Intron", "dn"))
EVENT_KEY_COLUMNS = ("chr", "strand", "exonStart_0base", "exonEnd",
                     "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
ADEQUATE_LABEL = "adequate"
TOOLS_DIR = Path(__file__).resolve().parent

target_keys = lib.target_keys
cluster_rows = lib.cluster_rows
ClusterDrawer = lib.ClusterDrawer


def settings(**overrides):
    """Run settings; paths are required inputs, statistics default to the locked lab values."""
    cfg = SimpleNamespace(counts_root=None, released_root=None, summary_a_root=None,
                          alias_table=None, out_root=None, event_dirs={}, underpowered=(),
                          positive_control=None, control_arms=(), seed=SEED,
                          stage1_perms=STAGE1_PERMS, stage2_perms=STAGE2_PERMS,
                          refine_threshold=REFINE_THRESHOLD, chunk=CHUNK)
    for key, value in overrides.items():
        if not hasattr(cfg, key):
            raise TypeError("unknown setting " + key)
        setattr(cfg, key, value)
    return cfg


def power_label(arm: str, underpowered, n_included: int, n_skipped: int) -> str:
    """The v1 power label; the underpowered arms are a named decision passed by the caller."""
    if arm in set(underpowered):
        return ("underpowered (n changed included/skipped = {}/{}) - landscape only, "
                "not for RBP ranking".format(n_included, n_skipped))
    return ADEQUATE_LABEL


def stage_rng(cfg, stage: int, direction: str):
    return lib.stage_rng(cfg.seed, stage, direction)


class SubsetModel:
    """calibrate_ranksum.MotifModel (row unit) restricted to a subset of rows.

    The arithmetic is identical to MotifModel; tests/test_unit_sensitivity.py proves that a
    full-row instance reproduces MotifModel attribute for attribute.
    """

    def __init__(self, npz_path: Path, direction: str, fg_rows=None, bg_rows=None):
        with np.load(npz_path, allow_pickle=False) as data:
            fg = io.group_csr(data, direction)
            bg = io.group_csr(data, "bg")
            self.region_index = np.asarray(data["region_of_position"])
            self.position = np.asarray(data["position"])
        if fg_rows is not None:
            fg = fg[np.asarray(fg_rows, dtype=np.int64)]
        if bg_rows is not None:
            bg = bg[np.asarray(bg_rows, dtype=np.int64)]
        self.n1 = int(fg.shape[0])
        self.n0 = int(bg.shape[0])
        self.n_total = self.n1 + self.n0
        self.n_windows = int(fg.shape[1])
        pooled = sparse.vstack([fg, bg], format="csr")
        kmax = int(pooled.data.max(initial=0))
        tied = io.csr_histograms(pooled, kmax).astype(np.float64)
        cumulative_below = np.cumsum(tied, axis=0) - tied
        midrank = cumulative_below + (tied + 1.0) / 2.0
        self.mu = self.n1 * self.n0 / 2.0
        tie_term = (tied ** 3 - tied).sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            variance = self.n1 * self.n0 / 12.0 * (
                (self.n_total + 1.0) - tie_term / (self.n_total * (self.n_total - 1.0)))
        self.sigma = np.sqrt(np.where(variance > 0, variance, np.nan))
        self.valid = np.isfinite(self.sigma)
        self.baseline = self.n1 * midrank[0]
        increment = midrank[pooled.data.astype(np.int64),
                            pooled.indices.astype(np.int64)] - midrank[0, pooled.indices]
        self.delta = sparse.csr_matrix(
            (increment.astype(np.float64), pooled.indices, pooled.indptr), shape=pooled.shape)
        self.observed_z = self.z(self.baseline + np.asarray(
            self.delta[: self.n1].sum(axis=0)).ravel())
        self.fg_carrying = np.asarray((fg > 0).sum(axis=0)).ravel()
        self.bg_carrying = np.asarray((bg > 0).sum(axis=0)).ravel()
        self.fg_total = np.asarray(fg.sum(axis=0)).ravel()
        self.bg_total = np.asarray(bg.sum(axis=0)).ravel()

    z = calib.MotifModel.z
    window_masks = calib.MotifModel.window_masks
    permutation_z = calib.MotifModel.permutation_z


def exon_ids_for_arm(cfg, arm: str, verify_all_motifs: bool = False):
    """8-field exon identifiers per group, in the row order of the count archives."""
    try:
        ids, motifs = lib.exon_ids_for_counts_dir(Path(cfg.counts_root) / arm, verify_all_motifs)
    except ValueError as exc:
        raise SystemExit("{}: {}".format(arm, exc))
    return ids, motifs, (len(motifs) if verify_all_motifs else 1)


def read_event_table(path: Path):
    """8-column key and |IncLevelDifference| per row, in file order."""
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {c: header.index(c) for c in EVENT_KEY_COLUMNS}
        dpsi_column = header.index("IncLevelDifference")
        keys, dpsi = [], []
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            keys.append(":".join(fields[index[c]] for c in EVENT_KEY_COLUMNS))
            try:
                dpsi.append(abs(float(fields[dpsi_column])))
            except ValueError:
                dpsi.append(math.nan)
    return np.asarray(keys, dtype="U"), np.asarray(dpsi, dtype=np.float64)


def strip_chr(keys: np.ndarray) -> np.ndarray:
    """Chromosome-agnostic key: engine-ready coord files may have renamed chr prefixes."""
    return np.asarray([k[3:] if k.startswith("chr") else k for k in keys], dtype="U")


def dpsi_for_group(cfg, arm: str, group: str, ids: np.ndarray):
    """|IncLevelDifference| aligned to the archive row order, or NaN if unrecoverable."""
    if arm not in cfg.event_dirs:
        return np.full(ids.shape[0], math.nan), "missing:no --event-dir for " + arm
    path = Path(cfg.event_dirs[arm]) / ("events_" + group + ".tsv")
    if not path.exists():
        return np.full(ids.shape[0], math.nan), "missing:" + str(path)
    keys, dpsi = read_event_table(path)
    if keys.shape[0] != ids.shape[0]:
        return np.full(ids.shape[0], math.nan), "row_count_mismatch:" + str(path)
    if np.array_equal(keys, ids):
        return dpsi, str(path)
    renamed = int((keys != ids).sum())
    if np.array_equal(strip_chr(keys), strip_chr(ids)):
        return dpsi, ("{} (chromosome-prefix normalised match, {} rows renamed by "
                      "the engine-ready prep)".format(path, renamed))
    return np.full(ids.shape[0], math.nan), "key_sequence_mismatch:" + str(path)


def dedup_index(keys: np.ndarray, dpsi: np.ndarray) -> np.ndarray:
    """One row per target exon: largest |dPSI|, ties and NaN broken by file order."""
    best = {}
    score = np.where(np.isfinite(dpsi), dpsi, -np.inf)
    for row, key in enumerate(keys):
        current = best.get(key)
        if current is None or score[row] > score[current]:
            best[key] = row
    return np.sort(np.asarray(sorted(best.values()), dtype=np.int64))


def make_drawer(fg_keys, bg_keys):
    try:
        return ClusterDrawer(fg_keys, bg_keys)
    except ValueError as exc:
        raise SystemExit(str(exc))


CONDENSED_COLUMNS = [
    "arm", "treatment", "power_label", "n_changed_included", "n_changed_skipped",
    "RBP", "rbp_table_names", "direction", "direction_label", "pooled_region", "plot",
    "selected_motif_key", "selected_native_region", "native_ranksum_p", "calibrated_p",
    "calibrated_q", "calib_perms_used", "calib_stage", "native_argmin_position",
    "n_fg_exons", "n_bg_exons", "n_fg_carrying", "n_bg_carrying", "fg_proportion",
    "bg_proportion", "enrichment_ratio", "fg_mean_count", "bg_mean_count", "count_ratio",
    "n_motifs_total", "n_motifs_calib_q_lt_0.05",
]


def run_treatment(cfg, arm: str, treatment: str, emit, dedup_groups=io.GROUPS, label=None):
    counts_dir = Path(cfg.counts_root) / arm
    out_dir = Path(cfg.out_root) / arm
    out_dir.mkdir(parents=True, exist_ok=True)
    alias = calib.load_alias(Path(cfg.alias_table))
    ids, motifs, _ = exon_ids_for_arm(cfg, arm)
    keys = {g: target_keys(ids[g]) for g in io.GROUPS}

    label = label or treatment
    keep, drawers, dpsi_source = None, None, {}
    if treatment == "B":
        if arm not in cfg.event_dirs and any(g in dedup_groups for g in io.GROUPS):
            raise SystemExit("treatment B needs --event-dir {}=<dir with events_<group>.tsv>: the "
                             "representative row is chosen by |IncLevelDifference|".format(arm))
        keep = {}
        for g in io.GROUPS:
            if g not in dedup_groups:
                keep[g] = np.arange(ids[g].shape[0], dtype=np.int64)
                dpsi_source[g] = "not deduplicated"
                emit("  {} {}: {} rows kept as-is (group excluded from dedup)".format(
                    arm, g, ids[g].shape[0]))
                continue
            dpsi, source = dpsi_for_group(cfg, arm, g, ids[g])
            dpsi_source[g] = source
            keep[g] = dedup_index(keys[g], dpsi)
            emit("  {} {}: {} rows -> {} exons (dPSI source {})".format(
                arm, g, ids[g].shape[0], keep[g].shape[0], source))
    else:
        drawers = {d: make_drawer(keys[d], keys["bg"]) for d in ("up", "dn")}
        for d in ("up", "dn"):
            emit("  {} {}: {} rows in {} clusters; pooled pool {} clusters; sizes {}".format(
                arm, d, drawers[d].n1, drawers[d].n_clusters_fg,
                drawers[d].n_clusters_pooled, drawers[d].need))

    roots = {d: io.read_root_table(Path(cfg.released_root) / arm / ("pVal." + d + ".vs.bg.RNAmap.txt"))
             for d in ("up", "dn")}

    stage1, meta, native, selection_cache = {}, {}, {}, {}
    started = time.perf_counter()
    for index, motif in enumerate(motifs, 1):
        table_name = motif.split(".", 1)[0]
        rbp = alias.get(table_name, table_name)
        for direction in ("up", "dn"):
            if treatment == "B":
                model = SubsetModel(counts_dir / (motif + ".counts.npz"), direction,
                                    keep[direction], keep["bg"])
            else:
                model = calib.MotifModel(counts_dir / (motif + ".counts.npz"), direction)
            masks = model.window_masks()
            cache_key = (direction, model.n_total, model.n1)
            if cache_key not in selection_cache:
                rng = stage_rng(cfg, 1, direction)
                selection_cache[cache_key] = (
                    calib.draw_selection(rng, model.n_total, model.n1, cfg.stage1_perms)
                    if treatment == "B" else drawers[direction].draw(rng, cfg.stage1_perms))
            stage1[(motif, direction)] = calib.calibrate(
                model, masks, selection_cache[cache_key], cfg.chunk)
            native_p = io.p_from_z(model.observed_z)
            argmin, effect, region_p = {}, {}, {}
            for region_id, region in enumerate(io.REGIONS):
                window_ids = np.flatnonzero(model.region_index == region_id)
                local = native_p[window_ids]
                best = int(window_ids[int(np.argmin(local))])
                argmin[region] = best
                effect[region] = calib.effect_at(model, best)
                region_p[region] = float(local.min())
                if treatment == "C":
                    released = roots[direction][motif][region]
                    if released > 0 and region_p[region] > 0:
                        if abs(math.log10(region_p[region]) - math.log10(released)) > 1e-6:
                            raise SystemExit("treatment C native p differs from the released "
                                             "table: {} {} {}".format(motif, direction, region))
            for pool, members in io.POOL_TO_REGIONS.items():
                best_region = min(members, key=lambda r: region_p[r])
                argmin[pool] = argmin[best_region]
                effect[pool] = effect[best_region]
                region_p[pool] = region_p[best_region]
            native[(motif, direction)] = region_p
            meta[(motif, direction)] = {"n1": model.n1, "n0": model.n0, "rbp": rbp,
                                        "table_name": table_name, "argmin": argmin,
                                        "effect": effect}
            del model
        if index % 25 == 0 or index == len(motifs):
            emit("  stage 1 [{} {}]: {}/{} motifs, {:.1f}s".format(
                arm, treatment, index, len(motifs), time.perf_counter() - started))
    stage1_seconds = time.perf_counter() - started

    qualifying = sorted(
        ((min(v["calibrated_p"] for v in stage1[k].values()
              if math.isfinite(v["calibrated_p"])), k)
         for k in stage1
         if any(math.isfinite(v["calibrated_p"]) and v["calibrated_p"] <= cfg.refine_threshold
                for v in stage1[k].values())),
        key=lambda x: (x[0], x[1]))
    selected = [k for _, k in qualifying]
    emit("stage 2 [{} {}]: {} of {} motif-direction pairs promoted".format(
        arm, treatment, len(selected), len(stage1)))

    stage = {k: 1 for k in stage1}
    final = {k: dict(v) for k, v in stage1.items()}
    stage2_started = time.perf_counter()
    for done, (motif, direction) in enumerate(selected, 1):
        if treatment == "B":
            model = SubsetModel(counts_dir / (motif + ".counts.npz"), direction,
                                keep[direction], keep["bg"])
        else:
            model = calib.MotifModel(counts_dir / (motif + ".counts.npz"), direction)
        masks = model.window_masks()
        rng = stage_rng(cfg, 2, direction)
        observed = calib.statistic_maxima(model.observed_z, masks)
        exceed = {name: 0 for name in masks}
        remaining = cfg.stage2_perms
        while remaining > 0:
            batch = min(cfg.chunk, remaining)
            selection = (calib.draw_selection(rng, model.n_total, model.n1, batch)
                         if treatment == "B" else drawers[direction].draw(rng, batch))
            maxima = calib.statistic_maxima(model.permutation_z(selection), masks)
            for name in masks:
                if math.isfinite(observed[name]):
                    values = maxima[name]
                    exceed[name] += int(np.sum(values[np.isfinite(values)] >= observed[name]))
            remaining -= batch
        for name in masks:
            if not math.isfinite(observed[name]):
                continue
            if abs(observed[name] - final[(motif, direction)][name]["observed_max_z"]) > 1e-9:
                raise SystemExit("observed statistic changed during refinement")
            final[(motif, direction)][name] = dict(final[(motif, direction)][name])
            final[(motif, direction)][name]["calibrated_p"] = (
                (1 + exceed[name]) / (1 + cfg.stage2_perms))
            final[(motif, direction)][name]["permutations"] = cfg.stage2_perms
        stage[(motif, direction)] = 2
        del model
        if done % 10 == 0 or done == len(selected):
            emit("  stage 2 [{} {}]: {}/{} pairs, {:.1f}s".format(
                arm, treatment, done, len(selected), time.perf_counter() - stage2_started))
    stage2_seconds = time.perf_counter() - stage2_started

    bh_keys = [(m, d, p) for m in motifs for d in ("up", "dn") for p in io.PLOT_POOLS]
    q_lookup = {k: float(q) for k, q in zip(bh_keys, calib.bh_adjust(
        [final[(m, d)][p]["calibrated_p"] for m, d, p in bh_keys]))}

    n_included = meta[(motifs[0], "up")]["n1"]
    n_skipped = meta[(motifs[0], "dn")]["n1"]
    power_columns = {"power_label": power_label(arm, cfg.underpowered, n_included, n_skipped),
                     "n_changed_included": n_included, "n_changed_skipped": n_skipped}

    condensed = []
    rbps = sorted({alias.get(m.split(".", 1)[0], m.split(".", 1)[0]) for m in motifs})
    for rbp in rbps:
        rbp_motifs = [m for m in motifs
                      if alias.get(m.split(".", 1)[0], m.split(".", 1)[0]) == rbp]
        table_names = sorted({m.split(".", 1)[0] for m in rbp_motifs})
        for direction in ("up", "dn"):
            for pool in io.POOL_TO_REGIONS:
                candidates = [(final[(m, direction)][pool]["calibrated_p"],
                               native[(m, direction)][pool], m) for m in rbp_motifs]
                finite = [c for c in candidates if math.isfinite(c[0])]
                motif = (min(finite, key=lambda c: (c[0], c[1], c[2]))[2] if finite
                         else rbp_motifs[0])
                best_region = min(io.POOL_TO_REGIONS[pool],
                                  key=lambda r: native[(motif, direction)][r])
                record = {
                    "arm": arm, "treatment": label, "RBP": rbp,
                    "rbp_table_names": ",".join(table_names), "direction": direction,
                    "direction_label": io.DIRECTION_LABEL[direction], "pooled_region": pool,
                    "plot": pool in io.PLOT_POOLS, "selected_motif_key": motif,
                    "selected_native_region": best_region,
                    "native_ranksum_p": native[(motif, direction)][pool],
                    "calibrated_p": final[(motif, direction)][pool]["calibrated_p"],
                    "calibrated_q": q_lookup.get((motif, direction, pool), math.nan),
                    "calib_perms_used": final[(motif, direction)][pool]["permutations"],
                    "calib_stage": stage[(motif, direction)],
                    "native_argmin_position": int(meta[(motif, direction)]["argmin"][pool]),
                    "n_fg_exons": meta[(motif, direction)]["n1"],
                    "n_bg_exons": meta[(motif, direction)]["n0"],
                    "n_motifs_total": len(rbp_motifs),
                    "n_motifs_calib_q_lt_0.05": sum(
                        1 for m in rbp_motifs
                        if math.isfinite(q_lookup.get((m, direction, pool), math.nan))
                        and q_lookup[(m, direction, pool)] < 0.05),
                }
                record.update(meta[(motif, direction)]["effect"][pool])
                record.update(power_columns)
                condensed.append(record)

    path = out_dir / ("treatment_" + label + "_condensed_per_rbp.tsv")
    calib.write_tsv(path, condensed, CONDENSED_COLUMNS)
    report = {"arm": arm, "treatment": label, "dedup_groups": list(dedup_groups),
              "seed": cfg.seed,
              "stage1_permutations": cfg.stage1_perms, "stage2_permutations": cfg.stage2_perms,
              "refine_threshold": cfg.refine_threshold, "pairs_total": len(stage1),
              "pairs_refined": len(selected), "bh_family_size": len(bh_keys),
              "n_changed_included": n_included, "n_changed_skipped": n_skipped,
              "stage1_wall_seconds": stage1_seconds, "stage2_wall_seconds": stage2_seconds,
              "dpsi_source": dpsi_source,
              "cluster_sizes": ({d: drawers[d].need for d in drawers} if drawers else None)}
    (out_dir / ("treatment_" + label + "_report.json")).write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    emit("wrote {} ({} rows); stage1 {:.1f}s stage2 {:.1f}s".format(
        path, len(condensed), stage1_seconds, stage2_seconds))
    return path


def duplication_table(cfg, arms, emit):
    rows = []
    for arm in arms:
        ids, motifs, checked = exon_ids_for_arm(cfg, arm, verify_all_motifs=True)
        keys = {g: target_keys(ids[g]) for g in io.GROUPS}
        sets = {g: set(keys[g].tolist()) for g in io.GROUPS}
        for g in io.GROUPS:
            clusters = cluster_rows(keys[g])
            sizes = np.asarray([c.shape[0] for c in clusters.values()])
            dpsi, source = dpsi_for_group(cfg, arm, g, ids[g])
            rows.append({
                "arm": arm, "group": g,
                "group_label": {"up": "changed INCLUDED", "dn": "changed SKIPPED",
                                "bg": "background"}[g],
                "n_rows": int(keys[g].shape[0]), "n_target_exons": int(sizes.shape[0]),
                "duplicate_row_fraction": float(1.0 - sizes.shape[0] / keys[g].shape[0]),
                "n_exons_with_duplicates": int((sizes > 1).sum()),
                "max_cluster_size": int(sizes.max()),
                "n_overlap_with_up": int(len(sets[g] & sets["up"])) if g != "up" else 0,
                "n_overlap_with_dn": int(len(sets[g] & sets["dn"])) if g != "dn" else 0,
                "n_overlap_with_bg": int(len(sets[g] & sets["bg"])) if g != "bg" else 0,
                "dpsi_recovered_rows": int(np.isfinite(dpsi).sum()),
                "dpsi_source": source,
                "motifs_with_identical_exon_axis": checked,
            })
        emit("  duplication: {} done ({} motifs share one exon axis)".format(arm, checked))
    columns = ["arm", "group", "group_label", "n_rows", "n_target_exons",
               "duplicate_row_fraction", "n_exons_with_duplicates", "max_cluster_size",
               "n_overlap_with_up", "n_overlap_with_dn", "n_overlap_with_bg",
               "dpsi_recovered_rows", "dpsi_source", "motifs_with_identical_exon_axis"]
    path = Path(cfg.out_root) / "duplication_table.tsv"
    calib.write_tsv(path, rows, columns)
    emit("wrote " + str(path))
    return rows


read_tsv = lib.read_tsv
as_float = lib.as_float


def rho(left, right):
    """Spearman over the RBPs where both treatments give a finite value."""
    pairs = [(a, b) for a, b in zip(left, right) if math.isfinite(a) and math.isfinite(b)]
    if len(pairs) < 3:
        return math.nan
    return float(spearmanr([p[0] for p in pairs], [p[1] for p in pairs]).statistic)


def panel_rank(records):
    """Rank within a panel by (q, calibrated p, native p, RBP); NA sorts last."""
    order = sorted(records, key=lambda r: (
        r["q"] if math.isfinite(r["q"]) else math.inf,
        r["p"] if math.isfinite(r["p"]) else math.inf,
        r["native"] if math.isfinite(r["native"]) else math.inf, r["RBP"]))
    return {r["RBP"]: i for i, r in enumerate(order, 1)}


def build_comparison(cfg, arms, emit):
    out_root = Path(cfg.out_root)
    data = {}
    for arm in arms:
        a_rows = read_tsv(Path(cfg.summary_a_root) / arm / "condensed_per_rbp.tsv")
        b_rows = read_tsv(out_root / arm / "treatment_B_condensed_per_rbp.tsv")
        c_rows = read_tsv(out_root / arm / "treatment_C_condensed_per_rbp.tsv")
        for label, rows in (("A", a_rows), ("B", b_rows), ("C", c_rows)):
            for row in rows:
                if row["pooled_region"] not in io.PLOT_POOLS:
                    continue
                data[(arm, row["direction"], row["pooled_region"], label, row["RBP"])] = {
                    "RBP": row["RBP"], "motif": row["selected_motif_key"],
                    "p": as_float(row["calibrated_p"]), "q": as_float(row["calibrated_q"]),
                    "native": as_float(row["native_ranksum_p"]),
                    "n_fg": row["n_fg_exons"], "stage": row["calib_stage"],
                }
    panels = sorted({(a, d, p) for (a, d, p, _, _) in data})
    ranks = {}
    for arm, direction, pool in panels:
        for label in ("A", "B", "C"):
            records = [v for (a, d, p, t, _), v in data.items()
                       if (a, d, p, t) == (arm, direction, pool, label)]
            ranks[(arm, direction, pool, label)] = panel_rank(records)

    comparison, panel_rows, lost_rows = [], [], []
    for arm, direction, pool in panels:
        rbps = sorted({r for (a, d, p, t, r) in data if (a, d, p) == (arm, direction, pool)})
        vectors = {label: [] for label in ("A", "B", "C")}
        native_vectors = {label: [] for label in ("A", "B", "C")}
        for rbp in rbps:
            row = {"arm": arm, "direction": direction,
                   "direction_label": io.DIRECTION_LABEL[direction],
                   "pooled_region": pool, "RBP": rbp}
            for label in ("A", "B", "C"):
                entry = data.get((arm, direction, pool, label, rbp))
                rank = ranks[(arm, direction, pool, label)].get(rbp)
                row[label + "_selected_motif"] = entry["motif"] if entry else "NA"
                row[label + "_native_ranksum_p"] = entry["native"] if entry else math.nan
                row[label + "_calibrated_p"] = entry["p"] if entry else math.nan
                row[label + "_calibrated_q"] = entry["q"] if entry else math.nan
                row[label + "_rank"] = rank if rank else None
                row[label + "_n_fg_exons"] = entry["n_fg"] if entry else "NA"
                vectors[label].append(entry["p"] if entry else math.nan)
                native_vectors[label].append(entry["native"] if entry else math.nan)
            for label in ("B", "C"):
                row["delta_rank_A_to_" + label] = (
                    row[label + "_rank"] - row["A_rank"]
                    if row["A_rank"] and row[label + "_rank"] else None)
            comparison.append(row)
            qa = row["A_calibrated_q"]
            for label in ("B", "C"):
                qx = row[label + "_calibrated_q"]
                if math.isfinite(qa) and qa < 0.05 and not (math.isfinite(qx) and qx < 0.05):
                    lost_rows.append({"arm": arm, "direction": direction,
                                      "direction_label": io.DIRECTION_LABEL[direction],
                                      "pooled_region": pool, "RBP": rbp,
                                      "comparison": "A_to_" + label,
                                      "A_calibrated_q": qa, "other_calibrated_q": qx,
                                      "A_calibrated_p": row["A_calibrated_p"],
                                      "other_calibrated_p": row[label + "_calibrated_p"]})
        summary = {"arm": arm, "direction": direction,
                   "direction_label": io.DIRECTION_LABEL[direction], "pooled_region": pool,
                   "n_RBPs": len(rbps)}
        for label in ("B", "C"):
            summary["spearman_calibrated_p_A_vs_" + label] = rho(
                vectors["A"], vectors[label])
            summary["spearman_native_p_A_vs_" + label] = rho(
                native_vectors["A"], native_vectors[label])
        for label in ("A", "B", "C"):
            summary["n_RBP_q_lt_0.05_" + label] = sum(
                1 for v in [data.get((arm, direction, pool, label, r)) for r in rbps]
                if v and math.isfinite(v["q"]) and v["q"] < 0.05)
        panel_rows.append(summary)
        emit("  panel {} {} {}: rho(A,B)={:.3f} rho(A,C)={:.3f} q<0.05 {}/{}/{}".format(
            arm, direction, pool, summary["spearman_calibrated_p_A_vs_B"],
            summary["spearman_calibrated_p_A_vs_C"], summary["n_RBP_q_lt_0.05_A"],
            summary["n_RBP_q_lt_0.05_B"], summary["n_RBP_q_lt_0.05_C"]))

    comparison_columns = ["arm", "direction", "direction_label", "pooled_region", "RBP"]
    for label in ("A", "B", "C"):
        comparison_columns += [label + s for s in
                               ("_selected_motif", "_native_ranksum_p", "_calibrated_p",
                                "_calibrated_q", "_rank", "_n_fg_exons")]
    comparison_columns += ["delta_rank_A_to_B", "delta_rank_A_to_C"]
    calib.write_tsv(out_root / "comparison_ABC.tsv", comparison, comparison_columns)

    panel_columns = ["arm", "direction", "direction_label", "pooled_region", "n_RBPs",
                     "spearman_calibrated_p_A_vs_B", "spearman_calibrated_p_A_vs_C",
                     "spearman_native_p_A_vs_B", "spearman_native_p_A_vs_C",
                     "n_RBP_q_lt_0.05_A", "n_RBP_q_lt_0.05_B", "n_RBP_q_lt_0.05_C"]
    lost_columns = ["arm", "direction", "direction_label", "pooled_region", "RBP",
                    "comparison", "A_calibrated_p", "other_calibrated_p",
                    "A_calibrated_q", "other_calibrated_q"]
    calib.write_tsv(out_root / "panel_summary_ABC.tsv", panel_rows, panel_columns)
    calib.write_tsv(out_root / "lost_significance_ABC.tsv", lost_rows, lost_columns)

    control, control_columns = [], []
    rbp = cfg.positive_control
    if rbp and cfg.control_arms:
        for arm in cfg.control_arms:
            for pool, direction in CONTROL_PANELS:
                row = {"arm": arm, "pooled_region": pool, "direction": direction,
                       "direction_label": io.DIRECTION_LABEL[direction],
                       "n_RBPs": len(ranks[(arm, direction, pool, "A")])}
                for label in ("A", "B", "C"):
                    entry = data.get((arm, direction, pool, label, rbp))
                    row[rbp + "_" + label + "_rank"] = ranks[(arm, direction, pool, label)].get(rbp)
                    row[rbp + "_" + label + "_selected_motif"] = entry["motif"] if entry else "NA"
                    row[rbp + "_" + label + "_native_ranksum_p"] = (
                        entry["native"] if entry else math.nan)
                    row[rbp + "_" + label + "_calibrated_p"] = entry["p"] if entry else math.nan
                    row[rbp + "_" + label + "_calibrated_q"] = entry["q"] if entry else math.nan
                    row[rbp + "_" + label + "_n_fg_exons"] = entry["n_fg"] if entry else "NA"
                control.append(row)
        control_columns = ["arm", "pooled_region", "direction", "direction_label", "n_RBPs"]
        for label in ("A", "B", "C"):
            control_columns += [rbp + "_" + label + s for s in
                                ("_rank", "_selected_motif", "_native_ranksum_p",
                                 "_calibrated_p", "_calibrated_q", "_n_fg_exons")]
        calib.write_tsv(out_root / (rbp.lower() + "_control_ABC.tsv"), control, control_columns)

    duplication_path = out_root / "duplication_table.tsv"
    duplication = read_tsv(duplication_path) if duplication_path.is_file() else []
    reports = {}
    for arm in arms:
        for label in ("B", "C"):
            path = out_root / arm / ("treatment_" + label + "_report.json")
            reports[(arm, label)] = json.loads(path.read_text(encoding="utf-8"))
    stage2 = sorted({r["stage2_permutations"] for r in reports.values()})
    families = sorted({r["bh_family_size"] for r in reports.values()})
    readme = [
        ("Question", "Does the rank-sum layer's RBP order and significance depend on whether "
                     "the statistical unit is the rMATS SE event row or the target exon?"),
        ("Treatment A", "One row per rMATS SE event, the released engine's unit. Read from "
                        "{}/<ARM>/condensed_per_rbp.tsv (row-unit summary of "
                        "tools/calibrate_ranksum.py), never recomputed here.".format(
                            cfg.summary_a_root)),
        ("Treatment B", "One row per target exon (chr, strand, exonStart, exonEnd). The kept "
                        "row is the one with the largest |IncLevelDifference|; ties and "
                        "unrecoverable dPSI fall back to first-in-file. Count matrices are "
                        "row-subset, so this is what the released engine would emit on "
                        "deduplicated inputs."),
        ("Treatment C", "Rows kept exactly as in A; only the permutation null changes. Labels "
                        "are assigned per target-exon cluster, all rows of an exon moving "
                        "together, with the observed foreground's exact cluster-size "
                        "composition drawn from the pooled cluster pool."),
        ("Flanking-exon caveat", "Duplicate rows of one target exon share identical counts in "
                                 "UpstreamIntron, TargetExon_5prime, TargetExon-3prime and "
                                 "DownstreamIntron, but differ in UpstreamExonIntron and "
                                 "DownstreamExonIntron, which are members of the pooled "
                                 "Upstream Intron and Downstream Intron regions. Under B the "
                                 "representative choice therefore does affect those two "
                                 "sub-regions. The Flanking Exon pool is excluded entirely: it "
                                 "is not plotted and carries no BH q in any treatment."),
        ("calibrated_p", "Westfall-Young min-P permutation p of the released rank-sum regional "
                         "minimum. Stage 1 = {} permutations; stage 2 = {} for pairs reaching "
                         "stage-1 p <= {}. Seed {} in all treatments.".format(
                             cfg.stage1_perms, ",".join(map(str, stage2)),
                             cfg.refine_threshold, cfg.seed)),
        ("calibrated_q", "Benjamini-Hochberg over motifs x 2 directions x 3 pooled regions "
                         "within one arm and one treatment (family size {} in the "
                         "treatment_<T>_report.json files). No correction across arms or "
                         "across treatments.".format(",".join(map(str, families)))),
        ("rank", "Rank within one panel (arm x direction x pooled region) over RBPs, ordered by "
                 "calibrated q, then calibrated p, then native p, then RBP name. RBPs at the "
                 "permutation floor 1/(stage-2 permutations + 1) tie, so rank differences at "
                 "the top are resolution limits, not reorderings."),
        ("Spearman", "Rank correlation of the per-RBP calibrated p across a panel. The native-p "
                     "correlation for A vs C is 1.0 by construction: C changes only the null, "
                     "never the observed statistic."),
        ("Species / assembly", "Human / GRCh38 (hg38) / GENCODE v49 unless the inputs say "
                               "otherwise; SE (skipped exon) events only."),
    ]
    sheets = {
        "README": (["field", "description"],
                   [{"field": k, "description": v} for k, v in readme]),
        "panel_summary": (panel_columns, panel_rows),
        "comparison_ABC": (comparison_columns, comparison),
    }
    if control:
        sheets[rbp.lower() + "_control"] = (control_columns, control)
    sheets["lost_significance"] = (lost_columns, lost_rows)
    if duplication:
        sheets["duplication_table"] = (list(duplication[0].keys()), duplication)
    calib.write_workbook(out_root / "comparison_ABC.xlsx", sheets)
    emit("wrote comparison_ABC.tsv ({} rows), panel_summary_ABC.tsv ({}), "
         "lost_significance_ABC.tsv ({}), control table ({}), comparison_ABC.xlsx"
         .format(len(comparison), len(panel_rows), len(lost_rows), len(control)))


def write_versions(cfg):
    with open(Path(cfg.out_root) / "versions.txt", "w", encoding="utf-8") as handle:
        handle.write("python\t{}\n".format(sys.version.split()[0]))
        handle.write("numpy\t{}\n".format(np.__version__))
        handle.write("scipy\t{}\n".format(scipy.__version__))
        handle.write("openpyxl\t{}\n".format(importlib.metadata.version("openpyxl")))
        handle.write("platform\t{}\n".format(platform.platform()))
        handle.write("script\t{}\n".format(Path(__file__).resolve()))
        handle.write("script_md5\t{}\n".format(io.md5(Path(__file__).resolve())))
        for name in ("calibrate_ranksum.py", "rmaps_calib_v2_lib.py", "rmaps_countdist_io.py"):
            handle.write("library_{}_md5\t{}\n".format(name[:-3], io.md5(TOOLS_DIR / name)))
        handle.write("seed\t{}\n".format(cfg.seed))
        handle.write("generated\t{}\n".format(time.strftime("%Y-%m-%dT%H:%M:%S")))
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            handle.write("{}\t{}\n".format(var, os.environ.get(var, "unset")))


def parse_event_dirs(values):
    out = {}
    for value in values or ():
        arm, sep, path = value.partition("=")
        if not sep or not arm or not path:
            raise SystemExit("--event-dir takes ARM=PATH, got " + value)
        out[arm] = path
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--mode", required=True,
                        choices=("duplication", "B", "C", "compare", "versions"))
    parser.add_argument("--arms", nargs="+", required=True)
    parser.add_argument("--out-root", required=True, help="output root; <ARM>/ subfolders")
    parser.add_argument("--counts-root", help="<ARM>/<motif>.counts.npz from tools/countdist_to_npz.py")
    parser.add_argument("--released-root", help="released run root holding <ARM>/pVal.*.vs.bg.RNAmap.txt")
    parser.add_argument("--alias-table", help="table_name -> HGNC symbol TSV")
    parser.add_argument("--summary-a-root",
                        help="compare: row-unit summary root (<ARM>/condensed_per_rbp.tsv)")
    parser.add_argument("--event-dir", action="append", default=[], metavar="ARM=PATH",
                        help="directory holding events_{up,dn,bg}.tsv with IncLevelDifference; "
                             "required for treatment B")
    parser.add_argument("--underpowered-arms", nargs="*", default=[])
    parser.add_argument("--positive-control", default=None, help="RBP for the control table")
    parser.add_argument("--control-arms", nargs="*", default=[])
    parser.add_argument("--dedup-groups", default="up,dn,bg",
                        help="groups deduplicated under treatment B")
    parser.add_argument("--label", default=None, help="output label; defaults to the mode")
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--permutations", type=int, default=STAGE1_PERMS)
    parser.add_argument("--refine-perms", type=int, default=STAGE2_PERMS)
    parser.add_argument("--refine-threshold", type=float, default=REFINE_THRESHOLD)
    args = parser.parse_args(argv)

    needs = {"duplication": ("counts_root",),
             "B": ("counts_root", "released_root", "alias_table"),
             "C": ("counts_root", "released_root", "alias_table"),
             "compare": ("summary_a_root",), "versions": ()}[args.mode]
    missing = ["--" + n.replace("_", "-") for n in needs if getattr(args, n) is None]
    if missing:
        parser.error("--mode {} requires {}".format(args.mode, ", ".join(missing)))
    cfg = settings(counts_root=args.counts_root, released_root=args.released_root,
                   summary_a_root=args.summary_a_root, alias_table=args.alias_table,
                   out_root=args.out_root, event_dirs=parse_event_dirs(args.event_dir),
                   underpowered=tuple(args.underpowered_arms),
                   positive_control=args.positive_control, control_arms=tuple(args.control_arms),
                   seed=args.seed, stage1_perms=args.permutations, stage2_perms=args.refine_perms,
                   refine_threshold=args.refine_threshold)

    random.seed(cfg.seed)
    np.random.seed(cfg.seed)
    Path(cfg.out_root).mkdir(parents=True, exist_ok=True)
    log = open(Path(cfg.out_root) / "command.log", "a", encoding="utf-8")

    def emit(message):
        print(message, flush=True)
        log.write(message + "\n")
        log.flush()

    started = time.perf_counter()
    emit("start={} cmd={}".format(time.strftime("%Y-%m-%dT%H:%M:%S"),
                                  " ".join([sys.executable] + sys.argv)))
    if args.mode == "duplication":
        duplication_table(cfg, args.arms, emit)
    elif args.mode in ("B", "C"):
        groups = tuple(g for g in args.dedup_groups.split(",") if g)
        for arm in args.arms:
            run_treatment(cfg, arm, args.mode, emit, dedup_groups=groups, label=args.label)
    elif args.mode == "compare":
        build_comparison(cfg, args.arms, emit)
    else:
        write_versions(cfg)
    emit("mode={} wall_seconds={:.1f} exit=0".format(args.mode, time.perf_counter() - started))
    log.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
