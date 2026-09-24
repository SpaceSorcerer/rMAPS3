"""Length-matched background sensitivity for the rMAPS3 rank-sum layer (SENSITIVITY, 2026-09-22).

Per arm and direction: bin target-exon length at the deciles of the changed set (unique target
exons), draw background target-exon CLUSTERS without replacement, stratified so the background
bin composition matches the changed set, then rerun the cluster-permutation calibration
(treatment C of tools/unit_sensitivity.py, the null of tools/calibrate_ranksum_v2.py) against
that matched background. --size-control draws the same number of clusters uniformly at random,
which separates background SIZE from exon LENGTH as the cause of any change.

A sensitivity, never the null. tools/build_region_lollipops_v4.py --length-root reads the
<ARM>/lengthmatched_report.json this tool writes. Ported from the 2026-09-22 lab script without
changing any arithmetic or draw stream; every lab path became a required argument.
Human / GRCh38 (hg38) / GENCODE v49 unless the inputs say otherwise; SE events only.
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

import numpy as np
import scipy
from scipy.stats import ks_2samp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import calibrate_ranksum as calib        # noqa: E402
import rmaps_countdist_io as io          # noqa: E402
import unit_sensitivity as unit          # noqa: E402

RESAMPLE_SEEDS = (20260922, 20260923)
SIZE_CONTROL_LABEL = "sizecontrol_random"
DECILES = np.arange(1, 10) / 10.0
QUANTS = (0.1, 0.25, 0.5, 0.75, 0.9)
TOOLS_DIR = Path(__file__).resolve().parent


def seed_labels(seeds):
    """First resample seed -> 'lengthmatched', the i-th further seed -> 'lengthmatched_seed<i+1>'."""
    return {s: ("lengthmatched" if i == 0 else "lengthmatched_seed{}".format(i + 1))
            for i, s in enumerate(seeds)}


def exon_lengths(ids: np.ndarray) -> np.ndarray:
    """exonEnd - exonStart(0-based) from the 8-field exon identifier."""
    return np.asarray([int(s.split(":")[3]) - int(s.split(":")[2]) for s in ids], dtype=np.int64)


def decile_edges(fg_keys: np.ndarray, fg_len: np.ndarray) -> np.ndarray:
    _, first = np.unique(fg_keys, return_index=True)
    return np.quantile(fg_len[first], DECILES)


def assign_bin(lengths: np.ndarray, edges: np.ndarray) -> np.ndarray:
    return np.searchsorted(edges, lengths, side="right")


def stratified_background(fg_keys, fg_len, bg_keys, bg_len, seed: int, direction: str,
                          size_control=False):
    """Rows of a background whose unique-exon length-bin composition matches the changed set.

    Whole target-exon clusters are drawn without replacement; the bin that is scarcest relative
    to its target share sets the matched size.
    """
    edges = decile_edges(fg_keys, fg_len)
    fg_u, fg_first = np.unique(fg_keys, return_index=True)
    fg_bins = assign_bin(fg_len[fg_first], edges)
    n_bins = edges.shape[0] + 1
    target = np.bincount(fg_bins, minlength=n_bins)

    clusters = unit.cluster_rows(bg_keys)
    c_keys = list(clusters.keys())
    c_len = np.asarray([bg_len[clusters[k][0]] for k in c_keys], dtype=np.int64)
    c_bins = assign_bin(c_len, edges)
    have = np.bincount(c_bins, minlength=n_bins)
    used = target > 0
    ratio = float(np.min(have[used] / target[used]))
    take = np.floor(ratio * target + 1e-9).astype(np.int64)
    take[~used] = 0

    rng = np.random.default_rng(np.random.SeedSequence([seed, io.DIRECTION_CODE[direction]]))
    chosen = []
    if size_control:
        n_bins = 0
        chosen = np.sort(rng.choice(len(c_keys), int(take.sum()), replace=False)).tolist()
    for b in range(n_bins):
        members = np.flatnonzero(c_bins == b)
        if take[b] == 0:
            continue
        pick = rng.choice(members.shape[0], int(take[b]), replace=False)
        chosen.extend(members[np.sort(pick)].tolist())
    rows = np.sort(np.concatenate([clusters[c_keys[i]] for i in chosen])).astype(np.int64)
    info = {"edges": edges.tolist(), "target_exons_per_bin": target.tolist(),
            "bg_exons_per_bin_before": have.tolist(), "bg_exons_per_bin_after": take.tolist(),
            "ratio_bg_to_changed": ratio, "n_bg_exons_before": len(c_keys),
            "n_bg_exons_after": int(take.sum()), "n_bg_rows_before": int(bg_keys.shape[0]),
            "n_bg_rows_after": int(rows.shape[0])}
    return rows, info


def quantiles(values: np.ndarray) -> list:
    return [float(v) for v in np.quantile(values, QUANTS)]


def unique_lengths(keys, lengths):
    _, first = np.unique(keys, return_index=True)
    return lengths[first]


def run_matched(cfg, arm: str, seed: int, label: str, emit, stage1_only=False,
                size_control=False):
    label = SIZE_CONTROL_LABEL if size_control else label
    counts_dir = Path(cfg.counts_root) / arm
    out_dir = Path(cfg.out_root) / arm
    out_dir.mkdir(parents=True, exist_ok=True)
    alias = calib.load_alias(Path(cfg.alias_table))
    ids, motifs, _ = unit.exon_ids_for_arm(cfg, arm)
    keys = {g: unit.target_keys(ids[g]) for g in io.GROUPS}
    lens = {g: exon_lengths(ids[g]) for g in io.GROUPS}

    bg_rows, drawers, match_info = {}, {}, {}
    for d in ("up", "dn"):
        bg_rows[d], info = stratified_background(keys[d], lens[d], keys["bg"], lens["bg"], seed, d,
                                                 size_control)
        fg_u = unique_lengths(keys[d], lens[d])
        bg_before = unique_lengths(keys["bg"], lens["bg"])
        bg_after = unique_lengths(keys["bg"][bg_rows[d]], lens["bg"][bg_rows[d]])
        info.update({
            "quantile_levels": list(QUANTS),
            "changed_length_quantiles": quantiles(fg_u),
            "bg_length_quantiles_before": quantiles(bg_before),
            "bg_length_quantiles_after": quantiles(bg_after),
            "ks_D_before": float(ks_2samp(fg_u, bg_before).statistic),
            "ks_D_after": float(ks_2samp(fg_u, bg_after).statistic),
        })
        match_info[d] = info
        drawers[d] = unit.make_drawer(keys[d], keys["bg"][bg_rows[d]])
        emit("  {} {} seed {}: bg exons {} -> {}, rows {} -> {}; median len changed {} bg {} -> {};"
             " KS D {:.3f} -> {:.3f}".format(
                 arm, d, seed, info["n_bg_exons_before"], info["n_bg_exons_after"],
                 info["n_bg_rows_before"], info["n_bg_rows_after"],
                 info["changed_length_quantiles"][2], info["bg_length_quantiles_before"][2],
                 info["bg_length_quantiles_after"][2], info["ks_D_before"], info["ks_D_after"]))

    def model_for(motif, direction):
        return unit.SubsetModel(counts_dir / (motif + ".counts.npz"), direction,
                                None, bg_rows[direction])

    stage1, meta, native, selection_cache = {}, {}, {}, {}
    started = time.perf_counter()
    for index, motif in enumerate(motifs, 1):
        table_name = motif.split(".", 1)[0]
        rbp = alias.get(table_name, table_name)
        for direction in ("up", "dn"):
            model = model_for(motif, direction)
            masks = model.window_masks()
            if direction not in selection_cache:
                selection_cache[direction] = drawers[direction].draw(
                    unit.stage_rng(cfg, 1, direction), cfg.stage1_perms)
            stage1[(motif, direction)] = calib.calibrate(
                model, masks, selection_cache[direction], cfg.chunk)
            native_p = io.p_from_z(model.observed_z)
            argmin, effect, region_p = {}, {}, {}
            for region_id, region in enumerate(io.REGIONS):
                window_ids = np.flatnonzero(model.region_index == region_id)
                local = native_p[window_ids]
                best = int(window_ids[int(np.argmin(local))])
                argmin[region] = best
                effect[region] = calib.effect_at(model, best)
                region_p[region] = float(local.min())
            for pool, members in io.POOL_TO_REGIONS.items():
                best_region = min(members, key=lambda r: region_p[r])
                argmin[pool] = argmin[best_region]
                effect[pool] = effect[best_region]
                region_p[pool] = region_p[best_region]
            native[(motif, direction)] = region_p
            meta[(motif, direction)] = {"n1": model.n1, "n0": model.n0, "rbp": rbp,
                                        "argmin": argmin, "effect": effect}
            del model
        if index % 25 == 0 or index == len(motifs):
            emit("  stage 1 [{} {}]: {}/{} motifs, {:.1f}s".format(
                arm, label, index, len(motifs), time.perf_counter() - started))
    stage1_seconds = time.perf_counter() - started

    selected = [] if stage1_only else [k for _, k in sorted(
        ((min(v["calibrated_p"] for v in stage1[k].values() if math.isfinite(v["calibrated_p"])), k)
         for k in stage1
         if any(math.isfinite(v["calibrated_p"]) and v["calibrated_p"] <= cfg.refine_threshold
                for v in stage1[k].values())), key=lambda x: (x[0], x[1]))]
    emit("stage 2 [{} {}]: {} of {} pairs promoted{}".format(
        arm, label, len(selected), len(stage1), " (STAGE 1 ONLY)" if stage1_only else ""))

    stage = {k: 1 for k in stage1}
    final = {k: dict(v) for k, v in stage1.items()}
    t2 = time.perf_counter()
    for done, (motif, direction) in enumerate(selected, 1):
        model = model_for(motif, direction)
        masks = model.window_masks()
        rng = unit.stage_rng(cfg, 2, direction)
        observed = calib.statistic_maxima(model.observed_z, masks)
        exceed = {name: 0 for name in masks}
        remaining = cfg.stage2_perms
        while remaining > 0:
            batch = min(cfg.chunk, remaining)
            maxima = calib.statistic_maxima(
                model.permutation_z(drawers[direction].draw(rng, batch)), masks)
            for name in masks:
                if math.isfinite(observed[name]):
                    v = maxima[name]
                    exceed[name] += int(np.sum(v[np.isfinite(v)] >= observed[name]))
            remaining -= batch
        for name in masks:
            if not math.isfinite(observed[name]):
                continue
            final[(motif, direction)][name] = dict(final[(motif, direction)][name])
            final[(motif, direction)][name]["calibrated_p"] = (
                (1 + exceed[name]) / (1 + cfg.stage2_perms))
            final[(motif, direction)][name]["permutations"] = cfg.stage2_perms
        stage[(motif, direction)] = 2
        del model
        if done % 10 == 0 or done == len(selected):
            emit("  stage 2 [{} {}]: {}/{} pairs, {:.1f}s".format(
                arm, label, done, len(selected), time.perf_counter() - t2))
    stage2_seconds = time.perf_counter() - t2

    bh_keys = [(m, d, p) for m in motifs for d in ("up", "dn") for p in io.PLOT_POOLS]
    q_lookup = {k: float(q) for k, q in zip(bh_keys, calib.bh_adjust(
        [final[(m, d)][p]["calibrated_p"] for m, d, p in bh_keys]))}
    n_inc, n_skp = meta[(motifs[0], "up")]["n1"], meta[(motifs[0], "dn")]["n1"]
    condensed = []
    for rbp in sorted({alias.get(m.split(".", 1)[0], m.split(".", 1)[0]) for m in motifs}):
        rbp_motifs = [m for m in motifs if alias.get(m.split(".", 1)[0], m.split(".", 1)[0]) == rbp]
        for direction in ("up", "dn"):
            for pool in io.POOL_TO_REGIONS:
                finite = [(final[(m, direction)][pool]["calibrated_p"], native[(m, direction)][pool], m)
                          for m in rbp_motifs
                          if math.isfinite(final[(m, direction)][pool]["calibrated_p"])]
                motif = min(finite)[2] if finite else rbp_motifs[0]
                record = {
                    "arm": arm, "treatment": label, "RBP": rbp,
                    "rbp_table_names": ",".join(sorted({m.split(".", 1)[0] for m in rbp_motifs})),
                    "direction": direction, "direction_label": io.DIRECTION_LABEL[direction],
                    "pooled_region": pool, "plot": pool in io.PLOT_POOLS,
                    "selected_motif_key": motif,
                    "selected_native_region": min(io.POOL_TO_REGIONS[pool],
                                                  key=lambda r: native[(motif, direction)][r]),
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
                    "power_label": unit.power_label(arm, cfg.underpowered, n_inc, n_skp),
                    "n_changed_included": n_inc, "n_changed_skipped": n_skp,
                }
                record.update(meta[(motif, direction)]["effect"][pool])
                condensed.append(record)
    path = out_dir / (label + "_condensed_per_rbp.tsv")
    calib.write_tsv(path, condensed, unit.CONDENSED_COLUMNS)
    report = {"arm": arm, "label": label, "resample_seed": seed, "permutation_seed": cfg.seed,
              "stage1_permutations": cfg.stage1_perms,
              "stage2_permutations": 0 if stage1_only else cfg.stage2_perms,
              "stage1_only": stage1_only, "pairs_total": len(stage1),
              "pairs_refined": len(selected), "bh_family_size": len(bh_keys),
              "stage1_wall_seconds": stage1_seconds, "stage2_wall_seconds": stage2_seconds,
              "length_match": match_info,
              "cluster_sizes": {d: drawers[d].need for d in drawers}}
    (out_dir / (label + "_report.json")).write_text(json.dumps(report, indent=2), encoding="utf-8")
    emit("wrote {} ({} rows); stage1 {:.1f}s stage2 {:.1f}s".format(
        path, len(condensed), stage1_seconds, stage2_seconds))
    return path


def load_panel(path: Path):
    out = {}
    for row in unit.read_tsv(path):
        if row["pooled_region"] in io.PLOT_POOLS:
            out[(row["direction"], row["pooled_region"], row["RBP"])] = {
                "RBP": row["RBP"], "p": unit.as_float(row["calibrated_p"]),
                "q": unit.as_float(row["calibrated_q"]),
                "native": unit.as_float(row["native_ranksum_p"])}
    return out


def compare(cfg, arms, seeds, emit):
    if len(seeds) != 2:
        raise SystemExit("compare needs exactly two resample seeds (M1 vs M2 = resampling noise)")
    labels = seed_labels(seeds)
    out_root = Path(cfg.out_root)
    rbp_control = cfg.positive_control
    control_panels = tuple(unit.CONTROL_PANELS) + (("Exon Body", "up"),)
    panel_rows, change_rows, control_rows, length_rows = [], [], [], []
    for arm in arms:
        data = {"C": load_panel(Path(cfg.c_root) / arm / "treatment_C_condensed_per_rbp.tsv")}
        for tag, seed in (("M1", seeds[0]), ("M2", seeds[1])):
            data[tag] = load_panel(out_root / arm / (labels[seed] + "_condensed_per_rbp.tsv"))
        for seed in seeds:
            rep = json.loads((out_root / arm / (labels[seed] + "_report.json")).read_text())
            for d, info in rep["length_match"].items():
                length_rows.append({
                    "arm": arm, "direction": d, "resample_seed": seed,
                    "n_bg_exons_before": info["n_bg_exons_before"],
                    "n_bg_exons_after": info["n_bg_exons_after"],
                    "n_bg_rows_before": info["n_bg_rows_before"],
                    "n_bg_rows_after": info["n_bg_rows_after"],
                    "changed_q10_q25_q50_q75_q90": ",".join(
                        "{:g}".format(v) for v in info["changed_length_quantiles"]),
                    "bg_before_q10_q25_q50_q75_q90": ",".join(
                        "{:g}".format(v) for v in info["bg_length_quantiles_before"]),
                    "bg_after_q10_q25_q50_q75_q90": ",".join(
                        "{:g}".format(v) for v in info["bg_length_quantiles_after"]),
                    "ks_D_before": info["ks_D_before"], "ks_D_after": info["ks_D_after"],
                    "stage1_only": rep["stage1_only"]})
        for d in ("up", "dn"):
            for pool in io.PLOT_POOLS:
                rbps = sorted({k[2] for t in data for k in data[t] if k[:2] == (d, pool)})
                vec = {t: [data[t].get((d, pool, r), {}).get("p", math.nan) for r in rbps]
                       for t in data}
                sig = {t: sorted(r for r in rbps
                                 if math.isfinite(data[t].get((d, pool, r), {}).get("q", math.nan))
                                 and data[t][(d, pool, r)]["q"] < 0.05) for t in data}
                ranks = {t: unit.panel_rank([data[t][(d, pool, r)] for r in rbps
                                             if (d, pool, r) in data[t]]) for t in data}
                row = {"arm": arm, "direction": d, "direction_label": io.DIRECTION_LABEL[d],
                       "pooled_region": pool, "n_RBPs": len(rbps),
                       "spearman_C_vs_M1": unit.rho(vec["C"], vec["M1"]),
                       "spearman_C_vs_M2": unit.rho(vec["C"], vec["M2"]),
                       "spearman_M1_vs_M2": unit.rho(vec["M1"], vec["M2"]),
                       "n_q_lt_0.05_C": len(sig["C"]), "n_q_lt_0.05_M1": len(sig["M1"]),
                       "n_q_lt_0.05_M2": len(sig["M2"])}
                for t in ("M1", "M2"):
                    row["gained_" + t] = ",".join(sorted(set(sig[t]) - set(sig["C"]))) or "none"
                    row["lost_" + t] = ",".join(sorted(set(sig["C"]) - set(sig[t]))) or "none"
                row["M1_vs_M2_discordant"] = ",".join(
                    sorted(set(sig["M1"]) ^ set(sig["M2"]))) or "none"
                panel_rows.append(row)
                for r in rbps:
                    change_rows.append({"arm": arm, "direction": d, "pooled_region": pool,
                                        "RBP": r, **{t + "_" + f: data[t].get((d, pool, r), {}).get(
                                            f, math.nan) for t in data for f in ("p", "q")},
                                        **{t + "_rank": ranks[t].get(r) for t in data}})
                if rbp_control and arm in cfg.control_arms and (pool, d) in control_panels:
                    control_rows.append({
                        "arm": arm, "pooled_region": pool, "direction": d, "n_RBPs": len(rbps),
                        **{rbp_control + "_" + t + "_rank": ranks[t].get(rbp_control) for t in data},
                        **{rbp_control + "_" + t + "_q": data[t].get((d, pool, rbp_control), {}).get(
                            "q", math.nan) for t in data}})
                emit("  {} {} {}: rho C/M1 {:.3f} C/M2 {:.3f} M1/M2 {:.3f}; q<0.05 C {} M1 {} M2 {}"
                     .format(arm, d, pool, row["spearman_C_vs_M1"], row["spearman_C_vs_M2"],
                             row["spearman_M1_vs_M2"], row["n_q_lt_0.05_C"],
                             row["n_q_lt_0.05_M1"], row["n_q_lt_0.05_M2"]))
    pcols = list(panel_rows[0].keys())
    calib.write_tsv(out_root / "panel_summary_C_vs_matched.tsv", panel_rows, pcols)
    calib.write_tsv(out_root / "comparison_C_vs_matched.tsv", change_rows, list(change_rows[0]))
    if control_rows:
        calib.write_tsv(out_root / (rbp_control.lower() + "_control_C_vs_matched.tsv"),
                        control_rows, list(control_rows[0]))
    calib.write_tsv(out_root / "length_match_summary.tsv", length_rows, list(length_rows[0]))
    families = sorted({json.loads((out_root / arm / (labels[s] + "_report.json")).read_text())
                       ["bh_family_size"] for arm in arms for s in seeds})
    readme = [
        ("C", "Treatment C from {}/<ARM>/treatment_C_condensed_per_rbp.tsv "
              "(tools/unit_sensitivity.py), full background, cluster-permutation null; read, not "
              "recomputed.".format(cfg.c_root)),
        ("M1 / M2", "Same cluster-permutation null (seed {}) against a length-matched background: "
                    "background target-exon clusters drawn without replacement, stratified on "
                    "deciles of the changed set's unique target-exon length; resample seeds "
                    "{} (M1) and {} (M2). M1 vs M2 is resampling noise.".format(
                        cfg.seed, seeds[0], seeds[1])),
        ("q", "BH over motifs x 2 directions x 3 plotted pools per arm (family size {}).".format(
            ",".join(map(str, families)))),
        ("Species", "Human / GRCh38 (hg38) / GENCODE v49 unless the inputs say otherwise; SE only."),
    ]
    sheets = {"README": (["field", "description"],
                         [{"field": k, "description": v} for k, v in readme]),
              "panel_summary": (pcols, panel_rows)}
    if control_rows:
        sheets[rbp_control.lower() + "_control"] = (list(control_rows[0]), control_rows)
    sheets["length_match"] = (list(length_rows[0]), length_rows)
    sheets["per_RBP"] = (list(change_rows[0]), change_rows)
    calib.write_workbook(out_root / "comparison_C_vs_matched.xlsx", sheets)
    emit("wrote comparison tables and comparison_C_vs_matched.xlsx")


def write_versions(cfg, seeds):
    with open(Path(cfg.out_root) / "versions.txt", "w", encoding="utf-8") as h:
        h.write("python\t{}\nnumpy\t{}\nscipy\t{}\nopenpyxl\t{}\nplatform\t{}\n".format(
            sys.version.split()[0], np.__version__, scipy.__version__,
            importlib.metadata.version("openpyxl"), platform.platform()))
        here = Path(__file__).resolve()
        h.write("script\t{}\nscript_md5\t{}\n".format(here, io.md5(here)))
        for name in ("unit_sensitivity.py", "calibrate_ranksum.py", "rmaps_calib_v2_lib.py",
                     "rmaps_countdist_io.py"):
            h.write("library_{}_md5\t{}\n".format(name[:-3], io.md5(TOOLS_DIR / name)))
        h.write("permutation_seed\t{}\nresample_seeds\t{}\n".format(
            cfg.seed, ",".join(map(str, seeds))))
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            h.write("{}\t{}\n".format(var, os.environ.get(var, "unset")))
        h.write("generated\t{}\n".format(time.strftime("%Y-%m-%dT%H:%M:%S")))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--mode", required=True, choices=("run", "compare", "versions"))
    parser.add_argument("--arms", nargs="+", required=True)
    parser.add_argument("--out-root", required=True, help="output root; <ARM>/ subfolders")
    parser.add_argument("--counts-root", help="run: <ARM>/<motif>.counts.npz archives")
    parser.add_argument("--alias-table", help="run: table_name -> HGNC symbol TSV")
    parser.add_argument("--c-root", help="compare: tools/unit_sensitivity.py output root "
                                         "(<ARM>/treatment_C_condensed_per_rbp.tsv)")
    parser.add_argument("--seeds", nargs="+", type=int, default=list(RESAMPLE_SEEDS))
    parser.add_argument("--stage1-only", action="store_true")
    parser.add_argument("--size-control", action="store_true",
                        help="uniform random background clusters, same count as the matched set")
    parser.add_argument("--underpowered-arms", nargs="*", default=[])
    parser.add_argument("--positive-control", default=None, help="RBP for the control table")
    parser.add_argument("--control-arms", nargs="*", default=[])
    parser.add_argument("--seed", type=int, default=unit.SEED, help="permutation seed")
    parser.add_argument("--permutations", type=int, default=unit.STAGE1_PERMS)
    parser.add_argument("--refine-perms", type=int, default=unit.STAGE2_PERMS)
    parser.add_argument("--refine-threshold", type=float, default=unit.REFINE_THRESHOLD)
    args = parser.parse_args(argv)

    needs = {"run": ("counts_root", "alias_table"), "compare": ("c_root",), "versions": ()}[args.mode]
    missing = ["--" + n.replace("_", "-") for n in needs if getattr(args, n) is None]
    if missing:
        parser.error("--mode {} requires {}".format(args.mode, ", ".join(missing)))
    cfg = unit.settings(counts_root=args.counts_root, alias_table=args.alias_table,
                        out_root=args.out_root, underpowered=tuple(args.underpowered_arms),
                        positive_control=args.positive_control,
                        control_arms=tuple(args.control_arms), seed=args.seed,
                        stage1_perms=args.permutations, stage2_perms=args.refine_perms,
                        refine_threshold=args.refine_threshold)
    cfg.c_root = args.c_root
    random.seed(cfg.seed)
    np.random.seed(cfg.seed)
    Path(cfg.out_root).mkdir(parents=True, exist_ok=True)
    log = open(Path(cfg.out_root) / "command.log", "a", encoding="utf-8")

    def emit(message):
        print(message, flush=True)
        log.write(message + "\n")
        log.flush()

    t0 = time.perf_counter()
    emit("start={} cmd={}".format(time.strftime("%Y-%m-%dT%H:%M:%S"),
                                  " ".join([sys.executable] + sys.argv)))
    labels = seed_labels(args.seeds)
    if args.mode == "run":
        for seed in args.seeds:
            for arm in args.arms:
                run_matched(cfg, arm, seed, labels[seed], emit, stage1_only=args.stage1_only,
                            size_control=args.size_control)
    elif args.mode == "compare":
        compare(cfg, args.arms, args.seeds, emit)
    else:
        write_versions(cfg, args.seeds)
    emit("mode={} wall_seconds={:.1f} exit=0".format(args.mode, time.perf_counter() - t0))
    log.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
