"""Cross-method rank comparison for the rMAPS3 SE motif maps.

Reads root tables only (no p-value is recomputed for any Fisher or engine
layer), pools the eight regions to the three RBP-RELI regions by the minimum of
the two sub-regions, groups motifs to RBPs by the HGNC alias table (best motif
per RBP), and compares ranks across every method whose root tables exist.

Re-runnable: methods whose engine runs have not landed yet are reported as
unavailable and the script can simply be run again.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import rankdata, spearmanr
from openpyxl import Workbook
from openpyxl.styles import Font

DIRECTIONS = ["up", "dn"]
DIRECTION_LABEL = {"up": "included", "dn": "skipped"}
POOL_TO_REGIONS = {
    "Upstream Intron": ("smallest_p_in_upstreamExonIntron", "smallest_p_in_upstreamIntron"),
    "Exon Body": ("smallest_p_in_targetExon-5prime", "smallest_p_in_targetExon-3prime"),
    "Downstream Intron": ("smallest_p_in_downstreamIntron", "smallest_p_in_downstreamExonIntron"),
}
POOLS = list(POOL_TO_REGIONS)

# name -> (path template, observation scale, closeness to the authors' tool, validity note)
METHODS = {
    "released_fisher": (
        r"{runs}\{arm}\released_b9a9dce\pVal.{direction}.vs.bg.RNAmap.txt",
        "motif hit counts in a 2x2 table",
        "the authors' default; the tool as released",
        "INVALID sampling model: an exon with k hits contributes k successes and "
        "n-k failures to a table whose margins are exon counts, so the 2x2 is not a "
        "count of independent Bernoulli exons and repetitive motifs are inflated.",
    ),
    "released_mannwhitney": (
        r"{stat}\released_mannwhitney\{arm}\pVal.{direction}.vs.bg.RNAmap.txt",
        "per-exon hit counts",
        "built-in --stat-method option of the released tool",
        "Correct sampling unit (one observation per eligible exon) but the engine "
        "calls scipy with the asymptotic normal approximation, which is not valid "
        "in this regime: counts are >99% zero and the groups are extremely "
        "unbalanced, so the tie-corrected variance collapses and tail p-values are "
        "orders of magnitude too small (see normal_approximation_gap.tsv). Tests "
        "stochastic dominance, not a rate.",
    ),
    "audited_fisher_binary": (
        r"{runs}\{arm}\fisher_c34776b\pVal.{direction}.vs.bg.RNAmap.txt",
        "binary exon carries motif / does not",
        "audited engine, same test family as the default",
        "Valid: exact Fisher on one Bernoulli observation per eligible exon; "
        "discards motif density entirely.",
    ),
    "audited_mannwhitney": (
        r"{stat}\audited_mannwhitney\{arm}\pVal.{direction}.vs.bg.RNAmap.txt",
        "binarized observations (WindowCounts.observations() emits 1s and 0s)",
        "built-in option of the audited engine",
        "Not count-aware and not calibrated: the audited engine binarizes before "
        "the rank test, so it scores the SAME 2x2 as audited_fisher_binary but "
        "through a normal approximation; the two differ only by the approximation "
        "error, which is what normal_approximation_gap.tsv measures.",
    ),
    "mw_counts": (
        r"{stat}\count_aware\{arm}\mw_counts_root_tables.{direction}.vs.bg.tsv",
        "per-exon hit counts reconstructed from stored hit spans",
        "downstream re-computation from the audited engine's positional archives",
        "Correct sampling unit and tie correction, and it reproduces the engines' "
        "rank test exactly; it therefore inherits the same normal-approximation "
        "miscalibration as released_mannwhitney. Usable for ORDER, not for the "
        "absolute p-value.",
    ),
    "poisson_rate": (
        r"{stat}\count_aware\{arm}\poisson_rate_root_tables.{direction}.vs.bg.tsv",
        "total hits per eligible exon (rate)",
        "downstream re-computation from the audited engine's positional archives",
        "Exact conditional on the total hit count, so it is calibrated in the sparse "
        "tail where the rank tests are not; it does assume hits are Poisson within a "
        "group, and extra-Poisson dispersion from repetitive motifs is not modelled, "
        "so it remains anticonservative for exactly those motifs.",
    ),
}
COUNT_AWARE = {"released_mannwhitney", "mw_counts", "poisson_rate"}


def read_alias(path):
    mapping = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("table_name") and row.get("hgnc_symbol"):
                mapping[row["table_name"]] = row["hgnc_symbol"]
    return mapping


def read_root_table(path):
    rows = {}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = row.get("RBP")
            if not key:
                continue
            values = {}
            for column in reader.fieldnames[1:]:
                text = (row.get(column) or "").strip()
                values[column] = (math.nan if text.upper() in {"", "NA", "NAN"} else float(text))
            rows[key] = values
    return rows


def method_path(method, arm, direction, runs_root, stat_root):
    template = METHODS[method][0]
    return Path(template.format(runs=runs_root, stat=stat_root, arm=arm, direction=direction))


def pool_values(table):
    """motif key -> {pool: min over the two sub-region columns}."""
    pooled = {}
    for key, values in table.items():
        entry = {}
        for pool, columns in POOL_TO_REGIONS.items():
            candidates = [values.get(c, math.nan) for c in columns]
            candidates = [c for c in candidates if isinstance(c, float) and math.isfinite(c)]
            entry[pool] = min(candidates) if candidates else math.nan
        pooled[key] = entry
    return pooled


def best_per_rbp(pooled, pool, alias):
    best = {}
    for key, entry in pooled.items():
        value = entry.get(pool, math.nan)
        table_name = key.split(".", 1)[0]
        rbp = alias.get(table_name, table_name)
        if not math.isfinite(value):
            continue
        if rbp not in best or value < best[rbp][0]:
            best[rbp] = (value, key)
    return best


def ranks_from_values(rbps, values):
    array = np.array([values.get(r, (math.nan,))[0] for r in rbps], dtype=float)
    filled = np.where(np.isfinite(array), array, 1.0)
    return array, rankdata(filled, method="average")


def top_n(rbps, ranks, n=10):
    order = np.argsort(ranks, kind="stable")
    return {rbps[i] for i in order[:n]}


def md5(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--arms", nargs="+", required=True)
    parser.add_argument("--runs-root", required=True,
                        help="root of the engine run directories, one per arm")
    parser.add_argument("--stat-root", required=True,
                        help="root holding released_mannwhitney/ and count_aware/")
    parser.add_argument("--alias", required=True,
                        help="HGNC alias table: table_name<TAB>hgnc_symbol")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    log_handle = open(out / "command.log", "a", encoding="utf-8")

    def log(message):
        print(message, flush=True)
        log_handle.write(message + "\n")
        log_handle.flush()

    started = time.time()
    log("start={} cmd={}".format(time.strftime("%Y-%m-%dT%H:%M:%S"), " ".join(sys.argv)))
    alias = read_alias(args.alias)
    log("alias rows={} from {}".format(len(alias), args.alias))

    availability = {}
    tables = {}
    for arm in args.arms:
        for method in METHODS:
            for direction in DIRECTIONS:
                path = method_path(method, arm, direction, args.runs_root, args.stat_root)
                ok = path.is_file() and path.stat().st_size > 0
                if ok:
                    table = read_root_table(path)
                    ok = len(table) > 1
                    if ok:
                        tables[(arm, method, direction)] = table
                availability[(arm, method, direction)] = (ok, str(path))
                log("available={} {} {} {} {}".format(int(ok), arm, method, direction, path))

    with open(out / "method_availability.tsv", "w", newline="\n", encoding="utf-8") as handle:
        handle.write("arm\tmethod\tdirection\tavailable\tpath\n")
        for (arm, method, direction), (ok, path) in sorted(availability.items()):
            handle.write("{}\t{}\t{}\t{}\t{}\n".format(arm, method, direction, int(ok), path))

    pooled = {k: pool_values(v) for k, v in tables.items()}

    # ---- long per-panel table -------------------------------------------------
    long_rows = []
    panel_best = {}
    for arm in args.arms:
        for direction in DIRECTIONS:
            for pool in POOLS:
                per_method = {}
                for method in METHODS:
                    key = (arm, method, direction)
                    if key not in pooled:
                        continue
                    per_method[method] = best_per_rbp(pooled[key], pool, alias)
                if not per_method:
                    continue
                rbps = sorted(set().union(*[set(v) for v in per_method.values()]))
                for method, values in per_method.items():
                    array, ranks = ranks_from_values(rbps, values)
                    panel_best[(arm, direction, pool, method)] = (rbps, array, ranks, values)
                    for rbp, value, rank in zip(rbps, array, ranks):
                        long_rows.append({
                            "arm": arm, "pooled_region": pool,
                            "direction": direction,
                            "direction_label": DIRECTION_LABEL[direction],
                            "rbp": rbp, "method": method,
                            "best_motif": values.get(rbp, (math.nan, "NA"))[1],
                            "p_value": value, "rank": rank,
                            "n_rbp_in_panel": len(rbps),
                        })
    long_path = out / "method_rank_comparison.tsv"
    with open(long_path, "w", newline="\n", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n",
                                fieldnames=["arm", "pooled_region", "direction", "direction_label",
                                            "rbp", "method", "best_motif", "p_value", "rank",
                                            "n_rbp_in_panel"])
        writer.writeheader()
        for row in long_rows:
            row = dict(row)
            row["p_value"] = "NA" if not math.isfinite(row["p_value"]) else repr(float(row["p_value"]))
            writer.writerow(row)
    log("wrote {} rows -> {}".format(len(long_rows), long_path))

    # ---- summary: per arm, per method pair ------------------------------------
    summary_rows = []
    methods_present = [m for m in METHODS
                       if any(k[3] == m for k in panel_best)]
    for arm in args.arms:
        for i, a in enumerate(methods_present):
            for b in methods_present[i + 1:]:
                rhos, overlaps, panels = [], [], 0
                for direction in DIRECTIONS:
                    for pool in POOLS:
                        ka = (arm, direction, pool, a)
                        kb = (arm, direction, pool, b)
                        if ka not in panel_best or kb not in panel_best:
                            continue
                        rbps_a, va, ra, _ = panel_best[ka]
                        rbps_b, vb, rb, _ = panel_best[kb]
                        shared = sorted(set(rbps_a) & set(rbps_b))
                        if len(shared) < 3:
                            continue
                        ia = [rbps_a.index(r) for r in shared]
                        ib = [rbps_b.index(r) for r in shared]
                        rho = spearmanr(ra[ia], rb[ib]).statistic
                        if np.isfinite(rho):
                            rhos.append(float(rho))
                        overlaps.append(len(top_n(rbps_a, ra) & top_n(rbps_b, rb)))
                        panels += 1
                if panels:
                    summary_rows.append({
                        "arm": arm, "method_a": a, "method_b": b, "panels_compared": panels,
                        "mean_spearman_rho": float(np.mean(rhos)) if rhos else math.nan,
                        "mean_top10_overlap": float(np.mean(overlaps)) if overlaps else math.nan,
                    })
    summary_path = out / "method_pair_summary.tsv"
    with open(summary_path, "w", newline="\n", encoding="utf-8") as handle:
        handle.write("arm\tmethod_a\tmethod_b\tpanels_compared\tmean_spearman_rho\tmean_top10_overlap\n")
        for row in summary_rows:
            handle.write("{arm}\t{method_a}\t{method_b}\t{panels_compared}\t{rho}\t{ov}\n".format(
                rho=("NA" if not math.isfinite(row["mean_spearman_rho"]) else "%.4f" % row["mean_spearman_rho"]),
                ov=("NA" if not math.isfinite(row["mean_top10_overlap"]) else "%.2f" % row["mean_top10_overlap"]),
                **row))
    log("wrote {} method-pair rows -> {}".format(len(summary_rows), summary_path))

    # ---- QKI positive control -------------------------------------------------
    control_panels = [("up", "Upstream Intron"), ("dn", "Downstream Intron")]
    qki_rows = []
    for arm in [a for a in args.arms if a.startswith("QKI_KO")]:
        for direction, pool in control_panels:
            for method in METHODS:
                key = (arm, direction, pool, method)
                if key not in panel_best:
                    qki_rows.append([arm, pool, DIRECTION_LABEL[direction], method,
                                     "NA", "NA", "NA", "not_available"])
                    continue
                rbps, array, ranks, values = panel_best[key]
                if "QKI" not in rbps:
                    qki_rows.append([arm, pool, DIRECTION_LABEL[direction], method,
                                     "NA", "NA", str(len(rbps)), "QKI_absent_from_panel"])
                    continue
                index = rbps.index("QKI")
                qki_rows.append([arm, pool, DIRECTION_LABEL[direction], method,
                                 "%g" % ranks[index],
                                 ("NA" if not math.isfinite(array[index]) else repr(float(array[index]))),
                                 str(len(rbps)), values["QKI"][1]])
    qki_path = out / "qki_positive_control.tsv"
    with open(qki_path, "w", newline="\n", encoding="utf-8") as handle:
        handle.write("arm\tpooled_region\tdirection\tmethod\tqki_rank\tqki_p\tn_rbp_in_panel\tbest_motif_or_reason\n")
        for row in qki_rows:
            handle.write("\t".join(row) + "\n")
    log("wrote {} rows -> {}".format(len(qki_rows), qki_path))

    # ---- normal-approximation calibration ------------------------------------
    # audited_mannwhitney and audited_fisher_binary score the SAME binary 2x2;
    # they differ only in exact hypergeometric vs normal approximation, so the
    # log10 gap between their root-table values isolates the approximation error.
    gap_rows = []
    for arm in args.arms:
        for direction in DIRECTIONS:
            exact_key = (arm, "audited_fisher_binary", direction)
            approx_key = (arm, "audited_mannwhitney", direction)
            if exact_key not in pooled or approx_key not in pooled:
                continue
            for pool in POOLS:
                gaps = []
                for motif in sorted(set(pooled[exact_key]) & set(pooled[approx_key])):
                    exact = pooled[exact_key][motif][pool]
                    approx = pooled[approx_key][motif][pool]
                    if (math.isfinite(exact) and math.isfinite(approx)
                            and exact > 0 and approx > 0):
                        gaps.append(math.log10(exact) - math.log10(approx))
                if gaps:
                    gap_rows.append([arm, pool, DIRECTION_LABEL[direction], str(len(gaps)),
                                     "%.2f" % float(np.median(gaps)),
                                     "%.2f" % float(np.max(gaps))])
    gap_path = out / "normal_approximation_gap.tsv"
    with open(gap_path, "w", newline="\n", encoding="utf-8") as handle:
        handle.write("arm\tpooled_region\tdirection\tn_motifs\t"
                     "median_log10_exact_minus_log10_approx\tmax_log10_exact_minus_log10_approx\n")
        for row in gap_rows:
            handle.write("\t".join(row) + "\n")
    log("wrote {} rows -> {}".format(len(gap_rows), gap_path))

    # ---- repetitive motif shift ----------------------------------------------
    motif_rows = []
    for arm in args.arms:
        density_path = Path(args.stat_root) / "count_aware" / arm / "motif_density.tsv"
        if not density_path.is_file():
            log("skip repetitive_motif_shift for {}: missing {}".format(arm, density_path))
            continue
        density = {}
        with open(density_path, newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                text = row["mean_hits_per_carrying_exon_window"]
                density[row["RBP"]] = math.nan if text == "NA" else float(text)
        motif_ranks = {}
        for method in METHODS:
            collected = {}
            for direction in DIRECTIONS:
                key = (arm, method, direction)
                if key not in pooled:
                    continue
                for pool in POOLS:
                    keys = sorted(pooled[key])
                    array = np.array([pooled[key][k][pool] for k in keys], dtype=float)
                    ranks = rankdata(np.where(np.isfinite(array), array, 1.0), method="average")
                    for motif, rank in zip(keys, ranks):
                        collected.setdefault(motif, []).append(rank)
            if collected:
                motif_ranks[method] = {m: float(np.mean(v)) for m, v in collected.items()}
        if "released_fisher" not in motif_ranks:
            log("skip repetitive_motif_shift for {}: released_fisher tables absent".format(arm))
            continue
        base = motif_ranks["released_fisher"]
        for motif in sorted(base):
            row = {"arm": arm, "motif": motif,
                   "mean_hits_per_carrying_exon_window": density.get(motif, math.nan),
                   "mean_rank_released_fisher": base[motif]}
            for method in motif_ranks:
                if method == "released_fisher":
                    continue
                value = motif_ranks[method].get(motif, math.nan)
                row["mean_rank_" + method] = value
                row["delta_rank_" + method] = value - base[motif]
            motif_rows.append(row)
    shift_path = out / "repetitive_motif_shift.tsv"
    if motif_rows:
        fields = ["arm", "motif", "mean_hits_per_carrying_exon_window", "mean_rank_released_fisher"]
        for method in METHODS:
            if method == "released_fisher":
                continue
            for prefix in ("mean_rank_", "delta_rank_"):
                name = prefix + method
                if any(name in r for r in motif_rows):
                    fields.append(name)
        with open(shift_path, "w", newline="\n", encoding="utf-8") as handle:
            handle.write("\t".join(fields) + "\n")
            for row in motif_rows:
                cells = []
                for field in fields:
                    value = row.get(field, math.nan)
                    if isinstance(value, float):
                        cells.append("NA" if not math.isfinite(value) else "%.4f" % value)
                    else:
                        cells.append(str(value))
                handle.write("\t".join(cells) + "\n")
        log("wrote {} rows -> {}".format(len(motif_rows), shift_path))

    # ---- workbook -------------------------------------------------------------
    workbook = Workbook()
    readme = workbook.active
    readme.title = "README"
    readme.append(["rMAPS3 SE motif map: cross-method rank comparison"])
    readme.append(["generated", time.strftime("%Y-%m-%dT%H:%M:%S")])
    readme.append(["scope", "human / GRCh38 (hg38) / GENCODE v49; SE events only"])
    readme.append([])
    readme.append(["method", "observation scale", "closeness to the authors' tool", "validity"])
    for method, (_, scale, closeness, validity) in METHODS.items():
        readme.append([method, scale, closeness, validity])
    readme.append([])
    readme.append(["panel", "pooled RBP-RELI region x direction; pooled p = min of the two sub-regions"])
    readme.append(["rbp value", "best (smallest) p across that RBP's motifs"])
    readme.append(["rank", "ascending rank of p within the panel; NA p ranked as p=1"])
    readme["A1"].font = Font(bold=True)

    sheet = workbook.create_sheet("Summary")
    sheet.append(["arm", "method_a", "method_b", "panels_compared",
                  "mean_spearman_rho", "mean_top10_overlap"])
    for row in summary_rows:
        sheet.append([row["arm"], row["method_a"], row["method_b"], row["panels_compared"],
                      None if not math.isfinite(row["mean_spearman_rho"]) else round(row["mean_spearman_rho"], 4),
                      None if not math.isfinite(row["mean_top10_overlap"]) else round(row["mean_top10_overlap"], 2)])

    control = workbook.create_sheet("QKI_positive_control")
    control.append(["arm", "pooled_region", "direction", "method", "qki_rank", "qki_p",
                    "n_rbp_in_panel", "best_motif_or_reason"])
    for row in qki_rows:
        control.append(row)

    for arm in args.arms:
        arm_sheet = workbook.create_sheet(arm[:31])
        header = ["pooled_region", "direction", "rbp"]
        present = [m for m in METHODS if any(k[0] == arm and k[3] == m for k in panel_best)]
        for method in present:
            header += ["p_" + method, "rank_" + method, "motif_" + method]
        arm_sheet.append(header)
        for direction in DIRECTIONS:
            for pool in POOLS:
                keys = [(arm, direction, pool, m) for m in present]
                available = [k for k in keys if k in panel_best]
                if not available:
                    continue
                rbps = sorted(set().union(*[set(panel_best[k][0]) for k in available]))
                for rbp in rbps:
                    row = [pool, DIRECTION_LABEL[direction], rbp]
                    for method in present:
                        key = (arm, direction, pool, method)
                        if key not in panel_best:
                            row += [None, None, None]
                            continue
                        panel_rbps, array, ranks, values = panel_best[key]
                        if rbp not in panel_rbps:
                            row += [None, None, None]
                            continue
                        index = panel_rbps.index(rbp)
                        row += [None if not math.isfinite(array[index]) else float(array[index]),
                                float(ranks[index]), values[rbp][1]]
                    arm_sheet.append(row)
        arm_sheet.auto_filter.ref = arm_sheet.dimensions
    for worksheet in workbook.worksheets:
        worksheet.freeze_panes = "A2"
    xlsx_path = out / "method_rank_comparison.xlsx"
    workbook.save(xlsx_path)
    log("wrote {}".format(xlsx_path))

    with open(out / "versions.txt", "w", encoding="utf-8") as handle:
        import openpyxl
        handle.write("python\t{}\n".format(sys.version.split()[0]))
        handle.write("numpy\t{}\n".format(np.__version__))
        handle.write("scipy\t{}\n".format(scipy.__version__))
        handle.write("openpyxl\t{}\n".format(openpyxl.__version__))
        handle.write("platform\t{}\n".format(platform.platform()))
        handle.write("generated\t{}\n".format(time.strftime("%Y-%m-%dT%H:%M:%S")))
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            handle.write("{}\t{}\n".format(var, os.environ.get(var, "unset")))

    with open(out / "md5sums.txt", "w", encoding="utf-8") as handle:
        for path in sorted(out.glob("*")):
            if path.is_file() and path.name != "md5sums.txt":
                handle.write("{}  {}\n".format(md5(path), path.name))
        for (arm, method, direction), (ok, path) in sorted(availability.items()):
            if ok:
                handle.write("{}  {}\n".format(md5(Path(path)), path))

    write_report(out, args, availability, summary_rows, qki_rows, motif_rows, gap_rows, log)
    log("wall_seconds={:.1f}".format(time.time() - started))
    log("exit=0")
    log_handle.close()


def write_report(out, args, availability, summary_rows, qki_rows, motif_rows, gap_rows, log):
    """Descriptive tables first, then one validity line per method, then a shortlist."""
    lines = ["# Statistical-method comparison for the rMAPS3 SE motif map", ""]
    lines.append("Generated {}. Human / GRCh38 (hg38) / GENCODE v49; SE events only. Arms: {}.".format(
        time.strftime("%Y-%m-%d %H:%M"), ", ".join(args.arms)))
    lines.append("No p-value was recomputed for any Fisher or engine layer; root tables were read as produced. "
                 "Panels are the three RBP-RELI pooled regions x two directions; per-RBP value is the best motif. "
                 "Full tables: `method_rank_comparison.xlsx/.tsv`, `method_pair_summary.tsv`, "
                 "`repetitive_motif_shift.tsv`, `qki_positive_control.tsv`, `method_availability.tsv`.")
    lines.append("")
    lines.append("## Availability")
    lines.append("")
    for arm in args.arms:
        missing = sorted({m for m in METHODS for d in DIRECTIONS
                          if not availability.get((arm, m, d), (False, ""))[0]})
        lines.append("- {}: {}".format(arm, "all six methods present" if not missing
                                       else "missing " + ", ".join(missing)))
    lines.append("")
    lines.append("## QKI positive control (rank; p in parentheses)")
    lines.append("")
    control_methods = [m for m in METHODS if any(r[3] == m and r[4] != "NA" for r in qki_rows)]
    lines.append("| arm x panel | " + " | ".join(control_methods) + " |")
    lines.append("|---" * (len(control_methods) + 1) + "|")
    seen = []
    for row in qki_rows:
        if (row[0], row[1], row[2]) not in seen:
            seen.append((row[0], row[1], row[2]))
    for arm, pool, direction in seen:
        cells = []
        for method in control_methods:
            match = [r for r in qki_rows if r[0] == arm and r[1] == pool
                     and r[2] == direction and r[3] == method]
            if not match or match[0][4] == "NA":
                cells.append("n/a")
            else:
                cells.append("{} ({:.2g})".format(match[0][4], float(match[0][5])))
        lines.append("| {} {} x {} | ".format(arm, pool, direction) + " | ".join(cells) + " |")
    lines.append("")
    lines.append("## Agreement with released_fisher (mean over the six panels per arm)")
    lines.append("")
    lines.append("| method | " + " | ".join(args.arms) + " |")
    lines.append("|---" * (len(args.arms) + 1) + "|")
    for method in METHODS:
        if method == "released_fisher":
            continue
        cells = []
        for arm in args.arms:
            match = [r for r in summary_rows if r["arm"] == arm
                     and "released_fisher" in (r["method_a"], r["method_b"])
                     and method in (r["method_a"], r["method_b"])]
            if not match or not math.isfinite(match[0]["mean_spearman_rho"]):
                cells.append("n/a")
            else:
                cells.append("rho {:.3f}, top10 {:.1f}".format(
                    match[0]["mean_spearman_rho"], match[0]["mean_top10_overlap"]))
        lines.append("| {} | ".format(method) + " | ".join(cells) + " |")
    lines.append("")
    lines.append("All other method pairs are in `method_pair_summary.tsv`; no value is averaged across methods.")
    lines.append("")
    lines.append("## Motifs that fall furthest under the authors' own count-aware option")
    lines.append("")
    field = "delta_rank_released_mannwhitney"
    if motif_rows and any(field in r for r in motif_rows):
        lines.append("| arm | motif | hits per carrying exon-window | mean rank released_fisher | delta rank |")
        lines.append("|---|---|---|---|---|")
        for arm in args.arms:
            subset = [r for r in motif_rows if r["arm"] == arm and math.isfinite(r.get(field, math.nan))]
            subset.sort(key=lambda r: -r[field])
            for row in subset[:3]:
                density = row["mean_hits_per_carrying_exon_window"]
                lines.append("| {} | {} | {} | {:.1f} | +{:.1f} |".format(
                    arm, row["motif"],
                    "NA" if not math.isfinite(density) else "%.3f" % density,
                    row["mean_rank_released_fisher"], row[field]))
        lines.append("")
        lines.append("Same table for every other method, all motifs, in `repetitive_motif_shift.tsv`.")
    else:
        lines.append("Not available: released_mannwhitney root tables or motif density absent for all arms.")
    lines.append("")
    lines.append("## What the rank test ranks, and why its p-values go so low")
    lines.append("")
    probe = Path(args.stat_root) / "count_aware" / "rank_unit_probe.tsv"
    if probe.is_file():
        with open(probe, newline="", encoding="utf-8") as handle:
            probes = list(csv.DictReader(handle, delimiter="\t"))
        lines.append("The ranked unit is ONE observation per eligible exon per window, not per hit and not "
                     "per exon-position pair: materialising the per-exon 0/1 vector explicitly and handing it "
                     "to `scipy.stats.mannwhitneyu` reproduces the engine root-table value exactly "
                     "(`count_aware/rank_unit_probe.tsv`). The extreme p-values are the tie-corrected normal "
                     "approximation, not a larger sample: with >99% of exons at zero the tie correction "
                     "collapses the variance, so the usual untied bound on |z| does not apply.")
        lines.append("")
        lines.append("| probe | changed / bg exons | with motif | rank sum (per exon) | scipy MWU on the explicit vector | exact hypergeometric, same 2x2 |")
        lines.append("|---|---|---|---|---|---|")
        for row in probes:
            lines.append("| {} {} {} x {} | {} / {} | {} / {} | {:.3g} | {:.3g} | {:.3g} |".format(
                row["arm"], row["motif"], row["region"], row["direction"],
                row["n_changed_exons"], row["n_bg_exons"],
                row["changed_with_motif"], row["bg_with_motif"],
                float(row["p_binary_rank_sum_per_exon"]),
                float(row["p_scipy_mwu_explicit_per_exon_vector"]),
                float(row["p_exact_hypergeometric_same_2x2"])))
    else:
        lines.append("Not available: {} absent.".format(probe))
    lines.append("")
    lines.append("## Calibration of the rank tests (same 2x2, exact vs normal approximation)")
    lines.append("")
    if gap_rows:
        lines.append("`audited_mannwhitney` and `audited_fisher_binary` score the same binary 2x2 and differ "
                     "only by the normal approximation, so log10(exact) - log10(approximate) is the "
                     "approximation error. Positive values mean the rank test reports a smaller p than the "
                     "exact test on identical data.")
        lines.append("")
        lines.append("| arm | worst panel | n motifs | median log10 gap there | max log10 gap there |")
        lines.append("|---|---|---|---|---|")
        for arm in args.arms:
            subset = [r for r in gap_rows if r[0] == arm]
            if not subset:
                continue
            worst = max(subset, key=lambda r: float(r[5]))
            lines.append("| {} | {} x {} | {} | {} | {} |".format(*worst))
        lines.append("")
        lines.append("All twelve panels are in `normal_approximation_gap.tsv`.")
    else:
        lines.append("Not computable: audited_mannwhitney or audited_fisher_binary root tables absent.")
    lines.append("")
    lines.append("## Validity, closeness to the authors' tool, and effect on the order")
    lines.append("")
    for method, (_, scale, closeness, validity) in METHODS.items():
        per_arm = []
        for row in summary_rows:
            if "released_fisher" in (row["method_a"], row["method_b"]) and method in (row["method_a"], row["method_b"]):
                if math.isfinite(row["mean_spearman_rho"]):
                    per_arm.append("{} {:.2f}".format(row["arm"], row["mean_spearman_rho"]))
        order = ("reference order" if method == "released_fisher"
                 else ("rho vs released_fisher: " + "; ".join(per_arm) if per_arm
                       else "not comparable (tables absent)"))
        lines.append("- **{}** — tests {}; {}. {} Order: {}.".format(
            method, scale, closeness, validity, order))
    lines.append("")
    lines.append("## Shortlist: valid AND closest to the authors' tool")
    lines.append("")
    lines.append("1. `audited_fisher_binary` — the only layer that is both a correct sampling model and exactly "
                 "calibrated; same test family as the tool's default, but count-blind by construction.")
    lines.append("2. `poisson_rate` — the only count-aware layer that stays calibrated in this sparse tail "
                 "(exact conditional on the total), with an interpretable rate ratio; downstream of the "
                 "authors' tool, and still anticonservative for over-dispersed repetitive motifs.")
    lines.append("3. `released_mannwhitney` — closest to the authors' tool (a built-in `--stat-method`) and the "
                 "correct sampling unit, but its p-values are not usable as p-values here; usable for ORDER only.")
    lines.append("4. `mw_counts` — the same statistic as 3, recomputed downstream; a reproducibility check on 3 "
                 "and the layer that isolates count-awareness from engine differences. Same calibration caveat.")
    lines.append("")
    lines.append("`released_fisher` is excluded on its sampling model. `audited_mannwhitney` is excluded because "
                 "it binarizes before ranking: it is a normal-approximate version of `audited_fisher_binary`, "
                 "not a count-aware test, so it cannot answer the question that motivated the comparison. "
                 "No recommendation is made beyond this shortlist.")
    path = out / "REPORT.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log("wrote {} ({} lines)".format(path, len(lines)))


if __name__ == "__main__":
    main()
