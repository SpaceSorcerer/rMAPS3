"""Count-aware per-window statistics for rMAPS3 SE motif maps.

Reads the audited engine's sparse positional archives (schema 2), reconstructs
the per-exon motif hit COUNT in every window from the stored hit spans, and
computes two count-aware window tests per motif x region x direction:

  mw_counts    one-sided Mann-Whitney U (changed > background) on per-exon counts
  poisson_rate one-sided two-sample Poisson rate test (conditional binomial)
               on total hits given eligible exons; effect size = rate ratio

Regional minima over windows are taken as the tool does and written in the
tool's 9-column root-table layout, one file per direction.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import norm, binom

REGIONS = (
    "UpstreamExon_3prime",
    "UpstreamExonIntron",
    "UpstreamIntron",
    "TargetExon_5prime",
    "TargetExon-3prime",
    "DownstreamIntron",
    "DownstreamExonIntron",
    "DownstreamExon_5prime",
)
ROOT_COLUMNS = {
    "UpstreamExon_3prime": "smallest_p_in_upstreamExon-3prime",
    "UpstreamExonIntron": "smallest_p_in_upstreamExonIntron",
    "UpstreamIntron": "smallest_p_in_upstreamIntron",
    "TargetExon_5prime": "smallest_p_in_targetExon-5prime",
    "TargetExon-3prime": "smallest_p_in_targetExon-3prime",
    "DownstreamIntron": "smallest_p_in_downstreamIntron",
    "DownstreamExonIntron": "smallest_p_in_downstreamExonIntron",
    "DownstreamExon_5prime": "smallest_p_in_downstreamExon-5prime",
}
DIRECTION_CODE = {"up": 0, "dn": 1}
BG_LABEL = 2


def load_counts(path):
    """Per-exon hit COUNT per window, plus eligibility and set labels."""
    with np.load(path, allow_pickle=False) as data:
        if int(np.asarray(data["schema_version"]).item()) != 2:
            raise ValueError("Unsupported sparse schema_version in " + str(path))
        window, step, length = (int(data[k]) for k in ("window", "step", "region_length"))
        positions = np.arange(0, length - window + 1, step)
        labels = data["set_label"].copy()
        lo, hi = data["elig_lo"], data["elig_hi"]
        eligible = ((positions[None, :] >= lo[:, None])
                    & (positions[None, :] < hi[:, None])
                    & (lo[:, None] >= 0))
        counts = np.zeros(eligible.shape, dtype=np.int16)
        for event, start, end in zip(data["event_index"], data["hit_start"], data["hit_end"]):
            first = max(0, (int(start) - window) // step + 1)
            last = min(len(positions), (int(end) - 1) // step + 1)
            counts[int(event), first:last] += 1
    counts[~eligible] = 0
    return counts, eligible, labels


def group_histograms(counts, eligible, rows, kmax):
    """Per-window frequency of each count value 0..kmax within `rows`."""
    sub_counts = counts[rows]
    sub_elig = eligible[rows]
    hist = np.empty((kmax + 1, counts.shape[1]), dtype=np.int64)
    hist[0] = (sub_elig & (sub_counts == 0)).sum(axis=0)
    for k in range(1, kmax + 1):
        hist[k] = (sub_counts == k).sum(axis=0)
    return hist


def mw_counts_pvalues(hist1, hist0):
    """One-sided (greater) Mann-Whitney U from per-window count histograms.

    Normal approximation with tie correction and continuity correction; this is
    scipy.stats.mannwhitneyu(method='asymptotic', alternative='greater').
    """
    hist1 = np.asarray(hist1, dtype=np.float64)
    hist0 = np.asarray(hist0, dtype=np.float64)
    n1 = hist1.sum(axis=0)
    n0 = hist0.sum(axis=0)
    tied = hist1 + hist0
    cum_below = np.cumsum(tied, axis=0) - tied
    midrank = cum_below + (tied + 1.0) / 2.0
    rank_sum = (hist1 * midrank).sum(axis=0)
    u1 = rank_sum - n1 * (n1 + 1.0) / 2.0
    n = n1 + n0
    mu = n1 * n0 / 2.0
    tie_term = (tied ** 3 - tied).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        var = n1 * n0 / 12.0 * ((n + 1.0) - tie_term / (n * (n - 1.0)))
        z = (u1 - mu - 0.5) / np.sqrt(var)
        p = norm.sf(z)
    p = np.where((n1 > 0) & (n0 > 0) & (var > 0), p, np.nan)
    return p


def poisson_rate_pvalues(hist1, hist0):
    """Conditional-binomial two-sample Poisson rate test (greater) + rate ratio."""
    hist1 = np.asarray(hist1, dtype=np.float64)
    hist0 = np.asarray(hist0, dtype=np.float64)
    values = np.arange(hist1.shape[0], dtype=np.float64)[:, None]
    n1 = hist1.sum(axis=0)
    n0 = hist0.sum(axis=0)
    x1 = (hist1 * values).sum(axis=0)
    x0 = (hist0 * values).sum(axis=0)
    total = x1 + x0
    with np.errstate(invalid="ignore", divide="ignore"):
        prob = n1 / (n1 + n0)
        p = binom.sf(x1 - 1, total, prob)
        ratio = (x1 / n1) / (x0 / n0)
    valid = (n1 > 0) & (n0 > 0) & (total > 0)
    return np.where(valid, p, np.nan), np.where(valid, ratio, np.nan)


def min_over_windows(p, extra=None):
    p = np.asarray(p, dtype=np.float64)
    if not np.any(np.isfinite(p)):
        return float("nan"), float("nan")
    index = int(np.nanargmin(p))
    return float(p[index]), (float("nan") if extra is None else float(extra[index]))


def write_root_table(path, keys, table):
    header = ["RBP"] + [ROOT_COLUMNS[r] for r in REGIONS]
    with open(path, "w", newline="\n", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for key in keys:
            row = [key]
            for region in REGIONS:
                value = table.get((key, region), float("nan"))
                row.append("NA" if not np.isfinite(value) else repr(value))
            handle.write("\t".join(row) + "\n")


def md5(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def process_arm(arm, runs_root, out_root, engine_run, log):
    positional = Path(runs_root) / arm / engine_run / "positional"
    if not positional.is_dir():
        raise FileNotFoundError("missing positional directory: " + str(positional))
    archives = sorted(positional.glob("*.hits.npz"))
    if not archives:
        raise FileNotFoundError("no sparse archives in " + str(positional))
    out_dir = Path(out_root) / arm
    out_dir.mkdir(parents=True, exist_ok=True)

    tables = {d: {m: {} for m in ("mw_counts", "poisson_rate")} for d in DIRECTION_CODE}
    ratios = {d: {} for d in DIRECTION_CODE}
    density = {}
    keys = []
    start = time.time()
    for index, archive in enumerate(archives, 1):
        stem = archive.name[: -len(".hits.npz")]
        key, region = stem.rsplit(".", 1)
        if region not in ROOT_COLUMNS:
            raise ValueError("unexpected region in " + archive.name)
        if key not in density:
            density[key] = [0, 0]
            keys.append(key)
        counts, eligible, labels = load_counts(archive)
        density[key][0] += int(counts[eligible].sum())
        density[key][1] += int((counts[eligible] >= 1).sum())
        kmax = int(counts.max())
        bg_rows = np.flatnonzero(labels == BG_LABEL)
        hist_bg = group_histograms(counts, eligible, bg_rows, kmax)
        for direction, code in DIRECTION_CODE.items():
            fg_rows = np.flatnonzero(labels == code)
            hist_fg = group_histograms(counts, eligible, fg_rows, kmax)
            p_mw = mw_counts_pvalues(hist_fg, hist_bg)
            p_pois, ratio = poisson_rate_pvalues(hist_fg, hist_bg)
            tables[direction]["mw_counts"][(key, region)] = min_over_windows(p_mw)[0]
            best_p, best_ratio = min_over_windows(p_pois, ratio)
            tables[direction]["poisson_rate"][(key, region)] = best_p
            ratios[direction][(key, region)] = best_ratio
        if index % 200 == 0 or index == len(archives):
            log("  {}: {}/{} archives, {:.1f}s".format(arm, index, len(archives), time.time() - start))

    for direction in DIRECTION_CODE:
        for method in ("mw_counts", "poisson_rate"):
            write_root_table(out_dir / "{}_root_tables.{}.vs.bg.tsv".format(method, direction),
                             keys, tables[direction][method])
        write_root_table(out_dir / "poisson_rate_ratio.{}.vs.bg.tsv".format(direction),
                         keys, ratios[direction])
    with open(out_dir / "motif_density.tsv", "w", newline="\n", encoding="utf-8") as handle:
        handle.write("RBP\ttotal_hits\tcarrying_exon_windows\tmean_hits_per_carrying_exon_window\n")
        for key in keys:
            total, carrying = density[key]
            mean = "NA" if not carrying else repr(total / carrying)
            handle.write("{}\t{}\t{}\t{}\n".format(key, total, carrying, mean))
    elapsed = time.time() - start
    log("  {}: wrote {} in {:.1f}s".format(arm, out_dir, elapsed))
    return len(keys), elapsed


def probe_rank_unit(runs_root, out_root, probes, engine_run, log):
    """Show what the engines' rank test actually ranks, at named windows.

    For each probe the per-exon 0/1 vector is materialised explicitly and handed
    to scipy.stats.mannwhitneyu; if that reproduces the engine's root-table value
    the ranked unit is one observation per eligible exon per window, and any gap
    to the exact hypergeometric on the same 2x2 is approximation error, not a
    different sampling unit.
    """
    from scipy.stats import hypergeom, mannwhitneyu

    path = Path(out_root) / "rank_unit_probe.tsv"
    with open(path, "w", newline="\n", encoding="utf-8") as handle:
        handle.write("arm\tmotif\tregion\tdirection\twindow\tn_changed_exons\tn_bg_exons\t"
                     "changed_with_motif\tbg_with_motif\tmax_hits_per_exon\t"
                     "p_binary_rank_sum_per_exon\tp_scipy_mwu_explicit_per_exon_vector\t"
                     "p_exact_hypergeometric_same_2x2\tp_rank_sum_on_counts\n")
        for arm, motif, region, direction in probes:
            archive = Path(runs_root) / arm / engine_run / "positional" / (
                "{}.{}.hits.npz".format(motif, region))
            if not archive.is_file():
                log("probe skipped, missing " + str(archive))
                continue
            counts, eligible, labels = load_counts(archive)
            fg_rows = np.flatnonzero(labels == DIRECTION_CODE[direction])
            bg_rows = np.flatnonzero(labels == BG_LABEL)
            binary = (counts > 0).astype(np.int16)
            b1 = group_histograms(binary, eligible, fg_rows, 1)
            b0 = group_histograms(binary, eligible, bg_rows, 1)
            p_binary = mw_counts_pvalues(b1, b0)
            kmax = int(counts.max())
            p_counts = mw_counts_pvalues(group_histograms(counts, eligible, fg_rows, kmax),
                                         group_histograms(counts, eligible, bg_rows, kmax))
            window = int(np.nanargmin(p_binary))
            a, c = int(b1[1][window]), int(b0[1][window])
            n1, n0 = int(b1[:, window].sum()), int(b0[:, window].sum())
            explicit_fg = np.concatenate([np.ones(a), np.zeros(n1 - a)])
            explicit_bg = np.concatenate([np.ones(c), np.zeros(n0 - c)])
            p_scipy = float(mannwhitneyu(explicit_fg, explicit_bg,
                                         alternative="greater", method="asymptotic").pvalue)
            p_exact = float(hypergeom.sf(a - 1, n1 + n0, a + c, n1))
            handle.write("\t".join([arm, motif, region, direction, str(window), str(n1), str(n0),
                                    str(a), str(c), str(kmax),
                                    repr(float(p_binary[window])), repr(p_scipy), repr(p_exact),
                                    repr(float(p_counts[window]))]) + "\n")
            log("probe {} {} {} {}: per-exon rank sum {:.3g} vs scipy {:.3g} vs exact {:.3g}".format(
                arm, motif, region, direction, float(p_binary[window]), p_scipy, p_exact))
    log("wrote " + str(path))
    return path


def parse_probe(text):
    fields = text.split(":")
    if len(fields) != 4 or fields[3] not in DIRECTION_CODE:
        raise ValueError("--probe must be ARM:MOTIF:REGION:{up|dn}, got " + text)
    return tuple(fields)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--arms", nargs="+", required=True)
    parser.add_argument("--runs-root", required=True,
                        help="root holding <arm>/<engine-run>/positional/*.hits.npz")
    parser.add_argument("--out-root", required=True)
    parser.add_argument("--engine-run", default="fisher_c34776b",
                        help="run directory under <runs-root>/<arm> that holds positional/")
    parser.add_argument("--probe", action="append", default=[],
                        metavar="ARM:MOTIF:REGION:DIRECTION",
                        help="repeatable; writes rank_unit_probe.tsv")
    parser.add_argument("--probe-only", action="store_true",
                        help="write rank_unit_probe.tsv and exit")
    args = parser.parse_args()

    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    handle = open(out_root / "command.log", "a", encoding="utf-8")

    def log(message):
        print(message, flush=True)
        handle.write(message + "\n")
        handle.flush()

    log("start={} cmd={}".format(time.strftime("%Y-%m-%dT%H:%M:%S"), " ".join(sys.argv)))
    probes = [parse_probe(text) for text in args.probe]
    if probes:
        probe_rank_unit(args.runs_root, out_root, probes, args.engine_run, log)
    if args.probe_only:
        log("exit=0")
        handle.close()
        return
    for arm in args.arms:
        motifs, seconds = process_arm(arm, args.runs_root, out_root, args.engine_run, log)
        for direction in ("up", "dn"):
            for method in ("mw_counts", "poisson_rate"):
                path = out_root / arm / "{}_root_tables.{}.vs.bg.tsv".format(method, direction)
                log("md5 {}  {}".format(md5(path), path))
        log("arm={} motifs={} wall_seconds={:.1f}".format(arm, motifs, seconds))
    with open(out_root / "versions.txt", "w", encoding="utf-8") as vh:
        vh.write("python\t{}\n".format(sys.version.split()[0]))
        vh.write("numpy\t{}\n".format(np.__version__))
        vh.write("scipy\t{}\n".format(scipy.__version__))
        vh.write("platform\t{}\n".format(platform.platform()))
        vh.write("generated\t{}\n".format(time.strftime("%Y-%m-%dT%H:%M:%S")))
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            vh.write("{}\t{}\n".format(var, os.environ.get(var, "unset")))
    log("exit=0")
    handle.close()


if __name__ == "__main__":
    main()
