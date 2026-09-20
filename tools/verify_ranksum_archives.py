"""Recompute the released rMAPS3 mannwhitney p-values from the npz archives only.

Checks, per arm:
  root tables - every value of pVal.{up,dn}.vs.bg.RNAmap.txt (126 motifs x 8
                sub-regions x 2 directions) against the regional minimum of the
                recomputed per-position p-values;
  per-position - a seeded random sample of per-motif pVal temp-file values.

The recomputation uses only the archives written by countdist_to_npz.py and the
tie-corrected normal approximation with continuity correction, which is what
scipy.stats.mannwhitneyu(alternative='greater') applies in this regime.
"""

from __future__ import annotations

import argparse
import math
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rmaps_countdist_io as io  # noqa: E402

SEED = 149
ROOT_TOLERANCE = 1e-6


def log10_gap(a: float, b: float) -> float:
    """|log10 a - log10 b|, with both-zero treated as identical."""
    if a == b:
        return 0.0
    if a <= 0.0 or b <= 0.0:
        return math.inf
    return abs(math.log10(a) - math.log10(b))


def motif_pvalues(npz_path: Path):
    """Per-position p for both directions, plus the position axis."""
    with np.load(npz_path, allow_pickle=False) as data:
        region_index = np.asarray(data["region_of_position"])
        position = np.asarray(data["position"])
        bg = io.group_csr(data, "bg")
        groups = {d: io.group_csr(data, d) for d in ("up", "dn")}
        kmax = int(max([bg.data.max(initial=0)]
                       + [g.data.max(initial=0) for g in groups.values()]))
        hist_bg = io.csr_histograms(bg, kmax)
        out = {}
        for direction, matrix in groups.items():
            hist_fg = io.csr_histograms(matrix, kmax)
            _, _, _, z, _ = io.rank_statistics(hist_fg, hist_bg)
            out[direction] = io.p_from_z(z)
    return region_index, position, out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", required=True)
    parser.add_argument("--released-root", required=True,
                        help="root holding <arm>/pVal.{up,dn}.vs.bg.RNAmap.txt")
    parser.add_argument("--counts-root", required=True,
                        help="root holding <arm>/*.counts.npz from countdist_to_npz.py")
    parser.add_argument("--positional-sample", type=int, default=200)
    args = parser.parse_args(argv)

    arm_dir = Path(args.released_root) / args.arm
    out_dir = Path(args.counts_root) / args.arm
    log = open(out_dir / "command.log", "a", encoding="utf-8")

    def emit(message):
        print(message, flush=True)
        log.write(message + "\n")
        log.flush()

    emit("start={} cmd={}".format(time.strftime("%Y-%m-%dT%H:%M:%S"),
                                  " ".join([sys.executable] + sys.argv)))
    roots = {d: io.read_root_table(arm_dir / ("pVal." + d + ".vs.bg.RNAmap.txt"))
             for d in ("up", "dn")}
    motifs = sorted(p.name[: -len(".counts.npz")] for p in out_dir.glob("*.counts.npz"))
    missing = [m for d in ("up", "dn") for m in roots[d] if m not in set(motifs)]
    if missing:
        raise SystemExit("root table motifs absent from archives: " + ", ".join(sorted(set(missing))))

    rng = np.random.default_rng(SEED)
    recomputed = {}
    root_rows = []
    worst_root = (0.0, None)
    exact_root = 0
    started = time.perf_counter()
    for index, motif in enumerate(motifs, 1):
        region_index, position, pvalues = motif_pvalues(out_dir / (motif + ".counts.npz"))
        recomputed[motif] = (region_index, position, pvalues)
        for direction in ("up", "dn"):
            for region_id, region in enumerate(io.REGIONS):
                mask = region_index == region_id
                mine = float(pvalues[direction][mask].min())
                released = roots[direction][motif][region]
                gap = log10_gap(mine, released)
                exact_root += int(mine == released)
                if gap > worst_root[0]:
                    worst_root = (gap, (motif, direction, region, mine, released))
                root_rows.append((motif, direction, region, released, mine, gap))
        if index % 25 == 0 or index == len(motifs):
            emit("  root check {}/{} motifs, {:.1f}s".format(
                index, len(motifs), time.perf_counter() - started))

    n_root = len(root_rows)
    root_pass = worst_root[0] <= ROOT_TOLERANCE

    sample_rows = []
    worst_sample = (0.0, None)
    picks = rng.choice(len(motifs) * 2, size=min(args.positional_sample, len(motifs) * 2),
                       replace=False)
    for pick in picks:
        motif = motifs[int(pick) // 2]
        direction = ("up", "dn")[int(pick) % 2]
        regions, positions, released = io.read_positional_p(
            arm_dir / "temp" / (motif + ".pVal." + direction + ".vs.bg.txt"))
        region_index, position, pvalues = recomputed[motif]
        order = np.array([io.REGIONS.index(r) for r in regions], dtype=np.int8)
        if not (np.array_equal(order, region_index) and np.array_equal(positions, position)):
            raise SystemExit("positional axis mismatch for " + motif + " " + direction)
        slot = int(rng.integers(0, released.size))
        mine = float(pvalues[direction][slot])
        gap = log10_gap(mine, float(released[slot]))
        if gap > worst_sample[0]:
            worst_sample = (gap, (motif, direction, io.REGIONS[order[slot]], int(positions[slot])))
        sample_rows.append((motif, direction, io.REGIONS[order[slot]], int(positions[slot]),
                            float(released[slot]), mine, gap))
    sample_pass = worst_sample[0] <= ROOT_TOLERANCE

    probe = []
    probe_motif = "QKI.ACTAAC[ACG]" if "QKI.ACTAAC[ACG]" in set(motifs) else motifs[0]
    with np.load(out_dir / (probe_motif + ".counts.npz"), allow_pickle=False) as data:
        bg = io.group_csr(data, "bg").toarray()
        fg = io.group_csr(data, "dn").toarray()
    region_index, position, pvalues = recomputed[probe_motif]
    window = int(np.argmin(pvalues["dn"]))
    for label, kwargs in (("scipy default", {}),
                          ("use_continuity=False", {"use_continuity": False})):
        value = float(mannwhitneyu(fg[:, window], bg[:, window],
                                   alternative="greater", **kwargs).pvalue)
        probe.append((label, value, log10_gap(value, float(pvalues["dn"][window]))))

    with open(out_dir / "detail_root_comparison.tsv", "w", newline="\n", encoding="utf-8") as handle:
        handle.write("motif\tdirection\tregion\treleased_p\trecomputed_p\tabs_log10_gap\n")
        for row in root_rows:
            handle.write("\t".join(repr(v) if isinstance(v, float) else str(v) for v in row) + "\n")
    with open(out_dir / "detail_positional_sample.tsv", "w", newline="\n", encoding="utf-8") as handle:
        handle.write("motif\tdirection\tregion\tposition\treleased_p\trecomputed_p\tabs_log10_gap\n")
        for row in sample_rows:
            handle.write("\t".join(repr(v) if isinstance(v, float) else str(v) for v in row) + "\n")

    lines = [
        "# Stage A verification - " + args.arm,
        "",
        "Recomputed from `" + str(out_dir) + "` archives only; released tables read, never written.",
        "Statistic: tie-corrected normal approximation with continuity correction, one-sided "
        "changed > background, non-finite replaced by 1.0 (the released engine convention).",
        "",
        "| check | n compared | max abs log10 p gap | exact bit-for-bit | verdict |",
        "|---|---|---|---|---|",
        "| root tables (126 motifs x 8 sub-regions x 2 directions) | {} | {:.3e} | {} | {} |".format(
            n_root, worst_root[0], exact_root, "PASS" if root_pass else "FAIL"),
        "| per-position sample from temp pVal files | {} | {:.3e} | {} | {} |".format(
            len(sample_rows), worst_sample[0],
            sum(1 for r in sample_rows if r[4] == r[5]), "PASS" if sample_pass else "FAIL"),
        "",
        "Tolerance: max abs log10 p gap <= {:g}.".format(ROOT_TOLERANCE),
        "",
        "## Continuity setting of the released engine",
        "",
        "Probe: {} direction dn, window index {} of the 1200-position axis.".format(
            probe_motif, window),
        "",
        "| scipy call | p | abs log10 gap to recomputed |",
        "|---|---|---|",
    ]
    for label, value, gap in probe:
        lines.append("| mannwhitneyu(alternative='greater', {}) | {:.6e} | {:.3e} |".format(
            label, value, gap))
    lines += [
        "",
        "The released engine calls scipy with its defaults, so use_continuity=True; the "
        "recomputation matches that call and not the uncorrected one.",
        "",
        "Worst root-table case: " + (
            "none (all exact)" if worst_root[1] is None else
            "{} {} {} recomputed {:.6e} vs released {:.6e}".format(*worst_root[1])),
        "",
        "Detail: `detail_root_comparison.tsv`, `detail_positional_sample.tsv`.",
        "Generated " + time.strftime("%Y-%m-%dT%H:%M:%S") + "; seed " + str(SEED) + ".",
    ]
    (out_dir / "VERIFY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    with open(out_dir / "versions.txt", "a", encoding="utf-8") as handle:
        handle.write("verify_script\t{}\n".format(Path(__file__).resolve()))
        handle.write("verify_script_md5\t{}\n".format(io.md5(Path(__file__).resolve())))
        handle.write("verify_generated\t{}\n".format(time.strftime("%Y-%m-%dT%H:%M:%S")))
        handle.write("verify_numpy\t{}\n".format(np.__version__))
        handle.write("verify_scipy\t{}\n".format(scipy.__version__))
        handle.write("verify_platform\t{}\n".format(platform.platform()))
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            handle.write("verify_{}\t{}\n".format(var, os.environ.get(var, "unset")))

    emit("root: n={} max_log10_gap={:.3e} exact={} verdict={}".format(
        n_root, worst_root[0], exact_root, "PASS" if root_pass else "FAIL"))
    emit("positional sample: n={} max_log10_gap={:.3e} verdict={}".format(
        len(sample_rows), worst_sample[0], "PASS" if sample_pass else "FAIL"))
    emit("exit={}".format(0 if (root_pass and sample_pass) else 1))
    log.close()
    return 0 if (root_pass and sample_pass) else 1


if __name__ == "__main__":
    raise SystemExit(main())
