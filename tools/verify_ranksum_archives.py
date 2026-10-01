"""Recompute the released rMAPS3 mannwhitney p-values from the npz archives only, exactly.

Checks, per arm, with exact float equality (no tolerance; the archives reproduce the engine bit for bit):
  motif set    - the archives and both root tables carry the same motif keys;
  root tables  - every value of pVal.{up,dn}.vs.bg.RNAmap.txt (motifs x 8 sub-regions x 2 directions)
                 against the regional minimum of the recomputed per-position p-values;
  per-position - every value of every temp/<motif>.pVal.{up,dn}.vs.bg.txt, on the same position axis;
  countDist    - every temp/<motif>.countDist.{up,dn,bg}.txt still hashes to the md5 recorded by
                 countdist_to_npz.py in conversion_manifest.tsv when its archive was written.

Only when every check passes does it write verified_temporaries.tsv (path, bytes, md5, kind): the
exact list of temporaries the archives now stand in for. tools/rmaps3_skill_run.py deletes those
files and nothing else, and only after this script exits 0. Any failure writes no list.

The recomputation uses the tie-corrected normal approximation with continuity correction, which is
what scipy.stats.mannwhitneyu(alternative='greater') applies in this regime.
"""

from __future__ import annotations

import argparse
import csv
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

MANIFEST_NAME = "verified_temporaries.tsv"
MANIFEST_COLUMNS = ("path", "bytes", "md5", "kind")


def log10_gap(a: float, b: float) -> float:
    """|log10 a - log10 b|, recorded for mismatches only; the verdict is exact equality."""
    if a == b:
        return 0.0
    if not (a > 0.0 and b > 0.0):
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


def read_conversion_manifest(path: Path) -> dict:
    """{motif: {group: countDist md5}} as recorded by countdist_to_npz.py."""
    with open(path, "r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {r["motif"]: {g: r[g + "_source_md5"] for g in io.GROUPS} for r in rows}


def read_manifest(path: Path):
    """The verified list as [(path, bytes, md5, kind)]: the wrapper's only deletion source."""
    with open(path, "r", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if tuple(header) != MANIFEST_COLUMNS:
            raise ValueError("unexpected verified-temporaries header in " + str(path))
        out = []
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            out.append((Path(fields[0]), int(fields[1]), fields[2], fields[3]))
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", required=True)
    parser.add_argument("--released-root", required=True,
                        help="root holding <arm>/pVal.{up,dn}.vs.bg.RNAmap.txt and <arm>/temp/")
    parser.add_argument("--counts-root", required=True,
                        help="root holding <arm>/*.counts.npz and conversion_manifest.tsv from countdist_to_npz.py")
    args = parser.parse_args(argv)

    arm_dir = Path(args.released_root) / args.arm
    temp_dir = arm_dir / "temp"
    out_dir = Path(args.counts_root) / args.arm
    manifest_path = out_dir / MANIFEST_NAME
    if manifest_path.exists():
        manifest_path.unlink()  # a stale list must never outlive a new verification
    log = open(out_dir / "command.log", "a", encoding="utf-8")

    def emit(message):
        print(message, flush=True)
        log.write(message + "\n")
        log.flush()

    emit("start={} cmd={}".format(time.strftime("%Y-%m-%dT%H:%M:%S"), " ".join(
        [sys.executable] + (sys.argv if argv is None else [__file__] + list(argv)))))
    failures = []
    roots = {d: io.read_root_table(arm_dir / ("pVal." + d + ".vs.bg.RNAmap.txt")) for d in ("up", "dn")}
    motifs = sorted(p.name[: -len(".counts.npz")] for p in out_dir.glob("*.counts.npz"))
    if not motifs:
        failures.append("no *.counts.npz archives in " + str(out_dir))
    for d in ("up", "dn"):
        if set(roots[d]) != set(motifs):
            failures.append("motif set differs between the archives and the {} root table: only in the root "
                            "table {}, only in the archives {}".format(
                                d, sorted(set(roots[d]) - set(motifs))[:5], sorted(set(motifs) - set(roots[d]))[:5]))
    conversion_path = out_dir / "conversion_manifest.tsv"
    conversion = read_conversion_manifest(conversion_path) if conversion_path.is_file() else {}
    if not conversion:
        failures.append("conversion_manifest.tsv absent or empty; countDist provenance cannot be checked")

    root_rows, position_rows, verified = [], [], []
    n_root = exact_root = n_position = exact_position = 0
    started = time.perf_counter()
    for index, motif in enumerate(motifs, 1):
        region_index, position, pvalues = motif_pvalues(out_dir / (motif + ".counts.npz"))
        motif_ok = True
        positional = []
        for direction in ("up", "dn"):
            for region_id, region in enumerate(io.REGIONS):
                mine = float(pvalues[direction][region_index == region_id].min())
                released = roots[direction].get(motif, {}).get(region, math.nan)
                n_root += 1
                if mine == released:
                    exact_root += 1
                else:
                    motif_ok = False
                    root_rows.append((motif, direction, region, released, mine, log10_gap(mine, released)))
            path = temp_dir / (motif + ".pVal." + direction + ".vs.bg.txt")
            if not path.is_file():
                failures.append("positional table absent: " + str(path))
                motif_ok = False
                continue
            regions, positions, released = io.read_positional_p(path)
            order = np.array([io.REGIONS.index(r) for r in regions], dtype=np.int8)
            if not (np.array_equal(order, region_index) and np.array_equal(positions, position)):
                failures.append("positional axis mismatch: " + str(path))
                motif_ok = False
                continue
            mine = np.asarray(pvalues[direction], dtype=np.float64)
            equal = mine == released
            n_position += int(equal.size)
            exact_position += int(equal.sum())
            for slot in np.flatnonzero(~equal)[:20]:
                position_rows.append((motif, direction, io.REGIONS[order[slot]], int(positions[slot]),
                                      float(released[slot]), float(mine[slot]),
                                      log10_gap(float(mine[slot]), float(released[slot]))))
            if equal.all():
                positional.append((path, "positional_pvalue"))
            else:
                motif_ok = False
        recorded = conversion.get(motif)
        sources = []
        for group in io.GROUPS:
            source = temp_dir / (motif + ".countDist." + group + ".txt")
            if recorded is None or not source.is_file():
                failures.append("countDist absent or not in conversion_manifest.tsv: " + str(source))
                motif_ok = False
            elif io.md5(source) != recorded[group]:
                failures.append("countDist changed since conversion: " + str(source))
                motif_ok = False
            else:
                sources.append((source, "countDist"))
        if motif_ok:
            verified += positional + sources
        if index % 25 == 0 or index == len(motifs):
            emit("  exact check {}/{} motifs, {:.1f}s".format(index, len(motifs), time.perf_counter() - started))

    root_pass = n_root > 0 and exact_root == n_root
    position_pass = n_position > 0 and exact_position == n_position
    passed = root_pass and position_pass and not failures

    probe = []
    if motifs:
        probe_motif = "QKI.ACTAAC[ACG]" if "QKI.ACTAAC[ACG]" in set(motifs) else motifs[0]
        with np.load(out_dir / (probe_motif + ".counts.npz"), allow_pickle=False) as data:
            bg = io.group_csr(data, "bg").toarray()
            fg = io.group_csr(data, "dn").toarray()
        _, _, pvalues = motif_pvalues(out_dir / (probe_motif + ".counts.npz"))
        window = int(np.argmin(pvalues["dn"]))
        for label, kwargs in (("scipy default", {}), ("use_continuity=False", {"use_continuity": False})):
            value = float(mannwhitneyu(fg[:, window], bg[:, window], alternative="greater", **kwargs).pvalue)
            probe.append((label, value, log10_gap(value, float(pvalues["dn"][window]))))

    with open(out_dir / "detail_root_comparison.tsv", "w", newline="\n", encoding="utf-8") as handle:
        handle.write("motif\tdirection\tregion\treleased_p\trecomputed_p\tabs_log10_gap\n")
        for row in root_rows:
            handle.write("\t".join(repr(v) if isinstance(v, float) else str(v) for v in row) + "\n")
    with open(out_dir / "detail_positional_mismatches.tsv", "w", newline="\n", encoding="utf-8") as handle:
        handle.write("motif\tdirection\tregion\tposition\treleased_p\trecomputed_p\tabs_log10_gap\n")
        for row in position_rows:
            handle.write("\t".join(repr(v) if isinstance(v, float) else str(v) for v in row) + "\n")

    lines = [
        "# Archive verification (exact) - " + args.arm, "",
        "Recomputed from `" + str(out_dir) + "` archives only; released tables read, never written.",
        "Statistic: tie-corrected normal approximation with continuity correction, one-sided "
        "changed > background, non-finite replaced by 1.0 (the released engine convention).",
        "Criterion: exact float equality for every value; no tolerance.", "",
        "| check | n compared | n exactly equal | verdict |", "|---|---|---|---|",
        "| root tables ({} motifs x 8 sub-regions x 2 directions) | {} | {} | {} |".format(
            len(motifs), n_root, exact_root, "PASS" if root_pass else "FAIL"),
        "| every per-position value of the temp pVal tables | {} | {} | {} |".format(
            n_position, exact_position, "PASS" if position_pass else "FAIL"),
        "| motif set; countDist md5 against conversion_manifest.tsv | - | - | {} |".format(
            "PASS" if not failures else "FAIL"), "",
    ]
    lines += ["- " + f for f in failures[:50]]
    lines += ["", "## Continuity setting of the released engine", "",
              "| scipy call | p | abs log10 gap to recomputed |", "|---|---|---|"]
    for label, value, gap in probe:
        lines.append("| mannwhitneyu(alternative='greater', {}) | {:.6e} | {:.3e} |".format(label, value, gap))
    lines += ["", "Deletion list: " + ("`{}` ({} files).".format(MANIFEST_NAME, len(verified)) if passed else
                                      "none written; verification did not pass, so nothing may be deleted."),
              "Detail: `detail_root_comparison.tsv`, `detail_positional_mismatches.tsv` (mismatches only).",
              "Generated " + time.strftime("%Y-%m-%dT%H:%M:%S") + "."]
    (out_dir / "VERIFY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    if passed:
        partial = manifest_path.with_name(MANIFEST_NAME + ".partial")
        with open(partial, "w", newline="\n", encoding="utf-8") as handle:
            handle.write("\t".join(MANIFEST_COLUMNS) + "\n")
            for path, kind in verified:
                handle.write("{}\t{}\t{}\t{}\n".format(path.resolve(), path.stat().st_size, io.md5(path), kind))
        os.replace(partial, manifest_path)

    with open(out_dir / "versions.txt", "a", encoding="utf-8") as handle:
        handle.write("verify_script\t{}\n".format(Path(__file__).resolve()))
        handle.write("verify_script_md5\t{}\n".format(io.md5(Path(__file__).resolve())))
        handle.write("verify_generated\t{}\n".format(time.strftime("%Y-%m-%dT%H:%M:%S")))
        handle.write("verify_numpy\t{}\n".format(np.__version__))
        handle.write("verify_scipy\t{}\n".format(scipy.__version__))
        handle.write("verify_platform\t{}\n".format(platform.platform()))
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            handle.write("verify_{}\t{}\n".format(var, os.environ.get(var, "unset")))

    emit("root: n={} exact={} verdict={}".format(n_root, exact_root, "PASS" if root_pass else "FAIL"))
    emit("positional: n={} exact={} verdict={}".format(n_position, exact_position,
                                                       "PASS" if position_pass else "FAIL"))
    for failure in failures[:10]:
        emit("FAIL " + failure)
    emit("verified temporaries listed: {}".format(len(verified) if passed else 0))
    emit("exit={}".format(0 if passed else 1))
    log.close()
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
