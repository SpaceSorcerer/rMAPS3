"""Convert the released rMAPS3 countDist text tables of one arm to npz archives.

One archive per motif holds the per-exon x per-position motif hit COUNT matrix
of the three groups (up, dn, bg), the exon identifiers of each group in engine
order, and the position axis (region label plus within-region index).
"""

from __future__ import annotations

import argparse
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rmaps_countdist_io as io  # noqa: E402


def convert_motif(temp_dir: Path, fasta_dir: Path, motif: str, arm: str,
                  exon_axis: dict, out_path: Path) -> dict:
    store = {
        "schema_version": np.asarray(io.SCHEMA_VERSION),
        "arm": np.asarray(arm),
        "motif": np.asarray(motif),
        "region_name": np.asarray(io.REGIONS, dtype="U"),
    }
    stats = {"motif": motif}
    reference_axis = None
    for group in io.GROUPS:
        source = temp_dir / (motif + ".countDist." + group + ".txt")
        region_index, position, matrix = io.parse_countdist(source)
        if reference_axis is None:
            store["region_of_position"] = region_index
            store["position"] = position
            reference_axis = (region_index, position)
        elif not (np.array_equal(reference_axis[0], region_index)
                  and np.array_equal(reference_axis[1], position)):
            raise ValueError("position axis differs between groups for " + motif)
        ids = exon_axis[group]
        if matrix.shape[0] != ids.size:
            raise ValueError(
                "exon count mismatch for {} {}: countDist {} vs fasta {}".format(
                    motif, group, matrix.shape[0], ids.size))
        store[group + "_exon_id"] = ids
        io.pack_group(store, group, matrix)
        stats[group + "_n_exons"] = int(matrix.shape[0])
        stats[group + "_nnz"] = int(np.count_nonzero(matrix))
        stats[group + "_max"] = int(matrix.max())
        stats[group + "_zero_fraction"] = float(store[group + "_zero_fraction"])
        stats[group + "_source_md5"] = io.md5(source)
    np.savez_compressed(out_path, **store)
    stats["npz_bytes"] = out_path.stat().st_size
    return stats


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", required=True)
    parser.add_argument("--released-root", required=True,
                        help="root holding <arm>/temp/*.countDist.*.txt and <arm>/fasta/")
    parser.add_argument("--out-root", required=True)
    parser.add_argument("--limit", type=int, default=0, help="first N motifs only (smoke test)")
    args = parser.parse_args(argv)

    arm_dir = Path(args.released_root) / args.arm
    temp_dir = arm_dir / "temp"
    fasta_dir = arm_dir / "fasta"
    out_dir = Path(args.out_root) / args.arm
    out_dir.mkdir(parents=True, exist_ok=True)

    log_path = out_dir / "command.log"
    log = open(log_path, "a", encoding="utf-8")

    def emit(message):
        print(message, flush=True)
        log.write(message + "\n")
        log.flush()

    emit("start={} cmd={}".format(time.strftime("%Y-%m-%dT%H:%M:%S"),
                                  " ".join([sys.executable] + sys.argv)))
    exon_axis = {g: io.read_exon_axis(fasta_dir, g) for g in io.GROUPS}
    emit("exon axis: " + ", ".join("{}={}".format(g, exon_axis[g].size) for g in io.GROUPS))

    motifs = io.motif_keys(temp_dir)
    if args.limit:
        motifs = motifs[: args.limit]
    emit("motifs={}".format(len(motifs)))

    manifest = out_dir / "conversion_manifest.tsv"
    columns = ["motif", "npz_path", "npz_bytes"]
    for group in io.GROUPS:
        columns += [group + "_n_exons", group + "_nnz", group + "_max",
                    group + "_zero_fraction", group + "_source_md5"]
    started = time.perf_counter()
    with open(manifest, "w", newline="\n", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for index, motif in enumerate(motifs, 1):
            out_path = out_dir / (motif + ".counts.npz")
            stats = convert_motif(temp_dir, fasta_dir, motif, args.arm, exon_axis, out_path)
            stats["npz_path"] = str(out_path)
            handle.write("\t".join(str(stats[c]) for c in columns) + "\n")
            handle.flush()
            if index % 10 == 0 or index == len(motifs):
                emit("  {}/{} motifs, {:.1f}s".format(index, len(motifs),
                                                      time.perf_counter() - started))
    total_bytes = sum(p.stat().st_size for p in out_dir.glob("*.counts.npz"))
    emit("wrote {} archives, {:.2f} GB, {:.1f}s".format(
        len(motifs), total_bytes / 1e9, time.perf_counter() - started))

    with open(out_dir / "versions.txt", "w", encoding="utf-8") as handle:
        handle.write("python\t{}\n".format(sys.version.split()[0]))
        handle.write("numpy\t{}\n".format(np.__version__))
        handle.write("scipy\t{}\n".format(scipy.__version__))
        handle.write("platform\t{}\n".format(platform.platform()))
        handle.write("script\t{}\n".format(Path(__file__).resolve()))
        handle.write("script_md5\t{}\n".format(io.md5(Path(__file__).resolve())))
        handle.write("module_md5\t{}\n".format(
            io.md5(Path(__file__).resolve().parent / "rmaps_countdist_io.py")))
        handle.write("generated\t{}\n".format(time.strftime("%Y-%m-%dT%H:%M:%S")))
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            handle.write("{}\t{}\n".format(var, os.environ.get(var, "unset")))
    emit("exit=0")
    log.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
