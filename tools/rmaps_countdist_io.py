"""Shared IO and rank-sum kernel for the released rMAPS3 mannwhitney layer.

The released engine (legacy/motifMapSE_MP.py) writes, per motif and group, a
countDist text table of per-exon motif hit COUNTS for every window of the eight
sub-regions, then calls scipy.stats.mannwhitneyu(changed, background,
alternative='greater') with scipy defaults (method='asymptotic' in this regime,
use_continuity=True) and replaces any non-finite p with 1.0.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import norm

SCHEMA_VERSION = 1

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
POOL_TO_REGIONS = {
    "Upstream Intron": ("UpstreamExonIntron", "UpstreamIntron"),
    "Exon Body": ("TargetExon_5prime", "TargetExon-3prime"),
    "Downstream Intron": ("DownstreamIntron", "DownstreamExonIntron"),
    "Flanking Exon": ("UpstreamExon_3prime", "DownstreamExon_5prime"),
}
REGION_TO_POOL = {r: p for p, rs in POOL_TO_REGIONS.items() for r in rs}
PLOT_POOLS = ("Upstream Intron", "Exon Body", "Downstream Intron")
DIRECTION_LABEL = {"up": "INCLUDED", "dn": "SKIPPED"}
DIRECTION_CODE = {"up": 0, "dn": 1}
GROUPS = ("up", "dn", "bg")
FASTA_FOR_REGION = {
    "UpstreamExon_3prime": "UpstreamExon",
    "UpstreamExonIntron": "UpstreamExonIntron",
    "UpstreamIntron": "UpstreamIntron",
    "TargetExon_5prime": "TargetExon",
    "TargetExon-3prime": "TargetExon",
    "DownstreamIntron": "DownstreamIntron",
    "DownstreamExonIntron": "DownstreamExonIntron",
    "DownstreamExon_5prime": "DownstreamExon",
}
SPARSE_ZERO_FRACTION = 0.90


def md5(path) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_countdist(path):
    """Return (region_of_position, position, matrix[n_exons, n_positions])."""
    regions, positions, rows = [], [], []
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header != ["Region", "position", "sum", "values"]:
            raise ValueError("unexpected countDist header in " + str(path))
        for line in handle:
            region, position, total, values = line.rstrip("\n").split("\t")
            vector = np.fromstring(values[1:-1], dtype=np.int16, sep=",")
            if int(total) != int(vector.sum()):
                raise ValueError("countDist sum mismatch in " + str(path))
            regions.append(region)
            positions.append(int(position))
            rows.append(vector)
    order = list(dict.fromkeys(regions))
    if tuple(order) != REGIONS:
        raise ValueError("unexpected region order in " + str(path))
    lengths = {len(v) for v in rows}
    if len(lengths) != 1:
        raise ValueError("ragged exon axis in " + str(path))
    region_index = np.array([REGIONS.index(r) for r in regions], dtype=np.int8)
    return region_index, np.asarray(positions, dtype=np.int32), np.vstack(rows).T


def read_exon_axis(fasta_dir: Path, group: str):
    """Exon identifiers for a group; verified identical across its seven fastas."""
    reference = None
    for name in sorted(set(FASTA_FOR_REGION.values())):
        ids = []
        with open(fasta_dir / (group + "." + name + ".fasta"), "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(">"):
                    ids.append(line[1:].strip())
        if reference is None:
            reference = ids
        elif ids != reference:
            raise ValueError("exon axis differs between fastas of group " + group + ": " + name)
    return np.asarray(reference, dtype="U")


def pack_group(store: dict, group: str, matrix: np.ndarray) -> None:
    zero_fraction = 1.0 - (np.count_nonzero(matrix) / matrix.size)
    store[group + "_shape"] = np.asarray(matrix.shape, dtype=np.int64)
    if zero_fraction > SPARSE_ZERO_FRACTION:
        csr = sparse.csr_matrix(matrix)
        store[group + "_format"] = np.asarray("csr")
        store[group + "_indptr"] = csr.indptr.astype(np.int64)
        store[group + "_indices"] = csr.indices.astype(np.int32)
        store[group + "_data"] = csr.data.astype(np.int16)
    else:
        store[group + "_format"] = np.asarray("dense")
        store[group + "_dense"] = matrix.astype(np.int16)
    store[group + "_zero_fraction"] = np.asarray(zero_fraction)


def unpack_group(data, group: str) -> np.ndarray:
    shape = tuple(int(x) for x in data[group + "_shape"])
    if str(data[group + "_format"]) == "csr":
        return sparse.csr_matrix(
            (data[group + "_data"], data[group + "_indices"], data[group + "_indptr"]),
            shape=shape,
        ).toarray().astype(np.int16)
    return np.asarray(data[group + "_dense"], dtype=np.int16)


def group_csr(data, group: str):
    """Counts of one group as a CSR matrix (exons x positions), never densified."""
    shape = tuple(int(x) for x in data[group + "_shape"])
    if str(data[group + "_format"]) == "csr":
        return sparse.csr_matrix(
            (np.asarray(data[group + "_data"], dtype=np.int16),
             np.asarray(data[group + "_indices"], dtype=np.int32),
             np.asarray(data[group + "_indptr"], dtype=np.int64)), shape=shape)
    return sparse.csr_matrix(np.asarray(data[group + "_dense"], dtype=np.int16))


def count_histograms(matrix: np.ndarray, kmax: int) -> np.ndarray:
    hist = np.empty((kmax + 1, matrix.shape[1]), dtype=np.int64)
    for k in range(kmax + 1):
        hist[k] = (matrix == k).sum(axis=0)
    return hist


def csr_histograms(csr, kmax: int) -> np.ndarray:
    """Per-window frequency of each count value 0..kmax, from CSR storage."""
    n_exons, n_windows = csr.shape
    hist = np.zeros((kmax + 1, n_windows), dtype=np.int64)
    columns = np.asarray(csr.indices, dtype=np.int64)
    values = np.asarray(csr.data, dtype=np.int64)
    for k in range(1, kmax + 1):
        hist[k] = np.bincount(columns[values == k], minlength=n_windows)
    hist[0] = n_exons - hist[1:].sum(axis=0)
    return hist


def rank_statistics(hist_fg: np.ndarray, hist_bg: np.ndarray):
    """Midranks, mu, sigma and observed z of the released one-sided rank test.

    Everything but the rank sum is invariant under permutation of the
    changed/background labels, because the pooled multiset of counts at a window
    is fixed by the permutation.
    """
    fg = np.asarray(hist_fg, dtype=np.float64)
    bg = np.asarray(hist_bg, dtype=np.float64)
    n1 = fg.sum(axis=0)
    n0 = bg.sum(axis=0)
    tied = fg + bg
    cumulative_below = np.cumsum(tied, axis=0) - tied
    midrank = cumulative_below + (tied + 1.0) / 2.0
    rank_sum = (fg * midrank).sum(axis=0)
    total = n1 + n0
    mu = n1 * n0 / 2.0
    tie_term = (tied ** 3 - tied).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        variance = n1 * n0 / 12.0 * ((total + 1.0) - tie_term / (total * (total - 1.0)))
    sigma = np.sqrt(np.where(variance > 0, variance, np.nan))
    z = z_from_rank_sum(rank_sum, n1, mu, sigma)
    return midrank, mu, sigma, z, n1


def z_from_rank_sum(rank_sum, n1, mu, sigma):
    """z of scipy.stats.mannwhitneyu(alternative='greater', use_continuity=True)."""
    u1 = rank_sum - n1 * (n1 + 1.0) / 2.0
    with np.errstate(invalid="ignore", divide="ignore"):
        return (u1 - mu - 0.5) / sigma


def p_from_z(z):
    """Released convention: a window with no usable variance reports p = 1.0."""
    p = norm.sf(z)
    return np.where(np.isfinite(p), p, 1.0)


def read_root_table(path: Path):
    table = {}
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        expected = ["RBP"] + [ROOT_COLUMNS[r] for r in REGIONS]
        if header != expected:
            raise ValueError("unexpected root table header in " + str(path))
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            table[fields[0]] = {r: float(v) for r, v in zip(REGIONS, fields[1:])}
    return table


def read_positional_p(path: Path):
    """Per-position released p, in file order (region-major, position in region)."""
    regions, positions, values = [], [], []
    with open(path, "r", encoding="utf-8") as handle:
        handle.readline()
        for line in handle:
            region, position, value = line.rstrip("\n").split("\t")
            regions.append(region)
            positions.append(int(position))
            values.append(float(value))
    return regions, np.asarray(positions, dtype=np.int32), np.asarray(values, dtype=np.float64)


def motif_keys(temp_dir: Path):
    suffix = ".countDist.bg.txt"
    return sorted(p.name[: -len(suffix)] for p in temp_dir.glob("*" + suffix))
