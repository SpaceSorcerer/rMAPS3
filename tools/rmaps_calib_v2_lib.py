"""Library for calibration v2 of the rMAPS3 rank-sum supplement (tools/calibrate_ranksum_v2.py).

Target-exon cluster permutation (rows sharing chr, strand, exonStart, exonEnd move together, at
the observed foreground's exact cluster-size composition), RBP-level combination of per-motif
statistics inside the permutation (min-P primary; max-z and mean-z sensitivity), and the
unique-k-mer motif family. Pure functions only; the arm runner lives in calibrate_ranksum_v2.py.
Ported from the 2026-09-22 lab scripts without changing any arithmetic or draw stream.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import calibrate_ranksum as calib  # noqa: E402
import rmaps_countdist_io as io  # noqa: E402


# ------------------------------------------------------------------ target-exon clusters
def target_keys(ids: np.ndarray) -> np.ndarray:
    """chr:strand:exonStart:exonEnd of each eight-field engine exon identifier."""
    return np.asarray([":".join(s.split(":")[:4]) for s in ids], dtype="U")


def cluster_rows(keys: np.ndarray):
    """{target exon key: row indices}, in first-appearance order."""
    groups = {}
    for row, key in enumerate(keys):
        groups.setdefault(key, []).append(row)
    return {k: np.asarray(v, dtype=np.int64) for k, v in groups.items()}


class ClusterDrawer:
    """Size-matched cluster-label permutation over the pooled row axis.

    Draws the observed foreground's exact cluster-size composition from the pooled cluster pool,
    so every draw carries exactly n1 rows and no target exon is ever split across the
    changed/background boundary.
    """

    def __init__(self, fg_keys: np.ndarray, bg_keys: np.ndarray):
        self.n1 = int(fg_keys.shape[0])
        self.n_total = self.n1 + int(bg_keys.shape[0])
        pooled_clusters = cluster_rows(np.concatenate([fg_keys, bg_keys]))
        self.n_clusters_pooled = len(pooled_clusters)
        by_size = {}
        for rows in pooled_clusters.values():
            by_size.setdefault(int(rows.shape[0]), []).append(rows)
        self.pool = {s: np.vstack(v).astype(np.int32) for s, v in by_size.items()}
        fg_clusters = cluster_rows(fg_keys)
        self.n_clusters_fg = len(fg_clusters)
        need = {}
        for rows in fg_clusters.values():
            need[int(rows.shape[0])] = need.get(int(rows.shape[0]), 0) + 1
        self.need = sorted(need.items())
        for size, count in self.need:
            have = 0 if size not in self.pool else int(self.pool[size].shape[0])
            if have < count:
                raise ValueError("pooled cluster pool too small at size {}: need {}, have {}"
                                 .format(size, count, have))
        self.max_cluster_fg = max(s for s, _ in self.need)
        self.overlap = len(set(fg_keys.tolist()) & set(bg_keys.tolist()))
        if self.overlap:
            raise ValueError("target exon present in both foreground and background")

    def draw(self, rng, batch: int) -> np.ndarray:
        out = np.empty((batch, self.n1), dtype=np.int32)
        for i in range(batch):
            start = 0
            for size, count in self.need:
                pool = self.pool[size]
                picked = pool[rng.choice(pool.shape[0], count, replace=False)].ravel()
                out[i, start: start + picked.shape[0]] = picked
                start += int(picked.shape[0])
            if start != self.n1:
                raise ValueError("cluster draw cardinality mismatch")
        return out


def stage_rng(seed: int, stage: int, direction: str):
    """Independent draw stream per (seed, stage, direction)."""
    return np.random.default_rng(np.random.SeedSequence(
        [seed, stage, io.DIRECTION_CODE[direction]]))


def exon_ids_for_counts_dir(counts_dir: Path, verify_all_motifs: bool = True):
    """Eight-field exon identifiers per group, in archive row order; identical across motifs."""
    counts_dir = Path(counts_dir)
    motifs = sorted(p.name[: -len(".counts.npz")] for p in counts_dir.glob("*.counts.npz"))
    if not motifs:
        raise ValueError("no *.counts.npz archives in " + str(counts_dir))
    reference = None
    for motif in motifs:
        with np.load(counts_dir / (motif + ".counts.npz"), allow_pickle=False) as data:
            ids = {g: np.asarray(data[g + "_exon_id"]) for g in io.GROUPS}
        if reference is None:
            reference = ids
        else:
            for g in io.GROUPS:
                if not np.array_equal(ids[g], reference[g]):
                    raise ValueError("exon axis differs between motifs: " + motif)
        if not verify_all_motifs:
            break
    return reference, motifs


# ------------------------------------------------------------------ motif families
def kmer_of(motif_key: str) -> str:
    return motif_key.split(".", 1)[1]


def table_name_of(motif_key: str) -> str:
    return motif_key.split(".", 1)[0]


def unique_motif_groups(motifs):
    """{k-mer: sorted carrier keys}; the representative is the first carrier."""
    groups = {}
    for key in sorted(motifs):
        groups.setdefault(kmer_of(key), []).append(key)
    return groups


def rbp_unique_kmers(motifs, alias):
    """{RBP: [unique k-mers in key order]} under the alias grouping."""
    out = {}
    for key in sorted(motifs):
        rbp = alias.get(table_name_of(key), table_name_of(key))
        kmers = out.setdefault(rbp, [])
        if kmer_of(key) not in kmers:
            kmers.append(kmer_of(key))
    return out


# ------------------------------------------------------------------ permutation statistics
def permutation_p(observed: float, null: np.ndarray) -> float:
    """(1 + #{null >= observed}) / (1 + B), the v1 convention."""
    if not math.isfinite(observed):
        return math.nan
    finite = np.isfinite(null)
    return (1 + int(np.sum(null[finite] >= observed))) / (1 + null.shape[0])


def rbp_combine(observed: np.ndarray, null: np.ndarray):
    """Max and mean over motifs, computed on observed and null columns alike.

    observed: (m,) pooled regional-max z of the RBP's motifs; null: (m, B) under B shared
    permutations. Returns (obs_max, null_max, obs_mean, null_mean, argmax_row).
    """
    keep = np.isfinite(observed)
    if not keep.any():
        empty = np.full(null.shape[1], np.nan)
        return math.nan, empty, math.nan, empty, None
    block = np.column_stack([observed[keep], null[keep]])
    if not np.all(np.isfinite(block)):
        raise ValueError("non-finite permutation maximum for a motif with finite observed z")
    maxima = block.max(axis=0)
    means = block.mean(axis=0)
    argmax_row = int(np.flatnonzero(keep)[int(np.argmax(observed[keep]))])
    return float(maxima[0]), maxima[1:], float(means[0]), means[1:], argmax_row


def exceed_counts(values: np.ndarray) -> np.ndarray:
    """For each entry, how many entries of the same vector are >= it (itself included)."""
    ordered = np.sort(values)
    return values.shape[0] - np.searchsorted(ordered, values, side="left")


def rbp_minp(observed: np.ndarray, null: np.ndarray):
    """Westfall-Young min-P over motifs, on the observed and null columns alike.

    In each column (observed + B permutations) each motif's statistic becomes its own tail count
    #{columns with z >= this z} out of B + 1; the RBP statistic per column is the minimum count
    over motifs. Calibrated p = #{columns with min count <= observed min count} / (B + 1), in
    integers. With one motif it equals the motif's own permutation p exactly.
    Returns (p, row of the motif with the smallest observed p).
    """
    keep = np.isfinite(observed)
    if not keep.any():
        return math.nan, None
    block = np.column_stack([observed[keep], null[keep]])
    if not np.all(np.isfinite(block)):
        raise ValueError("non-finite permutation maximum for a motif with finite observed z")
    counts = np.vstack([exceed_counts(row) for row in block])
    minimum = counts.min(axis=0)
    p = int(np.sum(minimum <= minimum[0])) / block.shape[1]
    best = int(np.lexsort((-observed[keep], counts[:, 0]))[0])
    return p, int(np.flatnonzero(keep)[best])


def calibrate_keep(model, masks, selection: np.ndarray, chunk: int):
    """calibrate_ranksum.calibrate plus the per-permutation maxima, for RBP-level reuse."""
    observed = calib.statistic_maxima(model.observed_z, masks)
    total = selection.shape[0]
    maxima_all = {name: np.empty(total, dtype=np.float64) for name in masks}
    for start in range(0, total, chunk):
        batch = selection[start: start + chunk]
        maxima = calib.statistic_maxima(model.permutation_z(batch), masks)
        for name in masks:
            maxima_all[name][start: start + batch.shape[0]] = maxima[name]
    result = {}
    for name in masks:
        values = maxima_all[name]
        finite = np.isfinite(values)
        result[name] = {
            "observed_max_z": float(observed[name]),
            "observed_min_p": float(io.p_from_z(np.asarray(observed[name])))
            if math.isfinite(observed[name]) else math.nan,
            "calibrated_p": permutation_p(float(observed[name]), values),
            "permutations": total,
            "usable_permutations": int(finite.sum()),
            "reason": "" if math.isfinite(observed[name]) else "no window with usable variance",
        }
    return result, maxima_all


def panel_rank(records, q_key, p_key, tie_key):
    """Rank by q, then p, then the tie key (smaller first), then RBP name."""
    def value(x):
        return x if (x is not None and math.isfinite(x)) else math.inf
    order = sorted(records, key=lambda r: (value(r[q_key]), value(r[p_key]),
                                           value(r[tie_key]), r["RBP"]))
    return {r["RBP"]: i for i, r in enumerate(order, 1)}


# ------------------------------------------------------------------ small IO
def read_tsv(path: Path):
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in handle]


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan
