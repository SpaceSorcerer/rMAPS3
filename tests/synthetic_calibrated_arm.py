"""A synthetic arm with both figure layers: released-layout root tables, count archives, a real
calibration v2.1 run at small B, and the figure builder's naming inputs. One motif carries no hit in
any event, so its calibrated cells are untestable (p NA) and the BH divisor is below the family."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import calibrate_ranksum as calib  # noqa: E402
import calibrate_ranksum_v2 as v2  # noqa: E402
import rmaps_countdist_io as io  # noqa: E402
import rmaps3_skill_run as skill  # noqa: E402
from test_skill_run import fork_motif_keys  # noqa: E402

ARM = "QKI_KO_T"
PER_REGION = 4
SIZES = {"up": [1] * 16 + [2] * 2, "dn": [1] * 14 + [2] * 2, "bg": [1] * 180 + [2] * 10}
OFFSETS = {"up": 10_000, "dn": 200_000, "bg": 1_000_000}
COUNTS = {"n_up": 20, "n_dn": 18, "n_bg": 200, "n_expr_unknown_in_fg": 0, "n_expr_unknown_in_bg": 2}
ALIAS = ROOT / "data" / "rbp_alias_hgnc_2026-09-17.tsv"
COMMIT = "abc1234"  # not the builder default, so propagation into every text is visible


def zero_motif(motifs):
    """A motif whose table name carries one motif only, so its RBP-level cells are untestable too."""
    alias = skill.load_alias(ALIAS)
    names = [alias.get(m.split(".", 1)[0], m.split(".", 1)[0]) for m in motifs]
    kmers = [m.split(".", 1)[1] for m in motifs]
    return next(m for m in motifs if names.count(alias.get(m.split(".", 1)[0], m.split(".", 1)[0])) == 1
                and kmers.count(m.split(".", 1)[1]) == 1 and not m.startswith(("QKI.", "motif_")))


def build(tmp_path: Path, b1=99, b2=499, commit=None):
    commit = commit or COMMIT
    motifs = sorted(fork_motif_keys())
    zero = zero_motif(motifs)
    released, counts_root, calibrated = tmp_path / "released", tmp_path / "counts", tmp_path / "summary"
    (released / ARM).mkdir(parents=True)
    (counts_root / ARM).mkdir(parents=True)
    rng = np.random.default_rng(149)
    n_windows = PER_REGION * len(io.REGIONS)
    region = np.repeat(np.arange(len(io.REGIONS), dtype=np.int8), PER_REGION)
    position = np.tile(np.arange(PER_REGION, dtype=np.int32), len(io.REGIONS))
    roots = {"up": {}, "dn": {}}
    by_kmer = {}
    for motif in motifs:
        store = {"schema_version": np.asarray(io.SCHEMA_VERSION), "arm": np.asarray(ARM),
                 "motif": np.asarray(motif), "region_name": np.asarray(io.REGIONS, dtype="U"),
                 "region_of_position": region, "position": position}
        kmer = motif.split(".", 1)[1]
        if kmer in by_kmer:  # keys sharing a k-mer get bit-identical counts, as the engine writes them
            store.update({k: v for k, v in by_kmer[kmer].items() if k != "motif"})
        for group, sizes in SIZES.items():
            if kmer in by_kmer:
                break
            rate = 0.0 if motif == zero else (0.9 if motif.startswith("QKI.") and group != "bg" else 0.2)
            rows, ids = [], []
            for c, size in enumerate(sizes):
                counts = rng.poisson(rate, n_windows).astype(np.int16)
                for r in range(size):
                    rows.append(counts)
                    ids.append("chr1:+:{}:{}:{}:1:2:3".format(OFFSETS[group] + 10 * c, OFFSETS[group] + 10 * c + 5, r))
            io.pack_group(store, group, np.vstack(rows))
            store[group + "_exon_id"] = np.asarray(ids, dtype="U")
        by_kmer.setdefault(kmer, dict(store))
        path = counts_root / ARM / (motif + ".counts.npz")
        np.savez_compressed(path, **store)
        for d in ("up", "dn"):
            model = calib.MotifModel(path, d)
            p = io.p_from_z(model.observed_z)
            roots[d][motif] = {r: float(p[model.region_index == i].min()) for i, r in enumerate(io.REGIONS)}
    for d in ("up", "dn"):
        with open(released / ARM / f"pVal.{d}.vs.bg.RNAmap.txt", "w", newline="\n", encoding="utf-8") as h:
            h.write("\t".join(["RBP"] + [io.ROOT_COLUMNS[r] for r in io.REGIONS]) + "\n")
            for motif in motifs:
                h.write("\t".join([motif] + [repr(roots[d][motif][r]) for r in io.REGIONS]) + "\n")
    (released / f"{ARM}_command.log").write_text(
        f"set={ARM} engine=released (synthetic @ {commit}) stat=mannwhitney\nexit=0 wall_s=1\n", encoding="utf-8")
    counts_json = tmp_path / "counts.json"
    counts_json.write_text(json.dumps(COUNTS), encoding="utf-8")
    common = ["--arm", ARM, "--counts-root", str(counts_root), "--released-root", str(released),
              "--alias-table", str(ALIAS), "--permutations", str(b1), "--refine-perms", str(b2),
              "--refine-threshold", "0.05"]
    rowunit = tmp_path / "summary_rowunit"
    assert calib.main(common + ["--out-root", str(rowunit), "--permutation-unit", "row"]) == 0
    assert v2.main(common + ["--out-root", str(calibrated), "--rowunit-root", str(rowunit)]) == 0
    alias = skill.load_alias(ALIAS)
    names = {m.split(".", 1)[0] for m in motifs if not m.startswith("motif_")}
    symbols = sorted({alias.get(n, n) for n in names} | set(alias.values()))
    gtf = tmp_path / "genes.gtf"
    gtf.write_text("".join(f'chr1\tT\tgene\t1\t2\t.\t+\t.\tgene_id "ENSG{i:011d}.1"; gene_name "{s}";\n'
                           for i, s in enumerate(symbols)), encoding="utf-8")
    (tmp_path / "splice.txt").write_text("SRSF1\nSNRPA\n", encoding="utf-8")
    (tmp_path / "broad.txt").write_text("PTBP1\nHNRNPC\n", encoding="utf-8")
    return {"released": released, "counts": counts_root, "calibrated": calibrated, "counts_json": counts_json,
            "gtf": gtf, "splice": tmp_path / "splice.txt", "broad": tmp_path / "broad.txt", "zero": zero,
            "motifs": motifs}


def draw(arm_inputs, figures: Path, extra=()):
    import build_region_lollipops_v4 as lol
    lol.main(["--arms", ARM, "--out-root", str(figures), "--released-root", str(arm_inputs["released"]),
              "--calibrated-root", str(arm_inputs["calibrated"]), "--counts-json", f"{ARM}={arm_inputs['counts_json']}",
              "--alias-table", str(ALIAS), "--gtf", str(arm_inputs["gtf"]),
              "--spliceosome-list", str(arm_inputs["splice"]), "--broad-binders-list", str(arm_inputs["broad"]),
              "--released-commit", COMMIT, "--gate-rule", f"{ARM}=A", "--arm-label", f"{ARM}=QKI knockout test arm",
              "--positive-control", "QKI", *extra])
    return figures
