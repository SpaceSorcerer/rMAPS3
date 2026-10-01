"""Exact archive verification and manifest-restricted temp deletion (tools/verify_ranksum_archives.py,
tools/rmaps3_skill_run.py convert_and_verify)."""
import argparse
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import rmaps3_skill_run as skill  # noqa: E402
import rmaps_countdist_io as io  # noqa: E402
import verify_ranksum_archives as verify  # noqa: E402

ARM = "SYN_A"
MOTIFS = ["QKI.ACTAAC[ACG]", "OTHER.TTTT"]
PER_REGION = 3
SIZES = {"up": 9, "dn": 8, "bg": 40}


def released_run(tmp_path):
    """A released-engine layout: root tables, temp countDist + positional pVal tables, fasta headers."""
    engine_root = tmp_path / "engine"
    arm_dir = engine_root / ARM
    temp, fasta = arm_dir / "temp", arm_dir / "fasta"
    temp.mkdir(parents=True)
    fasta.mkdir()
    rng = np.random.default_rng(149)
    n_windows = PER_REGION * len(io.REGIONS)
    region_index = np.repeat(np.arange(len(io.REGIONS)), PER_REGION)
    for group, n in SIZES.items():
        ids = ["chr1:+:{}:{}:1:2:3:4".format(1000 * (1 + list(SIZES).index(group)) + i, 5000 + i) for i in range(n)]
        for name in sorted(set(io.FASTA_FOR_REGION.values())):
            (fasta / f"{group}.{name}.fasta").write_text("".join(f">{i}\nACGT\n" for i in ids), encoding="utf-8")
    roots = {"up": {}, "dn": {}}
    for motif in MOTIFS:
        matrices = {g: rng.poisson(0.5, size=(n, n_windows)).astype(np.int16) for g, n in SIZES.items()}
        for group, matrix in matrices.items():
            with open(temp / f"{motif}.countDist.{group}.txt", "w", newline="\n", encoding="utf-8") as h:
                h.write("Region\tposition\tsum\tvalues\n")
                for w in range(n_windows):
                    column = matrix[:, w]
                    h.write("{}\t{}\t{}\t[{}]\n".format(io.REGIONS[region_index[w]], w % PER_REGION,
                                                        int(column.sum()), ",".join(map(str, column))))
        kmax = int(max(m.max() for m in matrices.values()))
        for d in ("up", "dn"):
            _, _, _, z, _ = io.rank_statistics(io.count_histograms(matrices[d], kmax),
                                               io.count_histograms(matrices["bg"], kmax))
            p = io.p_from_z(z)
            with open(temp / f"{motif}.pVal.{d}.vs.bg.txt", "w", newline="\n", encoding="utf-8") as h:
                h.write("Region\tposition\tpValue\n")
                for w in range(n_windows):
                    h.write("{}\t{}\t{!r}\n".format(io.REGIONS[region_index[w]], w % PER_REGION, float(p[w])))
            roots[d][motif] = {r: float(p[region_index == i].min()) for i, r in enumerate(io.REGIONS)}
    write_roots(arm_dir, roots)
    return engine_root, arm_dir, roots


def write_roots(arm_dir, roots):
    for d in ("up", "dn"):
        with open(arm_dir / f"pVal.{d}.vs.bg.RNAmap.txt", "w", newline="\n", encoding="utf-8") as h:
            h.write("\t".join(["RBP"] + [io.ROOT_COLUMNS[r] for r in io.REGIONS]) + "\n")
            for motif in MOTIFS:
                h.write("\t".join([motif] + [repr(roots[d][motif][r]) for r in io.REGIONS]) + "\n")


def wrapper_args():
    return argparse.Namespace(arm=ARM, engine="released", stat_method="mannwhitney", keep_temp=False)


def run_convert(tmp_path, engine_root):
    out = tmp_path / "out"
    (out / "logs").mkdir(parents=True)
    return out, skill.convert_and_verify(wrapper_args(), engine_root / ARM, out, lambda m: None, {})


def listed_files(arm_dir):
    return sorted(p.name for p in (arm_dir / "temp").iterdir()
                  if ".countDist." in p.name or ".pVal." in p.name)


def test_exact_verification_deletes_only_the_listed_files(tmp_path):
    engine_root, arm_dir, _ = released_run(tmp_path)
    stray = arm_dir / "temp" / "STRAY.pVal.up.vs.bg.txt"      # matches the old name filter, never verified
    stray.write_text("Region\tposition\tpValue\n", encoding="utf-8")
    other = arm_dir / "temp" / "engine_scratch.txt"
    other.write_text("x", encoding="utf-8")
    out, (counts_root, result) = run_convert(tmp_path, engine_root)
    assert result["verified"] and result["deleted_files"] == len(MOTIFS) * 5
    manifest = verify.read_manifest(counts_root / ARM / verify.MANIFEST_NAME)
    assert len(manifest) == len(MOTIFS) * 5
    assert {kind for *_, kind in manifest} == {"countDist", "positional_pvalue"}
    assert stray.is_file(), "an unlisted temp file must survive"
    assert other.is_file()
    assert listed_files(arm_dir) == ["STRAY.pVal.up.vs.bg.txt"]
    report = (counts_root / ARM / "VERIFY.md").read_text(encoding="utf-8")
    assert "exact float equality" in report and "| {} | {} | PASS |".format(
        len(MOTIFS) * 2 * len(io.REGIONS), len(MOTIFS) * 2 * len(io.REGIONS)) in report


def test_one_ulp_root_mismatch_blocks_every_deletion(tmp_path):
    engine_root, arm_dir, roots = released_run(tmp_path)
    value = roots["dn"][MOTIFS[1]]["DownstreamIntron"]
    roots["dn"][MOTIFS[1]]["DownstreamIntron"] = float(np.nextafter(value, 1.0))
    write_roots(arm_dir, roots)
    before = listed_files(arm_dir)
    with pytest.raises(RuntimeError, match="archive verification FAILED"):
        run_convert(tmp_path, engine_root)
    assert listed_files(arm_dir) == before
    assert not (tmp_path / "out" / "counts" / ARM / verify.MANIFEST_NAME).exists()


def test_one_ulp_positional_mismatch_blocks_every_deletion(tmp_path):
    engine_root, arm_dir, _ = released_run(tmp_path)
    path = arm_dir / "temp" / f"{MOTIFS[0]}.pVal.up.vs.bg.txt"
    lines = path.read_text(encoding="utf-8").splitlines()
    region, position, value = lines[5].split("\t")
    lines[5] = "\t".join([region, position, repr(float(np.nextafter(float(value), 2.0)))])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    before = listed_files(arm_dir)
    with pytest.raises(RuntimeError, match="archive verification FAILED"):
        run_convert(tmp_path, engine_root)
    assert listed_files(arm_dir) == before


def test_a_missing_positional_table_makes_verification_partial_and_blocks_deletion(tmp_path):
    engine_root, arm_dir, _ = released_run(tmp_path)
    (arm_dir / "temp" / f"{MOTIFS[1]}.pVal.dn.vs.bg.txt").rename(arm_dir / "moved_out.txt")
    before = listed_files(arm_dir)
    with pytest.raises(RuntimeError, match="archive verification FAILED"):
        run_convert(tmp_path, engine_root)
    assert listed_files(arm_dir) == before


def test_a_file_changed_after_verification_is_not_deleted(tmp_path):
    engine_root, arm_dir, _ = released_run(tmp_path)
    counts_root = tmp_path / "counts"
    import countdist_to_npz
    assert countdist_to_npz.main(["--arm", ARM, "--released-root", str(engine_root),
                                  "--out-root", str(counts_root)]) == 0
    assert verify.main(["--arm", ARM, "--released-root", str(engine_root), "--counts-root", str(counts_root)]) == 0
    target = arm_dir / "temp" / f"{MOTIFS[0]}.countDist.bg.txt"
    target.write_text(target.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    before = listed_files(arm_dir)
    with pytest.raises(RuntimeError, match="changed since verification; nothing deleted"):
        skill.delete_verified_temporaries(counts_root / ARM / verify.MANIFEST_NAME, arm_dir / "temp",
                                          tmp_path / "inventory.tsv")
    assert listed_files(arm_dir) == before


def test_verification_is_skipped_and_nothing_deleted_off_the_rank_sum_layer(tmp_path):
    engine_root, arm_dir, _ = released_run(tmp_path)
    before = listed_files(arm_dir)
    out = tmp_path / "out"
    (out / "logs").mkdir(parents=True)
    args = argparse.Namespace(arm=ARM, engine="released", stat_method="fisher", keep_temp=False)
    _, result = skill.convert_and_verify(args, engine_root / ARM, out, lambda m: None, {})
    assert not result["verified"] and result["deleted_files"] == 0
    assert listed_files(arm_dir) == before
