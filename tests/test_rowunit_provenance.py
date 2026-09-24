"""Row-unit sensitivity columns: exact axes, recorded unit and source md5 (review finding 10)."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import calibrate_ranksum as calib  # noqa: E402
import calibrate_ranksum_v2 as v2  # noqa: E402
from test_two_stage_fdr import null_arm  # noqa: E402


def arm_with_rowunit(tmp_path):
    args = null_arm(tmp_path, "R_A", np.random.default_rng(7))
    args.rowunit_root = str(tmp_path / "rowunit")
    assert calib.main(["--arm", args.arm, "--counts-root", args.counts_root, "--released-root", args.released_root,
                       "--out-root", args.rowunit_root, "--alias-table", args.alias_table, "--permutations", "49",
                       "--refine-perms", "99", "--permutation-unit", "row"]) == 0
    return args, Path(args.rowunit_root) / args.arm


def rewrite(path, edit):
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(edit(lines)) + "\n", encoding="utf-8")


def run(args):
    return v2.run_arm(args, lambda m: None)


def test_matching_rowunit_is_carried_with_its_unit_and_md5_record(tmp_path):
    args, _ = arm_with_rowunit(tmp_path)
    _, per_motif, _, _, report, motifs = run(args)
    assert report["rowunit_permutation_unit"] == "row"
    assert report["rowunit_archives_md5_matched"] == len(motifs)
    assert set(report["rowunit_tables_md5"]) == {"per_motif_regions.tsv", "condensed_per_rbp.tsv"}
    assert all(np.isfinite(r["calibrated_p_rowunit"]) for r in per_motif)


def test_a_renamed_key_with_the_same_cardinality_is_refused(tmp_path):
    args, base = arm_with_rowunit(tmp_path)
    def rename(lines):
        assert "\tRA.AAAC\t" in lines[1]
        return [lines[0], lines[1].replace("\tRA.AAAC\t", "\tRZ.AAAC\t", 1)] + lines[2:]
    rewrite(base / "per_motif_regions.tsv", rename)
    with pytest.raises(ValueError, match="exactly this run's"):
        run(args)


def test_a_missing_permutation_unit_is_refused(tmp_path):
    args, base = arm_with_rowunit(tmp_path)

    def drop_unit(lines):
        header = lines[0].split("\t")
        i = header.index("permutation_unit")
        return ["\t".join(f for j, f in enumerate(line.split("\t")) if j != i) for line in lines]
    rewrite(base / "per_motif_regions.tsv", drop_unit)
    with pytest.raises(ValueError, match="permutation_unit 'row' on every row"):
        run(args)


def test_a_non_row_unit_in_the_report_is_refused(tmp_path):
    args, base = arm_with_rowunit(tmp_path)
    report = json.loads((base / "refinement_report.json").read_text(encoding="utf-8"))
    report["permutation_unit"] = "cluster"
    (base / "refinement_report.json").write_text(json.dumps(report), encoding="utf-8")
    with pytest.raises(ValueError, match="not 'row'"):
        run(args)


def test_missing_or_different_source_md5_is_refused(tmp_path):
    args, base = arm_with_rowunit(tmp_path)
    rewrite(base / "input_md5s.tsv", lambda lines: [l if ".counts.npz" not in l or "RB.CCCA" not in l
                                                    else l.rsplit("\t", 1)[0] + "\t" + "0" * 32 for l in lines])
    with pytest.raises(ValueError, match="does not record this run's count archives"):
        run(args)
    (base / "input_md5s.tsv").unlink()
    with pytest.raises(ValueError, match="provenance incomplete"):
        run(args)
