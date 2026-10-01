import csv
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("builder", ROOT / "tools" / "build_event_sets.py")
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


def event(identifier, start, fdr="0.01", delta="0.10", coverage="10,10"):
    return {"ID": identifier, "GeneID": "ENSG000001.1", "geneSymbol": "TEST",
            "chr": "chr1", "strand": "+", "exonStart_0base": str(start),
            "exonEnd": str(start + 100), "upstreamES": str(start - 500),
            "upstreamEE": str(start - 400), "downstreamES": str(start + 500),
            "downstreamEE": str(start + 600), "FDR": fdr,
            "IncLevelDifference": delta, "IncLevel1": "0.6,0.6" if float(delta) > 0 else "0.4,0.4",
            "IncLevel2": "0.5,0.5", "IJC_SAMPLE_1": coverage, "SJC_SAMPLE_1": "0,0",
            "IJC_SAMPLE_2": "10,10", "SJC_SAMPLE_2": "0,0"}


def write_rmats(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_frozen_rule_a_boundaries_and_collision_exclusion(tmp_path):
    rows = [event("up", 1000), event("dn", 3000, delta="-0.1"),
            event("fdr_boundary", 5000, fdr="0.05"), event("low_dpsi", 7000, delta="0.099"),
            event("low_coverage", 9000, coverage="10,9"), event("bg", 11000, fdr="0.5"),
            event("collision", 1000, fdr="0.9")]
    source = tmp_path / "SE.MATS.JC.txt"
    write_rmats(source, rows)
    sets, counts, audit = builder.build_sets(source, {"base_mean_floor": 50, "expr_unknown": "retain"}, {})
    assert [row["ID"] for row in sets["up"]] == ["up"]
    assert [row["ID"] for row in sets["dn"]] == ["dn"]
    assert [row["ID"] for row in sets["bg"]] == ["bg"]
    assert counts["bg_collisions_removed"] == 1
    assert counts["n_expr_unknown_in_fg"] == 2
    assert len(audit) == 7


@pytest.mark.parametrize("mean,policy,expected", [(50, "retain", 0), (50.01, "retain", 1),
                                               (None, "retain", 1), (None, "exclude", 0)])
def test_expression_floor_is_strict_and_unknown_policy_explicit(tmp_path, mean, policy, expected):
    source = tmp_path / "SE.MATS.JC.txt"
    write_rmats(source, [event("up", 1000)])
    expression = {} if mean is None else {"ENSG000001": mean}
    sets, _, _ = builder.build_sets(source, {"base_mean_floor": 50, "expr_unknown": policy}, expression)
    assert len(sets["up"]) == expected


def test_foreground_wrong_orientation_rejected(tmp_path):
    source = tmp_path / "SE.MATS.JC.txt"
    row = event("up", 1000)
    row["IncLevel1"] = "0.2,0.2"
    write_rmats(source, [row])
    with pytest.raises(ValueError, match="orientation"):
        builder.build_sets(source, {}, {})


def test_background_window_collision_ignores_strand_and_exon_identity(tmp_path):
    source = tmp_path / "SE.MATS.JC.txt"
    up = event("up", 1000)
    bg = event("bg_collision", 1200, fdr="0.9")
    bg.update(strand="-", upstreamES="500", upstreamEE="600", downstreamES="1500", downstreamEE="1600")
    write_rmats(source, [up, bg, event("bg_kept", 5000, fdr="0.5")])
    sets, counts, _ = builder.build_sets(source, {}, {})
    assert [row["ID"] for row in sets["bg"]] == ["bg_kept"]
    assert counts["bg_collisions_removed"] == 1


def test_repeated_rmats_id_columns_must_agree(tmp_path):
    source = tmp_path / "SE.MATS.JC.txt"
    row = event("up", 1000)
    source.write_text("\t".join([*row, "ID"]) + "\n" + "\t".join([*row.values(), "different"]) + "\n")
    with pytest.raises(ValueError, match="Repeated rMATS ID columns disagree"):
        builder.build_sets(source, {}, {})
