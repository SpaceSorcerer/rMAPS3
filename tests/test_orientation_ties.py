"""Orientation guard: exact-zero IncLevelDifference ties must not count as disagreements.

Drop-in addition for E:\\Claude\\rMAPS3\\tests\\test_event_set_builder.py. Follows that
file's existing `event` / `write_rmats` helpers and its importlib loading of the builder.

Set RMAPS3_BUILDER to point at a specific build_event_sets.py; otherwise the repo's own
tools/build_event_sets.py is used.
"""
import csv
import importlib.util
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BUILDER_PATH = Path(os.environ.get("RMAPS3_BUILDER", ROOT / "tools" / "build_event_sets.py"))
# the builder imports its siblings (rmaps3_lab_run) by module name
TOOLS = Path(os.environ.get("RMAPS3_TOOLS", BUILDER_PATH.parent))
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))
spec = importlib.util.spec_from_file_location("builder", BUILDER_PATH)
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


def event(identifier, start, fdr="0.01", delta="0.10", coverage="10,10",
          inc1=None, inc2="0.5,0.5"):
    if inc1 is None:
        inc1 = "0.6,0.6" if float(delta) > 0 else "0.4,0.4"
    return {"ID": identifier, "GeneID": "ENSG000001.1", "geneSymbol": "TEST",
            "chr": "chr1", "strand": "+", "exonStart_0base": str(start),
            "exonEnd": str(start + 100), "upstreamES": str(start - 500),
            "upstreamEE": str(start - 400), "downstreamES": str(start + 500),
            "downstreamEE": str(start + 600), "FDR": fdr,
            "IncLevelDifference": delta, "IncLevel1": inc1, "IncLevel2": inc2,
            "IJC_SAMPLE_1": coverage, "SJC_SAMPLE_1": "0,0",
            "IJC_SAMPLE_2": "10,10", "SJC_SAMPLE_2": "0,0"}


def write_rmats(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def tie(identifier, start, fdr="0.9"):
    """dPSI rounded to exactly 0.000 while the recomputed mean difference is nonzero.

    This is what rMATS emits for a large cohort whose true delta is below 0.0005.
    """
    return event(identifier, start, fdr=fdr, delta="0.000",
                 inc1="0.5001,0.5001", inc2="0.5,0.5")


def build(tmp_path, rows):
    source = tmp_path / "SE.MATS.JC.txt"
    write_rmats(source, rows)
    return builder.build_sets(source, {}, {})


def test_exact_zero_dpsi_ties_do_not_trip_the_orientation_guard(tmp_path):
    """Ties must be excluded from the denominator, not scored as disagreements.

    One real up event, one real down event, one clean background, and 60 ties. Before the
    fix the ties dominate the evaluable denominator (3 of 63 agree = 0.048) and the build
    aborts even though every informative row agrees.
    """
    rows = [event("up", 1000), event("dn", 3000, delta="-0.10"),
            event("bg", 5000, fdr="0.9")]
    rows += [tie(f"tie{i}", 10000 + 2000 * i) for i in range(60)]
    sets, counts, audit = build(tmp_path, rows)
    assert [r["ID"] for r in sets["up"]] == ["up"]
    assert [r["ID"] for r in sets["dn"]] == ["dn"]
    assert sum(1 for r in audit if not r["informative"]) == 60
    # every tie still reaches the background; exclusion is from the DENOMINATOR only
    assert len(sets["bg"]) == 61


def test_a_genuine_orientation_flip_is_still_refused(tmp_path):
    """The guard must keep firing when nonzero-dPSI rows really do disagree."""
    rows = [event("a", 1000, delta="0.10", inc1="0.2,0.2"),
            event("b", 3000, delta="0.20", inc1="0.1,0.1"),
            event("bg", 5000, fdr="0.9")]
    with pytest.raises(ValueError, match="orientation"):
        build(tmp_path, rows)


def test_ties_alone_cannot_satisfy_the_guard(tmp_path):
    """With no informative row at all the build must refuse, not pass vacuously."""
    rows = [tie(f"tie{i}", 1000 + 2000 * i, fdr="0.9") for i in range(20)]
    with pytest.raises(ValueError):
        build(tmp_path, rows)
