"""Orientation guard of tools/build_event_sets_lab_full.py: exact-zero ties are not disagreements.

Same test shape as tests/test_orientation_ties.py (the portable builder, fix 4e5d717), driven
through prepare(). The lab builder needs pandas; without it these tests skip. The fixture stops at
the DESeq2 table, which is read only after the orientation guard has passed.
"""
import csv
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("pandas")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import build_event_sets_lab_full as lab  # noqa: E402

TREATMENT, CONTROL = lab.EXPECTED["QKI_KO"]


def event(identifier, start, fdr="0.01", delta="0.10", inc1=None, inc2="0.5,0.5"):
    if inc1 is None:
        inc1 = "0.6,0.6" if float(delta) > 0 else "0.4,0.4"
    return {"ID": identifier, "GeneID": "ENSG000001.1", "geneSymbol": "TEST",
            "chr": "chr1", "strand": "+", "exonStart_0base": str(start),
            "exonEnd": str(start + 100), "upstreamES": str(start - 500),
            "upstreamEE": str(start - 400), "downstreamES": str(start + 500),
            "downstreamEE": str(start + 600), "FDR": fdr,
            "IncLevelDifference": delta, "IncLevel1": inc1, "IncLevel2": inc2,
            "IJC_SAMPLE_1": "10,10", "SJC_SAMPLE_1": "0,0",
            "IJC_SAMPLE_2": "10,10", "SJC_SAMPLE_2": "0,0"}


def tie(identifier, start, fdr="0.9"):
    """dPSI rounded to exactly 0.000 while the recomputed mean difference is nonzero."""
    return event(identifier, start, fdr=fdr, delta="0.000", inc1="0.5001,0.5001", inc2="0.5,0.5")


def run_prepare(tmp_path, rows):
    source = tmp_path / "SE.MATS.JC.txt"
    with source.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    for name, samples in (("b1.txt", TREATMENT), ("b2.txt", CONTROL)):
        (tmp_path / name).write_text(",".join(f"bams/{s}.markdup.sorted.bam" for s in samples))
    out = tmp_path / "out"
    out.mkdir()
    args = SimpleNamespace(arm="QKI_KO", b1=tmp_path / "b1.txt", b2=tmp_path / "b2.txt",
                           treatment_is="b1", rmats_se=source, outdir=out,
                           deseq2=tmp_path / "absent_deseq2.csv")
    try:
        lab.prepare(args)
    except FileNotFoundError as exc:  # guard passed; the fixture has no DESeq2 table
        assert "absent_deseq2.csv" in str(exc)
    return json.loads((out / "orientation.json").read_text(encoding="utf-8"))


def test_exact_zero_dpsi_ties_do_not_trip_the_orientation_guard(tmp_path):
    """60 ties beside 3 agreeing rows: the old evaluable denominator gave 3 / 63 and aborted."""
    rows = [event("up", 1000), event("dn", 3000, delta="-0.10"), event("bg", 5000, fdr="0.9")]
    rows += [tie(f"tie{i}", 10000 + 2000 * i) for i in range(60)]
    orientation = run_prepare(tmp_path, rows)
    assert orientation["evaluable"] == 63 and orientation["informative"] == 3
    assert orientation["agrees"] / orientation["evaluable"] < 0.95  # what the unfixed guard saw
    assert orientation["agreement"] == 1.0 == orientation["informative_agreement"]


def test_a_genuine_orientation_flip_is_still_refused(tmp_path):
    rows = [event("a", 1000, delta="0.10", inc1="0.2,0.2"),
            event("b", 3000, delta="0.20", inc1="0.1,0.1"),
            event("bg", 5000, fdr="0.9")]
    with pytest.raises(AssertionError, match="Orientation assertion failed"):
        run_prepare(tmp_path, rows)


def test_ties_alone_cannot_satisfy_the_guard(tmp_path):
    rows = [tie(f"tie{i}", 1000 + 2000 * i) for i in range(20)]
    with pytest.raises(AssertionError, match="Orientation assertion failed"):
        run_prepare(tmp_path, rows)
