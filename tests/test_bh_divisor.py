"""BH divisor actually used, and untestable cells, from calibration to figure (review finding 5)."""
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import false_discovery_control

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_region_lollipops_v4 as lol  # noqa: E402
import calibrate_ranksum as calib  # noqa: E402
import synthetic_calibrated_arm as syn  # noqa: E402


@pytest.fixture
def arm(calibrated_arm):
    return calibrated_arm


def test_bh_adjust_reports_the_divisor_it_used_and_leaves_nan_untouched():
    values = [0.01, math.nan, 0.04, 0.5, math.inf]
    q, divisor = calib.bh_adjust(values, return_divisor=True)
    assert divisor == 3
    assert np.isnan(q[1]) and np.isnan(q[4])
    np.testing.assert_allclose(q[[0, 2, 3]], false_discovery_control([0.01, 0.04, 0.5], method="bh"))
    assert isinstance(calib.bh_adjust(values), np.ndarray), "the one-argument form is unchanged"


def test_calibration_records_divisor_and_flags_untestable_cells(arm):
    base = arm["calibrated"] / syn.ARM
    report = json.loads((base / "refinement_report.json").read_text(encoding="utf-8"))
    untestable_motif = 2 * 3                      # one unique motif x 2 directions x 3 plotted pools
    assert report["motif_bh_divisor"] == report["motif_family_size"] - untestable_motif
    assert report["motif_untestable_cells"] == untestable_motif
    assert report["rbp_bh_divisor"]["minp"] == report["rbp_family_size"] - untestable_motif
    rows = calib_rows(base / "per_motif_regions.tsv")
    zero = [r for r in rows if r["motif_key"] == arm["zero"]]
    assert zero and all(r["untestable"] == "TRUE" and r["calibrated_q"] == "NA" for r in zero)
    assert all(r["untestable"] == "FALSE" for r in rows if r["motif_key"] != arm["zero"])
    condensed = calib_rows(base / "condensed_per_rbp.tsv")
    assert sum(r["rbp_untestable"] == "TRUE" and r["plot"] == "TRUE" for r in condensed) == untestable_motif
    readout = (base / "readout.md").read_text(encoding="utf-8")
    assert "BH divisor {} testable cells of the {}".format(report["motif_bh_divisor"],
                                                           report["motif_family_size"]) in readout
    assert "BH divisor {} of {}".format(report["rbp_bh_divisor"]["minp"], report["rbp_family_size"]) in readout
    import openpyxl
    book = openpyxl.load_workbook(base / f"{syn.ARM}_calibrated_ranksum_v2.xlsx", read_only=False)
    readme = " ".join(str(c.value) for row in book["README"].iter_rows() for c in row)
    assert "divisor {} = the testable cells".format(report["motif_bh_divisor"]) in readme
    header = [c.value for c in book["per_motif"][1]]
    assert "untestable" in header


def test_figures_legend_subtitle_sidecar_and_index_print_the_divisor_not_the_family(arm):
    base = arm["calibrated"] / syn.ARM
    report = json.loads((base / "refinement_report.json").read_text(encoding="utf-8"))
    figures = arm["figures"] / syn.ARM
    motif_div, rbp_div = report["motif_bh_divisor"], report["rbp_bh_divisor"]["minp"]
    by_motif = (figures / f"{syn.ARM}_SE_byMotif_calibrated_ranksum.svg").read_text(encoding="utf-8")
    by_rbp = (figures / f"{syn.ARM}_SE_byRBP_calibrated_ranksum.svg").read_text(encoding="utf-8")
    assert f"BH q over {motif_div} motif-level tests" in by_motif
    assert f"BH q over {rbp_div} RBP-level tests" in by_rbp
    assert f"BH q over {report['motif_family_size']} " not in by_motif
    assert "6 untestable" in by_motif and "6 untestable" in by_rbp
    sidecar = (figures / f"{syn.ARM}_figure_provenance_v43.md").read_text(encoding="utf-8")
    assert f"{motif_div} motif-level tests" in sidecar and f"{rbp_div} RBP-level tests" in sidecar
    index = (arm["figures"] / "index.html").read_text(encoding="utf-8")
    assert f"by-motif {motif_div} motif-level tests" in index and "726" not in index
    drawn = (arm["figures"] / "selection_audit_v43.tsv").read_text(encoding="utf-8")
    calibrated_rows = [line for line in drawn.splitlines() if "\tcalibrated_ranksum\t" in line]
    assert calibrated_rows and not any(f"\t{arm['zero']}\t" in line for line in calibrated_rows), \
        "an untestable cell must never be drawn on the supplement"


def test_divisor_that_disagrees_with_the_tables_is_refused(arm, tmp_path):
    import shutil
    copy = tmp_path / "summary"
    shutil.copytree(arm["calibrated"], copy)
    path = copy / syn.ARM / "refinement_report.json"
    report = json.loads(path.read_text(encoding="utf-8"))
    report["motif_bh_divisor"] = report["motif_family_size"]
    path.write_text(json.dumps(report), encoding="utf-8")
    mappings, _, _ = lol.load_naming(lol.NAMING_TABLE, lol.ESRP_TABLE, arm["gtf"], syn.ALIAS)
    with pytest.raises(ValueError, match="BH divisors in refinement_report.json disagree"):
        lol.calibrated_entries(syn.ARM, copy, mappings)


def calib_rows(path):
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in handle]
