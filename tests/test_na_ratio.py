"""An undefined ratio (background count 0 in a window) draws an open dot; it never crashes the build and never moves a
rank or a colour."""
import csv
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_region_lollipops_v4 as lol  # noqa: E402
import rmaps_countdist_io as io  # noqa: E402
import synthetic_calibrated_arm as syn  # noqa: E402

ZERO_REGION = io.REGIONS.index("UpstreamIntron")


def rows(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.fixture(scope="module")
def zero_bg_arm(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("zero_bg_arm")
    inputs = syn.build(tmp, zero_bg=("QKI.", ZERO_REGION))
    inputs["figures"] = syn.draw(inputs, tmp / "figures")
    return inputs


def test_the_calibration_summary_carries_an_na_ratio_for_the_zero_background_window(zero_bg_arm):
    table = rows(zero_bg_arm["calibrated"] / syn.ARM / "per_motif_regions.tsv")
    zero = [r for r in table if r["motif_key"].startswith("QKI.") and r["region"] == "UpstreamIntron"]
    assert zero and all(r["enrichment_ratio"] == "NA" and r["count_ratio"] == "NA" for r in zero)
    assert all(float(r["bg_mean_count"]) == 0 for r in zero)


def test_the_build_succeeds_and_the_open_dot_audit_lists_the_motif(zero_bg_arm):
    figures = zero_bg_arm["figures"]
    audit = rows(figures / "open_dot_audit_v43.tsv")
    qki = [r for r in audit if r["motif_key"].startswith("QKI.") and r["pooled_region"] == "Upstream Intron"]
    assert qki, audit[:3]
    assert {r["layer"] for r in qki} == {"released_ranksum_rawP", "calibrated_ranksum"}
    assert all(float(r["bg_mean_count"]) == 0 for r in qki)
    svg = (figures / syn.ARM / f"{syn.ARM}_SE_byMotif_released_ranksum_rawP.svg").read_text(encoding="utf-8")
    assert "open dot = ratio undefined (background count 0)" in svg


def test_ranking_and_colour_are_unaffected(zero_bg_arm):
    figures = zero_bg_arm["figures"]
    selection = rows(figures / "selection_audit_v43.tsv")
    panels = {}
    for r in selection:
        panels.setdefault((r["layer"], r["variant"], r["type"], r["direction"], r["pooled_region"]), []).append(r)
    for key, drawn in panels.items():
        ps = [float(r["p"]) for r in sorted(drawn, key=lambda r: int(r["rank"]))]
        assert ps == sorted(ps), key                       # rank order is the p order; an NA ratio breaks no tie
    first = {k: min(v, key=lambda r: int(r["rank"]))["label"] for k, v in panels.items()}
    assert first[("released_ranksum_rawP", "main", "byMotif", "INCLUDED", "Upstream Intron")] == "QKI"
    colours = rows(figures / "colour_audit_v43.tsv")
    assert colours and all(r["pass"] == "True" for r in colours)
    na = [r for r in selection if r["motif_key"].startswith("QKI.") and r["pooled_region"] == "Upstream Intron"
          and r["layer"] == "calibrated_ranksum"]
    assert na and all(r["enrichment_ratio"] == "NA" for r in na)


def test_ratio_helpers_order_undefined_after_defined():
    assert math.isnan(lol.ratio_value("NA")) and lol.ratio_value("2.5") == 2.5
    ordered = sorted([float("nan"), 1.0, 3.0], key=lol.ratio_order)
    assert ordered[:2] == [3.0, 1.0] and math.isnan(ordered[2])
