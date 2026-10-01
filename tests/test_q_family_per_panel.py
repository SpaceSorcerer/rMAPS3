"""--q-family per_panel: RBP-level BH within each direction x pooled-region panel; the default stays pooled."""
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_region_lollipops_v4 as lol  # noqa: E402
import synthetic_calibrated_arm as syn  # noqa: E402


def entries(arm, q_family):
    mappings, _, _ = lol.load_naming(lol.NAMING_TABLE, lol.ESRP_TABLE, arm["gtf"], syn.ALIAS)
    return lol.calibrated_entries(syn.ARM, arm["calibrated"], mappings, q_family)


def test_default_keeps_the_pooled_calibration_q(calibrated_arm):
    rows, _, report = entries(calibrated_arm, "pooled")
    assert "rbp_q_family" not in report
    assert all(e["_rbp_q"] == e["_rbp_q_minp"] for e in rows)


def test_per_panel_q_is_bh_within_each_panel_over_finite_p(calibrated_arm):
    level, _ = lol.load_rbp_level(calibrated_arm["calibrated"] / syn.ARM / "condensed_per_rbp.tsv", syn.ARM)
    pooled = {k: dict(v) for k, v in level.items()}
    divisors = lol.per_panel_rbp_q(level)
    panels = defaultdict(list)
    for k, v in level.items():
        panels[(k[1], k[2])].append(k)
    assert len(divisors) == len(panels) == 6
    for (direction, region), keys in panels.items():
        finite = [k for k in keys if math.isfinite(level[k]["_rbp_p"])]
        assert divisors[f"{direction}|{region}"] == len(finite)
        expected = false_discovery_control(np.array([level[k]["_rbp_p"] for k in finite]), method="bh")
        np.testing.assert_allclose([level[k]["_rbp_q"] for k in finite], expected)
        assert all(math.isnan(level[k]["_rbp_q"]) for k in keys if k not in finite), "untestable q stays NaN"
    assert all(level[k]["_rbp_q_minp"] == pooled[k]["_rbp_q_minp"] for k in level), "pooled q kept for sensitivity"


def test_per_panel_family_text_and_figure_wording(calibrated_arm, tmp_path):
    _, _, report = entries(calibrated_arm, "per_panel")
    assert report["rbp_q_family"] == "per_panel"
    assert "within each region × direction panel" in lol.family_text(report, True)
    assert lol.family_text(report, False) == lol.family_text({k: v for k, v in report.items()
                                                              if k != "rbp_q_family"}, False), "motif family unchanged"
    figures = syn.draw(calibrated_arm, tmp_path / "figures", extra=("--q-family", "per_panel"))
    svg = (figures / syn.ARM / f"{syn.ARM}_SE_byRBP_calibrated_ranksum.svg").read_text(encoding="utf-8")
    assert "RBP-level tests within each" in svg
    sidecar = (figures / syn.ARM / f"{syn.ARM}_figure_provenance_v43.md").read_text(encoding="utf-8")
    assert "--q-family per_panel" in sidecar
