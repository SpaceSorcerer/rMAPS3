"""Figure builder v4.3.2: the significance ramp is visibly graded (CIEDE2000 targets at the key ticks)."""
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import build_region_lollipops_v4 as lol  # noqa: E402

# Ramp ends seen in the lab build: permutation floor, arm q floors (with and without the 0.01 anchor), y caps.
RAMP_ENDS = [1 / 100_001, 0.003, 0.002, 0.00242, 0.006, 0.00726, 0.024, 0.029, 1e-15, 1e-20, 1e-30]


def test_ciede2000_matches_sharma_2005_reference_pairs():
    pairs = [((50, 2.6772, -79.7751), (50, 0, -82.7485), 2.0425),
             ((50, -1.3802, -84.2814), (50, 0, -82.7485), 1.0000),
             ((60.2574, -34.0099, 36.2677), (60.4626, -34.1751, 39.4387), 1.2644),
             ((22.7233, 20.0904, -46.6940), (23.0331, 14.9730, -42.5619), 2.0373),
             ((50, 2.5, 0), (73, 25, -18), 27.1492),
             ((90.8027, -2.0831, 1.4410), (91.1528, -1.6435, 0.0447), 1.4441),
             ((2.0776, 0.0795, -1.1350), (0.9033, -0.0636, -0.5514), 0.9082)]
    for lab1, lab2, expected in pairs:
        assert lol.delta_e2000_lab(lab1, lab2) == pytest.approx(expected, abs=1e-4)


@pytest.mark.parametrize("end", RAMP_ENDS)
def test_adjacent_key_ticks_differ_by_at_least_15_and_grey_vs_tint_by_20(end):
    rows = lol.contrast_rows(end)
    assert rows and all(r["pass"] for r in rows)
    for d in ("INCLUDED", "SKIPPED"):
        grey = [r for r in rows if r["direction"] == d and r["from"].startswith("grey")]
        ticks = [r for r in rows if r["direction"] == d and not r["from"].startswith("grey")]
        assert len(grey) == 1 and grey[0]["delta_e2000"] >= 20
        assert ticks and min(r["delta_e2000"] for r in ticks) >= 15


@pytest.mark.parametrize("end", RAMP_ENDS)
def test_dot_colours_at_the_key_ticks_are_the_anchor_colours(end):
    for value, t in lol.key_anchors(end):
        for d in ("INCLUDED", "SKIPPED"):
            s = value * (1 - 1e-12) if value == lol.SIG_ALPHA else value
            assert lol.significance_colour(d, s, end) == lol.ramp_hex(d, t)


def test_anchor_set_drops_ticks_too_close_to_the_end():
    assert [v for v, _ in lol.key_anchors(1e-20)] == [0.05, 0.01, 0.001, 1e-20]
    assert [v for v, _ in lol.key_anchors(0.003)] == [0.05, 0.01, 0.003]
    assert [v for v, _ in lol.key_anchors(0.0011)] == [0.05, 0.01, 0.0011]
    assert [v for v, _ in lol.key_anchors(0.009)] == [0.05, 0.009]


@pytest.mark.parametrize("direction", ["INCLUDED", "SKIPPED"])
def test_ramp_spans_a_pale_tint_to_a_darkened_hue(direction):
    L0 = lol.hex_to_oklch(lol.ramp_hex(direction, 0.0))[0]
    L1, C1, h1 = lol.hex_to_oklch(lol.ramp_hex(direction, 1.0))
    assert 0.89 <= L0 <= 0.93
    assert 0.39 <= L1 <= 0.46
    assert L1 < lol.hex_to_oklch(lol.DIRECTION_HUE[direction])[0]
    assert abs(math.remainder(h1 - lol.hex_to_oklch(lol.DIRECTION_HUE[direction])[2], math.tau)) < 0.05


@pytest.mark.parametrize("end", [0.003, 1e-20])
def test_depth_is_monotone_in_significance(end):
    grid = [10 ** -(1.31 + i * 0.05) for i in range(400)]
    depths = [lol.depth_fraction(s, end) for s in grid]
    assert depths == sorted(depths)
    assert depths[-1] == 1.0
