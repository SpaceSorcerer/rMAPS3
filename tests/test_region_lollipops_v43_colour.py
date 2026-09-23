"""Figure builder v4.3: dot colour = direction hue deepened by significance, grey at s >= 0.05."""
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import build_region_lollipops_v4 as lol  # noqa: E402

FLOOR = 1 / 100_001


def test_hues_are_the_rbp_reli_gold_and_blue():
    assert lol.DIRECTION_HUE == {"INCLUDED": "#E69F00", "SKIPPED": "#0072B2"}
    assert lol.NS_GREY == "#999999"


@pytest.mark.parametrize("direction", ["INCLUDED", "SKIPPED"])
def test_grey_at_and_above_alpha_full_hue_at_floor(direction):
    assert lol.significance_colour(direction, 0.05, FLOOR) == lol.NS_GREY
    assert lol.significance_colour(direction, 0.7, FLOOR) == lol.NS_GREY
    assert lol.significance_colour(direction, FLOOR, FLOOR) == lol.DIRECTION_HUE[direction]
    assert lol.significance_colour(direction, FLOOR / 10, FLOOR) == lol.DIRECTION_HUE[direction]


@pytest.mark.parametrize("direction", ["INCLUDED", "SKIPPED"])
def test_ramp_is_monotone_in_lightness_in_gamut_and_holds_the_hue(direction):
    monotone, in_gamut, tint = lol.ramp_check(direction)
    assert monotone and in_gamut
    hue = lol.hex_to_oklch(lol.DIRECTION_HUE[direction])[2]
    lightness = []
    for s in (0.0499, 0.01, 0.001, 1e-4, 2e-5):
        colour = lol.significance_colour(direction, s, FLOOR)
        assert colour != lol.NS_GREY
        L, _, h = lol.hex_to_oklch(colour)
        assert abs(math.remainder(h - hue, math.tau)) < 0.05
        lightness.append(L)
    assert lightness == sorted(lightness, reverse=True)


def test_main_layer_floor_is_the_y_cap_and_supplement_floor_is_the_permutation_floor():
    assert lol.colour_floor("released_ranksum_rawP", None, {"ymax": 20}) == pytest.approx(1e-20)
    assert lol.colour_floor("calibrated_ranksum", {"stage2_permutations": 100_000}, {"ymax": 6}) == FLOOR


def test_non_finite_significance_is_refused():
    with pytest.raises(ValueError):
        lol.significance_colour("INCLUDED", float("nan"), FLOOR)
