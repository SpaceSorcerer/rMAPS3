"""Synthetic tests for tools/motif_scores_ranksum.py. Seed 149 everywhere.

Ported from the 2026-09-21 lab tests; the part that read lab count archives is replaced by a
synthetic end-to-end run through tools/calibrate_ranksum.py (row unit), whose v1 summary this
tool augments. The port was also checked against the lab outputs on 2026-09-24 (see LESSONS.md).
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tests"))

import calibrate_ranksum as calib  # noqa: E402
import motif_scores_ranksum as ms  # noqa: E402
import rmaps_countdist_io as io  # noqa: E402
from test_unit_sensitivity import ARM, MOTIFS, build_arm  # noqa: E402

SEED = 149


def write_countdist(path: Path, matrix: np.ndarray, positions_per_region: int):
    """The released four-column countDist text format, as rmaps_countdist_io parses it."""
    with open(path, "w", newline="\n", encoding="utf-8") as handle:
        handle.write("Region\tposition\tsum\tvalues\n")
        column = 0
        for region in io.REGIONS:
            for position in range(positions_per_region):
                values = matrix[:, column]
                handle.write("{}\t{}\t{}\t[{}]\n".format(
                    region, position, int(values.sum()), ",".join(str(int(v)) for v in values)))
                column += 1


def test_score_reproduces_countdist_sums_exactly(tmp_path):
    rng = np.random.default_rng(SEED)
    n_exons, per_region = 7, 3
    matrix = rng.integers(0, 3, size=(n_exons, len(io.REGIONS) * per_region)).astype(np.int16)
    path = tmp_path / "TEST.MOTIF.countDist.up.txt"
    write_countdist(path, matrix, per_region)
    _region_index, _position, parsed = io.parse_countdist(path)
    assert np.array_equal(parsed, matrix)
    with open(path, "r", encoding="utf-8") as handle:
        handle.readline()
        file_sums = [int(line.split("\t")[2]) for line in handle]
    scores = ms.motif_scores(parsed, parsed)
    for window, total in enumerate(file_sums):
        assert float(scores["sum_hits_changed"][window]) == float(total)
        assert float(scores["motif_score_changed"][window]) == total / n_exons


def test_ratio_hand_check_on_tiny_matrix():
    """fg counts [2,1,0]; bg counts [0,0,1,0]; every value worked out by hand."""
    fg = np.array([[2], [1], [0]], dtype=np.int16)
    bg = np.array([[0], [0], [1], [0]], dtype=np.int16)
    scores = ms.motif_scores(fg, bg)
    assert int(scores["sum_hits_changed"][0]) == 3
    assert int(scores["sum_hits_background"][0]) == 1
    assert int(scores["n_changed"][0]) == 3
    assert int(scores["n_background"][0]) == 4
    assert scores["motif_score_changed"][0] == pytest.approx(1.0)
    assert scores["motif_score_background"][0] == pytest.approx(0.25)
    assert scores["motif_score_ratio"][0] == pytest.approx(4.0)
    assert ms.score_flag(scores, 0) == ""
    row = ms.scores_at(scores, 0)
    assert row["motif_score_ratio"] == pytest.approx(4.0)
    assert row["motif_score_flag"] == ""


def test_background_with_no_hits_gives_na_ratio_and_a_flag():
    fg = np.array([[1], [1], [0]], dtype=np.int16)
    bg = np.zeros((5, 1), dtype=np.int16)
    scores = ms.motif_scores(fg, bg)
    assert scores["motif_score_changed"][0] == pytest.approx(2 / 3)
    assert float(scores["motif_score_background"][0]) == 0.0
    assert math.isnan(float(scores["motif_score_ratio"][0]))
    assert ms.score_flag(scores, 0) == "background_zero_hits"
    assert ms.format_cell(float(scores["motif_score_ratio"][0])) == "NA"


def test_map_key_collapses_bracket_runs_like_the_released_engine():
    assert ms.map_key("9G8.[AT]GGAC[AG]A") == "9G8-_AT_GGAC_AG_A"
    assert ms.map_key("QKI.ACTAAC") == "QKI-ACTAAC"


def test_cli_end_to_end_on_a_v1_summary(tmp_path):
    cfg, _, _ = build_arm(tmp_path)
    summary = tmp_path / "summary"
    assert calib.main(["--arm", ARM, "--counts-root", str(cfg.counts_root), "--released-root",
                       str(cfg.released_root), "--out-root", str(summary), "--alias-table",
                       str(cfg.alias_table), "--permutations", "199", "--refine-perms", "499",
                       "--seed", str(SEED)]) == 0
    maps = Path(cfg.released_root) / ARM / "maps"
    maps.mkdir()
    for motif in MOTIFS:
        (maps / ("SE." + ms.map_key(motif) + ".png")).write_bytes(b"png")
    protected = {p.name: ms.md5(p) for p in (summary / ARM).iterdir()
                 if p.name in ("per_motif_regions.tsv", "condensed_per_rbp.tsv")}
    assert ms.main(["--arms", ARM, "--counts-root", str(cfg.counts_root), "--summary-root",
                    str(summary), "--released-root", str(cfg.released_root), "--report-root",
                    str(tmp_path / "report")]) == 0
    for name, digest in protected.items():
        assert ms.md5(summary / ARM / name) == digest
    verify, _ = ms.read_tsv(summary / ARM / "motif_scores_verification.tsv")
    assert len(verify) == len(MOTIFS) * 2 * len(io.REGIONS)
    assert all(r["score_equals_sum_over_n"] == "TRUE" for r in verify)
    v2, _ = ms.read_tsv(summary / ARM / "per_motif_regions_v2.tsv")
    for row in v2:
        if row["motif_score_flag"] == "" and row["count_ratio"] != "NA":
            assert float(row["motif_score_ratio"]) == pytest.approx(float(row["count_ratio"]))
    assert (summary / ARM / (ARM + "_calibrated_ranksum_v2.xlsx")).stat().st_size > 0
    assert "## Motif scores" in (summary / ARM / "readout.md").read_text(encoding="utf-8")
