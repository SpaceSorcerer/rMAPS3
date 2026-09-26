"""AUDIT F14 (CLIP): a window whose test is undefined is written NA with its reason; the CLIP map never crashes on it
and never writes p=1 for it. Brunner-Munzel is undefined wherever a group is constant, which on sparse CLIP peaks is
most windows of a small exon set (GitHub Actions run 36251515896 died on this with exit 247)."""
import csv
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

from rmaps_core.stat_utils import StatisticUnavailable, compute_locus_pvalue, locus_pvalue_or_na

ROOT = Path(__file__).resolve().parents[1]
HEADER = "chr\tstrand\texonStart_0base\texonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n"


def test_undefined_test_is_na_with_reason_and_other_failures_still_raise(monkeypatch):
    pvalue, reason = locus_pvalue_or_na([0, 0, 0], [0, 1, 2], "brunnermunzel_greater")
    assert math.isnan(pvalue) and "nonconstant observations" in reason
    with pytest.raises(StatisticUnavailable):
        compute_locus_pvalue([0, 0, 0], [0, 1, 2], "brunnermunzel_greater")
    pvalue, reason = locus_pvalue_or_na([0, 3, 1], [0, 1, 0, 0], "brunnermunzel_greater")
    assert 0.0 <= pvalue <= 1.0 and reason == ""
    monkeypatch.setattr("rmaps_core.stat_utils.stats.fisher_exact",
                        lambda *args, **kwargs: (0, float("nan")))
    with pytest.raises(ValueError, match="Invalid statistical p-value") as info:
        locus_pvalue_or_na([1, 0], [0, 1], "fisher")
    assert not isinstance(info.value, StatisticUnavailable)


def clip_inputs(tmp_path):
    exon = tmp_path / "exon"
    exon.mkdir()
    peaks = []
    position = 10_000
    for name, n in (("up", 3), ("dn", 3), ("bg", 20)):
        rows = []
        for i in range(n):
            s = position
            rows.append(f"chr1\t+\t{s}\t{s + 100}\t{s - 1000}\t{s - 900}\t{s + 1000}\t{s + 1100}\n")
            if (name == "up" and i == 0) or (name == "bg" and i % 2 == 0):
                peaks.append(f"chr1\t{s - 150}\t{s - 100}\tpeak\t{i + 1}\t+\n")
            position += 10_000
        (exon / f"{name}.coord.txt").write_text(HEADER + "".join(rows), encoding="utf-8")
    (tmp_path / "peaks.bed").write_text("".join(peaks), encoding="utf-8")
    return exon, tmp_path / "peaks.bed"


def run_clip(tmp_path, method):
    exon, peaks = clip_inputs(tmp_path)
    out = tmp_path / f"out_{method}"
    out.mkdir()
    env = dict(os.environ, RMAPS_STAT_METHOD=method)
    proc = subprocess.run([sys.executable, str(ROOT / "bin" / "RNA.map.noWiggle.SE.py"), str(exon), str(peaks),
                           "250", "50", "10", "1", "0.05", "TestRBP", str(out), "3", "3", "20", "0", "SE"],
                          capture_output=True, text=True, cwd=ROOT, env=env)
    assert proc.returncode == 0, (proc.stdout + proc.stderr)[-2000:]
    tables = {}
    for direction in ("up", "dn"):
        with open(out / f"pVal.{direction}.vs.bg.RNAmap.txt", newline="", encoding="utf-8") as handle:
            tables[direction] = list(csv.DictReader(handle, delimiter="\t"))
    return tables


def test_clip_brunnermunzel_writes_na_with_reason_instead_of_crashing(tmp_path):
    tables = run_clip(tmp_path, "brunnermunzel_greater")
    column = "brunnermunzel.greater.pVal"
    for rows in tables.values():
        assert rows and set(rows[0]) == {"Region", "position", column, "reason"}
        for row in rows:
            assert (row[column] == "NA") == bool(row["reason"])
            assert row[column] == "NA" or 0.0 <= float(row[column]) <= 1.0
    assert all(row[column] == "NA" and "nonconstant observations" in row["reason"] for row in tables["dn"])
    assert any(row[column] != "NA" for row in tables["up"])


def test_clip_fisher_has_no_unavailable_windows(tmp_path):
    tables = run_clip(tmp_path, "fisher")
    for rows in tables.values():
        assert rows and all(row["reason"] == "" and 0.0 <= float(row["fisher.exact.pVal"]) <= 1.0 for row in rows)
