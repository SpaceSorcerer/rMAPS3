"""AUDIT F5/F14: genome errors and unavailable statistics must remain visible."""
import pytest
from pyfaidx import Fasta

from rmaps_core.genome_access import fetch_seq, revcomp
from rmaps_core.stat_utils import compute_locus_pvalue


def test_fetch_edges_pad_before_reverse_complement(tmp_path):
    path = tmp_path / "tiny.fa"
    path.write_text(">chr1 description\nACGTAC\n", encoding="utf-8")
    with Fasta(str(path), as_raw=True) as genome:
        for start, end, expected in [(-2, 4, "NNACGT"), (4, 9, "ACNNN"),
                                      (-5, -2, "NNN"), (8, 11, "NNN")]:
            assert fetch_seq(genome, "+", "chr1", start, end) == expected
            assert fetch_seq(genome, "-", "chr1", start, end) == revcomp(expected)
        with pytest.raises(ValueError, match="missing:0-4"):
            fetch_seq(genome, "+", "missing", 0, 4)


@pytest.mark.parametrize("one,two,method", [
    ([], [0, 1], "fisher"),
    ([float("nan")], [0], "fisher"),
    ([0, 0], [0, 0], "brunnermunzel"),
])
def test_unavailable_statistics_raise(one, two, method):
    with pytest.raises(ValueError):
        compute_locus_pvalue(one, two, method)


def test_statistical_exception_is_not_p_one(monkeypatch):
    def broken(*args, **kwargs):
        raise RuntimeError("deliberate statistical failure")
    monkeypatch.setattr("rmaps_core.stat_utils.stats.fisher_exact", broken)
    with pytest.raises(ValueError, match="deliberate statistical failure"):
        compute_locus_pvalue([1, 0], [0, 1], "fisher")


def test_nonfinite_p_value_is_not_p_one(monkeypatch):
    monkeypatch.setattr("rmaps_core.stat_utils.stats.fisher_exact",
                        lambda *args, **kwargs: (0, float("nan")))
    with pytest.raises(ValueError, match="Invalid statistical p-value"):
        compute_locus_pvalue([1, 0], [0, 1], "fisher")
