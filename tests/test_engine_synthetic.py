"""Numerical SE regressions using two synthetic 5 kb chromosomes."""

# AUDIT F21: assert scientific outputs on a synthetic genome without external data.

import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
from pyfaidx import Fasta
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rmaps_core.genome_access import fetch_seq, revcomp
from rmaps_core.se_windows import (
    Event, RegionSequence, binary_density, binary_table, binary_windows, event_regions,
    overlapping_hits, read_event_sets,
)

HEADER = "chr\tstrand\texonStart\texonEnd\tfirstExonStart\tfirstExonEnd\tsecondExonStart\tsecondExonEnd\n"


@pytest.fixture
def synthetic(tmp_path):
    build = tmp_path / "genomes" / "synthetic"
    build.mkdir(parents=True)
    sequence = list("C" * 5000)
    for start, motif in [(580, "ATAT"), (605, "ATAT"), (980, "ATAT"),
                         (1004, "AAAA"), (1096, "AAAA"), (1105, "ATAT"),
                         (1480, "ATAT"), (1505, "ATAT"), (2004, "AAAA")]:
        sequence[start:start + len(motif)] = motif
    sequence = "".join(sequence)
    fasta = build / "synthetic.fa"
    fasta.write_text(f">chr1 1\n{sequence}\n>chr2 2\n{revcomp(sequence)}\n", encoding="ascii")
    genome = Fasta(str(fasta), as_raw=True, sequence_always_upper=True)
    yield genome, build.parent
    genome.close()


def event(start=1000, strand="+", chrom="chr1"):
    return Event(chrom, strand, start, start + 100, start - 500, start - 400,
                 start + 500, start + 600)


def write_set(path, rows):
    path.write_text(HEADER + "".join("\t".join(map(str, row)) + "\n" for row in rows), encoding="ascii")
    return path


def input_paths(tmp_path):
    return {
        "up": write_set(tmp_path / "up.tsv", [("chr1", "+", 1000, 1100, 500, 600, 1500, 1600)]),
        "dn": write_set(tmp_path / "dn.tsv", [("chr1", "+", 2000, 2100, 1500, 1600, 2500, 2600)]),
        "bg": write_set(tmp_path / "bg.tsv", [("chr1", "+", 3000, 3100, 2500, 2600, 3500, 3600)]),
    }


def test_known_hit_positions_and_half_open_overlap(synthetic):
    genome, _ = synthetic
    region = event_regions(genome, event(), intron=20, exon=10)[3]
    assert overlapping_hits("AA", region.sequence) == [(4, 6), (5, 7), (6, 8), (96, 98), (97, 99), (98, 100)]
    assert binary_windows(region, ["AA"], window=4) == [0, 1, 1, 1, 1, 1, 1, 1, 0, 0]
    assert binary_windows(region, ["AA"], window=4, step=3) == [0, 1, 1, 0]
    end_region = event_regions(genome, event(), intron=20, exon=10)[4]
    assert binary_windows(end_region, ["AA"], window=4) == [0, 0, 0, 1, 1, 1, 1, None, None, None]


def test_all_regions_have_reverse_complement_strand_symmetry(synthetic):
    genome, _ = synthetic
    plus = event_regions(genome, event(), intron=20, exon=10)
    minus_event = Event("chr2", "-", 3900, 4000, 3400, 3500, 4400, 4500)
    minus = event_regions(genome, minus_event, intron=20, exon=10)
    assert len(plus) == len(minus) == 8
    for positive, negative in zip(plus, minus):
        assert (positive.sequence, positive.origin, positive.length) == (negative.sequence, negative.origin, negative.length)
        assert binary_windows(positive, ["AA", "ATAT"], 4) == binary_windows(negative, ["AA", "ATAT"], 4)


def test_planted_enrichment_matches_hand_computed_fisher_table(synthetic):
    genome, _ = synthetic
    def hit(start):
        return binary_windows(event_regions(genome, event(start), 20, 10)[3], ["AA"], 4)[1]
    up = [hit(start) for start in (1000, 2000, 3000)]
    background = [hit(start) for start in (1000, 2500, 3000, 3500)]
    assert up == [1, 1, 0]
    assert background == [1, 0, 0, 0]
    table = binary_table(up + [None], background + [None])
    assert table == [[2, 1], [1, 3]]
    assert fisher_exact(table, alternative="greater").pvalue == pytest.approx(13 / 35)


def test_short_intron_and_exon_exclude_ineligible_windows(synthetic):
    genome, _ = synthetic
    short = Event("chr1", "+", 1000, 1003, 990, 997, 1006, 1015)
    regions = event_regions(genome, short, intron=10, exon=10)
    assert len(regions[1].sequence) == len(regions[2].sequence) == 3
    assert len(regions[5].sequence) == len(regions[6].sequence) == 3
    assert all(value is None for value in binary_windows(regions[1], ["C"], 4))
    assert all(value is None for value in binary_windows(regions[3], ["C"], 4))
    assert binary_table([1, None, 0], [None, 1]) == [[1, 1], [1, 0]]
    assert binary_density([1, None, 0]) == 0.5
    assert binary_density([None, 1]) == 1.0


def test_fetch_error_names_event(synthetic, monkeypatch):
    genome, _ = synthetic
    def fail(*args, **kwargs):
        raise OSError("synthetic read failure")
    monkeypatch.setattr("rmaps_core.se_windows.fetch_seq", fail)
    with pytest.raises(ValueError, match="chr1:.*1000.*synthetic read failure"):
        event_regions(genome, event(), 20, 10)


def test_chromosome_edge_and_minus_padding(synthetic):
    genome, _ = synthetic
    edge = Event("chr1", "+", 4900, 4998, 4700, 4800, 5003, 5010)
    downstream = event_regions(genome, edge, intron=10, exon=10)[5]
    assert binary_windows(downstream, ["C"], 2) == [1] + [None] * 9
    assert fetch_seq(genome, "-", "chr1", -2, 2) == "GGNN"
    with pytest.raises((KeyError, ValueError), match="missing"):
        fetch_seq(genome, "+", "missing", 0, 3)


def test_overlapping_case_rna_and_variable_span_regex():
    assert overlapping_hits("AA", "aaaa") == [(0, 2), (1, 3), (2, 4)]
    assert overlapping_hits("AU[AT]", "atatta") == [(0, 3), (2, 5)]
    assert overlapping_hits("A+", "AAACA") == [(0, 3), (1, 3), (2, 3), (4, 5)]
    assert binary_windows(RegionSequence("AAAA", 0, 4), ["AA"], 2) == [1, 1, 1, None]


@pytest.mark.parametrize("bad", [
    "chr1\t+\t1000\t1100\t500\t600\t1500\t1600\n",
    HEADER + "chr1\t+\t1000\t1100\t500\t600\t1500\n",
    HEADER + "chr1\t+\tx\t1100\t500\t600\t1500\t1600\n",
    HEADER + "chr1\t+\t1100\t1000\t500\t600\t1500\t1600\n",
    HEADER + "chr1\t?\t1000\t1100\t500\t600\t1500\t1600\n",
    HEADER + "missing\t+\t1000\t1100\t500\t600\t1500\t1600\n",
    HEADER,
])
def test_invalid_coordinate_sets_raise(synthetic, tmp_path, bad):
    genome, _ = synthetic
    paths = input_paths(tmp_path)
    paths["up"].write_text(bad, encoding="ascii")
    with pytest.raises(ValueError):
        read_event_sets(paths, genome)


def test_duplicate_and_between_set_overlap_validation(synthetic, tmp_path):
    genome, _ = synthetic
    paths = input_paths(tmp_path)
    up_text = paths["up"].read_text()
    paths["up"].write_text(up_text + up_text.splitlines(keepends=True)[1])
    with pytest.raises(ValueError, match="(?i)duplicate"):
        read_event_sets(paths, genome)
    paths["up"].write_text(up_text)
    paths["dn"].write_text(up_text)
    with pytest.raises(ValueError, match="(?i)overlap"):
        read_event_sets(paths, genome)
    with pytest.warns(UserWarning, match="overlapping events"):
        sets = read_event_sets(paths, genome, allow_overlap=True)
    assert [len(sets[key]) for key in ("up", "dn", "bg")] == [1, 1, 1]


def test_cli_retains_numeric_outputs_and_manifest(synthetic, tmp_path):
    _, fasta_root = synthetic
    paths = input_paths(tmp_path)
    motifs = tmp_path / "motifs.tsv"
    motifs.write_text("Protein_name\tregularExpression\nTEST\tAA\n", encoding="ascii")
    output = tmp_path / "out"
    command = [sys.executable, str(ROOT / "cli.py"), "motif-map", "se",
               "--known-motifs", str(motifs), "--motifs", "NA",
               "--fasta-root", str(fasta_root), "--genome", "synthetic",
               "--output", str(output), "--rMATS", "NA", "--miso", "NA",
               "--up", str(paths["up"]), "--down", str(paths["dn"]),
               "--background", str(paths["bg"]), "--intron", "20", "--exon", "10",
               "--window", "4", "--step", "2", "--workers", "1",
               "--fisher-alternative", "two-sided"]
    env = dict(os.environ, RMAPS_FORCE_MOTIF_FALLBACK="1")
    result = subprocess.run(command, capture_output=True, text=True, env=env, timeout=40, cwd=ROOT)
    assert result.returncode == 0, result.stdout + result.stderr
    matrices = sorted((output / "positional").glob("*.hits.tsv"))
    assert len(matrices) == 8
    for matrix in matrices:
        with matrix.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 3
        assert [row["set"] for row in rows] == ["up", "dn", "bg"]
        positions = [key for key in rows[0] if key.lstrip("-").isdigit()]
        assert positions
        assert all(row[key] in {"0", "1", "NA"} for row in rows for key in positions)
        if "TargetExon_5prime" in matrix.name:
            assert [row["2"] for row in rows] == ["1", "1", "0"]
        if "TargetExon-3prime" in matrix.name:
            assert [row["8"] for row in rows] == ["NA", "NA", "NA"]
    assert list((output / "temp").glob("*.pVal.*.txt"))
    manifest = json.loads((output / "run_manifest.json").read_text())
    assert manifest["statistical_method"] == "fisher"
    assert manifest["python"]["executable"]
    assert manifest["packages"]
    assert manifest["git_revision"]
    recorded = {entry["path"] for entry in manifest["outputs"]}
    for entry in manifest["outputs"]:
        artifact = output / entry["path"]
        assert artifact.stat().st_size == entry["size_bytes"]
        assert hashlib.sha256(artifact.read_bytes()).hexdigest() == entry["sha256"]
    assert all(matrix.relative_to(output).as_posix() in recorded for matrix in matrices)
    refusal = subprocess.run(command, capture_output=True, text=True, env=env, timeout=10, cwd=ROOT)
    assert refusal.returncode != 0
    assert "overwrite" in (refusal.stdout + refusal.stderr).lower()
    stale = output / "temp" / "STALE.fake.pVal.up.vs.bg.txt"
    stale.write_text("must not enter summaries\n")
    repeat = subprocess.run(command + ["--overwrite", "--delete-temp"], capture_output=True,
                            text=True, env=env, timeout=40, cwd=ROOT)
    assert repeat.returncode == 0, repeat.stdout + repeat.stderr
    assert not (output / "temp").exists()
    assert len(list((output / "positional").glob("*.hits.tsv"))) == 8
    assert "STALE" not in (output / "pVal.up.vs.bg.RNAmap.txt").read_text()
