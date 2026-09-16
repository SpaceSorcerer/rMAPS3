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
import numpy as np
from pyfaidx import Fasta
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rmaps_core.genome_access import fetch_seq, revcomp
from rmaps_core.se_windows import (
    Event, RegionSequence, binary_density, binary_table, binary_windows, event_regions,
    overlapping_hits, read_event_sets, REGION_NAMES,
)

HEADER = "chr\tstrand\texonStart\texonEnd\tfirstExonStart\tfirstExonEnd\tsecondExonStart\tsecondExonEnd\n"


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("synthetic_genome")
    build = tmp_path / "genomes" / "synthetic"
    build.mkdir(parents=True)
    sequence = list("C" * 5000)
    for start, motif in [(580, "ATAT"), (605, "ATAT"), (980, "ATAT"),
                         (1004, "AAAA"), (1096, "AAAA"), (1105, "GTGT"),
                         (1480, "GTGT"), (1505, "ATAT"), (2004, "AAAA")]:
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
    assert overlapping_hits("AA", region.sequence) == [(4, 6), (5, 7), (6, 8)]
    assert binary_windows(region, ["AA"], window=4) == [0, 1, 1, 1, 1, 1, 1]
    assert binary_windows(region, ["AA"], window=4, step=3) == [0, 1, 1]
    end_region = event_regions(genome, event(), intron=20, exon=10)[4]
    assert binary_windows(end_region, ["AA"], window=4) == [0, 0, 0, 1, 1, 1, 1]


def test_all_regions_have_reverse_complement_strand_symmetry(synthetic):
    genome, _ = synthetic
    plus = event_regions(genome, event(), intron=20, exon=10)
    minus_event = Event("chr2", "-", 3900, 4000, 3400, 3500, 4400, 4500)
    minus = event_regions(genome, minus_event, intron=20, exon=10)
    # AUDIT R6: upstream/downstream introns have deliberately different sequence.
    assert plus[1].sequence != plus[5].sequence
    assert plus[2].sequence != plus[6].sequence
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
    assert binary_windows(downstream, ["C"], 2) == [1] + [None] * 8
    assert fetch_seq(genome, "-", "chr1", -2, 2) == "GGNN"
    with pytest.raises((KeyError, ValueError), match="missing"):
        fetch_seq(genome, "+", "missing", 0, 3)


def test_overlapping_case_rna_and_variable_span_regex():
    assert overlapping_hits("AA", "aaaa") == [(0, 2), (1, 3), (2, 4)]
    assert overlapping_hits("AU[AT]", "atatta") == [(0, 3), (2, 5)]
    assert overlapping_hits("A+", "AAACA") == [(0, 3), (1, 3), (2, 3), (4, 5)]
    assert binary_windows(RegionSequence("AAAA", 0, 4), ["AA"], 2) == [1, 1, 1]


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


# AUDIT R1: a 50-nt exon region with W=50 has exactly one complete window.
@pytest.mark.parametrize("length,expected", [(250, 201), (50, 1)])
def test_default_region_complete_window_counts(length, expected):
    assert len(binary_windows(RegionSequence("C" * length, 0, length), ["AA"], 50)) == expected


def test_long_introns_fetch_only_region_plus_motif_padding(synthetic, monkeypatch):
    # AUDIT R4: a whole-intron fetch would violate this bounded-query assertion.
    genome, _ = synthetic
    requested = []
    def tracked_fetch(genome, strand, chrom, start, end):
        requested.append(end - start)
        return fetch_seq(genome, strand, chrom, start, end)
    monkeypatch.setattr("rmaps_core.se_windows.fetch_seq", tracked_fetch)
    long_event = Event("chr1", "+", 2000, 2100, 0, 100, 4500, 4600)
    regions = event_regions(genome, long_event, intron=20, exon=10, padding=3)
    assert len(regions) == len(requested) == 8
    assert all(size <= region.length + 6 for size, region in zip(requested, regions))
    assert max(requested) < 100


def test_sparse_boundary_overlaps_and_unavailable_event(tmp_path):
    # AUDIT R4: preserve hits crossing either region edge, with half-open overlap.
    from rmaps_core.positional_io import iter_hits, load_hits
    path = tmp_path / "edge.hits.npz"
    np.savez_compressed(path, event_index=np.array([0, 0, 1], dtype=np.int32),
                        hit_start=np.array([-2, 9, -1], dtype=np.int16),
                        hit_end=np.array([1, 12, 11], dtype=np.int16),
                        elig_lo=np.array([0, -1]), elig_hi=np.array([7, -1]),
                        set_label=np.array([0, 2], dtype=np.uint8),
                        window=4, step=2, region_length=10, schema_version=2)
    hits, eligible, labels = load_hits(path)
    np.testing.assert_array_equal(hits, [[1, 0, 0, 1], [0, 0, 0, 0]])
    np.testing.assert_array_equal(eligible, [[1, 1, 1, 1], [0, 0, 0, 0]])
    np.testing.assert_array_equal(labels, [0, 2])
    assert list(iter_hits(path)) == [(0, -2, 1), (0, 9, 12), (1, -1, 11)]


def test_sparse_stream_crosses_compression_chunk_boundary(tmp_path):
    # AUDIT R4: check alignment when the compressed columns span multiple chunks.
    from rmaps_core.positional_io import iter_hits
    count = 65539
    path = tmp_path / "chunked.hits.npz"
    indices = np.arange(count, dtype=np.int32)
    starts = (indices % 201).astype(np.int16)
    np.savez_compressed(path, event_index=indices, hit_start=starts, hit_end=starts + 4)
    observed_count = 0
    for expected_index, (event_index, lo, hi) in enumerate(iter_hits(path)):
        assert event_index == expected_index
        assert lo == expected_index % 201 and hi == lo + 4
        observed_count += 1
    assert observed_count == count


MOTIFS = ("TEST.AA", "OTHER.ATAT")
FASTA_REGIONS = ("UpstreamExon", "UpstreamExonIntron", "UpstreamIntron",
                 "TargetExon", "DownstreamIntron", "DownstreamExonIntron", "DownstreamExon")


def expected_inventory():
    # AUDIT R7: independent, exhaustive contract; never derive expected paths from manifest.
    paths = {"log.motifMap.txt", "pVal.up.vs.bg.RNAmap.txt", "pVal.dn.vs.bg.RNAmap.txt"}
    paths |= {f"exon/{label}.coord.txt" for label in ("up", "dn", "bg")}
    paths |= {f"fasta/{label}.{region}.fasta" for label in ("up", "dn", "bg") for region in FASTA_REGIONS}
    paths |= {f"temp/sequence_region_{index}.npy" for index in range(8)}
    paths.add("temp/sequence_metadata.npz")
    for motif in MOTIFS:
        paths |= {f"positional/{motif}.{region}.hits.npz" for region in REGION_NAMES}
        paths |= {f"temp/{motif}.countDist.{label}.txt" for label in ("up", "dn", "bg")}
        paths |= {f"temp/{motif}.pVal.{label}.vs.bg.txt" for label in ("up", "dn")}
        paths.add(f"temp/{motif}.txt")
        paths |= {f"maps/SE.{motif.replace('.', '-')}.{extension}" for extension in ("pdf", "png")}
    return paths


def run_cli(command):
    result = subprocess.run(command, capture_output=True, text=True,
                            env=dict(os.environ, RMAPS_FORCE_MOTIF_FALLBACK="1"),
                            timeout=60, cwd=ROOT)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "falling back to serial" not in (result.stdout + result.stderr)
    return result


@pytest.fixture(scope="module")
def cli_runs(synthetic, tmp_path_factory):
    # AUDIT R5/R6: real CLI, two motifs, two workers, both Fisher alternatives.
    _, fasta_root = synthetic
    work = tmp_path_factory.mktemp("cli_numeric")
    paths = input_paths(work)
    write_set(paths["bg"], [("chr1", "+", 3000, 3100, 2500, 2997, 3500, 3600)])
    motifs = work / "motifs.tsv"
    motifs.write_text("Protein_name\tregularExpression\nTEST\tAA\nOTHER\tATAT\n", encoding="ascii")
    base = [sys.executable, str(ROOT / "cli.py"), "motif-map", "se",
            "--known-motifs", str(motifs), "--motifs", "NA",
            "--fasta-root", str(fasta_root), "--genome", "synthetic",
            "--rMATS", "NA", "--miso", "NA", "--up", str(paths["up"]),
            "--down", str(paths["dn"]), "--background", str(paths["bg"]),
            "--intron", "20", "--exon", "10", "--window", "4", "--step", "2"]
    runs = {}
    for alternative, workers in (("greater", 2), ("two-sided", 1)):
        output = work / alternative
        command = base + ["--output", str(output), "--workers", str(workers),
                          "--fisher-alternative", alternative]
        run_cli(command)
        runs[alternative] = (output, command)
    return runs


def read_tsv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        return reader.fieldnames, rows


def test_cli_every_pvalue_roots_alternative_and_na(cli_runs):
    from rmaps_core.positional_io import load_hits
    observed_alternatives = {}
    na_count = 0
    for alternative, (output, _) in cli_runs.items():
        observed_alternatives[alternative] = []
        for direction, label in (("up", 0), ("dn", 1)):
            root_header, root_rows = read_tsv(output / f"pVal.{direction}.vs.bg.RNAmap.txt")
            assert len(root_header) == 9
            roots = {row["RBP"]: row for row in root_rows}
            assert set(roots) == set(MOTIFS)
            for motif in MOTIFS:
                header, rows = read_tsv(output / "temp" / f"{motif}.pVal.{direction}.vs.bg.txt")
                assert header == ["Region", "position", "fisher.exact.pVal", "reason"]
                expected_pairs = set()
                for region_index, region in enumerate(REGION_NAMES):
                    length = 10 if region_index in (0, 3, 4, 7) else 20
                    starts = list(range(0, length - 4 + 1, 2))
                    hits, eligible, labels = load_hits(output / "positional" / f"{motif}.{region}.hits.npz")
                    assert hits.shape == eligible.shape == (3, len(starts))
                    assert hits.dtype == eligible.dtype == np.bool_
                    np.testing.assert_array_equal(labels, [0, 1, 2])
                    region_rows = [row for row in rows if row["Region"] == region]
                    assert [int(row["position"]) for row in region_rows] == starts
                    values = []
                    for col, row in enumerate(region_rows):
                        expected_pairs.add((region, str(starts[col])))
                        counts = []
                        for selected_label in (label, 2):
                            selected = (labels == selected_label) & eligible[:, col]
                            positives = int(hits[selected, col].sum())
                            counts.append([positives, int(selected.sum()) - positives])
                        if any(sum(count) == 0 for count in counts):
                            assert row["fisher.exact.pVal"] == "NA"
                            assert row["reason"] == "no eligible events in one or both sets"
                            na_count += 1
                        else:
                            expected = fisher_exact(counts, alternative=alternative).pvalue
                            assert float(row["fisher.exact.pVal"]) == pytest.approx(expected, abs=1e-14)
                            assert row["reason"] == ""
                            values.append(expected)
                            observed_alternatives[alternative].append(float(row["fisher.exact.pVal"]))
                    root = roots[motif][root_header[region_index + 1]]
                    assert (root == "NA") if not values else float(root) == pytest.approx(min(values))
                assert {(row["Region"], row["position"]) for row in rows} == expected_pairs
                assert len(rows) == len(expected_pairs)
    assert na_count > 0
    assert observed_alternatives["greater"] != observed_alternatives["two-sided"]


def test_cli_sparse_schema_stream_and_parallel_agreement(cli_runs):
    from rmaps_core.positional_io import iter_hits, load_hits
    parallel = cli_runs["greater"][0]
    serial = cli_runs["two-sided"][0]
    for motif in MOTIFS:
        for region in REGION_NAMES:
            name = f"{motif}.{region}.hits.npz"
            path = parallel / "positional" / name
            with np.load(path, allow_pickle=False) as data:
                assert set(data.files) == {"event_index", "hit_start", "hit_end", "elig_lo", "elig_hi",
                                           "set_label", "window", "step", "region_length", "schema_version"}
                assert data["event_index"].dtype == np.int32
                assert data["hit_start"].dtype == data["hit_end"].dtype == np.int16
                expected_stream = list(zip(data["event_index"].tolist(), data["hit_start"].tolist(),
                                           data["hit_end"].tolist()))
                assert list(iter_hits(path)) == expected_stream
                dense, eligible, labels = load_hits(path)
                starts = list(range(0, int(data["region_length"]) - 4 + 1, 2))
                for event_index in range(3):
                    for col, start in enumerate(starts):
                        expected_eligibility = data["elig_lo"][event_index] >= 0 and data["elig_lo"][event_index] <= start < data["elig_hi"][event_index]
                        assert eligible[event_index, col] == expected_eligibility
                        expected_hit = expected_eligibility and any(i == event_index and lo < start + 4 and hi > start
                                                                    for i, lo, hi in expected_stream)
                        assert dense[event_index, col] == expected_hit
            for first, second in zip(load_hits(path), load_hits(serial / "positional" / name)):
                np.testing.assert_array_equal(first, second)


def test_cli_manifest_matches_explicit_inventory(cli_runs):
    for output, _ in cli_runs.values():
        manifest = json.loads((output / "run_manifest.json").read_text())
        assert manifest["schema_version"] == 2
        assert manifest["statistical_method"] == "fisher"
        assert manifest["python"]["executable"] and manifest["packages"] and manifest["git_revision"]
        entries = manifest["outputs"]
        assert len(entries) == len(expected_inventory())
        assert {entry["path"] for entry in entries} == expected_inventory()
        assert {path.relative_to(output).as_posix() for path in output.rglob("*") if path.is_file()} == expected_inventory() | {"run_manifest.json"}
        for entry in entries:
            artifact = output / entry["path"]
            assert entry["status"] == "retained"
            assert artifact.stat().st_size == entry["size_bytes"]
            assert hashlib.sha256(artifact.read_bytes()).hexdigest() == entry["sha256"]


def test_cli_delete_temp_preserves_unowned_and_hashes_summary_inputs(cli_runs, tmp_path):
    # AUDIT R2/R3/R8: deletion/archive must use the registered inventory only.
    source, base_command = cli_runs["greater"]
    output = tmp_path / "ownership"
    command = list(base_command)
    command[command.index("--output") + 1] = str(output)
    run_cli(command)
    before = json.loads((output / "run_manifest.json").read_text())
    stale = output / "temp" / "STALE.fake.pVal.up.vs.bg.txt"
    stale.write_text("user-owned; preserve these bytes\n")
    refusal = subprocess.run(command, capture_output=True, text=True, cwd=ROOT, timeout=15)
    assert refusal.returncode != 0
    assert "overwrite" in (refusal.stdout + refusal.stderr).lower()
    run_cli(command + ["--overwrite", "--delete-temp"])
    assert stale.read_text() == "user-owned; preserve these bytes\n"
    archives = list(output.glob("_previous_*"))
    assert len(archives) == 1
    archive = archives[0]
    assert {path.relative_to(archive).as_posix() for path in archive.rglob("*") if path.is_file()} == expected_inventory() | {"run_manifest.json"}
    for entry in before["outputs"]:
        assert hashlib.sha256((archive / entry["path"]).read_bytes()).hexdigest() == entry["sha256"]
    current = json.loads((output / "run_manifest.json").read_text())
    assert {entry["path"] for entry in current["outputs"]} == expected_inventory()
    for entry in current["outputs"]:
        deleted = entry["path"].startswith("temp/")
        assert entry["status"] == ("deleted" if deleted else "retained")
        assert (output / entry["path"]).exists() == (not deleted)
        if ".pVal." in entry["path"]:
            assert entry["sha256"] == hashlib.sha256((archive / entry["path"]).read_bytes()).hexdigest()
    assert "STALE" not in (output / "pVal.up.vs.bg.RNAmap.txt").read_text()


def test_cli_exon_window_forwarding(cli_runs, tmp_path):
    # AUDIT R1: an explicit exon window changes only exon window counts.
    from rmaps_core.positional_io import load_hits
    _, base_command = cli_runs["two-sided"]
    output = tmp_path / "exon_window"
    command = list(base_command)
    command[command.index("--output") + 1] = str(output)
    run_cli(command + ["--exon-window", "10"])
    for index, region in enumerate(REGION_NAMES):
        path = output / "positional" / f"TEST.AA.{region}.hits.npz"
        hits, _, _ = load_hits(path)
        assert hits.shape == (3, 1 if index in (0, 3, 4, 7) else 9)
        _, rows = read_tsv(output / "temp" / "TEST.AA.pVal.up.vs.bg.txt")
        assert len([row for row in rows if row["Region"] == region]) == hits.shape[1]


def test_cli_xlsx_delete_temp_preserves_converted_input_hash(cli_runs, tmp_path):
    # AUDIT R3: wrapper provenance must retain engine-deleted inputs when adding XLSX conversion.
    from openpyxl import Workbook
    book = Workbook()
    sheet = book.active
    sheet.append(["ID", "GeneID", "geneSymbol"] + HEADER.strip().split("\t") +
                 ["FDR", "IncLevel1", "IncLevel2", "IncLevelDifference"])
    for index, (start, fdr, delta) in enumerate(((1000, 0.01, 0.2), (2000, 0.01, -0.2),
                                               (3000, 0.8, 0.0))):
        sheet.append([f"fixture_{index}", "synthetic", "synthetic", "chr1", "+",
                      start, start + 100, start - 500, start - 400, start + 500,
                      start + 600, fdr, "0.5", "0.5", delta])
    source = tmp_path / "input.xlsx"
    book.save(source)
    book.close()
    _, base = cli_runs["two-sided"]
    output = tmp_path / "xlsx_out"
    command = list(base)
    command[command.index("--rMATS") + 1] = str(source)
    command[command.index("--output") + 1] = str(output)
    run_cli(command)
    converted = output / "temp" / "input.from_xlsx.tsv"
    expected_hash = hashlib.sha256(converted.read_bytes()).hexdigest()
    run_cli(command + ["--overwrite", "--delete-temp"])
    manifest = json.loads((output / "run_manifest.json").read_text())
    entries = {entry["path"]: entry for entry in manifest["outputs"]}
    assert set(entries) == expected_inventory() | {"temp/input.from_xlsx.tsv"}
    assert entries["temp/input.from_xlsx.tsv"]["sha256"] == expected_hash
    assert entries["temp/input.from_xlsx.tsv"]["status"] == "deleted"
    assert not converted.exists()
    assert manifest["parameters"]["original_rmats"] == str(source)
    archive, = output.glob("_previous_*")
    for path, entry in entries.items():
        if path.startswith("temp/"):
            assert entry["status"] == "deleted"
            assert entry["sha256"] == hashlib.sha256((archive / path).read_bytes()).hexdigest()
