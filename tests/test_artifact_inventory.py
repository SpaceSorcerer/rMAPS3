"""The completeness inventory is derived from the writers (review 4 must-fix 1).

Every writer registers what it writes (tools/rmaps_artifacts.py). Each test below runs the wrapper end to end on a
synthetic genome, with the released engine (a shared clone of this repository at the released commit) or the audited
engine, and asserts: REQUIRED_ARTIFACTS for the run, expanded, == the files the writers registered == the files on
disk; every registered file is hashed in run_manifest.json; and `--verify` reports INCOMPLETE when any one of them is
lost or changed after the run.
"""
import hashlib
import json
import os
import random
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import rmaps3_skill_run as skill  # noqa: E402
from test_engine_synthetic import synthetic, input_paths  # noqa: E402,F401

ARM = "SYN_ARM"
HEADER = "chr\tstrand\texonStart\texonEnd\tfirstExonStart\tfirstExonEnd\tsecondExonStart\tsecondExonEnd\n"


# ------------------------------------------------------------------ fixtures
@pytest.fixture(scope="module")
def released_engine(tmp_path_factory):
    """A clean checkout of the released engine commit: a shared clone of this repository, nothing written to it."""
    clone = tmp_path_factory.mktemp("released_engine") / "engine"
    for command in (["git", "clone", "-q", "--shared", "--no-checkout", str(ROOT), str(clone)],
                    ["git", "-C", str(clone), "checkout", "-q", skill.RELEASED_ENGINE_COMMIT]):
        proc = subprocess.run(command, capture_output=True, text=True)
        assert proc.returncode == 0, proc.stderr
    return clone


@pytest.fixture(scope="module")
def arm_inputs(tmp_path_factory):
    """20 / 18 / 200 events on a random 400 kb chromosome whose build directory is named hg38 (rank_stability.py
    accepts only an hg38 audited run), three frequent 4-mer motifs, and the figure naming inputs."""
    base = tmp_path_factory.mktemp("arm_inputs")
    rng = random.Random(149)
    build = base / "genomes" / "hg38"
    build.mkdir(parents=True)
    length = 400_000
    seq = "".join(rng.choice("ACGT") for _ in range(length))
    (build / "hg38.fa").write_text(">chr1\n" + "\n".join(seq[i:i + 60] for i in range(0, length, 60)) + "\n")
    import pyfaidx
    pyfaidx.Faidx(str(build / "hg38.fa")).close()
    position = 2000
    for name, n in (("up", 20), ("dn", 18), ("bg", 200)):
        rows = []
        for i in range(n):
            s = position
            rows.append(f"chr1\t{'+-'[i % 2]}\t{s}\t{s + 100}\t{s - 800}\t{s - 700}\t{s + 900}\t{s + 1000}\n")
            position += 1500
        (base / f"{name}.txt").write_text(HEADER + "".join(rows))
    (base / "motifs.txt").write_text("Protein_name\tregularExpression\nQKI\tACTA[AC]\nPTBP1\tTCTT\nRBFOX2\tTGCA\n")
    (base / "alias.tsv").write_text("table_name\thgnc_symbol\tevidence\tambiguity_note\n"
                                    "PTBP1\tPTBP1\tsynthetic identity row\t\n")
    (base / "genes.gtf").write_text("".join(f'chr1\tT\tgene\t1\t2\t.\t+\t.\tgene_id "ENSG{i:011d}.1"; gene_name '
                                            f'"{s}";\n' for i, s in enumerate(("QKI", "PTBP1", "RBFOX2"))))
    (base / "splice.txt").write_text("SRSF1\n")
    (base / "broad.txt").write_text("PTBP1\n")
    (base / "counts.json").write_text(json.dumps({"n_up": 20, "n_dn": 18, "n_bg": 200, "n_expr_unknown_in_fg": 0,
                                                  "n_expr_unknown_in_bg": 0, "rule": "A"}))
    return base


def argv_for(base, out, *cell):
    return ["--arm", ARM, "--out", str(out), "--up", str(base / "up.txt"), "--dn", str(base / "dn.txt"),
            "--bg", str(base / "bg.txt"), "--gate-counts", str(base / "counts.json"),
            "--genome-root", str(base / "genomes"), "--genome", "hg38", "--alias-table", str(base / "alias.tsv"),
            "--known-motifs", str(base / "motifs.txt"), "--additional-motifs", "NA", "--gtf", str(base / "genes.gtf"),
            "--spliceosome-list", str(base / "splice.txt"), "--broad-binders-list", str(base / "broad.txt"),
            "--blas-threads", "1", "--arm-label", "Synthetic arm", *cell]


def run(argv):
    proc = subprocess.run([sys.executable, str(ROOT / "tools" / "rmaps3_skill_run.py"), *argv], cwd=ROOT,
                          capture_output=True, text=True, timeout=900,
                          env=dict(os.environ, RMAPS_FORCE_MOTIF_FALLBACK="1"))
    out = Path(argv[argv.index("--out") + 1])
    manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "complete" and proc.returncode == 0, (manifest.get("error"), manifest.get("missing"),
                                                                         proc.stderr[-2000:])
    return out, manifest, skill.build_parser().parse_args(argv)


def on_disk(out):
    return {p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file()}


def expected_files(args, manifest):
    motifs, registered = manifest["root_motifs"], manifest["registered"]
    return ({rel for _, rel in skill.required_paths(args, motifs, registered=registered)}
            | {rel for _, rel in skill.required_paths(args, motifs, final=True)} | {"run_manifest.json"})


def assert_inventory_is_the_writers_output(out, manifest, args, arm=ARM):
    registered, disk, expected = set(manifest["registered"]), on_disk(out), expected_files(args, manifest)
    assert registered == disk, ("written but unregistered", sorted(disk - registered)[:10],
                                "registered but absent", sorted(registered - disk)[:10])
    assert expected == registered, ("registered, not in REQUIRED_ARTIFACTS", sorted(registered - expected)[:10],
                                    "in REQUIRED_ARTIFACTS, written by no writer", sorted(expected - registered)[:10])
    assert manifest["registered_not_required"] == [] and manifest["unregistered_files"] == []
    hashed = {row["path"]: row["sha256"] for row in manifest["inventory"]}
    assert set(hashed) == registered - {"run_manifest.json"}
    for rel, digest in hashed.items():
        assert hashlib.sha256((out / rel).read_bytes()).hexdigest() == digest, rel
    trees = [r for r in skill.required_spec(args) if r[2] == "tree"]
    assert {t[1].replace("{arm}", arm) for t in trees} <= {f"engine/{arm}", "stability/.matplotlib"}


def assert_every_lost_file_is_caught(out, manifest):
    """Remove each registered file in turn: --verify must say INCOMPLETE (exit 4) and name it; restore it after."""
    assert skill.verify_run(out) == ("complete", 0, [])
    for rel in manifest["registered"]:
        if rel == "run_manifest.json":
            continue
        path = out / rel
        aside = path.with_name(path.name + ".lost")
        path.rename(aside)
        try:
            status, code, problems = skill.verify_run(out)
        finally:
            aside.rename(path)
        assert (status, code) == ("INCOMPLETE", 4), rel
        assert any(line.startswith(rel + " absent") for line in problems), (rel, problems)
    assert skill.verify_run(out) == ("complete", 0, [])


# ------------------------------------------------------------------ cells
@pytest.fixture(scope="module")
def audited_fisher_run(arm_inputs, tmp_path_factory):
    out = tmp_path_factory.mktemp("cells") / "quick_audited_fisher"
    return run(argv_for(arm_inputs, out, "--mode", "quick", "--engine", "audited", "--workers", "1",
                        "--stat-method", "fisher"))


@pytest.fixture(scope="module")
def full_run(arm_inputs, released_engine, audited_fisher_run, tmp_path_factory):
    out = tmp_path_factory.mktemp("cells") / "full"
    return run(argv_for(arm_inputs, out, "--mode", "full", "--engine", "released", "--engine-root",
                        str(released_engine), "--permutations", "99", "--refine-perms", "199",
                        "--positive-control", "QKI", "--positive-control-advisory",
                        "--rank-stability-run", str(audited_fisher_run[0] / "engine" / ARM)))


@pytest.fixture(scope="module")
def quick_released_run(arm_inputs, released_engine, tmp_path_factory):
    out = tmp_path_factory.mktemp("cells") / "quick_released"
    return run(argv_for(arm_inputs, out, "--mode", "quick", "--engine", "released", "--engine-root",
                        str(released_engine), "--positive-control", "QKI", "--positive-control-advisory"))


def test_quick_released_rank_sum_with_figures_and_control(quick_released_run):
    out, manifest, args = quick_released_run
    assert_inventory_is_the_writers_output(out, manifest, args)
    registered = set(manifest["registered"])
    for rel in (f"figures/motif_scores/{ARM}/per_motif_regions.tsv", f"figures/motif_scores/{ARM}/command.log",
                "figures/positive_control_audit.tsv", "figures/selection_audit.tsv", "figures/layout_report.json",
                f"counts/{ARM}/VERIFY.md", "logs/temp_deletion.log"):
        assert rel in registered, rel


def test_quick_released_fisher_and_no_figures_keep_temp(arm_inputs, released_engine, tmp_path):
    for name, cell in (("fisher", ("--stat-method", "fisher")), ("keep", ("--no-figures", "--keep-temp"))):
        out, manifest, args = run(argv_for(arm_inputs, tmp_path / name, "--mode", "quick", "--engine", "released",
                                           "--engine-root", str(released_engine), *cell))
        assert_inventory_is_the_writers_output(out, manifest, args)
        assert not [r for r in manifest["registered"] if r.startswith("figures/")]
        assert "logs/temp_deletion.log" not in manifest["registered"]


def test_quick_audited_fisher_and_rank_sum(audited_fisher_run, arm_inputs, tmp_path):
    out, manifest, args = audited_fisher_run
    assert_inventory_is_the_writers_output(out, manifest, args)
    out, manifest, args = run(argv_for(arm_inputs, tmp_path / "mw", "--mode", "quick", "--engine", "audited",
                                       "--workers", "1"))
    assert_inventory_is_the_writers_output(out, manifest, args)


def test_full_released_rank_sum_with_stability_and_control(full_run):
    out, manifest, args = full_run
    assert_inventory_is_the_writers_output(out, manifest, args)
    registered = set(manifest["registered"])
    for rel in (f"figures/{ARM}/{ARM}_rank_comparison.xlsx", f"figures/{ARM}/{ARM}_figure_provenance_v43.md",
                f"figures/{ARM}/command_v43.log", f"figures/{ARM}/versions_v43.txt", "figures/selection_audit_v43.tsv",
                "figures/colour_audit_v43.tsv", "figures/figures_manifest_v43.tsv", "figures/positive_control_audit.tsv",
                f"stability/{ARM}_rank_stability.xlsx", f"summary/{ARM}/rbp_level.tsv"):
        assert rel in registered, rel


def test_raw_rmats_input_registers_the_gate_record(synthetic, tmp_path):  # noqa: F811
    _, genome_root = synthetic
    rows = []
    for ident, start, fdr, delta in (("up", 1000, "0.01", "0.1"), ("dn", 2000, "0.01", "-0.1"),
                                     ("bg", 3000, "0.9", "0.0")):
        rows.append({"ID": ident, "GeneID": "ENSG000001.1", "geneSymbol": "TEST", "chr": "chr1", "strand": "+",
                     "exonStart_0base": str(start), "exonEnd": str(start + 100), "upstreamES": str(start - 500),
                     "upstreamEE": str(start - 400), "downstreamES": str(start + 500),
                     "downstreamEE": str(start + 600), "FDR": fdr, "IncLevelDifference": delta,
                     "IncLevel1": {"up": "0.6,0.6", "dn": "0.4,0.4", "bg": "0.5,0.5"}[ident], "IncLevel2": "0.5,0.5", "IJC_SAMPLE_1": "10,10",
                     "SJC_SAMPLE_1": "10,10", "IJC_SAMPLE_2": "10,10", "SJC_SAMPLE_2": "10,10"})
    table = tmp_path / "SE.MATS.JC.txt"
    table.write_text("\t".join(rows[0]) + "\n" + "".join("\t".join(r.values()) + "\n" for r in rows))
    spec = tmp_path / "gates.json"
    spec.write_text(json.dumps({"fdr": 0.05, "min_abs_dpsi": 0.1, "min_jc_per_sample": 10, "bg_fdr_min": 0.5}))
    motifs = tmp_path / "motifs.tsv"
    motifs.write_text("Protein_name\tregularExpression\nTEST\tAA\nOTHER\tATAT\n")
    out, manifest, args = run(["--mode", "quick", "--arm", "SYNTH", "--out", str(tmp_path / "run"), "--rmats-se",
                               str(table), "--filter", str(spec), "--genome-root", str(genome_root), "--genome",
                               "synthetic", "--engine", "audited", "--stat-method", "fisher", "--known-motifs",
                               str(motifs), "--additional-motifs", "NA", "--window", "4", "--step", "2", "--intron",
                               "20", "--exon", "10", "--workers", "1", "--blas-threads", "1"])
    assert_inventory_is_the_writers_output(out, manifest, args, arm="SYNTH")
    assert {"event_set_config.json", "event_sets/counts.json", "logs/build_event_sets.stdout.log"} <= set(
        manifest["registered"])


# ------------------------------------------------------------------ the table cannot drift from the writers
def test_a_writer_output_missing_from_the_table_is_recorded_and_a_phantom_row_is_incomplete(arm_inputs, tmp_path,
                                                                                        monkeypatch):
    base = skill.REQUIRED_ARTIFACTS["quick"]["audited"]["fisher"]
    argv = argv_for(arm_inputs, tmp_path / "dropped", "--mode", "quick", "--engine", "audited", "--workers", "1",
                    "--stat-method", "fisher")
    monkeypatch.setattr(sys, "orig_argv", [sys.executable, "rmaps3_skill_run.py"], raising=False)
    monkeypatch.setenv("RMAPS_FORCE_MOTIF_FALLBACK", "1")
    monkeypatch.setitem(skill.REQUIRED_ARTIFACTS["quick"]["audited"], "fisher",
                        [r for r in base if r[1] != "quick_summary.xlsx"])
    assert skill.main(argv) == 0
    manifest = json.loads((tmp_path / "dropped" / "run_manifest.json").read_text())
    assert manifest["registered_not_required"] == ["quick_summary.xlsx"]
    assert "quick_summary.xlsx" in {row["path"] for row in manifest["inventory"]}, "still hashed"
    monkeypatch.setitem(skill.REQUIRED_ARTIFACTS["quick"]["audited"], "fisher",
                        base + [("phantom", "figures/never_written.tsv", None, None)])
    argv[argv.index("--out") + 1] = str(tmp_path / "phantom")
    assert skill.main(argv) == 4
    manifest = json.loads((tmp_path / "phantom" / "run_manifest.json").read_text())
    assert manifest["status"] == "INCOMPLETE" and manifest["missing"][0].startswith("phantom: 1 of 1 absent")


# ------------------------------------------------------------------ --verify after the run
def test_losing_any_registered_file_of_a_quick_run_makes_verify_incomplete(quick_released_run):
    out, manifest, _ = quick_released_run
    assert_every_lost_file_is_caught(out, manifest)


def test_losing_any_registered_file_of_a_full_run_makes_verify_incomplete(full_run):
    out, manifest, _ = full_run
    assert_every_lost_file_is_caught(out, manifest)


def test_verify_cli_reports_a_changed_file_and_exits_4(quick_released_run, tmp_path):
    import shutil
    out = tmp_path / "copy"
    shutil.copytree(quick_released_run[0], out)
    cli = [sys.executable, str(ROOT / "tools" / "rmaps3_skill_run.py"), "--verify", str(out)]
    proc = subprocess.run(cli, capture_output=True, text=True)
    assert proc.returncode == 0 and proc.stdout.startswith("VERIFY complete (exit 0)"), proc.stdout + proc.stderr
    target = out / "figures" / "positive_control_audit.tsv"
    target.write_text(target.read_text(encoding="utf-8") + "tampered\n", encoding="utf-8")
    proc = subprocess.run(cli, capture_output=True, text=True)
    assert proc.returncode == 4 and proc.stdout.startswith("VERIFY INCOMPLETE (exit 4)")
    assert "figures/positive_control_audit.tsv changed since the run" in proc.stdout
    (out / "run_manifest.json").unlink()
    assert skill.verify_run(out)[:2] == ("INCOMPLETE", 4)


def test_builder_manifest_lists_exactly_the_files_the_builder_wrote(calibrated_arm):
    import csv
    figures = calibrated_arm["figures"]
    with open(figures / "figures_manifest_v43.tsv", encoding="utf-8", newline="") as handle:
        listed = {Path(r["file"]).resolve() for r in csv.DictReader(handle, delimiter="\t") if r["status"] == "built"}
    disk = {p.resolve() for p in figures.rglob("*") if p.is_file() and p.name != "figures_manifest_v43.tsv"}
    assert listed == disk
