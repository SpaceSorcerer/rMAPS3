"""Quick mode may not exit 0 as complete without its promised archives and figures (review finding 6)."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import rmaps3_skill_run as skill  # noqa: E402
from test_engine_synthetic import synthetic, input_paths  # noqa: E402,F401
from test_skill_run import run_skill  # noqa: E402


def quick_args(tmp_path, genome_root, out, *extra):
    paths = input_paths(tmp_path)
    motifs = tmp_path / "motifs.tsv"
    motifs.write_text("Protein_name\tregularExpression\nTEST\tAA\nOTHER\tATAT\n")
    return ["--mode", "quick", "--arm", "SYNTH", "--out", str(out), "--up", str(paths["up"]),
            "--dn", str(paths["dn"]), "--bg", str(paths["bg"]), "--genome-root", str(genome_root),
            "--genome", "synthetic", "--engine", "audited", "--engine-root", str(ROOT), "--stat-method", "fisher",
            "--known-motifs", str(motifs), "--additional-motifs", "NA", "--window", "4", "--step", "2",
            "--intron", "20", "--exon", "10", "--workers", "1", "--blas-threads", "1", *extra]


def test_quick_run_without_archives_or_figures_is_incomplete_and_exits_4(synthetic, tmp_path):  # noqa: F811
    _, genome_root = synthetic
    out = tmp_path / "run"
    result = run_skill(quick_args(tmp_path, genome_root, out), tmp_path)
    assert result.returncode == 4, result.stdout + result.stderr
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["status"] == "INCOMPLETE"
    assert any(m.startswith("verified count archives") for m in manifest["missing"])
    assert any(m.startswith("main-layer figures") for m in manifest["missing"])
    assert "INCOMPLETE" in (out / "index.html").read_text(encoding="utf-8")


def test_no_figures_removes_only_the_figure_promise(synthetic, tmp_path):  # noqa: F811
    _, genome_root = synthetic
    out = tmp_path / "run"
    result = run_skill(quick_args(tmp_path, genome_root, out, "--no-figures"), tmp_path)
    assert result.returncode == 4
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert [m.split(":")[0] for m in manifest["missing"]] == ["verified count archives"]


def test_missing_deliverables_is_empty_for_a_verified_run_with_figures():
    args = skill.build_parser().parse_args(["--mode", "quick", "--arm", "S_A", "--out", "o"])
    assert skill.missing_deliverables(args, {"verified": True}, {"figures": [Path("f.png")]}) == []
    assert skill.missing_deliverables(args, {"verified": True}, {"skipped": "reason"}) == [
        "main-layer figures: reason"]
