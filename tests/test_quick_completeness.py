"""Mode-aware completeness of tools/rmaps3_skill_run.py (review 1 finding 6; review 2 regression).

Required artefacts per engine x statistic (beyond both root tables and quick_summary.xlsx):
released + mannwhitney = verified count archives, plus main-layer figures in quick mode unless --no-figures;
released + fisher = converted count archives (a Fisher run cannot be verified);
audited = the engine's positional/*.hits.npz.
"""
import json
import sys
from pathlib import Path

import pytest

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


def parsed(*extra):
    return skill.build_parser().parse_args(["--mode", "quick", "--arm", "S_A", "--out", "o", *extra])


def test_required_artefacts_are_defined_per_engine_and_statistic():
    assert skill.required_artefacts(parsed()) == ["verified count archives (VERIFY.md)", "main-layer figures"]
    assert skill.required_artefacts(parsed("--no-figures")) == ["verified count archives (VERIFY.md)"]
    assert skill.required_artefacts(parsed("--stat-method", "fisher"))[0].startswith("converted count archives")
    for stat in ("fisher", "mannwhitney"):
        assert skill.required_artefacts(parsed("--engine", "audited", "--stat-method", stat)) == [
            "audited positional archives (engine positional/*.hits.npz)"]


def test_released_rank_sum_needs_verified_archives_and_figures():
    args = parsed()
    assert skill.missing_deliverables(args, {"verified": True}, {"figures": [Path("f.png")]}) == []
    assert skill.missing_deliverables(args, {"verified": True}, {"skipped": "reason"}) == [
        "main-layer figures: reason"]
    missing = skill.missing_deliverables(args, {"converted": True, "verified": False, "reason": "r"}, None)
    assert [m.split(":")[0] for m in missing] == ["verified count archives", "main-layer figures"]


def test_released_fisher_needs_converted_archives_only():
    args = parsed("--stat-method", "fisher")
    # verification cannot run on a Fisher run and figures draw the rank-sum only: neither is missing
    assert skill.missing_deliverables(args, {"converted": True, "verified": False}, {"skipped": "fisher"}) == []
    missing = skill.missing_deliverables(args, {"converted": False, "reason": "no countDist"}, None)
    assert missing == ["converted count archives: no countDist"]


def test_audited_needs_its_positional_archives(tmp_path):
    args = parsed("--engine", "audited", "--stat-method", "fisher")
    (tmp_path / "positional").mkdir()
    assert skill.missing_deliverables(args, {"verified": False}, None, tmp_path)[0].startswith(
        "audited positional archives")
    (tmp_path / "positional" / "M.AA.UpstreamIntron.hits.npz").write_bytes(b"x")
    assert skill.missing_deliverables(args, {"verified": False}, {"skipped": "audited"}, tmp_path) == []


def test_audited_fisher_quick_run_is_complete_without_allow_partial(synthetic, tmp_path):  # noqa: F811
    _, genome_root = synthetic
    out = tmp_path / "run"
    result = run_skill(quick_args(tmp_path, genome_root, out), tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["status"] == "complete" and manifest["missing"] == []
    assert manifest["required_artefacts"] == ["audited positional archives (engine positional/*.hits.npz)"]
    index = (out / "index.html").read_text(encoding="utf-8")
    assert "INCOMPLETE" not in index and "audited positional archives" in index


def test_audited_run_without_its_positional_archives_is_incomplete_and_exits_4(synthetic, tmp_path,  # noqa: F811
                                                                               monkeypatch):
    _, genome_root = synthetic
    out = tmp_path / "run"
    engine = skill.run_engine

    def engine_then_lose_positional(*a, **k):
        result = engine(*a, **k)
        positional = result[0] / "positional"
        positional.rename(positional.with_name("positional_moved_by_test"))
        return result
    monkeypatch.setattr(skill, "run_engine", engine_then_lose_positional)
    monkeypatch.setattr(sys, "orig_argv", [sys.executable, "rmaps3_skill_run.py"], raising=False)
    monkeypatch.setenv("RMAPS_FORCE_MOTIF_FALLBACK", "1")
    assert skill.main(quick_args(tmp_path, genome_root, out)) == 4
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["status"] == "INCOMPLETE"
    assert [m.split(":")[0] for m in manifest["missing"]] == ["audited positional archives"]
    assert "INCOMPLETE" in (out / "index.html").read_text(encoding="utf-8")


@pytest.mark.parametrize("stat", ["fisher", "mannwhitney"])
def test_figures_are_not_promised_where_they_cannot_be_drawn(stat):
    args = parsed("--engine", "audited", "--stat-method", stat)
    assert skill.figure_skip_reason(args, {"verified": False})
    assert "main-layer figures" not in skill.required_artefacts(args)
