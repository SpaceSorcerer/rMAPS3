"""Completeness of tools/rmaps3_skill_run.py = an exact, hashed inventory per mode (review 3 must-fix 1).

REQUIRED_ARTIFACTS[mode][engine][stat] lists every path a complete run must have. Status `complete` needs every
path to exist, be non-empty and carry its sha256 in run_manifest.json; one missing file, one empty file or a
partial archive conversion (fewer npz than root-table motifs) is not complete; --allow-partial gives
complete_partial, never complete.
"""
import hashlib
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

MOTIFS = ["QKI.ACUAAY", "RBFOX2.UGCAUG"]


def quick_args(tmp_path, genome_root, out, *extra):
    paths = input_paths(tmp_path)
    motifs = tmp_path / "motifs.tsv"
    motifs.write_text("Protein_name\tregularExpression\nTEST\tAA\nOTHER\tATAT\n")
    return ["--mode", "quick", "--arm", "SYNTH", "--out", str(out), "--up", str(paths["up"]),
            "--dn", str(paths["dn"]), "--bg", str(paths["bg"]), "--genome-root", str(genome_root),
            "--genome", "synthetic", "--engine", "audited", "--engine-root", str(ROOT), "--stat-method", "fisher",
            "--known-motifs", str(motifs), "--additional-motifs", "NA", "--window", "4", "--step", "2",
            "--intron", "20", "--exon", "10", "--workers", "1", "--blas-threads", "1", "--gate-rule", "A", *extra]


def parsed(mode="quick", *extra):
    return skill.build_parser().parse_args(["--mode", mode, "--arm", "S_A", "--out", "o", *extra])


CASES = {  # (mode, extra flags) -> one engine x statistic cell of REQUIRED_ARTIFACTS each
    "quick_released_ranksum": ("quick", ()),
    "quick_released_fisher": ("quick", ("--stat-method", "fisher")),
    "quick_audited": ("quick", ("--engine", "audited", "--stat-method", "fisher")),
    "full_released_ranksum": ("full", ()),
}


def populate(out, args):
    """A fake run directory holding every required path, each non-empty."""
    for final in (False, True):
        for _, rel in skill.required_paths(args, MOTIFS, final=final):
            path = out / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"x")


def status_of(out, args):
    rows, missing = skill.inventory(out, skill.required_paths(args, MOTIFS)
                                    + skill.required_paths(args, MOTIFS, final=True))
    return skill.run_status(missing, args.allow_partial, False), rows, missing


def test_the_table_covers_every_engine_and_statistic_quick_mode_accepts_and_full_mode_only_rank_sum():
    for engine in ("released", "audited"):
        for stat in ("fisher", "mannwhitney"):
            assert skill.REQUIRED_ARTIFACTS["quick"][engine][stat]
    assert list(skill.REQUIRED_ARTIFACTS["full"]) == ["released"]
    assert list(skill.REQUIRED_ARTIFACTS["full"]["released"]) == ["mannwhitney"]


def test_quick_inventory_names_archives_verification_figures_index_and_provenance():
    req = skill.required_artefacts(parsed())
    for needed in ("engine/S_A/pVal.up.vs.bg.RNAmap.txt", "quick_summary.xlsx", "counts/S_A/{motif}.counts.npz",
                   "counts/S_A/VERIFY.md", "figures/S_A_SE_byRBP_released_ranksum_rawP.png", "index.html",
                   "md5.txt", "versions.txt"):
        assert needed in req, needed
    assert not [r for r in skill.required_artefacts(parsed("quick", "--no-figures")) if r.startswith("figures/")]
    fisher = skill.required_artefacts(parsed("quick", "--stat-method", "fisher"))
    assert "counts/S_A/{motif}.counts.npz" in fisher and "counts/S_A/VERIFY.md" not in fisher
    audited = skill.required_artefacts(parsed("quick", "--engine", "audited"))
    assert "engine/S_A/positional/{motif}.{region}.hits.npz" in audited
    assert not [r for r in audited if r.startswith(("counts/", "figures/"))]


def test_full_inventory_adds_calibration_per_unit_rank_workbook_supplement_figures_and_readouts():
    req = skill.required_artefacts(parsed("full"))
    for needed in ("summary/S_A/S_A_calibrated_ranksum_v2.xlsx", "summary/S_A/rbp_level.tsv",
                   "summary/S_A/readout.md", "summary_rowunit/S_A/readout.md", "figures/index.html",
                   "figures/S_A/S_A_rank_comparison.xlsx",
                   "figures/S_A/S_A_SE_byMotif_calibrated_ranksum_noSpliceosome_noBroad.svg",
                   "counts/S_A/{motif}.counts.npz", "md5.txt"):
        assert needed in req, needed
    assert "stability/S_A_rank_stability.xlsx" not in req
    assert "stability/S_A_rank_stability.xlsx" in skill.required_artefacts(
        parsed("full", "--rank-stability-run", "x"))


@pytest.mark.parametrize("case", list(CASES))
def test_a_fully_populated_run_is_complete_and_every_path_is_hashed(case, tmp_path):
    mode, extra = CASES[case]
    args = parsed(mode, *extra)
    populate(tmp_path, args)
    (status, code), rows, missing = status_of(tmp_path, args)
    assert (status, code, missing) == ("complete", 0, [])
    expected = [rel for _, rel in skill.required_paths(args, MOTIFS) + skill.required_paths(args, MOTIFS, final=True)]
    assert [r["path"] for r in rows] == expected
    assert all(r["sha256"] == hashlib.sha256(b"x").hexdigest() for r in rows)


@pytest.mark.parametrize("case", list(CASES))
@pytest.mark.parametrize("damage", ["missing", "empty", "partial_archive"])
def test_one_missing_one_empty_or_a_partial_archive_set_is_not_complete(case, damage, tmp_path):
    mode, extra = CASES[case]
    args = parsed(mode, *extra)
    populate(tmp_path, args)
    required = [rel for _, rel in skill.required_paths(args, MOTIFS)]
    if damage == "partial_archive":   # one motif's archive(s) absent: fewer archives than root-table motifs
        victims = [r for r in required if MOTIFS[1] in r and (r.endswith(".counts.npz") or r.endswith(".hits.npz"))]
        assert victims
        for rel in victims:
            (tmp_path / rel).unlink()
    elif damage == "missing":
        (tmp_path / required[len(required) // 2]).unlink()
    else:
        (tmp_path / required[-1]).write_bytes(b"")
    (status, code), _, missing = status_of(tmp_path, args)
    assert (status, code) == ("INCOMPLETE", 4) and missing
    args.allow_partial = True
    assert status_of(tmp_path, args)[0] == ("complete_partial", 0)


def test_missing_lines_count_the_absent_share_of_a_partial_conversion(tmp_path):
    args = parsed("quick", "--stat-method", "fisher")
    populate(tmp_path, args)
    (tmp_path / "counts" / "S_A" / f"{MOTIFS[0]}.counts.npz").unlink()
    _, _, missing = status_of(tmp_path, args)
    assert missing == [f"count archive: 1 of 2 absent or empty (counts/S_A/{MOTIFS[0]}.counts.npz absent)"]


def test_audited_fisher_quick_run_is_complete_and_hashes_every_required_path(synthetic, tmp_path):  # noqa: F811
    _, genome_root = synthetic
    out = tmp_path / "run"
    result = run_skill(quick_args(tmp_path, genome_root, out), tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["status"] == "complete" and manifest["missing"] == []
    hashed = {row["path"]: row["sha256"] for row in manifest["inventory"]}
    assert {"md5.txt", "command.log", "index.html", "versions.txt"} <= set(hashed)
    assert len([p for p in hashed if p.endswith(".hits.npz")]) == 2 * 8   # 2 motifs x 8 sub-regions
    for rel, digest in hashed.items():
        assert hashlib.sha256((out / rel).read_bytes()).hexdigest() == digest, rel
    index = (out / "index.html").read_text(encoding="utf-8")
    assert "INCOMPLETE" not in index and "positional/{motif}.{region}.hits.npz" in index


def lose_one_positional(monkeypatch):
    engine = skill.run_engine

    def run(*a, **k):
        result = engine(*a, **k)
        next(iter(sorted((result[0] / "positional").glob("*.hits.npz")))).unlink()
        return result
    monkeypatch.setattr(skill, "run_engine", run)
    monkeypatch.setattr(sys, "orig_argv", [sys.executable, "rmaps3_skill_run.py"], raising=False)
    monkeypatch.setenv("RMAPS_FORCE_MOTIF_FALLBACK", "1")


def test_one_lost_positional_archive_makes_the_run_incomplete_exit_4(synthetic, tmp_path, monkeypatch):  # noqa: F811
    _, genome_root = synthetic
    lose_one_positional(monkeypatch)
    out = tmp_path / "run"
    assert skill.main(quick_args(tmp_path, genome_root, out)) == 4
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["status"] == "INCOMPLETE"
    assert [m.split(":")[0] for m in manifest["missing"]] == ["audited positional archive"]
    assert "1 of 16 absent" in manifest["missing"][0]
    assert "INCOMPLETE" in (out / "index.html").read_text(encoding="utf-8")


def test_allow_partial_turns_it_into_complete_partial_never_complete(synthetic, tmp_path, monkeypatch):  # noqa: F811
    _, genome_root = synthetic
    lose_one_positional(monkeypatch)
    out = tmp_path / "run"
    assert skill.main(quick_args(tmp_path, genome_root, out, "--allow-partial")) == 0
    assert json.loads((out / "run_manifest.json").read_text())["status"] == "complete_partial"


@pytest.mark.parametrize("stat", ["fisher", "mannwhitney"])
def test_figures_are_not_promised_where_they_cannot_be_drawn(stat):
    args = parsed("quick", "--engine", "audited", "--stat-method", stat)
    assert skill.figure_skip_reason(args, {"verified": False})
    assert not [r for r in skill.required_artefacts(args) if r.startswith("figures/")]
