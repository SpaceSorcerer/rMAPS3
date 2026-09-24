"""Full-mode figure positive controls reach the wrapper's status and exit code (review 5 must-fix 1).

The v4 builder runs as a subprocess in full mode and writes one positive_control_audit.tsv row per layer x variant x
kind x panel. The wrapper reads that file back: a failed row on either layer gives complete_with_failed_positive_control,
exit 3, which --allow-partial never turns into 0; --positive-control-advisory keeps it advisory. The synthetic arm makes
PTBP1 rank 1 on every panel of both layers (PTBP1 is kept out of the broad-binder list, so both variants agree); the
raw-pass / calibrated-fail case rewrites the calibrated rows of the builder's audit as a failing builder would write them.
"""
import csv
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import rmaps3_skill_run as skill  # noqa: E402
from test_artifact_inventory import ARM, argv_for, arm_inputs, released_engine  # noqa: E402,F401

LAYERS = ("released_ranksum_rawP", "calibrated_ranksum")


def fail_layer_in_audit(audit: Path, layer: str) -> None:
    with open(audit, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, rows = reader.fieldnames, list(reader)
    assert {r["layer"] for r in rows} == set(LAYERS)
    for row in rows:
        if row["layer"] == layer:
            row.update({"first": "RBFOX2", "pass": "False"})
    with open(audit, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def full_run(arm_inputs, released_engine, tmp_path, monkeypatch, fail=None, *extra):  # noqa: F811
    real = skill.run_step

    def step(name, command, log, env, cwd, out):
        wall = real(name, command, log, env, cwd, out)
        if fail and name.startswith("figures_v"):
            fail_layer_in_audit(out / "figures" / "positive_control_audit.tsv", fail)
        return wall
    monkeypatch.setattr(skill, "run_step", step)
    monkeypatch.setattr(sys, "orig_argv", [sys.executable, "rmaps3_skill_run.py"], raising=False)
    monkeypatch.setenv("RMAPS_FORCE_MOTIF_FALLBACK", "1")
    broad = tmp_path / "broad_no_ptbp1.txt"
    broad.write_text("SRSF1\n")
    out = tmp_path / "run"
    code = skill.main(argv_for(arm_inputs, out, "--mode", "full", "--engine", "released", "--engine-root",
                               str(released_engine), "--permutations", "99", "--refine-perms", "199",
                               "--positive-control", "PTBP1", "--broad-binders-list", str(broad), *extra))
    return code, out, json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))


def test_raw_pass_calibrated_fail_is_exit_3_and_names_the_layer(arm_inputs, released_engine, tmp_path,  # noqa: F811
                                                                monkeypatch):
    code, out, manifest = full_run(arm_inputs, released_engine, tmp_path, monkeypatch, "calibrated_ranksum")
    assert (code, manifest["status"]) == (3, "complete_with_failed_positive_control"), manifest.get("error")
    assert manifest["missing"] == [] and all(r["pass"] for r in manifest["positive_control"])
    controls = manifest["figure_positive_control"]
    assert {(r["layer"], r["variant"], r["kind"]) for r in controls} == {
        (layer, variant, kind) for layer in LAYERS for variant in skill.FIGURE_VARIANTS for kind in skill.FIGURE_KINDS}
    assert all(r["pass"] for r in controls if r["layer"] == "released_ranksum_rawP")
    assert not any(r["pass"] for r in controls if r["layer"] == "calibrated_ranksum")
    assert manifest["failed_positive_control_layers"] == ["calibrated_ranksum"]
    index = (out / "index.html").read_text(encoding="utf-8")
    assert "Figure positive control PTBP1 FAILED on layer(s): calibrated_ranksum" in index
    assert "FAILED on layer(s): calibrated_ranksum" in (out / "command.log").read_text(encoding="utf-8")


def test_both_layers_pass_is_complete_exit_0(arm_inputs, released_engine, tmp_path, monkeypatch):  # noqa: F811
    code, out, manifest = full_run(arm_inputs, released_engine, tmp_path, monkeypatch)
    assert (code, manifest["status"]) == (0, "complete"), (manifest.get("error"), manifest.get("missing"))
    controls = manifest["figure_positive_control"]
    assert {r["layer"] for r in controls} == set(LAYERS) and all(r["pass"] for r in controls)
    assert manifest["failed_positive_control_layers"] == []
    assert "passed on every drawn layer" in (out / "index.html").read_text(encoding="utf-8")


def test_allow_partial_never_hides_a_failed_calibrated_control(arm_inputs, released_engine, tmp_path,  # noqa: F811
                                                               monkeypatch):
    code, _, manifest = full_run(arm_inputs, released_engine, tmp_path, monkeypatch, "calibrated_ranksum",
                                 "--allow-partial")
    assert (code, manifest["status"]) == (3, "complete_with_failed_positive_control")


@pytest.mark.parametrize("damage", ["absent", "one_layer"])
def test_an_absent_or_one_layer_audit_counts_as_a_failed_control(damage, tmp_path):
    audit = tmp_path / "positive_control_audit.tsv"
    if damage == "one_layer":
        audit.write_text("arm\tlayer\tvariant\tkind\tpanel\tsymbol\tfirst\tpass\n"
                         "A\treleased_ranksum_rawP\tmain\tbyRBP\tINCLUDED × Upstream Intron\tQKI\tQKI\tTrue\n",
                         encoding="utf-8")
    rows = skill.read_figure_controls(audit, "QKI")
    failed = [r for r in rows if not r["pass"]]
    assert failed and {r["layer"] for r in failed} == (
        {"calibrated_ranksum"} if damage == "one_layer" else set(LAYERS))
