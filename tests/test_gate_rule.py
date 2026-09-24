"""Effective gate rule in tools/rmaps3_skill_run.py (review 1 finding 4; review 2 regression on *_B names).

The rule comes from --gate-rule or the gate record (--gate-counts 'rule'); the arm name is a label read only when
neither states one. Raw --rmats-se input is rule A: a *_B / *_Beffect arm needs --gate-rule A to run on it.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import build_region_lollipops_v4 as lol  # noqa: E402
import rmaps3_skill_run as skill  # noqa: E402

RAW = ("--rmats-se", "SE.MATS.JC.txt", "--filter", "gates.json")
PRESPLIT = ("--up", "u", "--dn", "d", "--bg", "b", "--gate-counts", "c.json")


def args(arm, *extra):
    return skill.build_parser().parse_args(["--mode", "quick", "--arm", arm, "--out", "o", *extra])


def record(tmp_path, rule):
    path = tmp_path / "counts.json"
    body = {"n_up": 1, "n_dn": 1, "n_bg": 1, "n_expr_unknown_in_fg": 0, "n_expr_unknown_in_bg": 0}
    if rule is not None:
        body["rule"] = rule
    path.write_text(json.dumps(body), encoding="utf-8")
    return ("--up", "u", "--dn", "d", "--bg", "b", "--gate-counts", str(path))


@pytest.mark.parametrize("arm", ["QKI_KO_B", "QKI_KO_Beffect"])
def test_raw_input_for_a_rule_b_name_is_refused_without_an_explicit_rule_a(arm):
    with pytest.raises(ValueError, match="Rule B sets are frozen concordant files: supply pre-split inputs"):
        skill.validate(args(arm, *RAW))
    parsed = args(arm, *RAW, "--gate-rule", "A")
    skill.validate(parsed)                                    # the suffix is only a name once rule A is stated
    assert skill.check_gate_rule(parsed, True) == ("A", "--gate-rule A")


def test_raw_input_is_rule_a_only_and_any_name_runs_under_it():
    skill.validate(args("QKI_KO_A", *RAW, "--gate-rule", "A"))
    skill.validate(args("QKI_KO_A", *RAW))
    skill.validate(args("QKI_KO_X", *RAW, "--gate-rule", "A"))
    assert skill.check_gate_rule(args("MYRUN", *RAW), True) == ("A", "portable rule-A builder (--rmats-se)")
    for rule in ("B", "Beffect"):
        with pytest.raises(ValueError, match="implements rule A only"):
            skill.validate(args("QKI_KO_A", *RAW, "--gate-rule", rule))


def test_presplit_rule_comes_from_gate_rule_then_record_then_name(tmp_path):
    skill.validate(args("QKI_KO_B", *PRESPLIT, "--gate-rule", "B"))
    parsed = args("QKI_KO_B", *PRESPLIT, "--gate-rule", "A")     # a name ending _B is not a rule
    skill.validate(parsed)
    assert skill.check_gate_rule(parsed, False) == ("A", "--gate-rule")
    parsed = args("QKI_KO_A", *record(tmp_path, "B"))              # the record beats the name
    skill.validate(parsed)
    assert skill.check_gate_rule(parsed, False)[0] == "B"
    assert skill.check_gate_rule(args("QKI_KO_B", *record(tmp_path, None)), False)[0] == "B"
    assert skill.check_gate_rule(args("SYNTH", *record(tmp_path, None)), False) == (None, None)


def test_gate_rule_that_disagrees_with_the_gate_record_is_refused(tmp_path):
    with pytest.raises(ValueError, match="disagrees with the gate record"):
        skill.validate(args("QKI_KO_B", *record(tmp_path, "B"), "--gate-rule", "A"))
    with pytest.raises(ValueError, match="not one of"):
        skill.validate(args("QKI_KO_B", *record(tmp_path, "C")))


def test_figures_need_a_stated_rule(tmp_path):
    with pytest.raises(ValueError, match="none is stated"):
        skill.validate(args("SYNTH", *record(tmp_path, None)))
    skill.validate(args("SYNTH", *record(tmp_path, None), "--no-figures"))
    skill.validate(args("SYNTH", *record(tmp_path, None), "--gate-rule", "A"))


def test_figure_footer_prints_the_effective_rule_not_the_name_suffix(capsys):
    counts = {"n_up": 1, "n_dn": 1, "n_bg": 1, "n_expr_unknown_in_fg": 0, "n_expr_unknown_in_bg": 0}
    texts = {"gate_text": None, "rule": "A"}
    assert "rule A" in " ".join(lol.gate_lines(texts, "QKI_KO_B", counts))
    assert "rule B" in " ".join(lol.gate_lines({"gate_text": None, "rule": None}, "QKI_KO_B", counts))
    with pytest.raises(SystemExit):
        lol.main(["--arms", "X_A", "Y_A", "--out-root", "o", "--released-root", "r", "--calibrated-root", "c",
                  "--gtf", "g", "--spliceosome-list", "s", "--broad-binders-list", "b", "--gate-rule", "A"])
    assert "give one arm per call" in capsys.readouterr().err


def test_filter_spec_asking_for_rule_b_is_refused(tmp_path):
    spec = tmp_path / "gates.json"
    spec.write_text('{"fdr": 0.05, "rule": "B"}', encoding="utf-8")
    parsed = skill.build_parser().parse_args(["--mode", "quick", "--arm", "X_A", "--out", str(tmp_path),
                                              "--rmats-se", "se.txt", "--filter", str(spec)])
    with pytest.raises(ValueError, match="implements rule A only"):
        skill.build_event_sets(parsed, tmp_path, lambda m: None, {})


def test_cli_refuses_raw_input_for_a_rule_b_arm_before_writing(tmp_path):
    out = tmp_path / "run"
    proc = subprocess.run([sys.executable, str(ROOT / "tools" / "rmaps3_skill_run.py"), "--mode", "quick",
                           "--arm", "QKI_KO_B", "--out", str(out), *RAW, "--genome-root", "g",
                           "--engine-root", "e", "--no-figures"],
                          cwd=ROOT, capture_output=True, text=True, timeout=300, env=dict(os.environ))
    assert proc.returncode == 2
    assert "supply pre-split inputs" in proc.stderr
    assert not out.exists()
