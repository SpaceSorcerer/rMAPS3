"""Gate rule in tools/rmaps3_skill_run.py and tools/build_region_lollipops_v4.py (review 3 must-fix 2).

The rule comes only from --gate-rule or the gate record ('rule' of --gate-counts / --counts-json). The arm name is
never read as a rule, in the wrapper or the figure builder. Raw --rmats-se input is rule A (the portable builder).
The effective rule is printed in versions.txt, the workbook README and the figure footer.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_region_lollipops_v4 as lol  # noqa: E402
import rmaps3_skill_run as skill  # noqa: E402

RAW = ("--rmats-se", "SE.MATS.JC.txt", "--filter", "gates.json")
PRESPLIT = ("--up", "u", "--dn", "d", "--bg", "b", "--gate-counts", "c.json")
COUNTS = {"n_up": 1, "n_dn": 1, "n_bg": 1, "n_expr_unknown_in_fg": 0, "n_expr_unknown_in_bg": 0}


def args(arm, *extra):
    return skill.build_parser().parse_args(["--mode", "quick", "--arm", arm, "--out", "o", "--arm-label", "test arm",
                                            *extra])


def record(tmp_path, rule):
    path = tmp_path / "counts.json"
    body = dict(COUNTS)
    if rule is not None:
        body["rule"] = rule
    path.write_text(json.dumps(body), encoding="utf-8")
    return ("--up", "u", "--dn", "d", "--bg", "b", "--gate-counts", str(path))


@pytest.mark.parametrize("arm", ["QKI_KO_B", "QKI_KO_Beffect"])
def test_a_rule_b_named_raw_arm_runs_under_rule_a(arm):
    parsed = args(arm, *RAW, "--gate-rule", "A")
    skill.validate(parsed)
    assert skill.check_gate_rule(parsed, True) == ("A", "--gate-rule A")
    parsed = args(arm, *RAW)                                  # the name is a label; raw input is rule A
    skill.validate(parsed)
    assert skill.check_gate_rule(parsed, True) == ("A", "portable rule-A builder (--rmats-se)")


def test_raw_input_refuses_a_stated_rule_b():
    for rule in ("B", "Beffect"):
        with pytest.raises(ValueError, match="implements rule A only"):
            skill.validate(args("QKI_KO_A", *RAW, "--gate-rule", rule))


def test_presplit_rule_comes_from_gate_rule_or_record_never_the_name(tmp_path):
    parsed = args("QKI_KO_B", *PRESPLIT, "--gate-rule", "A")
    skill.validate(parsed)
    assert skill.check_gate_rule(parsed, False) == ("A", "--gate-rule")
    parsed = args("QKI_KO_A", *record(tmp_path, "B"))
    skill.validate(parsed)
    assert skill.check_gate_rule(parsed, False)[0] == "B"
    for arm in ("QKI_KO_B", "QKI_KO_Beffect", "QKI_KO_A", "SYNTH"):
        with pytest.raises(ValueError, match="no event-set rule is stated.*never read as a rule"):
            skill.check_gate_rule(args(arm, *record(tmp_path, None)), False)
        with pytest.raises(ValueError, match="no event-set rule is stated"):
            skill.validate(args(arm, *record(tmp_path, None), "--no-figures"))


def test_gate_rule_that_disagrees_with_the_gate_record_is_refused(tmp_path):
    with pytest.raises(ValueError, match="disagrees with the gate record"):
        skill.validate(args("QKI_KO_B", *record(tmp_path, "B"), "--gate-rule", "A"))
    with pytest.raises(ValueError, match="not one of"):
        skill.validate(args("QKI_KO_B", *record(tmp_path, "C")))


def test_no_suffix_logic_is_left_in_either_tool():
    wrapper = (ROOT / "tools" / "rmaps3_skill_run.py").read_text(encoding="utf-8")
    builder = (ROOT / "tools" / "build_region_lollipops_v4.py").read_text(encoding="utf-8")
    assert "RULE_B_SUFFIXES" not in wrapper and "def arm_rule" not in wrapper
    assert "arm.rsplit('_', 1)" not in builder and 'arm.rsplit("_", 1)' not in builder


def test_builder_footer_prints_the_stated_rule_and_refuses_none():
    texts = {"gate_text": None, "rule": "A"}
    assert "rule A" in " ".join(lol.gate_lines(texts, "QKI_KO_B", COUNTS))
    with pytest.raises(ValueError, match="no event-set rule stated for QKI_KO_B"):
        lol.gate_lines({"gate_text": None, "rule": None}, "QKI_KO_B", COUNTS)
    with pytest.raises(ValueError, match="no event-set rule stated"):
        lol.gate_lines({"gate_text": "n={n_up} rule {rule}", "rule": None}, "X_B", COUNTS)
    assert lol.gate_lines({"gate_text": "jc10 n={n_up}", "rule": None}, "X_B", COUNTS) == ["jc10 n=1"]


def test_builder_rule_from_flag_or_counts_record(tmp_path):
    path = tmp_path / "counts.json"
    path.write_text(json.dumps(dict(COUNTS, rule="Beffect")), encoding="utf-8")
    assert lol.arm_gate_rule("QKI_KO_A", None, str(path)) == "Beffect"
    assert lol.arm_gate_rule("QKI_KO_A", "Beffect", str(path)) == "Beffect"
    with pytest.raises(ValueError, match="disagrees with the gate record"):
        lol.arm_gate_rule("QKI_KO_A", "A", str(path))
    path.write_text(json.dumps(COUNTS), encoding="utf-8")
    assert lol.arm_gate_rule("QKI_KO_B", None, str(path)) is None
    with pytest.raises(ValueError, match="pass --counts-json QKI_KO_B=<path>"):
        lol.gate_counts("QKI_KO_B", None)


def test_builder_gate_rule_flag_is_per_arm(capsys):
    base = ["--out-root", "o", "--released-root", "r", "--calibrated-root", "c", "--gtf", "g",
            "--spliceosome-list", "s", "--broad-binders-list", "b"]
    with pytest.raises(SystemExit):
        lol.main(["--arms", "X_A", "Y_A", *base, "--gate-rule", "A"])
    assert "give one arm per call or ARM=RULE" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        lol.main(["--arms", "X_A", *base, "--gate-rule", "Z_A=A"])
    assert "expected [ARM=]RULE" in capsys.readouterr().err


def test_filter_spec_asking_for_rule_b_is_refused(tmp_path):
    spec = tmp_path / "gates.json"
    spec.write_text('{"fdr": 0.05, "rule": "B"}', encoding="utf-8")
    parsed = skill.build_parser().parse_args(["--mode", "quick", "--arm", "X_A", "--out", str(tmp_path),
                                              "--rmats-se", "se.txt", "--filter", str(spec)])
    with pytest.raises(ValueError, match="implements rule A only"):
        skill.build_event_sets(parsed, tmp_path, lambda m: None, {})


def test_cli_refuses_presplit_input_with_no_rule_before_writing(tmp_path):
    out = tmp_path / "run"
    proc = subprocess.run([sys.executable, str(ROOT / "tools" / "rmaps3_skill_run.py"), "--mode", "quick",
                           "--arm", "QKI_KO_B", "--out", str(out), "--up", "u", "--dn", "d", "--bg", "b",
                           "--genome-root", "g", "--engine-root", "e", "--no-figures"],
                          cwd=ROOT, capture_output=True, text=True, timeout=300, env=dict(os.environ))
    assert proc.returncode == 2
    assert "no event-set rule is stated for pre-split input: pass --gate-rule" in proc.stderr
    assert not out.exists()


def test_the_rule_reaches_versions_readme_and_event_counts(synthetic, tmp_path, monkeypatch):  # noqa: F811
    from test_quick_completeness import quick_args
    import openpyxl
    _, genome_root = synthetic
    monkeypatch.setattr(sys, "orig_argv", [sys.executable, "rmaps3_skill_run.py"], raising=False)
    monkeypatch.setenv("RMAPS_FORCE_MOTIF_FALLBACK", "1")
    out = tmp_path / "run"
    argv = quick_args(tmp_path, genome_root, out)
    argv[argv.index("--gate-rule") + 1] = "Beffect"
    argv[argv.index("--arm") + 1] = "SYNTH_A"                 # a rule-looking suffix that must be ignored
    assert skill.main(argv) == 0
    versions = (out / "versions.txt").read_text(encoding="utf-8")
    assert "gate_rule\tBeffect\t--gate-rule" in versions
    book = openpyxl.load_workbook(out / "quick_summary.xlsx")
    readme = {r[0]: r[1] for r in book["README"].iter_rows(min_row=2, values_only=True)}
    assert readme["Event-set rule"].startswith("rule Beffect (from --gate-rule")
    assert json.loads((out / "event_counts.json").read_text())["rule"] == "Beffect"


from test_engine_synthetic import synthetic  # noqa: E402,F401
