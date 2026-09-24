"""Effective gate rule versus the arm name in tools/rmaps3_skill_run.py (review finding 4)."""
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import rmaps3_skill_run as skill  # noqa: E402

RAW = ("--rmats-se", "SE.MATS.JC.txt", "--filter", "gates.json")
PRESPLIT = ("--up", "u", "--dn", "d", "--bg", "b", "--gate-counts", "c.json")


def args(arm, *extra):
    return skill.build_parser().parse_args(["--mode", "quick", "--arm", arm, "--out", "o", *extra])


@pytest.mark.parametrize("arm", ["QKI_KO_B", "QKI_KO_Beffect"])
def test_raw_rmats_input_for_a_rule_b_arm_is_refused(arm):
    with pytest.raises(ValueError, match="Rule B sets are frozen concordant files: supply pre-split inputs"):
        skill.validate(args(arm, *RAW))
    with pytest.raises(ValueError):
        skill.validate(args(arm, *RAW, "--gate-rule", "A"))


def test_raw_rmats_input_needs_rule_a_when_a_rule_is_stated():
    skill.validate(args("QKI_KO_A", *RAW, "--gate-rule", "A"))
    skill.validate(args("QKI_KO_A", *RAW))
    with pytest.raises(ValueError, match="implements rule A only"):
        skill.validate(args("QKI_KO_A", *RAW, "--gate-rule", "B"))
    with pytest.raises(ValueError, match="disagrees with the arm name"):
        skill.validate(args("QKI_KO_X", *RAW, "--gate-rule", "A"))


def test_presplit_rule_b_is_accepted_and_must_match_its_arm():
    skill.validate(args("QKI_KO_B", *PRESPLIT, "--gate-rule", "B"))
    skill.validate(args("QKI_KO_B", *PRESPLIT))
    with pytest.raises(ValueError, match="disagrees with the arm name"):
        skill.validate(args("QKI_KO_B", *PRESPLIT, "--gate-rule", "A"))


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
