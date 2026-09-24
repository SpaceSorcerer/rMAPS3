"""Release documentation contract: cited documents are tracked, CI runs the lab suite, the report states HEAD's count."""
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs" / "CONSOLIDATION_REPORT_2026-09-24.md"
SUITE = ["tests", "--ignore=tests/legacy", "--ignore=tests/test_motif.py"]


def tracked(relative: str) -> bool:
    proc = subprocess.run(["git", "-C", str(ROOT), "ls-files", "--error-unmatch", relative],
                          capture_output=True, text=True)
    return proc.returncode == 0


def test_consolidation_report_cited_by_lessons_is_tracked():
    text = (ROOT / "LESSONS.md").read_text(encoding="utf-8")
    cited = re.findall(r"`([^`]*CONSOLIDATION_REPORT_2026-09-24\.md)`", text)
    assert cited, "LESSONS.md must cite the consolidation report"
    for path in cited:
        assert tracked(path), f"LESSONS.md cites {path}, which a clone does not contain"


def test_ci_lab_job_installs_the_lock_compiles_tools_and_runs_the_suite_as_steps():
    yaml = pytest.importorskip("yaml")          # PyYAML is pinned in requirements-lock.txt
    ci = yaml.safe_load((ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8"))
    steps = ci["jobs"]["lab-fork"]["steps"]
    runs = [(i, " ".join(str(s.get("run", "")).split())) for i, s in enumerate(steps)]
    setup = [s for s in steps if str(s.get("uses", "")).startswith("actions/setup-python")]
    assert setup and str(setup[0]["with"]["python-version"]) == "3.11"
    install = [i for i, r in runs if "pip install -r requirements-lock.txt" in r]
    compile_ = [i for i, r in runs if r == "python -m compileall -q tools"]
    suite = [i for i, r in runs if r.startswith("python -m pytest -q " + " ".join(SUITE))]
    assert install and compile_ and suite, runs
    assert install[0] < compile_[0] < suite[0]
    assert steps[suite[0]]["env"]["RMAPS_FORCE_MOTIF_FALLBACK"] == "1"


def unmet_lock_pins():
    """Pinned distributions this interpreter lacks or holds at another version."""
    import importlib.metadata as metadata
    unmet = []
    for line in (ROOT / "requirements-lock.txt").read_text(encoding="utf-8").splitlines():
        if line.strip() and not line.startswith("#"):
            name, version = line.strip().split("==")
            try:
                if metadata.version(name) != version:
                    unmet.append(line)
            except metadata.PackageNotFoundError:
                unmet.append(line)
    return unmet


def test_consolidation_report_states_the_test_count_pytest_collects_at_head():
    """The count is that of the lock environment: a pandas-less interpreter skips a module and collects fewer."""
    unmet = unmet_lock_pins()
    if unmet:
        pytest.skip("counts only in an environment that satisfies requirements-lock.txt; unmet: " + ", ".join(unmet))
    text = REPORT.read_text(encoding="utf-8")
    stated = re.search(r"collects \*\*(\d+) tests\*\*", text)
    assert stated, "the report must state 'collects **N tests**' for the current HEAD"
    proc = subprocess.run([sys.executable, "-m", "pytest", "--collect-only", "-q", "-p", "no:cacheprovider", *SUITE],
                          cwd=ROOT, capture_output=True, text=True, timeout=600)
    collected = re.search(r"(\d+) tests? collected", proc.stdout)
    assert collected, proc.stdout[-2000:] + proc.stderr[-2000:]
    assert int(stated.group(1)) == int(collected.group(1))


# ------------------------------------------------------------------ v1.0.1c: skills and docs match the wrapper
SKILLS = {name: Path.home() / ".claude" / "skills" / name / "SKILL.md" for name in ("rmaps3-quick", "rmaps3-full")}
GATE_SENTENCE = ("The event-set rule comes only from `--gate-rule` or the gate record's `rule` field, never from the "
                 "arm name; the wrapper prints it in `versions.txt`, the workbook README and the figure footer.")
TWO_STAGE_CLAUSE = "reports its stage-2 p only when its own stage-1 p"
RETIRED_PHRASES = (
    "arm-name suffix", "arm suffix", "text after the last `_`", "<name>_<rule>", "needs an arm named",
    "hard-codes the reli_v121", "is also promoted when any of its motifs", "exactly valid for p",
    "exact only at or below", "e:/rmaps_venv/", "e:\rmaps_venv\\", "group-triggered promotion is used",
)


def lab_docs():
    docs = {"README.md": ROOT / "README.md", "tools/README.md": ROOT / "tools" / "README.md"}
    docs.update({k: v for k, v in SKILLS.items() if v.is_file()})
    return {k: v.read_text(encoding="utf-8").replace("\r\n", "\n") for k, v in docs.items()}


def test_docs_and_skills_state_the_gate_rule_contract_and_no_retired_phrase():
    for name, text in lab_docs().items():
        assert GATE_SENTENCE in " ".join(text.split()), name
        low = text.lower()
        hits = [p for p in RETIRED_PHRASES if p in low]
        assert not hits, f"{name} still says {hits}"


def test_skills_run_venv2_state_4_3_3_the_two_stage_contract_statuses_and_exit_codes():
    present = {k: v for k, v in SKILLS.items() if v.is_file()}
    if not present:
        pytest.skip("the lab skills live outside the repository; absent in this checkout")
    for name, path in present.items():
        text = " ".join(path.read_text(encoding="utf-8").split())
        for needed in ("E:/rmaps_venv2/Scripts/python.exe", "4.3.3", "--gate-rule", "`complete` (exit 0)",
                       "`INCOMPLETE` (exit 4)", "`complete_partial` (exit 0)", "(exit 3)", "`failed` (exit 1)"):
            assert needed in text, (name, needed)
    full = " ".join(SKILLS["rmaps3-full"].read_text(encoding="utf-8").split()) if "rmaps3-full" in present else ""
    if full:
        assert TWO_STAGE_CLAUSE in full and "super-uniform at every α" in full


def test_two_stage_note_calls_tied_permutation_p_conservative_not_uniform():
    text = (ROOT / "docs" / "two_stage_validity.md").read_text(encoding="utf-8")
    assert "conservative (super-uniform), not uniform" in text
    assert "is then uniform, with ties counted against it" not in text
    assert TWO_STAGE_CLAUSE in text
