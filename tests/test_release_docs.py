"""Release documentation contract: every in-repo document the lessons cite is tracked by git."""
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


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


def test_ci_runs_the_lab_suite_and_compiles_tools():
    ci = (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    assert "python -m compileall -q tools" in ci
    assert "python -m pytest -q tests --ignore=tests/legacy" in ci
    assert "pip install -r requirements-lock.txt" in ci
