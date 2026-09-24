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
