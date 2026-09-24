"""Run status precedence and the md5 record of tools/rmaps3_skill_run.py (review 2 regressions)."""
import hashlib
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import rmaps3_skill_run as skill  # noqa: E402
from test_engine_synthetic import synthetic  # noqa: E402,F401
from test_quick_completeness import quick_args  # noqa: E402


@pytest.mark.parametrize("missing, allow_partial, control_failed, failed, expected", [
    ([], False, False, False, ("complete", 0)),
    (["x"], True, False, False, ("complete_partial", 0)),
    (["x"], False, False, False, ("INCOMPLETE", 4)),
    ([], False, True, False, ("complete_with_failed_positive_control", 3)),
    (["x"], True, True, False, ("complete_with_failed_positive_control", 3)),   # --allow-partial never hides it
    (["x"], False, True, False, ("INCOMPLETE", 4)),
    (["x"], True, True, True, ("failed", 1)),
])
def test_status_precedence_table(missing, allow_partial, control_failed, failed, expected):
    assert skill.run_status(missing, allow_partial, control_failed, failed) == expected


def test_precedence_table_is_ordered_most_severe_first():
    assert [s for s, _, _ in skill.STATUS_PRECEDENCE] == [
        "failed", "INCOMPLETE", "complete_with_failed_positive_control", "complete_partial", "complete"]


def in_process(monkeypatch):
    monkeypatch.setattr(sys, "orig_argv", [sys.executable, "rmaps3_skill_run.py"], raising=False)
    monkeypatch.setenv("RMAPS_FORCE_MOTIF_FALLBACK", "1")


def test_allow_partial_keeps_a_failed_positive_control(synthetic, tmp_path, monkeypatch):  # noqa: F811
    _, genome_root = synthetic
    engine = skill.run_engine

    def engine_then_lose_positional(*a, **k):
        result = engine(*a, **k)
        (result[0] / "positional").rename(result[0] / "positional_moved_by_test")
        return result
    monkeypatch.setattr(skill, "run_engine", engine_then_lose_positional)
    in_process(monkeypatch)
    out = tmp_path / "run"
    code = skill.main(quick_args(tmp_path, genome_root, out, "--allow-partial", "--positive-control", "NOT_AN_RBP"))
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert (manifest["status"], code) == ("complete_with_failed_positive_control", 3)
    assert manifest["missing"]                                   # the missing artefact is still recorded


def test_md5_record_matches_every_file_as_left_including_the_closed_command_log(synthetic, tmp_path,  # noqa: F811
                                                                              monkeypatch):
    _, genome_root = synthetic
    in_process(monkeypatch)
    out = tmp_path / "run"
    assert skill.main(quick_args(tmp_path, genome_root, out)) == 0
    log = (out / "command.log").read_text(encoding="utf-8").splitlines()
    assert " exit=0 " in log[-1] + " "
    recorded = {}
    for line in (out / "md5.txt").read_text(encoding="utf-8").splitlines()[1:]:
        digest, size, path = line.split("\t", 2)
        recorded[Path(path)] = (digest, int(size))
    assert out / "command.log" in recorded
    for path, (digest, size) in recorded.items():
        data = path.read_bytes()
        assert (hashlib.md5(data).hexdigest(), len(data)) == (digest, size), path
    # run_manifest.json is written last and carries the sha256 of md5.txt; it cannot carry its own hash
    on_disk = {p for p in out.rglob("*") if p.is_file() and p.name not in ("md5.txt", "run_manifest.json")}
    assert on_disk == set(recorded)
    hashed = {row["path"]: row["sha256"] for row in json.loads((out / "run_manifest.json").read_text())["inventory"]}
    assert hashed["md5.txt"] == hashlib.sha256((out / "md5.txt").read_bytes()).hexdigest()
    assert hashed["command.log"] == hashlib.sha256((out / "command.log").read_bytes()).hexdigest()
