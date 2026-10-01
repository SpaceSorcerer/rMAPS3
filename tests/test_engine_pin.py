"""Released-engine pin, reference hashes and the pinned environment (review finding 7)."""
import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import rmaps3_skill_run as skill  # noqa: E402


def git(repo, *command):
    proc = subprocess.run(["git", "-C", str(repo), "-c", "user.name=test", "-c", "user.email=test@example.org",
                           "-c", "commit.gpgsign=false", *command], capture_output=True, text=True, check=True)
    return proc.stdout.strip()


@pytest.fixture
def engine_repo(tmp_path):
    repo = tmp_path / "engine"
    repo.mkdir()
    git(repo, "init", "-q")
    (repo / "cli.py").write_text("print('engine')\n", encoding="utf-8")
    git(repo, "add", "cli.py")
    git(repo, "commit", "-q", "-m", "first")
    first = git(repo, "rev-parse", "HEAD")
    (repo / "cli.py").write_text("print('engine v2')\n", encoding="utf-8")
    git(repo, "commit", "-q", "-am", "second")
    return repo, first, git(repo, "rev-parse", "HEAD")


def pin(repo, commit, engine="released"):
    return argparse.Namespace(engine=engine, engine_root=str(repo), engine_commit=commit)


def test_default_pin_is_the_released_engine_commit():
    args = skill.build_parser().parse_args(["--mode", "quick", "--arm", "A", "--out", "o"])
    assert args.engine_commit == "b9a9dce"


def test_clean_checkout_at_the_pinned_commit_passes(engine_repo):
    repo, _, second = engine_repo
    assert skill.check_engine_checkout(pin(repo, second[:7])) == second


def test_checkout_at_another_commit_is_refused(engine_repo):
    repo, first, _ = engine_repo
    with pytest.raises(ValueError, match="not at the expected commit"):
        skill.check_engine_checkout(pin(repo, first[:7]))
    with pytest.raises(ValueError, match="not at the expected commit"):
        skill.check_engine_checkout(pin(repo, "0000000"))


def test_dirty_or_untracked_checkout_is_refused(engine_repo):
    repo, _, second = engine_repo
    (repo / "cli.py").write_text("print('patched')\n", encoding="utf-8")
    with pytest.raises(ValueError, match="uncommitted changes"):
        skill.check_engine_checkout(pin(repo, second))
    git(repo, "checkout", "--", "cli.py")
    (repo / "extra_module.py").write_text("x = 1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="untracked files"):
        skill.check_engine_checkout(pin(repo, second))


def test_non_git_engine_root_is_refused_and_audited_engine_is_not_pinned(tmp_path):
    with pytest.raises(ValueError, match="not a git checkout"):
        skill.check_engine_checkout(pin(tmp_path, "b9a9dce"))
    assert skill.check_engine_checkout(pin(tmp_path, "b9a9dce", engine="audited")) is None


def test_cli_refuses_a_mismatched_engine_before_writing(engine_repo, tmp_path):
    repo, first, _ = engine_repo
    out = tmp_path / "run"
    proc = subprocess.run([sys.executable, str(ROOT / "tools" / "rmaps3_skill_run.py"), "--mode", "quick",
                           "--arm", "S_A", "--out", str(out), "--up", "u", "--dn", "d", "--bg", "b",
                           "--genome-root", "g", "--engine-root", str(repo), "--engine-commit", first[:7],
                           "--no-figures", "--gate-rule", "A"], cwd=ROOT, capture_output=True, text=True, timeout=300,
                          env=dict(os.environ))
    assert proc.returncode == 2
    assert "not at the expected commit" in proc.stderr
    assert not out.exists()


def test_reference_hashes_cover_every_given_reference_list_and_gate_input(tmp_path):
    files = {}
    for name in ("g.fa", "g.fa.fai", "known.txt", "esrp.txt", "alias.tsv", "genes.gtf", "splice.txt",
                 "broad.txt", "counts.json", "up.txt", "dn.txt", "bg.txt"):
        files[name] = tmp_path / name
        files[name].write_text(name, encoding="utf-8")
    args = skill.build_parser().parse_args([
        "--mode", "quick", "--arm", "S_A", "--out", "o", "--up", str(files["up.txt"]), "--dn", str(files["dn.txt"]),
        "--bg", str(files["bg.txt"]), "--alias-table", str(files["alias.tsv"]), "--gtf", str(files["genes.gtf"]),
        "--spliceosome-list", str(files["splice.txt"]), "--broad-binders-list", str(files["broad.txt"]),
        "--gate-counts", str(files["counts.json"])])
    hashes = skill.reference_hashes(args, files["g.fa"], files["known.txt"], files["esrp.txt"])
    labels = {label for label, _, _ in hashes}
    assert labels == {"genome_fasta", "genome_fai", "known_motifs", "additional_motifs", "alias_table", "gtf",
                      "spliceosome_list", "broad_binders_list", "gate_counts", "up", "dn", "bg"}
    assert all(re.fullmatch(r"[0-9a-f]{64}", digest) for _, _, digest in hashes)


MODULE_TO_DIST = {"yaml": "pyyaml", "PIL": "pillow"}


def canonical(name):
    return re.sub(r"[-_.]+", "-", name).lower()


def lock_pins():
    """{canonical distribution: version}; every non-comment line must be an exact name==version pin."""
    pins = {}
    for line in (ROOT / "requirements-lock.txt").read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        match = re.fullmatch(r"([A-Za-z0-9][A-Za-z0-9_.\-]*)==([A-Za-z0-9_.+!\-]+)", line.strip())
        assert match, f"not an exact pin: {line!r}"
        pins[canonical(match.group(1))] = match.group(2)
    return pins


def chain_imports():
    """Top-level third-party modules imported by cli.py, rmaps_core/, tools/, webui/ and tests/ (not legacy)."""
    import ast
    files = [ROOT / "cli.py"] + [p for d in ("rmaps_core", "tools", "webui", "tests") for p in (ROOT / d).rglob("*.py")
                                 if "legacy" not in p.relative_to(ROOT).parts]
    local = {p.stem for p in files} | {"rmaps_core", "webui", "tests", "tools"}
    found = set()
    for path in files:
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8-sig"))):
            names = ([a.name for a in node.names] if isinstance(node, ast.Import) else
                     [node.module] if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module else [])
            found |= {n.split(".")[0] for n in names}
    return {m for m in found if m not in sys.stdlib_module_names and m not in local}


def test_lock_is_an_exact_pin_set_covering_every_third_party_import_of_the_chain():
    pins = lock_pins()
    imports = chain_imports()
    assert {"numpy", "scipy", "pandas", "yaml", "pyfaidx", "typer", "flask"} <= imports   # the scan sees the chain
    for module in sorted(imports):
        assert canonical(MODULE_TO_DIST.get(module, module)) in pins, module
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "## Environment" in readme and "requirements-lock.txt" in readme


def test_every_pinned_distribution_installed_here_is_at_its_pinned_version():
    import importlib.metadata as metadata
    drift = {}
    for name, version in lock_pins().items():
        try:
            installed = metadata.version(name)
        except metadata.PackageNotFoundError:
            continue
        if installed != version:
            drift[name] = (installed, version)
    assert not drift, drift
