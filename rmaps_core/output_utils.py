from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import Iterable

from rmaps_core.path_utils import repo_root


# AUDIT F8: check before creating any preparatory files in the output directory.
def ensure_output_directory(output: Path | str, overwrite: bool = False) -> Path:
    output = Path(output).resolve()
    if output.exists() and not output.is_dir():
        raise ValueError(f"Output is not a directory: {output}")
    if output.exists() and any(output.iterdir()) and not overwrite:
        raise ValueError(f"Output directory is not empty: {output}; use --overwrite to reuse it")
    output.mkdir(parents=True, exist_ok=True)
    return output


# AUDIT F8: provenance includes only files explicitly written by the current run.
def write_run_manifest(
    output: Path | str, output_files: Iterable[Path | str], parameters: dict
) -> Path:
    output = Path(output).resolve()
    manifest_path = output / "run_manifest.json"
    outputs = []
    seen = set()
    for item in output_files:
        path = Path(item)
        path = path.resolve() if path.is_absolute() else (output / path).resolve()
        relative = path.relative_to(output).as_posix()
        if path == manifest_path or relative in seen:
            continue
        seen.add(relative)
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        outputs.append({"path": relative, "sha256": digest.hexdigest(), "size_bytes": path.stat().st_size})
    packages = {}
    for name in ("numpy", "scipy", "pyx", "pypdfium2", "Pillow", "pyfaidx", "typer"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    try:
        git_revision = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=repo_root(), capture_output=True, text=True, check=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        git_revision = None
    payload = {
        "schema_version": 1,
        "parameters": parameters,
        "statistical_method": parameters.get("stat_method", "fisher"),
        "permutations": parameters.get("stat_permutations"),
        "seed": parameters.get("stat_seed"),
        "python": {"version": platform.python_version(), "executable": sys.executable},
        "packages": packages,
        "git_revision": git_revision,
        "outputs": sorted(outputs, key=lambda entry: entry["path"]),
        "manifest_self_hash": "excluded: a manifest cannot include its own SHA256",
    }
    manifest_path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    return manifest_path
