from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from rmaps_core.path_utils import repo_root


# AUDIT R8: archive only paths registered by the previous run; preserve unowned files.
def ensure_output_directory(output: Path | str, overwrite: bool = False) -> Path:
    output = Path(output).resolve()
    if output.exists() and not output.is_dir():
        raise ValueError(f"Output is not a directory: {output}")
    if output.exists() and any(output.iterdir()) and not overwrite:
        raise ValueError(f"Output directory is not empty: {output}; use --overwrite to reuse it")
    output.mkdir(parents=True, exist_ok=True)
    manifest = output / "run_manifest.json"
    if manifest.exists():
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        registered = []
        for entry in payload["outputs"]:
            # A previously deleted path may now belong to the user.
            if entry.get("status", "retained") == "deleted":
                continue
            path = (output / entry["path"]).resolve()
            relative = path.relative_to(output)
            if not relative.parts or relative.parts[0].startswith("_previous_"):
                raise ValueError(f"Invalid previous-run output path: {entry['path']}")
            if path.exists():
                if not path.is_file():
                    raise ValueError(f"Registered output is not a file: {path}")
                registered.append((path, relative))
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
        archive = output / f"_previous_{stamp}"
        archive.mkdir()
        for path, relative in registered:
            destination = archive / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            path.rename(destination)
        manifest.rename(archive / manifest.name)
    return output


# AUDIT F8: provenance includes only files explicitly written by the current run.
def write_run_manifest(
    output: Path | str, output_files: Iterable[Path | str], parameters: dict,
    *, delete_files: Iterable[Path | str] = (), preserved_entries: Iterable[dict] = (),
) -> Path:
    output = Path(output).resolve()
    manifest_path = output / "run_manifest.json"
    # AUDIT R3: preserve hashed deleted inputs when the wrapper extends the inventory.
    outputs = [dict(entry) for entry in preserved_entries]
    seen = {entry["path"] for entry in outputs}
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
        outputs.append({"path": relative, "sha256": digest.hexdigest(), "size_bytes": path.stat().st_size, "status": "retained"})
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
        "schema_version": 2,  # AUDIT R9: sparse positions and retained/deleted provenance.
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
    # AUDIT R2/R3: publish all hashes before deletion, then delete only registered files.
    manifest_path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    entries = {entry["path"]: entry for entry in outputs}
    pending = {}
    for item in delete_files:
        path = Path(item)
        path = path.resolve() if path.is_absolute() else (output / path).resolve()
        relative = path.relative_to(output).as_posix()
        if relative not in entries or entries[relative]["status"] != "retained":
            raise ValueError(f"Refusing to delete an unregistered output: {path}")
        pending[relative] = path
    if pending:
        try:
            for relative, path in pending.items():
                path.unlink()
                entries[relative]["status"] = "deleted"
        finally:
            manifest_path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    return manifest_path
