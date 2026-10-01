from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from rmaps_core.path_utils import repo_root


# AUDIT R8: archive only paths registered by the previous run; preserve unowned files.
def ensure_output_directory(output: Path | str, overwrite: bool = False,
                            *, target_paths: Iterable[Path | str] = ()) -> Path:
    output = Path(output).resolve()
    if output.exists() and not output.is_dir():
        raise ValueError(f"Output is not a directory: {output}")
    nonempty = output.exists() and any(output.iterdir())
    if nonempty and not overwrite:
        raise ValueError(f"Output directory is not empty: {output}; use --overwrite to reuse it")
    manifest = output / "run_manifest.json"
    # AUDIT S1: unknown ownership cannot be made safe by --overwrite.
    if nonempty and not manifest.is_file():
        raise ValueError(f"Refusing non-empty output directory without run_manifest.json: {output}")
    registered = []
    owned = set()
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
            if path == manifest:
                raise ValueError("The manifest must not register itself as an output")
            owned.add(path)
            if path.exists():
                if not path.is_file():
                    raise ValueError(f"Registered output is not a file: {path}")
                registered.append((path, relative))
    # AUDIT S1: validate every planned write and its parents before mkdir/archive.
    for item in target_paths:
        path = Path(item)
        path = path.resolve() if path.is_absolute() else (output / path).resolve()
        relative = path.relative_to(output)
        if not relative.parts or relative.parts[0].startswith("_previous_"):
            raise ValueError(f"Invalid output target: {item}")
        if path.exists() and path not in owned and path != manifest:
            raise ValueError(f"Refusing unregistered output target: {path}")
        if path.exists() and not path.is_file():
            raise ValueError(f"Output target is not a file: {path}")
        for parent in path.parents:
            if parent == output:
                break
            if parent.exists() and not parent.is_dir():
                raise ValueError(f"Output target parent is not a directory: {parent}")
    output.mkdir(parents=True, exist_ok=True)
    if manifest.exists():
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
        archive = output / f"_previous_{stamp}"
        archive.mkdir()
        for path, relative in registered:
            destination = archive / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            path.rename(destination)
        manifest.rename(archive / manifest.name)
    return output


def se_output_targets(known_motifs: Path | str, motifs: Path | str = "NA",
                      rmats: Path | str = "NA", miso: Path | str = "NA") -> set[str]:
    # AUDIT S1: enumerate possible SE writes before the first output mutation.
    from rmaps_core.se_windows import REGION_NAMES

    targets = {"run_manifest.json", "log.motifMap.txt", "maps/jj.png",
               "pVal.up.vs.bg.RNAmap.txt", "pVal.dn.vs.bg.RNAmap.txt"}
    fasta_regions = ("UpstreamExon", "UpstreamExonIntron", "UpstreamIntron",
                     "TargetExon", "DownstreamIntron", "DownstreamExonIntron",
                     "DownstreamExon")
    for label in ("up", "dn", "bg"):
        targets.add(f"exon/{label}.coord.txt")
        targets.update(f"fasta/{label}.{region}.fasta" for region in fasta_regions)
    targets.update(f"temp/sequence_region_{index}.npy" for index in range(8))
    targets.add("temp/sequence_metadata.npz")
    if str(miso) != "NA":
        targets.add("temp/converted.rMATS.se.txt")
    if str(rmats) != "NA" and Path(rmats).suffix.lower() == ".xlsx":
        targets.add(f"temp/{Path(rmats).stem}.from_xlsx.tsv")
    for source in (known_motifs, motifs):
        if str(source) == "NA":
            continue
        with Path(source).open() as handle:
            next(handle)
            for line in handle:
                if not line.strip():
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 2:
                    raise ValueError(f"Invalid motif row in {source}: expected name and expression")
                name, pattern = fields[:2]
                motif_id = f"{name}.{pattern}"
                targets.add(f"temp/{motif_id}.txt")
                targets.update(f"temp/{motif_id}.countDist.{label}.txt" for label in ("up", "dn", "bg"))
                targets.update(f"temp/{motif_id}.pVal.{label}.vs.bg.txt" for label in ("up", "dn"))
                targets.update(f"positional/{motif_id}.{region}.hits.npz" for region in REGION_NAMES)
                safe_name = re.sub(r'[^A-Za-z0-9._-]+', '_', f"{name}-{pattern}")
                targets.update((f"maps/SE.{safe_name}.pdf", f"maps/SE.{safe_name}.png",
                                f"maps/{safe_name}.png"))
    return targets


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
