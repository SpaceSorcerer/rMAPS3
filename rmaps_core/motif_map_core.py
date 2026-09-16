from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import sys
import subprocess
import shutil
import os
import json

from rmaps_core.input_utils import maybe_prepare_rmats_input
from rmaps_core.path_utils import build_subprocess_env, repo_root, resolve_user_path
from rmaps_core.stat_utils import normalize_stat_method
from rmaps_core.output_utils import ensure_output_directory, write_run_manifest


PYTHON = sys.executable
REPO_ROOT = repo_root()


@dataclass(frozen=True)
class EventSpec:
    """
    Minimal shared description for each event type.

    Right now this encodes:
    - which motifMap script to run
    - which miso2rMATS converter to use (for the standalone convert CLI)
    """

    name: str
    script: str
    miso_converter: str


EVENT_SPECS: Dict[str, EventSpec] = {
    "se": EventSpec(
        name="SE",
        script="legacy/motifMapSE_MP.py",
        miso_converter="bin/miso2rMATS.SE.pl",
    ),
    "a3ss": EventSpec(
        name="A3SS",
        script="legacy/motifMapA3SS_MP.py",
        miso_converter="bin/miso2rMATS.A3SS.pl",
    ),
    "a5ss": EventSpec(
        name="A5SS",
        script="legacy/motifMapA5SS_MP.py",
        miso_converter="bin/miso2rMATS.A5SS.pl",
    ),
    "ri": EventSpec(
        name="RI",
        script="legacy/motifMapRI_MP.py",
        miso_converter="bin/miso2rMATS.RI.pl",
    ),
    "mxe": EventSpec(
        name="MXE",
        script="legacy/motifMapMXE_MP.py",
        miso_converter="bin/miso2rMATS.MXE.pl",
    ),
}


def event_script(event: str) -> Path:
    """
    Return the path to the motifMap script for the given event (se, a3ss, a5ss, ri, mxe).
    """
    key = event.lower()
    if key not in EVENT_SPECS:
        raise ValueError(f"Unsupported event type: {event}")
    return REPO_ROOT / EVENT_SPECS[key].script


def miso_converter_script(event: str) -> Path:
    """
    Return the path to the Perl miso2rMATS converter for the given event.
    """
    key = event.lower()
    if key not in EVENT_SPECS:
        raise ValueError(f"Unsupported event type: {event}")
    return REPO_ROOT / EVENT_SPECS[key].miso_converter


def run_subprocess(
    cmd: list[str],
    env_overrides: dict[str, str] | None = None,
    base_cwd: Path | None = None,
) -> int:
    """
    Shared helper for launching child processes from this project.
    """
    env = build_subprocess_env(REPO_ROOT, env_overrides)
    result = subprocess.run(cmd, cwd=base_cwd or Path.cwd(), env=env)
    return result.returncode


def run_motif_map(
    event: str,
    known_motifs: Path,
    motifs: str,
    fasta_root: Path,
    genome: str,
    output: Path,
    rmats: str,
    miso: str,
    up: str,
    down: str,
    background: str,
    label: str,
    intron: int,
    exon: int,
    window: int,
    step: int,
    sig_fdr: float,
    sig_delta_psi: float,
    separate: bool,
    stat_method: str = "fisher",
    stat_permutations: int | None = None,
    stat_seed: int | None = None,
    keep_temp: bool = False,
    base_cwd: Path | None = None,
    delete_temp: bool = False,
    overwrite: bool = False,
    allow_overlap: bool = False,
    fisher_alternative: str = "greater",
    workers: int | None = None,
) -> int:
    """
    Build and run the legacy motifMap* script for a given event type.

    This keeps all event-specific wiring in one place so the Typer CLI
    can call a single entrypoint per event.
    """
    script_path = event_script(event)
    stat_method = normalize_stat_method(stat_method)
    base_cwd = base_cwd or Path.cwd()
    known_motifs = resolve_user_path(known_motifs, base_cwd)
    motifs = resolve_user_path(motifs, base_cwd)
    fasta_root = resolve_user_path(fasta_root, base_cwd)
    output = Path(resolve_user_path(output, base_cwd))
    rmats = resolve_user_path(rmats, base_cwd)
    miso = resolve_user_path(miso, base_cwd)
    up = resolve_user_path(up, base_cwd)
    down = resolve_user_path(down, base_cwd)
    background = resolve_user_path(background, base_cwd)
    # AUDIT F8: reject reuse before XLSX preparation can modify an existing run.
    if event.lower() == "se":
        ensure_output_directory(output, overwrite=overwrite)
    original_rmats = rmats
    rmats = maybe_prepare_rmats_input(rmats, output)
    cmd: list[str] = [
        PYTHON,
        str(script_path),
        "-k",
        known_motifs,
        "-m",
        motifs,
        "--fasta-root",
        fasta_root,
        "-g",
        genome,
        "-o",
        str(output),
        "-r",
        rmats,
        "-mi",
        miso,
        "-u",
        up,
        "-d",
        down,
        "-b",
        background,
        "--label",
        label,
        "--intron",
        str(intron),
        "--exon",
        str(exon),
        "--window",
        str(window),
        "--step",
        str(step),
        "--sigFDR",
        str(sig_fdr),
        "--sigDeltaPSI",
        str(sig_delta_psi),
    ]
    if separate:
        cmd.append("--separate")

    if event.lower() == "se":
        # AUDIT F15: pass the explicit Fisher tail to the SE engine.
        if fisher_alternative not in ("greater", "two-sided"):
            raise ValueError("fisher_alternative must be greater or two-sided")
        cmd.extend(["--fisher-alternative", fisher_alternative])
        # AUDIT F20: bound SE concurrency, including single-CPU systems.
        workers = workers if workers is not None else min(4, max(1, (os.cpu_count() or 1) - 1))
        if workers < 1:
            raise ValueError("workers must be at least 1")
        cmd.extend(["--workers", str(workers)])
        # AUDIT F13: cross-set overlap remains an explicit opt-in.
        if allow_overlap:
            cmd.append("--allow-overlap")
        # AUDIT F8: prepared XLSX files are authorized by the parent preflight.
        if overwrite or rmats != original_rmats:
            cmd.append("--overwrite")
        # AUDIT F7: retain positional tables unless deletion is requested.
        if delete_temp:
            cmd.append("--delete-temp")

    env_overrides = {"RMAPS_STAT_METHOD": stat_method}
    if stat_permutations is not None:
        env_overrides["RMAPS_STAT_PERMUTATIONS"] = str(stat_permutations)
    if stat_seed is not None:
        env_overrides["RMAPS_STAT_SEED"] = str(stat_seed)

    code = run_subprocess(cmd, env_overrides, base_cwd=base_cwd)

    # AUDIT F7: SE owns explicit cleanup; other event types retain prior behavior.
    if code == 0 and event.lower() != "se" and not keep_temp:
        shutil.rmtree(output / "temp", ignore_errors=True)

    # AUDIT F8: include the wrapper-created XLSX conversion in this run's manifest.
    manifest = output / "run_manifest.json"
    if code == 0 and event.lower() == "se" and rmats != original_rmats and manifest.exists():
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        files = [entry["path"] for entry in payload["outputs"]]
        if Path(rmats).exists():
            files.append(rmats)
        parameters = payload["parameters"]
        parameters["original_rmats"] = original_rmats
        parameters["overwrite"] = overwrite
        write_run_manifest(output, files, parameters)

    return code

