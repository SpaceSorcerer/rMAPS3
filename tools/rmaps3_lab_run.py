#!/usr/bin/env python3
"""Run configured SE pre-splitting, chromosome normalization, mapping and calibration."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
COORD_HEADER = ["chr", "strand", "exonStart", "exonEnd", "firstExonStart",
                "firstExonEnd", "secondExonStart", "secondExonEnd"]


def load_config(path):
    text = path.read_text(encoding="utf-8-sig")
    if path.suffix.lower() in {".yaml", ".yml"}:
        try:
            import yaml
        except ImportError as exc:
            raise ValueError("YAML config requires already installed PyYAML; use a .json config instead") from exc
        config = yaml.safe_load(text)
    elif path.suffix.lower() == ".json":
        config = json.loads(text)
    else:
        raise ValueError("Config extension must be .json, .yaml or .yml")
    if not isinstance(config, dict):
        raise ValueError("Config must be an object")
    return config


def resolve_path(value, base):
    path = Path(value).expanduser()
    return (base / path).resolve() if not path.is_absolute() else path.resolve()


def fingerprint(path):
    path = Path(path).resolve()
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            md5.update(block)
            sha.update(block)
    return {"path": str(path), "size_bytes": path.stat().st_size,
            "md5": md5.hexdigest(), "sha256": sha.hexdigest()}


def versions():
    result = {"python": sys.version, "executable": sys.executable, "platform": platform.platform()}
    for package in ("numpy", "scipy", "openpyxl", "pyfaidx", "PyYAML"):
        try:
            result[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            result[package] = "not installed"
    return result


def git_state():
    def query(*args):
        proc = subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True)
        return proc.stdout.strip() if proc.returncode == 0 else "unavailable: " + proc.stderr.strip()
    return {"revision": query("rev-parse", "HEAD"), "status": query("status", "--porcelain")}


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def fresh_directory(path):
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError(f"Refusing non-empty output directory: {path}")
    path.mkdir(parents=True, exist_ok=True)


def positive_int(value, name):
    if isinstance(value, bool) or int(value) != value or int(value) < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def normalize_coordinates(sources, fai, destination):
    names = set()
    with fai.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                name = line.split("\t", 1)[0]
                if name in names:
                    raise ValueError(f"Duplicate FASTA index name: {name}")
                names.add(name)
    if not names:
        raise ValueError(f"Empty FASTA index: {fai}")
    counts, paths = {}, {}
    for group, source in sources.items():
        rows, renamed = [], {}
        with source.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            if next(reader, None) != COORD_HEADER:
                raise ValueError(f"Expected exact eight-column coordinate header: {source}")
            for number, row in enumerate(reader, 2):
                if len(row) != 8:
                    raise ValueError(f"Expected eight columns: {source}:{number}")
                original = row[0]
                if original not in names:
                    alias = original[3:] if original.startswith("chr") else "chr" + original
                    if alias not in names:
                        raise ValueError(f"Chromosome {original!r} absent from FASTA index: {source}:{number}")
                    row[0] = alias
                    renamed[original] = {"target": alias, "count": renamed.get(original, {}).get("count", 0) + 1}
                rows.append(row)
        if not rows:
            raise ValueError(f"Empty event set: {source}")
        target = destination / f"{group}.coord.txt"
        with target.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(COORD_HEADER)
            writer.writerows(rows)
        paths[group] = target
        counts[group] = {"n_events": len(rows), "renamed": renamed, "output": fingerprint(target)}
    return paths, counts


def run(config_path, perms=None):
    config_path = config_path.resolve()
    config = load_config(config_path)
    base = config_path.parent
    arm = str(config.get("arm", "run"))
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", arm) or ".." in arm:
        raise ValueError("arm must be a filename-safe name without '..' or path separators")
    out = resolve_path(config["output_root"], base)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError(f"Refusing non-empty output directory: {out}")
    inputs, genome, motifs = (config[name] for name in ("inputs", "genome", "motifs"))
    engine, summary = config.get("engine", {}), config.get("summary", {})
    if engine.get("stat_method", "fisher") != "fisher" or engine.get("fisher_alternative", "greater") != "greater":
        raise ValueError("Calibrated summarizer supports only stat_method=fisher with fisher_alternative=greater")
    workers = positive_int(engine.get("workers", 4), "engine.workers")
    summary_workers = positive_int(summary.get("workers", 1), "summary.workers")
    threads = positive_int(config.get("blas_threads", 4), "blas_threads")
    permutations = positive_int(perms if perms is not None else summary.get("perms", 2000), "summary.perms")
    sizes = {key: positive_int(engine.get(key, default), "engine." + key)
             for key, default in (("intron", 250), ("exon", 50), ("window", 50), ("step", 1))}
    if "exon_window" in engine:
        sizes["exon-window"] = positive_int(engine["exon_window"], "engine.exon_window")
    if not isinstance(genome["build"], str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", genome["build"]) or ".." in genome["build"]:
        raise ValueError("genome.build must be a single directory name")
    genome_root = resolve_path(genome["root"], base)
    fasta = genome_root / genome["build"] / (genome["build"] + ".fa")
    fai = Path(str(fasta) + ".fai")
    known = resolve_path(motifs["known"], base)
    additional = resolve_path(motifs["additional"], base) if motifs.get("additional") else None
    source_paths = {key: resolve_path(value, base) for key, value in inputs.items() if key in {"up", "dn", "bg", "rmats_se", "deseq2"}}
    presplit = all(key in source_paths for key in ("up", "dn", "bg"))
    if presplit == ("rmats_se" in source_paths) or (not presplit and any(key in source_paths for key in ("up", "dn", "bg"))):
        raise ValueError("Provide either all inputs.up/dn/bg or inputs.rmats_se")
    if not presplit and config.get("gates", {}).get("base_mean_floor") is not None and "deseq2" not in source_paths:
        raise ValueError("base_mean_floor requires inputs.deseq2; missing table is not an expr_unknown policy")
    initial = [fingerprint(p) for p in [config_path, fasta, fai, known, *source_paths.values(), *([additional] if additional else [])]]
    fresh_directory(out)
    env = os.environ.copy()
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[key] = str(threads)
    if engine.get("force_motif_fallback", True):
        env["RMAPS_FORCE_MOTIF_FALLBACK"] = "1"
    else:
        env.pop("RMAPS_FORCE_MOTIF_FALLBACK", None)
    environment = {key: value for key, value in env.items() if key.endswith("NUM_THREADS") or key in {"VECLIB_MAXIMUM_THREADS", "RMAPS_FORCE_MOTIF_FALLBACK"}}
    manifest = {"schema_version": 1, "status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
                "config": fingerprint(config_path), "effective_config": config, "perms_override": perms,
                "git": git_state(), "versions": versions(), "inputs": initial, "environment": environment,
                "scripts": [fingerprint(ROOT / path) for path in ("tools/rmaps3_lab_run.py", "tools/build_event_sets.py", "tools/summarize_rmaps_regions.py", "cli.py")],
                "steps": []}
    version_text = "\n".join(f"{key}={value}" for key, value in manifest["versions"].items()) + "\n"
    (out / "versions.txt").write_text(version_text, encoding="utf-8")
    (out / "command.log").write_text(subprocess.list2cmdline(sys.orig_argv) + "\n", encoding="utf-8")
    def save():
        write_json(out / "run_manifest.json", manifest)
    def step(name, command, paths):
        record = {"name": name, "status": "running", "command": command, "config_sha256": manifest["config"]["sha256"],
                  "git": manifest["git"], "versions": manifest["versions"], "inputs": [fingerprint(p) for p in paths]}
        manifest["steps"].append(record)
        save()
        with (out / (name + ".stdout.log")).open("w", encoding="utf-8") as stdout, (out / (name + ".stderr.log")).open("w", encoding="utf-8") as stderr:
            proc = subprocess.run(command, cwd=ROOT, env=env, stdout=stdout, stderr=stderr)
        record.update(returncode=proc.returncode, status="complete" if proc.returncode == 0 else "failed")
        save()
        if proc.returncode:
            raise RuntimeError(f"{name} failed (exit {proc.returncode}); see {out / (name + '.stderr.log')}")
    save()
    try:
        if presplit:
            sources = {key: source_paths[key] for key in ("up", "dn", "bg")}
        else:
            build_out = out / "event_sets"
            step("build_event_sets", [sys.executable, str(ROOT / "tools/build_event_sets.py"), "--config", str(config_path), "--out", str(build_out)], [config_path, *source_paths.values()])
            sources = {key: build_out / f"{key}.coord.txt" for key in ("up", "dn", "bg")}
        mapped = out / "inputs"
        fresh_directory(mapped)
        paths, counts = normalize_coordinates(sources, fai, mapped)
        mapping = {"name": "scaffold_mapping", "status": "complete", "config_sha256": manifest["config"]["sha256"],
                   "git": manifest["git"], "versions": manifest["versions"],
                   "rule": "Exact name preferred; otherwise add/remove literal chr prefix against FASTA index; no rows dropped",
                   "inputs": [fingerprint(p) for p in [*sources.values(), fai]], "sets": counts}
        write_json(mapped / "run_manifest.json", mapping)
        (mapped / "versions.txt").write_text(version_text, encoding="utf-8")
        (mapped / "command.log").write_text(subprocess.list2cmdline(sys.orig_argv) + "\n", encoding="utf-8")
        manifest["steps"].append(mapping)
        save()
        engine_out = out / "engine"
        command = [sys.executable, str(ROOT / "cli.py"), "motif-map", "se", "--known-motifs", str(known),
                   "--motifs", str(additional) if additional else "NA", "--fasta-root", str(genome_root), "--genome", genome["build"],
                   "--rMATS", "NA", "--miso", "NA", "--up", str(paths["up"]), "--down", str(paths["dn"]),
                   "--background", str(paths["bg"]), "--output", str(engine_out), "--stat-method", "fisher",
                   "--fisher-alternative", "greater", "--workers", str(workers), "--keep-temp"]
        for key, value in sizes.items():
            command.extend(["--" + key, str(value)])
        step("motif_map", command, [*paths.values(), fasta, fai, known, *([additional] if additional else [])])
        if not (engine_out / "run_manifest.json").is_file() or not list((engine_out / "positional").glob("*.hits.npz")):
            raise RuntimeError("Engine exit was zero but manifest or retained positional data is missing")
        (engine_out / "versions.txt").write_text(version_text, encoding="utf-8")
        (engine_out / "command.log").write_text(subprocess.list2cmdline(command) + "\n", encoding="utf-8")
        command = [sys.executable, str(ROOT / "tools/summarize_rmaps_regions.py"), "--run", str(engine_out),
                   "--out", str(out / "summary"), "--arm", arm, "--perms", str(permutations),
                   "--seed", str(summary.get("seed", 149)), "--workers", str(summary_workers),
                   "--chunk-size", str(positive_int(summary.get("chunk_size", 500), "summary.chunk_size"))]
        summary_inputs = [engine_out / "run_manifest.json", *sorted(engine_out.glob("pVal.*.RNAmap.txt")), *sorted((engine_out / "positional").glob("*.hits.npz")),
                          *sorted((engine_out / "temp").glob("*.txt"))]
        step("summarize", command, summary_inputs)
        for name in ("condensed_per_rbp.tsv", "per_motif_regions.tsv"):
            target = out / "summary" / name
            if not target.is_file() or target.stat().st_size == 0:
                raise RuntimeError(f"Summary artifact missing or empty: {target}")
        for original in initial:
            if fingerprint(original["path"])["sha256"] != original["sha256"]:
                raise RuntimeError(f"Input changed during run: {original['path']}")
        manifest["status"] = "complete"
        manifest["outputs"] = [fingerprint(p) for p in sorted(out.rglob("*")) if p.is_file() and p != out / "run_manifest.json"]
    except Exception as exc:
        manifest.update(status="failed", error=str(exc))
        raise
    finally:
        manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
        save()
    print(f"COMPLETE {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--perms", type=int, help="Override summary permutation count (use 20 for smoke tests)")
    args = parser.parse_args()
    try:
        run(args.config, args.perms)
    except (ValueError, KeyError, OSError, RuntimeError) as exc:
        parser.exit(2, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
