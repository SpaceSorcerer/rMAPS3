#!/usr/bin/env python3
"""Portable rule-A SE gate builder; VAST-concordant rule B requires pre-split files.

Adapted from the lab build_event_sets.py (2026-09-16), SHA256
b58e1a79b2b2b917d344350d95fad9d929b306a1f305d6c495eb3728997aa191.
Frozen-logic attribution: reli_v121/code/build_reli_foregrounds_v2.py,
source MD5 32be07fd0e6e2edff5f077384f45b352. region_windows is copied
verbatim; rule-A selection, ID-only expression and exact collision semantics
are translated from pandas into stdlib operations. The lab-specific VAST,
ledger regression and sensitivity-grid driver are not copied.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import re
import subprocess
import sys
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rmaps3_lab_run import (COORD_HEADER, fingerprint, fresh_directory, git_state,
                            load_config, resolve_path, versions, write_json)

COORD_COLS = ["chr", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"]
MISSING = {"", "NA", "NAN", "N/A", "NULL", ".", "-"}


def number(value, label):
    value = str(value).strip().strip('"')
    if value.upper() in MISSING:
        return math.nan
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Nonfinite {label}: {value}")
    return result


def rows(path, delimiter, repeated_id=False):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader, None)
        if not header:
            raise ValueError(f"Empty table: {path}")
        duplicates = {name for name in header if header.count(name) > 1}
        if duplicates and not (repeated_id and duplicates == {"ID"} and header.count("ID") == 2):
            raise ValueError(f"Duplicate columns: {path}: {duplicates}")
        ids = [i for i, value in enumerate(header) if value == "ID"]
        for line, values in enumerate(reader, 2):
            if len(values) != len(header):
                raise ValueError(f"Wrong column count: {path}:{line}")
            if len(ids) == 2 and values[ids[0]] != values[ids[1]]:
                raise ValueError(f"Repeated rMATS ID columns disagree: {path}:{line}")
            yield dict(zip(header, values))


def gene_id(value):
    return re.sub(r"\.\d+$", "", str(value).strip().strip('"'))


def expression_values(path, delimiter=","):
    result = {}
    for row in rows(path, delimiter):
        gid = gene_id(row["gene_id"])
        value = number(row["baseMean"], "DESeq2 baseMean")
        if not re.fullmatch(r"ENSG\d+", gid) or value < 0:
            raise ValueError(f"Invalid human DESeq2 gene ID/baseMean: {gid}")
        if gid in result and not (value == result[gid] or math.isnan(value) and math.isnan(result[gid])):
            raise ValueError(f"Conflicting duplicate DESeq2 gene: {gid}")
        result[gid] = value
    return result


def region_windows(row):
    """Recovered 2026-07-08 window arithmetic; see module docstring for the proof."""
    left = (row.chr, int(row.upstreamEE), int(row.exonStart_0base))
    exon = (row.chr, int(row.exonStart_0base), int(row.exonEnd))
    right = (row.chr, int(row.exonEnd), int(row.downstreamES))
    merged = (row.chr, int(row.upstreamEE), int(row.downstreamES))
    up, dn = (left, right) if row.strand == "+" else (right, left)
    return {"UPintr": up, "AltEX": exon, "DNintr": dn, "merged": merged}


def build_sets(source, gates, expression):
    defaults = {"fdr": 0.05, "abs_dpsi": 0.10, "min_sample_coverage": 10,
                "base_mean_floor": None, "expr_unknown": "retain", "background_fdr": 0.5,
                "treatment_is": "b1", "rule": "A"}
    unknown = set(gates) - set(defaults)
    if unknown:
        raise ValueError(f"Unknown gate keys: {sorted(unknown)}")
    gates = {**defaults, **gates}
    if gates["rule"] != "A":
        raise ValueError("Portable builder supports rule A; supply pre-split coordinates for rule B")
    if gates["expr_unknown"] not in {"retain", "exclude", "error"} or gates["treatment_is"] not in {"b1", "b2"}:
        raise ValueError("Invalid expr_unknown or treatment_is policy")
    for key in ("fdr", "abs_dpsi", "background_fdr", "min_sample_coverage", "base_mean_floor"):
        if gates[key] is not None:
            value = float(gates[key])
            if not math.isfinite(value) or value < 0 or key in {"fdr", "abs_dpsi", "background_fdr"} and value > 1:
                raise ValueError(f"Invalid gate {key}: {value}")
            gates[key] = value
    required = set(COORD_COLS + ["ID", "GeneID", "FDR", "IncLevelDifference", "IncLevel1", "IncLevel2",
                               "IJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_1", "SJC_SAMPLE_2"])
    records, audit, seen, sample_lengths = [], [], {}, {}
    for line, row in enumerate(rows(source, "\t", repeated_id=True), 2):
        if required - row.keys():
            raise ValueError(f"Missing rMATS columns: {sorted(required - row.keys())}")
        if not row["ID"].strip():
            raise ValueError("Empty rMATS event ID")
        if row["ID"] in seen and row != seen[row["ID"]]:
            raise ValueError(f"Conflicting duplicate rMATS ID: {row['ID']}")
        seen[row["ID"]] = row.copy()
        for col in COORD_COLS[2:]:
            value = number(row[col], col)
            if not math.isfinite(value) or value < 0 or not value.is_integer():
                raise ValueError(f"Invalid rMATS coordinate: {col}")
            row[col] = int(value)
        if row["strand"] not in {"+", "-"} or not re.fullmatch(r"[^\s:]+", row["chr"]):
            raise ValueError("Invalid rMATS chromosome/strand")
        if not (row["upstreamES"] < row["upstreamEE"] <= row["exonStart_0base"] < row["exonEnd"] <= row["downstreamES"] < row["downstreamEE"]):
            raise ValueError(f"Invalid SE geometry: {row['ID']}")
        fdr, dpsi = number(row["FDR"], "FDR"), number(row["IncLevelDifference"], "dPSI")
        if not (math.isnan(fdr) or 0 <= fdr <= 1) or not (math.isnan(dpsi) or -1 <= dpsi <= 1):
            raise ValueError("FDR/dPSI outside allowed range")
        covered, means = True, []
        for group in (1, 2):
            vectors = [[number(x, col) for x in row[col].split(",")]
                       for col in (f"IJC_SAMPLE_{group}", f"SJC_SAMPLE_{group}", f"IncLevel{group}")]
            length = len(vectors[0])
            if any(len(vec) != length for vec in vectors) or sample_lengths.get(group, length) != length:
                raise ValueError(f"Sample vector length mismatch: {row['ID']}")
            sample_lengths[group] = length
            ijc, sjc, inc = vectors
            if any(math.isfinite(x) and (x < 0 or not x.is_integer()) for vec in (ijc, sjc) for x in vec):
                raise ValueError(f"Invalid junction count: {row['ID']}")
            if any(math.isfinite(x) and not 0 <= x <= 1 for x in inc):
                raise ValueError(f"Invalid inclusion level: {row['ID']}")
            covered = covered and all(i + s >= gates["min_sample_coverage"] for i, s in zip(ijc, sjc))
            finite = [x for x in inc if math.isfinite(x)]
            means.append(sum(finite) / len(finite) if finite else math.nan)
        delta = (means[0] - means[1]) * (1 if gates["treatment_is"] == "b1" else -1)
        evaluable = math.isfinite(delta) and math.isfinite(dpsi)
        agrees = evaluable and ((delta > 0) - (delta < 0)) == ((dpsi > 0) - (dpsi < 0))
        gid = gene_id(row["GeneID"])
        if not re.fullmatch(r"ENSG\d+", gid):
            raise ValueError(f"Invalid human rMATS GeneID: {gid}")
        bm = expression.get(gid, math.nan)
        unknown_expr = math.isnan(bm)
        floor = gates["base_mean_floor"]
        if floor is not None and unknown_expr and gates["expr_unknown"] == "error":
            raise ValueError(f"expr_unknown: {gid}")
        eligible = floor is None or (gates["expr_unknown"] == "retain" if unknown_expr else bm > floor)
        fg = fdr < gates["fdr"] and abs(dpsi) >= gates["abs_dpsi"] and covered and eligible
        bg = fdr >= gates["background_fdr"] and covered and eligible
        if fg and not agrees:
            raise ValueError(f"Selected foreground orientation disagreement: {row['ID']}; no automatic sign flip")
        records.append((row, fg, bg, dpsi, unknown_expr))
        audit.append({"source_row": line, "ID": row["ID"], "orientation_evaluable": evaluable,
                      "orientation_agrees": agrees, "informative": evaluable and dpsi != 0,
                      "coverage_pass": covered, "expr_status": "expr_unknown" if unknown_expr else "pass" if eligible else "fail",
                      "baseMean": None if unknown_expr else bm, "foreground": fg, "background_before_collision": bg})
    if not records:
        raise ValueError("Empty rMATS table")
    evaluable = [r for r in audit if r["orientation_evaluable"]]
    informative = [r for r in audit if r["informative"]]
    if not evaluable or not informative or sum(r["orientation_agrees"] for r in evaluable) / len(evaluable) <= .95 or sum(r["orientation_agrees"] for r in informative) / len(informative) <= .95:
        raise ValueError("Global orientation agreement must exceed 95% for evaluable and nonzero-dPSI events")
    foreground = [r for r, fg, _, _, _ in records if fg]
    ids = {r["ID"] for r in foreground}
    exons = {(r["chr"], r["exonStart_0base"], r["exonEnd"]) for r in foreground}
    windows = {w for r in foreground for w in region_windows(SimpleNamespace(**r)).values()}
    sets, collisions, unknown_fg, unknown_bg = {"up": [], "dn": [], "bg": []}, 0, 0, 0
    for row, fg, bg, dpsi, unknown_expr in records:
        collision = row["ID"] in ids or (row["chr"], row["exonStart_0base"], row["exonEnd"]) in exons or any(w in windows for w in region_windows(SimpleNamespace(**row)).values())
        collisions += int(bg and collision)
        if fg:
            if dpsi != 0:
                sets["up" if dpsi > 0 else "dn"].append(row)
            unknown_fg += int(unknown_expr)
        if bg and not collision:
            sets["bg"].append(row)
            unknown_bg += int(unknown_expr)
    counts = {"n_input": len(records), "bg_collisions_removed": collisions,
              "n_expr_unknown_in_fg": unknown_fg, "n_expr_unknown_in_bg": unknown_bg,
              **{"n_" + key: len(value) for key, value in sets.items()}}
    return sets, counts, audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    try:
        config_path = args.config.resolve()
        config = load_config(config_path)
        base = config_path.parent
        inputs, gates = config["inputs"], config.get("gates", {})
        source = resolve_path(inputs["rmats_se"], base)
        if gates.get("base_mean_floor") is not None and not inputs.get("deseq2"):
            raise ValueError("base_mean_floor requires inputs.deseq2")
        de = resolve_path(inputs["deseq2"], base) if inputs.get("deseq2") else None
        delimiter = inputs.get("deseq2_delimiter", "," if de and de.suffix.lower() == ".csv" else "\t")
        expression = expression_values(de, delimiter) if de else {}
        sets, counts, audit = build_sets(source, gates, expression)
        for group, rows_ in sets.items():
            if not rows_:
                raise ValueError(f"Gate produced empty {group} set; cannot run SE enrichment")
        out = args.out.resolve()
        fresh_directory(out)
        for group, rows_ in sets.items():
            with (out / f"{group}.coord.txt").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(COORD_HEADER)
                writer.writerows([[row[key] for key in COORD_COLS] for row in rows_])
        write_json(out / "counts.json", counts)
        write_json(out / "event_audit.json", audit)
        version = versions()
        command = subprocess.list2cmdline(sys.orig_argv)
        (out / "versions.txt").write_text("\n".join(f"{k}={v}" for k, v in version.items()) + "\n", encoding="utf-8")
        (out / "command.log").write_text(command + "\n", encoding="utf-8")
        write_json(out / "run_manifest.json", {"schema_version": 1, "status": "complete", "config": fingerprint(config_path),
                   "git": git_state(), "versions": version, "command": command, "gates": gates,
                   "input_orientation": "Treatment-control sign verified; no sign flip", "builder": fingerprint(Path(__file__)),
                   "inputs": [fingerprint(p) for p in [source, *([de] if de else [])]], "counts": counts,
                   "outputs": [fingerprint(p) for p in out.iterdir() if p.is_file()]})
        print("BUILT " + str(out))
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(2, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
