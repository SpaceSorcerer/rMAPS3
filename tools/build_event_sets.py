#!/usr/bin/env python3
# Homo sapiens / GRCh38 (hg38) / GENCODE v49; count-only SE gate sensitivity.
# Functions/constants below copied verbatim from:
# F:\Publication Work\00_PAPER_FIGURES\Fig3_rMAPS_and_RBP-RELI\02_RESUME_2026-09-08\reli_v121\code\build_reli_foregrounds_v2.py
# source MD5: 32be07fd0e6e2edff5f077384f45b352
from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import platform
import random
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

import numpy as np
import pandas as pd


VERSION = "2.0.0-draft"
REGIONS = ("UPintr", "AltEX", "DNintr", "merged")
MISSING = {"", "NA", "NAN", "N/A", "NULL", ".", "-"}
EXPECTED = {
    "MIAT_KD": (["WT_MIAT_KD_1", "WT_MIAT_KD_2", "WT_MIAT_KD_3"],
                ["WT_NT_1", "WT_NT_2", "WT_NT_3"]),
    "MIAT_OE": (["MIAT_dox1", "MIAT_dox2", "MIAT_dox3"],
                ["NOG_dox1", "NOG_dox2", "NOG_dox3"]),
    "QKI_KO": (["UDQKI1", "UDQKI3"], ["UDWT1", "UDWT3"]),
}
# Explicit aliases based on the inspected OE header and ledger sample names;
# no fuzzy matching or selecting the first six columns of its nine-library table.
OE_ALIASES = {
    "MIAT_dox1": "NKX-MIAT3_dox1_S7_R1_001",
    "MIAT_dox2": "NKX-MIAT3_dox2_S8_R1_001",
    "MIAT_dox3": "NKX-MIAT3_dox3_S9_R1_001",
    "NOG_dox1": "NKX-no-gRNA_dox1_S1_R1_001",
    "NOG_dox2": "NKX-no-gRNA_dox2_S2_R1_001",
    "NOG_dox3": "NKX-no-gRNA_dox3_S3_R1_001",
}
COORD = re.compile(r"^(chr[^:]+):(\d+)-(\d+)$")
GENE_ID = re.compile(r"^ENSG\d+$")


def require(condition, message):
    """Runtime assertion, deliberately not disabled by python -O."""
    if not bool(condition):
        raise AssertionError(message)


def clean(value):
    return str(value).strip().strip('"')


def gene_id(value):
    return re.sub(r"\.\d+$", "", clean(value))


def number(value, context):
    text = clean(value)
    if text.upper() in MISSING:
        return math.nan
    try:
        result = float(text)
    except ValueError as exc:
        raise AssertionError(f"Invalid numeric value at {context}: {text!r}") from exc
    require(math.isfinite(result), f"Nonfinite value at {context}: {text!r}")
    return result


def read_table(path, delimiter, rmats=False):
    # rMATS really has two columns called ID. Validate their equality instead of
    # accepting pandas' automatic ID.1 renaming without checking the source.
    with path.open(encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle, delimiter=delimiter))
    duplicates = {x for x in header if header.count(x) > 1}
    require(not duplicates or (rmats and duplicates == {"ID"} and header.count("ID") == 2),
            f"Duplicate column names in {path}: {sorted(duplicates)}")
    table = pd.read_csv(path, sep=delimiter, dtype=str, keep_default_na=False,
                        encoding="utf-8-sig")
    if rmats and "ID" in duplicates:
        require((table["ID"] == table["ID.1"]).all(), "rMATS repeated ID columns disagree")
        table = table.drop(columns="ID.1")
    return table


def columns(table, names, context):
    missing = sorted(set(names) - set(table.columns))
    require(not missing, f"Missing columns in {context}: {missing}")


def manifest(path):
    entries = path.read_text(encoding="utf-8-sig").strip().split(",")
    result = []
    for entry in entries:
        name = entry.strip().replace("\\", "/").rsplit("/", 1)[-1]
        require(name.endswith(".markdup.sorted.bam"),
                f"Unrecognized BAM suffix in manifest {path}: {name!r}")
        result.append(name.removesuffix(".markdup.sorted.bam"))
    require(len(result) == len(set(result)), f"Duplicate samples in {path}")
    return result


def vectors(table, column, samples, counts=False):
    result = []
    for event, value in zip(table["ID"], table[column]):
        tokens = value.split(",")
        require(len(tokens) == len(samples),
                f"Sample vector length mismatch: event {event}, {column}, "
                f"expected {len(samples)}, got {len(tokens)}")
        nums = [number(x, f"{event}/{column}") for x in tokens]
        for x in nums:
            if math.isfinite(x):
                require(x >= 0 and (x.is_integer() if counts else x <= 1),
                        f"Out-of-domain value at {event}/{column}: {x}")
        result.append(nums)
    return pd.DataFrame(result, columns=samples, index=table.index, dtype=float)


def expression(table, de):
    columns(de, ["gene_id", "gene_name", "baseMean"], "DESeq2")
    by_id = defaultdict(list)
    by_symbol = defaultdict(set)
    for row in de.itertuples(index=False):
        gid = gene_id(row.gene_id)
        require(bool(GENE_ID.fullmatch(gid)), f"Invalid human DESeq2 gene ID: {gid!r}")
        value = number(row.baseMean, f"DESeq2/{gid}/baseMean")
        require(math.isnan(value) or value >= 0, f"Negative baseMean: {gid}")
        by_id[gid].append(value)
        symbol = clean(row.gene_name)
        if symbol and symbol.upper() not in MISSING:
            by_symbol[symbol].add(gid)
    resolved = {}
    for gid, vals in by_id.items():
        unique = {"NA" if math.isnan(x) else x for x in vals}
        require(len(unique) == 1, f"Conflicting duplicate DESeq2 rows after version stripping: {gid}")
        resolved[gid] = vals[0]
    records = []
    for row in table.itertuples(index=False):
        gid, symbol = gene_id(row.GeneID), clean(row.geneSymbol)
        require(bool(GENE_ID.fullmatch(gid)), f"Invalid human rMATS GeneID: {gid!r}")
        fallback = by_symbol.get(symbol, set())
        route, chosen = "ensembl", gid
        if gid not in resolved:
            require(len(fallback) <= 1, f"Ambiguous expression symbol fallback: {symbol!r}")
            if fallback:
                chosen, route = next(iter(fallback)), "symbol_fallback"
            else:
                # Both explicit reconciliation attempts are required before absence.
                require(gid not in resolved and not fallback,
                        f"Unreconciled expression absence: {gid}/{symbol}")
                chosen, route = gid, "absent_after_id_and_symbol"
        value = resolved.get(chosen, math.nan)
        status = "expr_unknown" if math.isnan(value) else ("pass" if value > 50 else "fail")
        reason = "missing_baseMean" if chosen in resolved and math.isnan(value) else route
        records.append((gid, symbol, chosen, route, value, status, reason))
    return pd.DataFrame(records, index=table.index, columns=[
        "gene_id", "gene_symbol", "expression_gene_id", "expr_mapping", "baseMean",
        "expr_status", "expr_reason"]), by_symbol


def vast_geometry(row):
    match = COORD.fullmatch(row["COORD"])
    require(match is not None, f"Unsupported HsaEX COORD: {row['EVENT']}/{row['COORD']}")
    chrom, first, last = match.group(1), int(match.group(2)), int(match.group(3))
    length = number(row["LENGTH"], f"{row['EVENT']}/LENGTH")
    require(first >= 1 and last >= first and length == last - first + 1,
            f"VAST one-based closed coordinate/LENGTH assertion failed: {row['EVENT']}")
    # HsaEX FullCO contains upstream junction(s), cassette exon, downstream
    # junction(s) in transcript order; '+' within a field joins alternative sites,
    # and MUST NOT be read as a plus-strand annotation.
    full = row["FullCO"]
    require(full.startswith(chrom + ":"), f"FullCO chromosome mismatch: {row['EVENT']}")
    parts = full[len(chrom) + 1:].split(",")
    require(len(parts) == 3, f"Unsupported HsaEX FullCO grammar: {row['EVENT']}/{full}")
    # The cassette-exon field may carry ALTERNATIVE donor/acceptor sites appended with
    # '+' (24,507 / 223,565 HsaEX rows in the MIAT-KD table), so it is NOT required to
    # equal "first-last". COORD stays authoritative for the exon; the flanking fields
    # are used only to orient the event.
    require(all(re.fullmatch(r"\d+(?:[-+]\d+)*", p) for p in parts),
            f"Unsupported HsaEX FullCO sites: {row['EVENT']}/{full}")
    left = [int(x) for x in re.split(r"[-+]", parts[0])]
    right = [int(x) for x in re.split(r"[-+]", parts[2])]
    plus = max(left) < first and min(right) > last
    minus = min(left) > last and max(right) < first
    # Strand is left UNRESOLVED (never guessed) when the flanks do not bracket the exon
    # monotonically; such events cannot satisfy the freeze's "same strand" requirement
    # and are counted in VAST_AUDIT.tsv / VAST_PARSE.json rather than dropped silently.
    strand = "+" if plus and not minus else ("-" if minus and not plus else None)
    return chrom, first - 1, last, strand


def vast_diff(path, treatment_group, control_group):
    """Per-contrast `vast-tools diff` tab: GENE EVENT <treat> <ctrl> E[dPsi] MV[dPsi]_at_0.95."""
    table = read_table(path, "\t")
    header = list(table.columns)
    require(header[:2] == ["GENE", "EVENT"] and header[4:] == ["E[dPsi]", "MV[dPsi]_at_0.95"],
            f"Unsupported vast-tools diff schema in {path}: {header}")
    require(header[2] == treatment_group and header[3] == control_group,
            f"vast-tools diff group order is not treatment,control in {path}: {header[2:4]}")
    table = table.loc[table["EVENT"].str.startswith("HsaEX")].drop_duplicates()
    require(table["EVENT"].is_unique, f"Conflicting duplicate HsaEX rows in {path}")
    out = {}
    for row in table.to_dict("records"):
        event = row["EVENT"]
        treat = number(row[treatment_group], f"{event}/{treatment_group}")
        ctrl = number(row[control_group], f"{event}/{control_group}")
        dpsi = number(row["E[dPsi]"], f"{event}/E[dPsi]")
        mv = number(row["MV[dPsi]_at_0.95"], f"{event}/MV")
        require(math.isnan(dpsi) or -1 <= dpsi <= 1, f"E[dPsi] outside [-1,1]: {event}")
        require(math.isnan(mv) or mv >= 0, f"Negative MV[dPsi]: {event}")
        # Orientation of the diff tab is asserted, never assumed: E[dPsi] must equal
        # treatment-group PSI minus control-group PSI to the tab's own precision.
        if all(map(math.isfinite, (treat, ctrl, dpsi))):
            # Tolerance is one unit in the tab's last printed decimal (6 dp): the two
            # group PSIs and E[dPsi] are rounded independently. Measured max residual
            # on the QKI-KO tab = 1.000e-06 over 135,096 HsaEX rows. A swapped group
            # pair would give a residual of 2*|dPSI|, which this still catches.
            require(abs((treat - ctrl) - dpsi) <= 2e-6,
                    f"vast-tools diff E[dPsi] != treatment-control: {event}")
        out[event] = (dpsi, mv)
    return out


def vast_index(vast, diff, samples, treatment, symbols, args):
    mapping = {s: OE_ALIASES[s] if args.arm == "MIAT_OE" else s for s in samples}
    columns(vast, ["GENE", "EVENT", "COORD", "LENGTH", "FullCO"] +
            [c for s in samples for c in (mapping[s], mapping[s] + "-Q")], "VAST inclusion")
    selected = vast.loc[vast["EVENT"].str.startswith("HsaEX")].copy()
    require(not selected.empty, "No VAST HsaEX records")
    relevant = ["GENE", "EVENT", "COORD", "LENGTH", "FullCO"] + [
        c for s in samples for c in (mapping[s], mapping[s] + "-Q")]
    selected = selected[relevant].drop_duplicates()
    require(selected["EVENT"].is_unique, "Conflicting duplicate VAST HsaEX event IDs")
    index = defaultdict(list)
    records = []
    parse = defaultdict(int)
    accepted_q = {"VLOW", "LOW", "OK", "SOK"}
    for row in selected.to_dict("records"):
        chrom, start, end, strand = vast_geometry(row)
        parse["hsaex_rows"] += 1
        raw_gene = clean(row["GENE"])
        mapped_gene = gene_id(raw_gene)
        gids = symbols.get(raw_gene, set())
        # Symbol matching is allowed only if the symbol maps to one gene in the
        # combined rMATS/DESeq2 namespace. Ambiguous symbols never become matches.
        resolved_gene = mapped_gene if GENE_ID.fullmatch(mapped_gene) else (
            next(iter(gids)) if len(gids) == 1 else None)
        psi = {s: number(row[mapping[s]], f"{row['EVENT']}/{mapping[s]}") for s in samples}
        require(all(math.isnan(v) or 0 <= v <= 100 for v in psi.values()),
                f"VAST PSI outside [0,100]: {row['EVENT']}")
        cov = {}
        for s in samples:
            token = row[mapping[s] + "-Q"].split(",")[0]
            require(token in accepted_q | {"N", "NA", "", "-"},
                    f"Unrecognized VAST primary coverage code: {row['EVENT']}/{s}/{token}")
            cov[s] = token in accepted_q and math.isfinite(psi[s])
        dpsi, mv = diff.get(row["EVENT"], (math.nan, math.nan))
        in_diff = row["EVENT"] in diff
        parse["absent_from_diff_tab"] += 0 if in_diff else 1
        parse["strand_unresolved"] += 0 if strand else 1
        parse["gene_unresolved"] += 0 if resolved_gene else 1
        record = dict(event=row["EVENT"], chrom=chrom, start=start, end=end,
                      strand=strand or "", gene=resolved_gene, raw_gene=raw_gene,
                      coord_raw=row["COORD"], fullco_raw=row["FullCO"],
                      in_diff_tab=in_diff, dpsi=dpsi, mv=mv,
                      coverage=all(cov.values()), range_pass=bool(mv > 0),
                      sample_cov=cov)
        records.append(record)
        if resolved_gene is not None and strand is not None:
            index[(resolved_gene, chrom, strand)].append(record)
    for key in index:
        index[key].sort(key=lambda v: (v["start"], v["end"], v["event"]))
    starts = {k: [v["start"] for v in vals] for k, vals in index.items()}
    (args.outdir / "VAST_PARSE.json").write_text(json.dumps(dict(parse), indent=2), encoding="utf-8")
    print("VAST_PARSE " + json.dumps(dict(parse)))
    return index, starts, mapping, records


def match_vast(row, index, starts):
    key = (row.expression_gene_id, row.chr, row.strand)
    vals = index.get(key, [])
    positions = starts.get(key, [])
    lo = bisect.bisect_left(positions, row.exonStart_0base - 5)
    hi = bisect.bisect_right(positions, row.exonStart_0base + 5)
    matches = [v for v in vals[lo:hi] if abs(v["end"] - row.exonEnd) <= 5]
    # Geometric choice precedes all significance/coverage tests. Never choose
    # the candidate with the largest dPSI or a passing gate to manufacture B.
    matches.sort(key=lambda v: (abs(v["start"] - row.exonStart_0base) +
                               abs(v["end"] - row.exonEnd),
                               max(abs(v["start"] - row.exonStart_0base),
                                   abs(v["end"] - row.exonEnd)), v["event"]))
    return matches


def region_windows(row):
    """Recovered 2026-07-08 window arithmetic; see module docstring for the proof."""
    left = (row.chr, int(row.upstreamEE), int(row.exonStart_0base))
    exon = (row.chr, int(row.exonStart_0base), int(row.exonEnd))
    right = (row.chr, int(row.exonEnd), int(row.downstreamES))
    merged = (row.chr, int(row.upstreamEE), int(row.downstreamES))
    up, dn = (left, right) if row.strand == "+" else (right, left)
    return {"UPintr": up, "AltEX": exon, "DNintr": dn, "merged": merged}


def exclude_foreground(ledger, fg, bg, harmonise_vast=False):
    foreground = ledger.loc[fg]
    ids = set(foreground["ID"])
    # Ignore strand/gene in exact coordinate disjointness: RELI BED4 is unstranded.
    exons = set(zip(foreground.chr, foreground.exonStart_0base, foreground.exonEnd))
    vast_ids = ({x for v in foreground.vast_candidate_ids for x in json.loads(v)}
                if harmonise_vast else set())
    windows = {w for r in foreground.itertuples() for w in region_windows(r).values()}
    collisions = []
    for row in ledger.itertuples():
        collision = (row.ID in ids or (row.chr, row.exonStart_0base, row.exonEnd) in exons or
                     bool(set(json.loads(row.vast_candidate_ids)) & vast_ids) or
                     any(w in windows for w in region_windows(row).values()))
        collisions.append(collision)
    collisions = pd.Series(collisions, index=ledger.index)
    remaining = bg & ~collisions
    back = ledger.loc[remaining]
    require(not ids & set(back.ID), "Foreground/background event-ID collision")
    require(not exons & set(zip(back.chr, back.exonStart_0base, back.exonEnd)),
            "Foreground/background harmonised exon-coordinate collision")
    require(not vast_ids & {x for v in back.vast_candidate_ids for x in json.loads(v)},
            "Foreground/background harmonised VAST-event collision")
    require(not windows & {w for r in back.itertuples() for w in region_windows(r).values()},
            "Foreground/background emitted-window collision")
    return remaining, bg & collisions



# Preparation adapted from same frozen source; expression is strictly Ensembl-ID-only.
def prepare(args):
    b1, b2 = manifest(args.b1), manifest(args.b2)
    expected_t, expected_c = EXPECTED[args.arm]
    treat, control = (b1, b2) if args.treatment_is == "b1" else (b2, b1)
    require(set(treat) == set(expected_t) and set(control) == set(expected_c),
            "Arm/sample/treatment manifest assertion failed; consult ledger. No automatic sign flip.")
    require(not set(b1) & set(b2), "A sample occurs in both groups")
    table = read_table(args.rmats_se, "\t", rmats=True)
    columns(table, ["ID", "GeneID", "geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd",
                    "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "FDR",
                    "IncLevelDifference", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "SJC_SAMPLE_1",
                    "IJC_SAMPLE_2", "SJC_SAMPLE_2"], "rMATS SE JC")
    require(not table.empty and table.ID.map(lambda x: bool(clean(x))).all(), "Empty rMATS table/event ID")
    distinct = table.drop_duplicates()
    require(distinct.ID.is_unique, "Conflicting duplicate rMATS event IDs; cannot pick a representative")
    table["source_row"] = np.arange(2, len(table) + 2)
    for c in ["exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"]:
        nums = [number(v, c) for v in table[c]]
        require(all(math.isfinite(v) and v >= 0 and v.is_integer() for v in nums),
                f"Invalid rMATS coordinates: {c}")
        table[c] = [int(v) for v in nums]
    require(table.strand.isin(["+", "-"]).all(), "Unknown rMATS strand")
    require(table.chr.str.fullmatch(r"chr[^\s:]+").all(), "Invalid rMATS chromosome")
    require(((table.upstreamES < table.upstreamEE) & (table.upstreamEE <= table.exonStart_0base) &
             (table.exonStart_0base < table.exonEnd) & (table.exonEnd <= table.downstreamES) &
             (table.downstreamES < table.downstreamEE)).all(), "Invalid rMATS SE geometry")
    ledger = table.copy()
    ledger["duplicate_id_rows"] = ledger.groupby("ID")["ID"].transform("size")
    ledger["rmats_fdr"] = [number(v, "FDR") for v in table.FDR]
    ledger["rmats_dpsi"] = [number(v, "IncLevelDifference") for v in table.IncLevelDifference]
    require((ledger.rmats_fdr.isna() | ledger.rmats_fdr.between(0, 1)).all(), "FDR outside [0,1]")
    require((ledger.rmats_dpsi.isna() | ledger.rmats_dpsi.between(-1, 1)).all(), "rMATS dPSI outside [-1,1]")
    means = []
    for group, samples in ((1, b1), (2, b2)):
        inc = vectors(table, f"IncLevel{group}", samples)
        means.append(inc.mean(axis=1, skipna=True))
        jc = vectors(table, f"IJC_SAMPLE_{group}", samples, True) + vectors(
            table, f"SJC_SAMPLE_{group}", samples, True)
        for sample in samples:
            ledger[f"coverage_count__{sample}"] = jc[sample]
            ledger[f"coverage_pass__{sample}"] = jc[sample] >= 10
    delta = means[0] - means[1] if args.treatment_is == "b1" else means[1] - means[0]
    ledger["mean_inc_treatment_minus_control"] = delta
    evaluable = delta.notna() & ledger.rmats_dpsi.notna()
    agrees = evaluable & (np.sign(delta) == np.sign(ledger.rmats_dpsi))
    informative = evaluable & (ledger.rmats_dpsi != 0)
    ledger["orientation_evaluable"] = evaluable
    ledger["orientation_agrees"] = agrees
    orientation = dict(total=len(ledger), evaluable=int(evaluable.sum()),
                       unavailable=int((~evaluable).sum()), informative=int(informative.sum()),
                       agrees=int(agrees.sum()),
                       agreement=float(agrees.sum() / evaluable.sum()) if evaluable.any() else 0,
                       informative_agreement=float((agrees & informative).sum() / informative.sum())
                       if informative.any() else 0)
    (args.outdir / "orientation.json").write_text(json.dumps(orientation, indent=2), encoding="utf-8")
    ledger[["source_row", "ID", "IncLevel1", "IncLevel2", "rmats_dpsi",
            "mean_inc_treatment_minus_control", "orientation_evaluable", "orientation_agrees"]].to_csv(
        args.outdir / "ORIENTATION_AUDIT.tsv", sep="\t", index=False, na_rep="NA")
    print("ORIENTATION " + json.dumps(orientation))
    require(orientation["agreement"] > .95 and orientation["informative_agreement"] > .95,
            "Orientation assertion failed: treatment-control versus IncLevelDifference sign must agree "
            "for >95% of evaluable events AND >95% of nonzero-dPSI evaluable events; "
            "inspect ORIENTATION_AUDIT.tsv, b1/b2 and --treatment-is. No sign flip was applied.")
    expr, symbols = expression_id_only(table, read_table(args.deseq2, ","))
    ledger = pd.concat([ledger, expr], axis=1)
    ledger["expr_eligible"] = ledger.expr_status != "fail"
    require(not ledger.frozen_expr_mapping.eq('ambiguous_symbol_fallback_unavailable').any(),
            'Ambiguous expression symbol fallback prevents frozen VAST gene reconciliation')
    ledger['vast_matching_gene_id'] = ledger.frozen_expression_gene_id
    ledger['vast_gene_mapping'] = ledger.frozen_expr_mapping.replace({'absent_after_id': 'absent_after_id_and_symbol'})
    for gid, symbol in zip(ledger.vast_matching_gene_id, ledger.gene_symbol):
        if symbol and symbol.upper() not in MISSING:
            symbols[symbol].add(gid)
    diff = vast_diff(args.vast_diff, args.vast_treatment_group, args.vast_control_group)
    index, starts, mapping, vast_records = vast_index(read_table(args.vast_inclusion, "\t"), diff,
                                                     b1 + b2, set(treat), symbols, args)
    (args.outdir / "sample_mapping.json").write_text(json.dumps(
        dict(b1=b1, b2=b2, treatment=treat, control=control, vast=mapping), indent=2), encoding="utf-8")
    # Keep the copied matcher and its frozen gene reconciliation separate from
    # the dispatch's strict Ensembl-ID expression gate.
    matching_rows = ledger[['chr', 'strand', 'exonStart_0base', 'exonEnd']].copy()
    matching_rows['expression_gene_id'] = ledger.vast_matching_gene_id
    matches = [match_vast(r, index, starts) for r in matching_rows.itertuples()]
    chosen = [m[0] if m else None for m in matches]
    ledger["vast_match"] = [bool(m) for m in matches]
    ledger["vast_candidate_count"] = [len(m) for m in matches]
    ledger["vast_candidate_ids"] = [json.dumps([v["event"] for v in m]) for m in matches]
    ledger["vast_event"] = [v["event"] if v else "" for v in chosen]
    ledger["vast_start_0base"] = [v["start"] if v else math.nan for v in chosen]
    ledger["vast_end"] = [v["end"] if v else math.nan for v in chosen]
    ledger["vast_strand"] = [v["strand"] if v else "" for v in chosen]
    ledger["vast_gene_id"] = [v["gene"] if v else "" for v in chosen]
    ledger["vast_start_distance_bp"] = (ledger.vast_start_0base - ledger.exonStart_0base).abs()
    ledger["vast_end_distance_bp"] = (ledger.vast_end - ledger.exonEnd).abs()
    matched = ledger.loc[ledger.vast_match]
    require(((matched.vast_start_distance_bp <= 5) & (matched.vast_end_distance_bp <= 5) &
             (matched.vast_strand == matched.strand) &
             (matched.vast_gene_id == matched.vast_matching_gene_id)).all(),
            "VAST match violates boundary tolerance, strand, or reconciled gene identity")
    ledger["vast_dpsi_value"] = [v["dpsi"] if v else math.nan for v in chosen]
    ledger["vast_mv"] = [v["mv"] if v else math.nan for v in chosen]
    ledger["vast_in_diff_tab"] = [bool(v["in_diff_tab"]) if v else False for v in chosen]
    ledger["vast_coverage"] = [v["coverage"] if v else False for v in chosen]
    # freeze:22 "|dPSI| >= 10" / freeze:23 "|dPSI| < 5" are percentage points; the
    # vast-tools diff tab is on the 0-1 scale, hence 0.10 / 0.05.
    ledger["vast_dpsi"] = ledger.vast_dpsi_value.abs() >= 0.10
    ledger["vast_range"] = [bool(v["range_pass"]) if v else False for v in chosen]
    ledger["vast_background_dpsi"] = ledger.vast_dpsi_value.abs() < 0.05
    ledger["concordant_direction"] = ((ledger.rmats_dpsi * ledger.vast_dpsi_value) > 0)
    for sample in b1 + b2:
        ledger[f"vast_coverage_pass__{sample}"] = [v["sample_cov"][sample] if v else False for v in chosen]
    return ledger, orientation


from types import SimpleNamespace
from itertools import product
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment

ROOT = Path(r'E:\rmaps_work')
FROZEN = Path(r'F:\Publication Work\00_PAPER_FIGURES\Fig3_rMAPS_and_RBP-RELI\02_RESUME_2026-09-08\reli_v121')
SOURCE = FROZEN / 'code' / 'build_reli_foregrounds_v2.py'
SOURCE_MD5 = '32be07fd0e6e2edff5f077384f45b352'
PACKAGES = Path(r'F:\rMAPS\inputs\rmaps_packages_2026-09-08')
FAI = Path(r'E:\references\rmaps_genomes\hg38\hg38.fa.fai')
COORD_COLS = ['chr', 'strand', 'exonStart_0base', 'exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE']
COORD_HEADER = ['chr', 'strand', 'exonStart', 'exonEnd', 'firstExonStart', 'firstExonEnd', 'secondExonStart', 'secondExonEnd']
INPUT_HASHES = {}
CREATED = set()


def tracked(path, role='input', expected=None):
    path = Path(path)
    require(path.is_absolute() and path.is_file(), f'Missing absolute input: {path}')
    if str(path) not in INPUT_HASHES:
        with path.open('rb') as h:
            md5 = hashlib.file_digest(h, 'md5').hexdigest()
        INPUT_HASHES[str(path)] = dict(path=str(path), md5=md5, bytes=path.stat().st_size, role=role)
    if expected:
        require(INPUT_HASHES[str(path)]['md5'] == expected.lower(), f'Input MD5 mismatch: {path}')
    return path


def write_text(path, content):
    path.write_text(content, encoding='utf-8')
    CREATED.add(str(path))


def expression_id_only(table, de):
    columns(de, ['gene_id', 'baseMean'], 'DESeq2')
    resolved, symbols = {}, defaultdict(set)
    for row in de.itertuples(index=False):
        gid, value = gene_id(row.gene_id), number(row.baseMean, 'DESeq2/baseMean')
        require(bool(GENE_ID.fullmatch(gid)), f'Invalid human DESeq2 gene ID: {gid!r}')
        require(math.isnan(value) or value >= 0, f'Negative baseMean: {gid}')
        if gid in resolved:
            previous = resolved[gid]
            require(value == previous or (math.isnan(value) and math.isnan(previous)), f'Conflicting duplicate DESeq2 rows after version stripping: {gid}')
        resolved[gid] = value
        symbol = clean(getattr(row, 'gene_name', ''))
        if symbol and symbol.upper() not in MISSING:
            symbols[symbol].add(gid)
    rows = []
    for row in table.itertuples(index=False):
        gid, symbol = gene_id(row.GeneID), clean(row.geneSymbol)
        require(bool(GENE_ID.fullmatch(gid)), f'Invalid human rMATS GeneID: {gid!r}')
        value = resolved.get(gid, math.nan)
        route = 'ensembl' if gid in resolved else 'absent_after_id'
        reason = 'missing_baseMean' if gid in resolved and math.isnan(value) else route
        fallback, frozen_gid, frozen_route = symbols.get(symbol, set()), gid, route
        if gid not in resolved:
            if len(fallback) == 1:
                frozen_gid, frozen_route = next(iter(fallback)), 'symbol_fallback'
            elif len(fallback) > 1:
                frozen_route = 'ambiguous_symbol_fallback_unavailable'
        rows.append(dict(gene_id=gid, gene_symbol=symbol, expression_gene_id=gid,
                         expr_mapping=route, baseMean=value, expr_reason=reason,
                         expr_status='expr_unknown' if math.isnan(value) else ('pass' if value > 50 else 'fail'),
                         frozen_expression_gene_id=frozen_gid, frozen_expr_mapping=frozen_route,
                         frozen_baseMean=resolved.get(frozen_gid, math.nan)))
    return pd.DataFrame(rows, index=table.index), symbols


def ledger_expected(path):
    text = path.read_text(encoding='utf-8')
    out = {}
    for line in text.splitlines():
        cells = [c.strip().replace('**', '').replace(',', '') for c in line.strip('|').split('|')]
        if len(cells) == 9 and cells[0] in EXPECTED and cells[1] in ('A', 'B'):
            out[(cells[0], cells[1])] = tuple(int(cells[i]) for i in (2, 3, 5, 7, 8))
    require(len(out) == 6, 'Could not parse six count rows from frozen ledger section 2')
    return out


def check_log_orientation(prov, arm):
    require(prov['treatment_is'] == 'b1', f'{arm}: frozen provenance does not say treatment=b1')
    b1, b2 = manifest(Path(prov['b1'])), manifest(Path(prov['b2']))
    require(b1 == EXPECTED[arm][0] and b2 == EXPECTED[arm][1], f'{arm}: sample-order mismatch')
    logs = []
    for label in ('A', 'B', 'ledgers'):
        p = tracked(FROZEN / 'foregrounds' / f'{arm}_{label}' / 'command.log', 'orientation evidence')
        s = p.read_text(encoding='utf-8-sig')
        require(re.search(r'--treatment-is\s+b1(?:\s|$)', s), f'{arm}: command.log missing treatment-is b1: {p}')
        for key in ('rmats_se', 'b1', 'b2', 'vast_diff', 'deseq2'):
            require(prov[key].lower() in s.lower(), f'{arm}: {key} disagrees with command.log: {p}')
        logs.append(str(p))
    evidence_path = tracked(ROOT / 'logs' / 'orientation_evidence.json', 'parent original-command verification')
    evidence = [x for x in json.loads(evidence_path.read_text(encoding='utf-8-sig')) if x['arm'] == arm]
    require(len(evidence) == 1 and evidence[0]['treatment_is'] == 'b1', f'{arm}: missing parent original command evidence')
    original = tracked(evidence[0]['rmats_command_log'], 'original rMATS command log', evidence[0]['command_md5'])
    logs.append(str(original))
    return dict(status='PASS', b1=b1, b2=b2, treatment='b1', control='b2', logs=logs, original_command_verification=evidence[0]['result'])


def source_metadata(arm):
    p = tracked(FROZEN / 'foregrounds' / f'{arm}_ledgers' / 'PROVENANCE.json')
    prov = json.loads(p.read_text(encoding='utf-8-sig'))
    require(prov['arm'] == arm, 'Frozen arm provenance mismatch')
    for key, record in prov['input_fingerprints'].items():
        tracked(record['path'], f'{arm} frozen {key}', record['md5'])
    for label in ('A', 'B', 'ledgers'):
        folder = FROZEN / 'foregrounds' / f'{arm}_{label}'
        for name in ('versions.txt', 'MANIFEST.md5', 'PROVENANCE.json'):
            path = folder / name
            if path.is_file():
                tracked(path, 'frozen metadata')
    rmats = PACKAGES / ('MIAT_KD_STRANDFIX_CANDIDATE' if arm == 'MIAT_KD' else arm) / 'SE.MATS.JC.txt'
    tracked(rmats, f'{arm} rMATS byte copy', prov['input_fingerprints']['rmats_se']['md5'])
    return prov, rmats


def expr_gate(values, floor):
    if floor == 'none':
        return pd.Series(True, index=values.index)
    return values.isna() | (values >= 10 if floor == 'bm10' else values > 50)


def select_events(ledger, coverage, floor, bg_fdr, rule, samples, vast_conf='mv'):
    require(vast_conf in ('mv', 'effect-only'), f'Unknown VAST confidence mode: {vast_conf}')
    f = {'n_input': len(ledger)}
    fg = ledger.rmats_fdr < .05
    f['fg_after_fdr'] = int(fg.sum())
    fg &= ledger.rmats_dpsi.abs() >= .10
    f['fg_after_dpsi'] = int(fg.sum())
    bg = ledger.rmats_fdr >= float(bg_fdr)
    f['bg_after_fdr'] = int(bg.sum())
    cutoff = 10 if coverage == 'persample10' else 5
    for i, sample in enumerate(samples, 1):
        cov = ledger[f'coverage_count__{sample}'] >= cutoff
        fg &= cov
        bg &= cov
        f[f'fg_after_coverage_sample{i}'] = int(fg.sum())
        f[f'bg_after_coverage_sample{i}'] = int(bg.sum())
    f['fg_after_coverage'], f['bg_after_coverage'] = int(fg.sum()), int(bg.sum())
    expr = expr_gate(ledger.baseMean, floor)
    fg &= expr
    bg &= expr
    f['fg_after_expression'], f['bg_after_expression'] = int(fg.sum()), int(bg.sum())
    if rule == 'B':
        gates = ('vast_match', 'vast_coverage', 'vast_dpsi', 'vast_range', 'concordant_direction')
        for gate in gates:
            if gate == 'vast_range' and vast_conf == 'effect-only':
                continue
            fg &= ledger[gate]
            f['fg_after_' + gate] = int(fg.sum())
        for gate in ('vast_match', 'vast_coverage', 'vast_background_dpsi'):
            bg &= ledger[gate]
            f['bg_after_' + gate] = int(bg.sum())
    f['bg_before_fg_removal'] = int(bg.sum())
    bg, removed = exclude_foreground(ledger, fg, bg, harmonise_vast=rule == 'B')
    f['bg_collisions_removed'] = int(removed.sum())
    f['bg_after_fg_removal'] = int(bg.sum())
    require(not (fg & bg).any(), 'Foreground/background source-row intersection')
    require(ledger.loc[fg, 'orientation_agrees'].all(), 'Selected foreground orientation disagreement')
    unknown = ledger.baseMean.isna()
    f.update(n_up=int((fg & (ledger.rmats_dpsi > 0)).sum()),
             n_dn=int((fg & (ledger.rmats_dpsi < 0)).sum()), n_bg=int(bg.sum()),
             n_expr_unknown_in_fg=int((fg & unknown).sum()), n_expr_unknown_in_bg=int((bg & unknown).sum()))
    return fg, bg, f


def fai_chromosome(chrom, chroms):
    """Keep exact FAI names; otherwise accept literal chr-prefix removal only."""
    if chrom in chroms:
        return chrom
    candidate = chrom[3:] if chrom.startswith('chr') else chrom
    require(candidate in chroms, f'Chromosome cannot be mapped to FAI: {chrom}')
    return candidate


def write_combination(dest, raw, ledger, fg, bg, row, command, versions, chroms=None):
    dest.mkdir(parents=True, exist_ok=True)
    floor = row['floor']
    annotated = ledger.copy()
    annotated['expr_eligible'] = expr_gate(ledger.baseMean, floor)
    annotated['expr_status'] = np.where(ledger.baseMean.isna(), 'expr_unknown', np.where(annotated.expr_eligible, 'pass', 'fail'))
    extra = [c for c in annotated.columns if c not in raw.columns]
    for name, mask in [('up', fg & (ledger.rmats_dpsi > 0)), ('dn', fg & (ledger.rmats_dpsi < 0)), ('bg', bg)]:
        coords = ledger.loc[mask, COORD_COLS].copy()
        coords.columns = COORD_HEADER
        if chroms is not None:
            original_chroms = coords['chr'].copy()
            coords['chr'] = original_chroms.map(lambda chrom: fai_chromosome(chrom, chroms))
            row[f'n_{name}_chromosome_names_mapped'] = int((coords['chr'] != original_chroms).sum())
            row[f'{name}_chromosome_name_mapping'] = dict(zip(
                original_chroms[coords['chr'] != original_chroms],
                coords.loc[coords['chr'] != original_chroms, 'chr']))
        p = dest / f'{name}.coord.txt'
        coords.to_csv(p, sep='\t', index=False, **({'lineterminator': '\n'} if chroms is not None else {}))
        CREATED.add(str(p))
        # ID.1 retains the second raw ID column; all raw values remain unmodified.
        events = pd.concat([raw.loc[mask], annotated.loc[mask, extra]], axis=1)
        p = dest / f'events_{name}.tsv'
        if chroms is not None:
            events['chr_fai'] = coords['chr']
        events.to_csv(p, sep='\t', index=False, na_rep='NA', quoting=csv.QUOTE_MINIMAL,
                      **({'lineterminator': '\n'} if chroms is not None else {}))
        CREATED.add(str(p))
        require(len(coords) == row['n_' + name], f'Output count mismatch {dest}/{name}')
    if chroms is not None:
        row['all_chromosomes_in_fai'] = True
        row['n_chromosome_names_mapped'] = sum(row[f'n_{name}_chromosome_names_mapped'] for name in ('up', 'dn', 'bg'))
    write_text(dest / 'counts.json', json.dumps(row, indent=2) + '\n')
    write_text(dest / 'command.log', command + '\n')
    write_text(dest / 'versions.txt', versions)


def write_xlsx(table, path):
    table.to_excel(path, index=False, sheet_name='gate_sensitivity')
    wb = openpyxl.load_workbook(path)
    ws = wb.active
    ws.freeze_panes = 'F2'
    ws.auto_filter.ref = ws.dimensions
    for row in ws:
        for cell in row:
            cell.font = Font(name='Arial', size=10)
    for cell in ws[1]:
        cell.font = Font(name='Arial', size=10, bold=True, color='FFFFFF')
        cell.fill = PatternFill('solid', fgColor='264653')
        cell.alignment = Alignment(wrap_text=True)
    ws.row_dimensions[1].height = 42
    for col in ws.columns:
        ws.column_dimensions[col[0].column_letter].width = min(34, max(15, len(str(col[0].value)) + 2))
    notes = wb.create_sheet('definitions')
    for row in [
        ('scope', 'Homo sapiens; GRCh38/hg38; GENCODE v49 input provenance; no GTF loaded'),
        ('foreground', 'FDR < 0.05 and absolute treatment-control dPSI >= 0.10'),
        ('persample10', 'IJC+SJC >= 10 in every replicate'),
        ('groupmin5', 'Minimum IJC+SJC across replicates >= 5 in each group, exactly as dispatched'),
        ('none/bm10/bm50', 'No floor / baseMean >= 10 / baseMean > 50; unknown retained'),
        ('expression mapping', 'Ensembl ID with version removed. No expression symbol fallback; missing IDs remain unknown.'),
        ('VAST gene matching', 'Frozen reconciled gene namespace retained independently in vast_matching_gene_id and vast_gene_mapping.'),
        ('rule B', 'Frozen geometric HsaEX matcher; VLOW+; abs E[dPsi]>=0.10; MV>0; concordant sign'),
        ('background B', 'FDR gate + coverage/expression + HsaEX match + VLOW+ + abs E[dPsi]<0.05'),
        ('collisions', 'Frozen ID, exon, four-region window and rule-B HsaEX candidate exclusions'),
        ('source', str(FROZEN / 'FOREGROUND_LEDGER.md')),
        ('counts', 'Static counts imported from complete source event tables; not an editable analysis model'),
        ('not_applicable', 'Blank funnel entries are inapplicable, not zero'),
    ]:
        notes.append(row)
    notes.column_dimensions['A'].width = 25
    notes.column_dimensions['B'].width = 110
    for row in notes:
        for cell in row:
            cell.font = Font(name='Arial', size=10)
            cell.alignment = Alignment(wrap_text=True, vertical='top')
    wb.save(path)
    check = pd.read_excel(path, sheet_name='gate_sensitivity', dtype={'bg_fdr': str})
    require(len(check) == len(table) and list(check.columns) == list(table.columns), 'XLSX shape differs from TSV')
    for c in ['n_up', 'n_dn', 'n_bg', 'n_expr_unknown_in_fg', 'n_expr_unknown_in_bg']:
        require(check[c].tolist() == table[c].tolist(), f'XLSX count mismatch: {c}')
    CREATED.add(str(path))


def regression(rows, expected, out, frozen_diagnostic, chromosome_failures):
    lines = ['# Regression versus reli_v121', '', '| arm | rule | generated up/dn/bg/unknown-fg/unknown-bg | ledger section 2 | on-disk up/dn/bg | frozen-expression diagnostic | result |', '|---|---|---|---|---|---|---|']
    exact = []
    for row in rows:
        if row['rule'] == 'Beffect':
            continue
        if (row['coverage'], row['floor'], row['bg_fdr']) != ('persample10', 'bm50', '0.5'):
            continue
        arm, rule = row['arm'], row['rule']
        counts = tuple(row[c] for c in ('n_up', 'n_dn', 'n_bg', 'n_expr_unknown_in_fg', 'n_expr_unknown_in_bg'))
        disk = []
        for name in ('INCL_merged.snp', 'SKIP_merged.snp', 'BG_full_merged.bed'):
            path = tracked(FROZEN / 'foregrounds' / f'{arm}_{rule}' / name, 'regression reference')
            with path.open(encoding='utf-8') as h:
                disk.append(sum(1 for line in h if line.strip()))
        status = 'EXACT MATCH' if counts == expected[(arm, rule)] and counts[:3] == tuple(disk) else 'MISMATCH'
        exact.append(status == 'EXACT MATCH')
        diagnostic = frozen_diagnostic[(arm, rule)]
        lines.append(f'| {arm} | {rule} | {counts} | {expected[(arm, rule)]} | {tuple(disk)} | {diagnostic} | {status} |')
    chr_status = ('FAILED: selected events have chromosome names absent from the hg38 FAI; affected inputs are NOT ENGINE READY. See BLOCKED_CHROMOSOMES.md and CHROMOSOME_FAILURES.tsv.' if chromosome_failures else 'PASS: all selected-event chromosome names are present in the hg38 FAI.')
    lines += ['', '- Foreground/background disjointness asserted for every computed combination, using frozen ID/exon/window/HsaEX collision semantics.', '- Chromosome assertion: ' + chr_status, '- All events are retained unchanged; coordinates remain 0-based, half-open. No chromosome renaming or filtering was applied.', '- The dispatch explicitly uses Ensembl-ID expression matching. Frozen code rescues absent IDs through unique symbols; this implementation retains those expression IDs as expr_unknown. VAST matching separately retains the frozen reconciled gene namespace (vast_matching_gene_id and vast_gene_mapping). Frozen baseMean is included for audit; count differences are reported without threshold changes.', '- Rule-B background applies its own foreground collision exclusion, reproducing the frozen implementation. It is not filtered by rule-A foreground windows first.', '- Ledger section 1 rule-B counts are cross-checked against section 2 at load.']
    write_text(out / 'REGRESSION_vs_reli_v121.md', '\n'.join(lines) + '\n')
    return bool(exact) and all(exact)


def arguments():
    p = argparse.ArgumentParser(description='Human hg38 SE event sets; all gate combinations by default.')
    p.add_argument('--arms', nargs='+', choices=list(EXPECTED), default=list(EXPECTED))
    p.add_argument('--coverage', nargs='+', choices=['persample10', 'groupmin5'], default=['persample10', 'groupmin5'])
    p.add_argument('--floor', nargs='+', choices=['none', 'bm10', 'bm50'], default=['none', 'bm10', 'bm50'])
    p.add_argument('--bg-fdr', nargs='+', choices=['0.5', '0.05'], default=['0.5', '0.05'])
    p.add_argument('--rule', nargs='+', choices=['A', 'B'], default=['A', 'B'])
    p.add_argument('--vast-conf', choices=['mv', 'effect-only'], default='mv',
                   help='Rule B confidence: frozen MV>0 (default), or effect-only ruleBeffect.')
    p.add_argument('--outdir', type=Path, default=ROOT / 'event_sets')
    p.add_argument('--logdir', type=Path, default=ROOT / 'logs')
    p.add_argument('--fai', type=Path, default=FAI)
    p.add_argument('--overwrite', action='store_true', help='Replace only this script\'s named output files in selected combinations.')
    return p.parse_args()


def main():
    args = arguments()
    for path in (args.outdir, args.logdir):
        require(path.is_absolute() and path.resolve().is_relative_to(ROOT.resolve()), f'Output must be under {ROOT}: {path}')
    if not args.overwrite:
        for name in ('gate_sensitivity.tsv', 'gate_sensitivity.xlsx', 'PROVENANCE.md', 'REGRESSION_vs_reli_v121.md', 'CHROMOSOME_FAILURES.tsv', 'BLOCKED_CHROMOSOMES.md', 'command.log', 'versions.txt'):
            require(not (args.outdir / name).exists(), f'Output exists; use --overwrite explicitly: {args.outdir / name}')
        for arm, coverage, floor, bg_fdr, rule in product(args.arms, args.coverage, args.floor, args.bg_fdr, args.rule):
            disk_rule = 'Beffect' if rule == 'B' and args.vast_conf == 'effect-only' else rule
            dest = args.outdir / arm / f'{coverage}_{floor}_bgfdr{bg_fdr}_rule{disk_rule}'
            require(not dest.exists() or not any(dest.iterdir()), f'Output exists; use --overwrite explicitly: {dest}')
    for path in (args.outdir, args.logdir):
        path.mkdir(parents=True, exist_ok=True)
    command = subprocess.list2cmdline(sys.orig_argv)
    versions = f'builder=1.0\npython={sys.version}\npandas={pd.__version__}\nnumpy={np.__version__}\nopenpyxl={openpyxl.__version__}\nplatform={platform.platform()}\nscope=Homo sapiens; GRCh38/hg38; GENCODE v49 input provenance; no GTF loaded\n'
    write_text(args.logdir / 'command.log', command + '\n')
    write_text(args.outdir / 'command.log', command + '\n')
    write_text(args.outdir / 'versions.txt', versions)
    tracked(SOURCE, 'copied source', SOURCE_MD5)
    metadata_record = ROOT / 'logs' / 'frozen_metadata_read.json'
    if metadata_record.is_file():
        tracked(metadata_record, 'parent frozen metadata read record')
    ledger_path = tracked(FROZEN / 'FOREGROUND_LEDGER.md', 'protocol/count reference')
    expected = ledger_expected(ledger_path)
    ledger_text = ledger_path.read_text(encoding='utf-8')
    require('**PASSED** all 3 arms' in ledger_text, 'Ledger orientation PASS statement missing')
    for arm in EXPECTED:
        pattern = rf'\| \*\*{arm}\*\* \| (\d+) \| (\d+) \| \*\*(\d+)\*\*'
        match = re.search(pattern, ledger_text)
        require(match and tuple(map(int, match.groups()[:2])) == expected[(arm, 'B')][:2], 'Ledger sections 1/2 disagree')
    fai = tracked(args.fai, 'human hg38 chromosome index')
    chroms = {line.split('\t')[0] for line in fai.read_text().splitlines() if line.strip()}
    orientations, expression_audits, rows, frozen_diagnostic = {}, {}, [], {}
    chromosome_failures, raw_chromosome_findings = [], {}
    skill_path = tracked(ROOT / 'logs' / 'skill_sources.json', 'live scientific source hash registry')
    for skill in json.loads(skill_path.read_text(encoding='utf-8-sig')):
        p = tracked(skill['path'], 'live skill source', skill['md5'])
        with p.open('rb') as handle:
            sha = hashlib.file_digest(handle, 'sha256').hexdigest()
        require(sha == skill['sha256'], f'Skill source changed during run: {p}')
        INPUT_HASHES[str(p)]['sha256'] = sha
    combinations = list(product(args.coverage, args.floor, args.bg_fdr, args.rule))
    for arm in args.arms:
        prov, rmats_path = source_metadata(arm)
        orientations[arm] = check_log_orientation(prov, arm)
        auditdir = args.logdir / arm
        auditdir.mkdir(parents=True, exist_ok=True)
        ns = SimpleNamespace(**{k: Path(prov[k]) for k in ('b1', 'b2', 'deseq2', 'vast_inclusion', 'vast_diff')},
                             arm=arm, rmats_se=rmats_path, outdir=auditdir, treatment_is='b1',
                             vast_treatment_group=prov['vast_treatment_group'], vast_control_group=prov['vast_control_group'])
        ledger, orientation = prepare(ns)
        orientations[arm]['arithmetic'] = orientation
        require(orientation['informative_agreement'] == 1.0, f'{arm}: nonzero dPSI direction disagreement')
        unsupported = ~ledger.chr.isin(chroms)
        raw_chromosome_findings[arm] = {str(k): int(v) for k, v in ledger.loc[unsupported, 'chr'].value_counts().items()}
        raw_chr_path = auditdir / 'RAW_CHROMOSOMES_ABSENT_FROM_FAI.tsv'
        ledger.loc[unsupported, ['source_row', 'ID', 'GeneID', 'chr', 'strand', 'exonStart_0base', 'exonEnd']].to_csv(raw_chr_path, sep='\t', index=False)
        CREATED.add(str(raw_chr_path))
        if args.vast_conf == 'effect-only':
            unsupported &= ~ledger.chr.map(lambda chrom: chrom[3:] if chrom.startswith('chr') else chrom).isin(chroms)
        diagnostic_dir = auditdir / 'frozen_expression_diagnostic'
        diagnostic_dir.mkdir(parents=True, exist_ok=True)
        frozen_ledger = ledger.copy()
        frozen_ledger['baseMean'] = ledger.frozen_baseMean
        for diagnostic_rule in ('A', 'B'):
            _, _, diagnostic_counts = select_events(frozen_ledger, 'persample10', 'bm50', '0.5', diagnostic_rule, EXPECTED[arm][0] + EXPECTED[arm][1])
            frozen_diagnostic[(arm, diagnostic_rule)] = tuple(diagnostic_counts[c] for c in ('n_up', 'n_dn', 'n_bg', 'n_expr_unknown_in_fg', 'n_expr_unknown_in_bg'))
        del frozen_ledger
        write_text(diagnostic_dir / 'DIAGNOSTIC_METHOD.md',
                   '# Frozen expression diagnostic\n\n- Reuses the prepared ledger and the identical frozen VAST matching namespace.\n'
                   '- Replaces baseMean with frozen_baseMean for every row; select_events recalculates expression eligibility and unknown counts.\n'
                   '- No second VAST parse, coordinate remapping, or gate alteration.\n')
        write_text(diagnostic_dir / 'command.log', command + '\n')
        write_text(diagnostic_dir / 'versions.txt', versions)
        for p in diagnostic_dir.iterdir():
            if p.is_file():
                CREATED.add(str(p))
        raw = pd.read_csv(rmats_path, sep='\t', dtype=str, keep_default_na=False, encoding='utf-8-sig')
        require(len(raw) == len(ledger), 'Raw/annotated row alignment mismatch')
        fallback = ledger.frozen_expr_mapping == 'symbol_fallback'
        expression_audits[arm] = dict(frozen_symbol_fallback_rows=int(fallback.sum()), id_only_unknown_rows=int(ledger.baseMean.isna().sum()),
                                     frozen_fallback_measured_rows=int((fallback & ledger.frozen_baseMean.notna()).sum()))
        p = auditdir / 'EXPRESSION_MAPPING_DIFFERENCES.tsv'
        ledger.loc[fallback, ['source_row', 'ID', 'gene_id', 'gene_symbol', 'baseMean', 'expr_status', 'frozen_expression_gene_id', 'frozen_baseMean']].to_csv(p, sep='\t', index=False, na_rep='NA')
        CREATED.add(str(p))
        write_text(auditdir / 'command.log', command + '\n')
        write_text(auditdir / 'versions.txt', versions)
        for coverage, floor, bg_fdr, rule in combinations:
            disk_rule = 'Beffect' if rule == 'B' and args.vast_conf == 'effect-only' else rule
            dest = args.outdir / arm / f'{coverage}_{floor}_bgfdr{bg_fdr}_rule{disk_rule}'
            require(args.overwrite or not dest.exists() or not any(dest.iterdir()), f'Output exists; use --overwrite explicitly: {dest}')
            fg, bg, counts = select_events(ledger, coverage, floor, bg_fdr, rule, EXPECTED[arm][0] + EXPECTED[arm][1], vast_conf=args.vast_conf)
            counts.update(arm=arm, coverage=coverage, floor=floor, bg_fdr=bg_fdr, rule=disk_rule,
                          fg_frozen_symbol_fallback_rows=int((fg & fallback).sum()), bg_frozen_symbol_fallback_rows=int((bg & fallback).sum()),
                          foreground_background_disjoint=True, all_chromosomes_in_fai=not bool(((fg | bg) & unsupported).any()))
            for population, mask in [('up', fg & (ledger.rmats_dpsi > 0)), ('dn', fg & (ledger.rmats_dpsi < 0)), ('bg', bg)]:
                missing_chr_counts = ledger.loc[mask & unsupported, 'chr'].value_counts()
                counts[f'n_{population}_chromosome_absent_from_fai'] = int(missing_chr_counts.sum())
                counts[f'{population}_chromosomes_absent_from_fai'] = json.dumps(sorted(missing_chr_counts.index.tolist()))
                for chrom, count in missing_chr_counts.items():
                    chromosome_failures.append(dict(arm=arm, coverage=coverage, floor=floor, bg_fdr=bg_fdr,
                                                    rule=disk_rule, population=population, chromosome=chrom, n_events=int(count),
                                                    candidate_fai_name=chrom[3:] if chrom.startswith('chr') and chrom[3:] in chroms else ''))
            rows.append(counts)
            write_combination(dest, raw, ledger, fg, bg, counts, command, versions,
                              chroms=chroms if args.vast_conf == 'effect-only' else None)
            print(json.dumps({k: counts[k] for k in ['arm', 'coverage', 'floor', 'bg_fdr', 'rule', 'n_up', 'n_dn', 'n_bg']}), flush=True)
        for p in auditdir.iterdir():
            if p.is_file():
                CREATED.add(str(p))
    table = pd.DataFrame(rows).sort_values(['arm', 'rule', 'coverage', 'floor', 'bg_fdr']).reset_index(drop=True)
    first = ['arm', 'coverage', 'floor', 'bg_fdr', 'rule', 'n_up', 'n_dn', 'n_bg', 'n_expr_unknown_in_fg', 'n_expr_unknown_in_bg']
    table = table[first + [c for c in table if c not in first]]
    p = args.outdir / 'gate_sensitivity.tsv'
    table.to_csv(p, sep='\t', index=False, na_rep='NA')
    CREATED.add(str(p))
    write_xlsx(table, args.outdir / 'gate_sensitivity.xlsx')
    chr_columns = ['arm', 'coverage', 'floor', 'bg_fdr', 'rule', 'population', 'chromosome', 'n_events', 'candidate_fai_name']
    chr_path = args.outdir / 'CHROMOSOME_FAILURES.tsv'
    pd.DataFrame(chromosome_failures, columns=chr_columns).to_csv(chr_path, sep='\t', index=False)
    CREATED.add(str(chr_path))
    blocked_combinations = [row for row in rows if not row['all_chromosomes_in_fai']]
    chr_report = ['# Chromosome validation', '',
                  '**BLOCKED: affected inputs are NOT ENGINE READY.**' if chromosome_failures else '**PASS: every selected event chromosome is present in the supplied hg38 FAI.**', '',
                  f'- Reference: `{fai}`.', f'- Affected combinations: {len(blocked_combinations)} of {len(rows)}.',
                  '- No event was dropped and no chromosome name or coordinate was changed.',
                  '- CHROMOSOME_FAILURES.tsv gives per-combination, population, chromosome event counts; these are repeated gate memberships, not unique-event totals.',
                  '- Each counts.json and sensitivity-table row carries all_chromosomes_in_fai plus absent-event counts and chromosome lists for up/dn/bg.',
                  '- Counts and expression/regression findings remain reviewable; chromosome failure prevents complete status.', '',
                  '## Raw input findings', '```json', json.dumps(raw_chromosome_findings, indent=2), '```']
    candidate_aliases = {chrom: chrom[3:] for finding in raw_chromosome_findings.values() for chrom in finding if chrom.startswith('chr') and chrom[3:] in chroms}
    chr_report += ['', '## Exact-name discrepancy',
                   '- The following raw chromosome strings are absent, but removal of the literal chr prefix produces names present in the supplied FAI. This is a candidate naming repair only; reference-sequence identity has not been established here.',
                   '- Any normalization requires explicit reviewed authorization. None was applied in this dispatch.',
                   '```json', json.dumps(candidate_aliases, indent=2), '```']
    if args.vast_conf == 'effect-only':
        chr_report = ['# Chromosome validation', '',
                      '- Effect-only exports retain exact FAI names; otherwise the literal chr prefix is removed only when the result exists in the FAI.',
                      '- Original event chromosome names remain in events_*.tsv; chr_fai records the exported name.',
                      '- Mapping counts and dictionaries are recorded per population in counts.json.',
                      f'- Unmappable selected-event records: {len(chromosome_failures)}.',
                      f'- Reference: `{fai}`.']
    write_text(args.outdir / 'BLOCKED_CHROMOSOMES.md', '\n'.join(chr_report) + '\n')
    if args.vast_conf == 'effect-only':
        exact = None
        write_text(args.outdir / 'REGRESSION_vs_reli_v121.md',
                   '# Regression versus reli_v121\n\nNOT APPLICABLE: effect-only ruleBeffect drops the MV>0 gate. '
                   'Frozen rule-B regression must be run separately with --vast-conf mv.\n')
    else:
        exact = regression(rows, expected, args.outdir, frozen_diagnostic, chromosome_failures)
    write_text(args.logdir / 'orientation_checks.json', json.dumps(orientations, indent=2) + '\n')
    write_text(args.logdir / 'expression_mapping_summary.json', json.dumps(expression_audits, indent=2) + '\n')
    for record in INPUT_HASHES.values():
        with Path(record['path']).open('rb') as handle:
            current = hashlib.file_digest(handle, 'md5').hexdigest()
        require(current == record['md5'], f"Input changed during run: {record['path']}")
    script = tracked(Path(__file__).resolve(), 'delivered script')
    md = ['# Provenance', '', '- Scope: Homo sapiens / GRCh38 (hg38) / GENCODE v49, inherited from verified input provenance; no FASTA sequence or GTF loaded.', f'- Chromosome reference: `{fai}`.', '- Git status: not queried; this workspace is not treated as a repository.', '- No enrichment, plotting, sampling, capping, or inferential statistics.', '- Expression matching: version-stripped Ensembl ID only; absent or missing-baseMean retained. Frozen symbol fallback recorded for comparison.', '- groupmin5 is the explicitly dispatched minimum over replicates >=5 in each group; it is not the historical group-mean filter described in the frozen ledger caveat.', '- Raw rMATS repeated ID is preserved under pandas column ID.1; both IDs are verified equal before selection.', '- Rule-B matching/collision functions copied verbatim from frozen source; no runtime source import.', '- VAST matching retains the frozen reconciled gene namespace independently of expression gating: vast_matching_gene_id and vast_gene_mapping. expression_gene_id/baseMean remain ID-only.', '- Frozen-expression diagnostic reuses identical VAST matches and substitutes frozen_baseMean for every row before expression gating.', f'- Exact invocation: `{command}`', f'- Generated gate combinations: {len(table)}.', f'- Canonical regression: {"EXACT MATCH" if exact else "MISMATCH or canonical gates not selected"}.', '', '## Versions', '```', versions.rstrip(), '```', '', '## Orientation', '```json', json.dumps(orientations, indent=2), '```', '', '## Expression mapping audit', '```json', json.dumps(expression_audits, indent=2), '```', '', '## Inputs and script', '| role | path | bytes | md5 |', '|---|---|---:|---|']
    for record in INPUT_HASHES.values():
        md.append(f"| {record['role']} | `{record['path']}` | {record['bytes']} | {record['md5']} |")
    md += ['', '## Chromosome readiness',
           f'- Chromosome assertion: {"FAILED; affected inputs NOT ENGINE READY" if chromosome_failures else "PASS"}.',
           f'- Affected combinations: {len(blocked_combinations)} of {len(rows)}.',
           '- All event rows, chromosome names and coordinates retained unchanged. See BLOCKED_CHROMOSOMES.md and CHROMOSOME_FAILURES.tsv.',
           '- Candidate FAI aliases are reported only when literal prefix removal gives an existing FAI name; no normalization applied.']
    if args.vast_conf == 'effect-only':
        md = [line for line in md if not any(token in line for token in (
            'Canonical regression:', 'All event rows, chromosome names', 'Candidate FAI aliases'))]
        md += ['', '## Effect-only rule',
               '- Rule A selection remains unchanged; ruleBeffect differs from ruleB only by omitting the VAST MV>0 foreground gate.',
               '- Frozen matching, orientation, background gates and foreground-collision exclusion are reused.',
               '- Canonical regression: NOT APPLICABLE to effect-only; frozen rule B must be checked separately.',
               '- Coordinate exports use LF line endings and names in the supplied FAI. Literal chr-prefix removal is allowed only for existing FAI names.',
               '- Raw chromosome names remain in event tables alongside chr_fai; per-population mapping counts are in counts.json.']
    write_text(args.outdir / 'PROVENANCE.md', '\n'.join(md) + '\n')
    write_text(args.outdir / 'INPUT_HASHES.json', json.dumps(list(INPUT_HASHES.values()), indent=2) + '\n')
    CREATED.add(str(script))
    manifest_path = args.outdir / 'DELIVERABLE_MANIFEST.txt'
    CREATED.add(str(manifest_path))
    write_text(manifest_path, '\n'.join(sorted(CREATED)) + '\n')
    print(f'BUILD_ARTIFACTS_WRITTEN combinations={len(table)} canonical_regression={"NOT_APPLICABLE_EFFECT_ONLY" if args.vast_conf == "effect-only" else ("EXACT_MATCH" if exact else "MISMATCH_OR_NOT_SELECTED")} chromosome_readiness={"BLOCKED" if chromosome_failures else "PASS"}', flush=True)
    try:
        require(not chromosome_failures, f'Chromosome assertion failed: {len(blocked_combinations)} combinations contain selected chromosome names absent from {fai}; affected inputs NOT ENGINE READY')
    except AssertionError as exc:
        print(f'AssertionError: {exc}', file=sys.stderr, flush=True)
        return 3
    return 0 if exact or args.vast_conf == 'effect-only' else 2


if __name__ == '__main__':
    sys.exit(main())
