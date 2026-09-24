---
name: rmaps3-quick
description: "Run rMAPS3 the simple way on any rMATS SE result and get the tool's own answer back as an explorable workbook plus the main-layer region lollipops. Invoke for 'run rMAPS', 'rMAPS3 motif map', 'RNA map for these splicing events', 'which RBP motifs are enriched around these exons', 'rMAPS on my SE.MATS.JC.txt', or 'redo the rMAPS run for <arm>'. Produces the authors' released-engine layer ONLY: raw p, no BH, no calibration, no lab statistics painted onto a third-party tool's output; its four figures draw that layer alone. For the calibrated supplement, its figures and rank stability use rmaps3-full instead. Human / GRCh38 (hg38) / GENCODE v49 by default; mouse must be named explicitly."
---

# rMAPS3 quick run — the authors' layer, one command

One wrapper, `tools/rmaps3_skill_run.py` in the lab fork `E:\Claude\rMAPS3` (branch `lab/miat-qki`).
It runs the collaborators' released engine exactly as the lab launchers do, proves the outputs
reproduce, and packages them with four main-layer figures. It never edits engine code and never
recomputes an engine number.

Companion: **`rmaps3-full`** adds the calibrated supplement, which carries the p and q to report
(target-exon cluster permutation, RBP-level min-P), its version 4.3.3 lollipops, rank stability and
the rank workbook. Depth on every footgun below: `E:\Claude\rMAPS3\LESSONS.md`; tool index:
`E:\Claude\rMAPS3\tools\README.md` (release tag `lab-v1.0-2026-09-24`).

## When to use

| Use this skill | Use something else |
|---|---|
| A new or repeated SE motif map for one arm | Non-SE event classes — this wrapper is SE only |
| "What does rMAPS itself say about these events?" | A calibrated p or the supplement figures → `rmaps3-full` |
| Reproducing an existing released run byte for byte | CLIP maps, MISO conversion → the engine CLI directly |

## Exact command

```bash
E:/rmaps_venv2/Scripts/python.exe E:/Claude/rMAPS3/tools/rmaps3_skill_run.py \
  --mode quick --arm <ARM> --out <OUT_DIR> \
  --up <up.coord.txt> --dn <dn.coord.txt> --bg <bg.coord.txt> \
  --genome-root E:/references/rmaps_genomes --genome hg38 \
  --engine released --engine-root E:/Claude/rMAPS3_upstream --stat-method mannwhitney \
  --gate-counts <gate record counts.json> --gate-rule <A|B|Beffect> \
  --gtf E:/references/gencode_v49/gencode.v49.primary_assembly.annotation.gtf \
  --spliceosome-list "F:/RNA-SEQ-ANALYSIS/RBP-RELI/RBP-RELI/reli/data/spliceosome_census.txt" \
  --broad-binders-list "F:/RNA-SEQ-ANALYSIS/RBP-RELI/RBP-RELI/reli/data/broad_binders.txt" \
  --positive-control QKI --arm-label "<title>"
```

`--arm-label` is the figure and index title and is required whenever figures are drawn; the wrapper
never derives a title from the arm name. For the seven dissertation arms copy the `arm_label` of the
arm's row in `E:\Claude\rMAPS3\data\arm_labels_dissertation.tsv`; for any other arm write a title
that states the event source and rule, e.g. "QKI knockout (rMATS-strict)" for a rule-A arm.

Since 2026-09-24 the wrapper carries no site paths: `--genome-root`, `--engine-root` (released
engine) and, whenever figures are drawn, `--gtf` / `--spliceosome-list` / `--broad-binders-list`
have no default, and a missing one is refused before any output exists. With `--no-figures` the
last three can be dropped. Run it with `E:\rmaps_venv2`, the interpreter that satisfies
`requirements-lock.txt`.

The event-set rule comes only from `--gate-rule` or the gate record's `rule` field, never from the arm name; the wrapper prints it in `versions.txt`, the workbook README and the figure footer. Pre-split input without either is refused (exit 2). Raw `--rmats-se`
input is rule A, because the portable builder implements rule A only; `--gate-rule B` or `Beffect`
with raw input is refused, whatever the arm is called.

A run is `complete` (exit 0) only when every artefact its mode, engine and statistic require (`REQUIRED_ARTIFACTS` in the wrapper) exists, is non-empty and has its sha256 in `run_manifest.json`; otherwise it is `INCOMPLETE` (exit 4), or `complete_partial` (exit 0) with `--allow-partial`. A failed positive control gives `complete_with_failed_positive_control` (exit 3), an exception `failed` (exit 1), a refused argument exit 2 before any output exists.

For a lab arm the gate record is `F:\rMAPS\event_sets\<ARM base>\<gate>_rule<R>\counts.json`,
e.g. `F:\rMAPS\event_sets\QKI_KO\persample10_bm50_bgfdr0.5_ruleB\counts.json` for `QKI_KO_B`.

From a raw rMATS table instead of pre-split files, swap the three set flags and `--gate-counts`
for the line below; the wrapper then writes its own gate record:

```bash
  --rmats-se <SE.MATS.JC.txt> --filter <gates.json>
```

`gates.json` keys: `fdr`, `min_abs_dpsi`, `min_jc_per_sample`, `bg_fdr_min`, `expr_table`,
`base_mean_floor`. Any other key is refused, not ignored.

## Inputs and locked defaults

| Flag | Default (lab value when required) | Note |
|---|---|---|
| `--genome-root` / `--genome` | required (lab: `E:\references\rmaps_genomes`) / `hg38` | FASTA read is `<root>/<build>/<build>.fa`; its `.fai` must exist |
| `--engine` / `--engine-root` | `released` / required for released (lab: `E:\Claude\rMAPS3_upstream`, the worktree at tag `upstream-base-2026-09-16`) | `released` = the collaborators' code at `b9a9dce`; `audited` = the fork's repaired engine, default root = the fork itself |
| `--stat-method` | `mannwhitney` | the published observational unit; `fisher` is the tool default and counts hits, not exons |
| `--window` / `--step` / `--intron` / `--exon` | 50 / 1 / 250 / 50 | locked; these are also the engine defaults |
| motif tables | the engine's own `knownMotifs.human.mouse.txt` + `ESRP.like.motif.txt` | 126 motif keys = known RBP motifs plus the twelve GU-rich ESRP-like hexamers |
| `--workers` | 1 | audited engine only; the released engine has no worker flag |
| `--alias-table` | `E:\Claude\rMAPS3\data\rbp_alias_hgnc_2026-09-17.tsv` | maps motif-table names to HGNC symbols |
| `--blas-threads` | 4 | exported to every `*_NUM_THREADS` variable before any child starts |
| `--gate-counts` | none | pre-split inputs only: the gate record `counts.json`. Its `n_up`/`n_dn`/`n_bg` must equal the input files or the run stops; its expression-unknown counts go into the figure footer. Required when figures are drawn from pre-split inputs |
| `--gate-rule` | none; required for pre-split input unless the gate record states `rule` | the event-set rule (`A`, `B`, `Beffect`). It must agree with the gate record's `rule` when both are given. Raw `--rmats-se` input is rule A |
| `--allow-partial` | off | a run missing a required artefact exits 0 as `complete_partial` instead of exiting 4 as `INCOMPLETE`; it never yields `complete` |
| `--no-figures` | off | skip the four figures; the workbook and index are otherwise unchanged |
| `--gtf` | required with figures (lab: `E:\references\gencode_v49\gencode.v49.primary_assembly.annotation.gtf`) | figure naming only: every motif-table name must resolve to a GENCODE v49 `gene_name`, directly or through the alias table |
| `--spliceosome-list` / `--broad-binders-list` | required with figures (lab: `F:\RNA-SEQ-ANALYSIS\RBP-RELI\RBP-RELI\reli\data\spliceosome_census.txt` / `broad_binders.txt`) | define the `_noSpliceosome_noBroad` figure variant; SRSF1 is always kept |
| `--top-n` / `--arm-label` / `--underpowered-arms` | 10 / required with figures / none | rows per panel; the figure and index title, from `data/arm_labels_dissertation.tsv` for a dissertation arm; the underpowered banner |

**Gates, when the wrapper pre-splits for you.** The lab default is the RBP-RELI `reli_v121` ledger
gate verbatim, and it is a composite whose components must each be named: a significance filter on
the rMATS FDR, a magnitude filter on `|dPSI|`, a per-sample junction-coverage floor on `IJC+SJC` in
**every** sample, and a DESeq2 `baseMean` expression floor; background takes the same coverage and
expression gates, a high-FDR floor instead of the significance filter, and has the foreground and
its overlapping windows removed. Read the live values off the ledger source
`F:\Publication Work\00_PAPER_FIGURES\Fig3_rMAPS_and_RBP-RELI\02_RESUME_2026-09-08\reli_v121\code\build_reli_foregrounds_v2.py`
and the gate tree `F:\rMAPS\event_sets\<ARM>\` — never from memory, and never let one shorthand
stand for the whole composite. Genes absent from the expression table are `expr_unknown` and are
retained.

## Species and assembly

Human / GRCh38 (hg38) / GENCODE v49; FASTA `E:\references\rmaps_genomes\hg38\hg38.fa`. A run with a
`mm*` build or a mouse `--species` prints a MOUSE REFERENCE IN PLAY banner on stderr and records it
in the README sheet. Never merge or compare a mouse arm with a human arm.

## Outputs and how to read them

| File | What it is |
|---|---|
| `index.html` | start here: event counts, links to the tool's own tables and maps, the positive control, the statistic caveat |
| `figures/<ARM>_SE_byRBP_released_ranksum_rawP.{png,svg}` | per pooled region × direction, the top RBPs by their best motif's released rank-sum p |
| `figures/<ARM>_SE_byMotif_released_ranksum_rawP.{png,svg}` | the same panels, one dot per motif |
| `figures/<ARM>_SE_{byRBP,byMotif}_released_ranksum_rawP_noSpliceosome_noBroad.{png,svg}` | the same two figures after dropping the lab core-spliceosome and broad-binder lists |
| `figures/positive_control_audit.tsv` | with `--positive-control`: rank 1 checked on every drawn control panel, both kinds, both variants |
| `figures/selection_audit.tsv`, `exclusion_audit.tsv`, `naming_audit.tsv`, `layout_report.json` | every drawn row with its motif, sub-region, p and ratio; what the exclusion variant dropped; the HGNC mapping; the layout audit, y-scale cap, truncated stems and dot-size key |
| `figures/motif_scores/<ARM>/per_motif_regions.tsv` | the engine's own motif scores per sub-region, read from the verified archives, released-layer columns only; dot size comes from here |
| `quick_summary.xlsx` | README sheet (inputs, gates, engine commit, statistic, caveats); `INCLUDED` and `SKIPPED` sheets, one row per motif with the smallest p in each of the 8 sub-regions, the 3 pooled-region minima plus Flanking Exon, and the engine's own motif score at the window that produced each minimum; `RBP_best_motif` sheet ranking each RBP by its best motif within a panel; `positive_control` sheet when `--positive-control` is given |
| `engine/<ARM>/pVal.{up,dn}.vs.bg.RNAmap.txt` | the tool's own root tables, unmodified |
| `engine/<ARM>/maps/` | the authors' per-motif RNA maps, unmodified, one PNG + PDF per motif |
| `counts/<ARM>/*.counts.npz` + `VERIFY.md` | per-motif count archives and the proof they reproduce the root tables |
| `command.log`, `versions.txt`, `md5.txt`, `logs/` | provenance; `versions.txt` carries the `gate_rule` line; `md5.txt` covers every file except itself and `run_manifest.json` |
| `run_manifest.json` | written last: `status`, `missing`, `required_artefacts`, `registered` (every file a writer wrote) and `inventory`, the sha256 of every registered file, `md5.txt` and `command.log` included |

**What the figures show.** They are drawn by the functions of `tools/build_region_lollipops_v4.py`,
figure version 4.3.3, unmodified. Stem = −log10 of the released rank-sum p, raw, smallest over the
50-nt windows of a pooled region. Colour = direction hue (included gold `#E69F00`, skipped blue
`#0072B2`) deepening with the same raw p from the 0.05 tint to the y-cap floor; grey `#999999` at
p ≥ 0.05. Dot size = the tool's own motif-score
ratio, changed ÷ background hits per event per window, at the window that gave the minimum.
Included rises, skipped hangs. One shared y-scale is capped at 1.25 × the largest value outside the
top RBP; stems above it are truncated and labelled. The footer gives n events, which are rMATS SE
rows. The legend says "use for RBP ORDER only".

**What they do not show.** No calibrated p, no q, no multiple-testing adjustment, and no binding
call. The panels match the main layer `rmaps3-full` draws for the same run: checked 2026-09-22 on
QKI_KO_B, where the selections, motif scores and size key were identical and the SVG text differed
only in the corrected expression-unknown footer count.

Pooled regions: Upstream Intron = min(UpstreamExonIntron, UpstreamIntron); Exon Body =
min(TargetExon_5prime, TargetExon-3prime); Downstream Intron = min(DownstreamIntron,
DownstreamExonIntron). Flanking Exon is kept as a column and is never plotted.

## Verification

1. **Archive verification is automatic and blocking.** Every root-table value is recomputed from the
   archives alone; a mismatch beyond `1e-6` in `|log10 p|` aborts the run. Read `VERIFY.md`.
2. **Reproduction check.** Re-running an existing arm must give md5-identical root tables. Proven
   2026-09-21 on QKI_KO_B against `E:\rmaps_stat\released_mannwhitney\QKI_KO_B\`.
3. **Positive control.** For a QKI-KO arm, `--positive-control QKI` expects QKI first in
   INCLUDED × Upstream Intron and SKIPPED × Downstream Intron. A failure exits 3 unless
   `--positive-control-advisory` is passed. The same rank-1 check runs on every drawn figure panel
   and feeds the same exit code. Do not expect QKI to lead a MIAT-KD or MIAT-OE arm; the data there
   are presented agnostically.
4. **Event counts** are printed per set and written to `event_counts.json`. Compare them with the
   arm's own gate record before trusting anything downstream.

## Footguns

- **The venv launcher stub shows 0 CPU.** `E:\rmaps_venv2\Scripts\python.exe` is a stub; the
  interpreter that burns CPU is `C:\Users\ambur\miniconda3\envs\codex_py\python.exe`. Filter
  `Win32_Process` by CommandLine or ExecutablePath, never by the venv path, or a healthy job looks
  like a hang.
- **`--workers 1`.** Multiprocessing spawn under this venv can die with
  `PermissionError: [WinError 5] DuplicateHandle`. Parallelise across arms as separate processes.
- **Disk.** One arm's `temp/` reached about 8 GB on a 39,607-event background before deletion. Stage
  on `E:`, never on `F:`: `F:` is an SMR drive and sustained writes collapse it.
- **Temporaries are deleted only after they verify**, and every deleted file is listed with its size
  and md5 in `logs/temp_deletion.log`. `--keep-temp` keeps them. A `fisher` run is never verified by
  this route, so its temporaries are always kept.
- **Scaffold names.** Input contigs are matched to the FASTA index, adding or removing a literal
  `chr` only where that resolves an existing key. An unresolved name fails the run. This is name
  normalisation, not a liftover.
- **Heredoc backslashes.** Inside a Bash-tool heredoc a doubled backslash collapses, so a Windows
  path in a non-raw Python string arrives mangled. Use raw strings or the Write tool.
- **Duplicate target exons.** One row is one rMATS SE event, not one distinct target exon; the same
  exon paired with different flanking exons contributes several rows. Report events and exons
  separately.
- **State the rule.** The arm name is a label only: `QKI_KO_B` on raw input runs as rule A, and a
  pre-split arm needs `--gate-rule` or a gate record with `rule`. Check the `gate_rule` line of
  `versions.txt` before shipping figures.
- **The default footer text is the dissertation gate.** Without `--gate-text` or `--gate-record` the
  footer prints the reli_v121 gate sentence with this run's counts and rule. For any other gate pass
  `--gate-text` or `--gate-record`; otherwise the footer is wrong.
- **Completeness is an inventory.** A partial archive conversion, one missing figure or one empty
  file makes the run `INCOMPLETE`; read `missing` in `run_manifest.json`. The inventory is the set of
  files the writers register, motif scores, positive-control audit, figure audits, provenance sidecars
  and rank workbook included. Before shipping a run that has been copied or moved, re-check it with
  `E:/rmaps_venv2/Scripts/python.exe E:/Claude/rMAPS3/tools/rmaps3_skill_run.py --verify <OUT_DIR>`:
  one lost or changed file gives `INCOMPLETE` (exit 4).
- **Figures are drawn only off the released rank-sum layer** with verified archives. An
  `--engine audited` or `--stat-method fisher` run logs `SKIP figures` with its reason.
- **Detached jobs.** Launch long runs with `Invoke-CimMethod Win32_Process Create`, not
  PowerShell `Start-Process`, or they die when the session restarts.

## What NOT to claim

- **Do not quote either statistic's p as a calibrated p.** `fisher` counts motif hits, not exons, so
  its table margins are not event counts. `mannwhitney` p is scipy's asymptotic normal approximation
  with tie and continuity correction, and with almost every eligible exon carrying no hit it is
  severely anti-conservative. Report the **RBP order**; say so explicitly.
- **No multiple-testing adjustment is applied here, and none should be added to this layer.** Never
  decorate a third-party tool's own output with our statistics; that is what `rmaps3-full`'s
  separate, labelled supplement is for.
- Do not call an enrichment "binding". Do not carry a count, a p or a rank out of this skill into a
  document without re-reading it from the run's own files.
