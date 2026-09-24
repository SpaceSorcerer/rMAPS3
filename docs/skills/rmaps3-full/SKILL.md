---
name: rmaps3-full
description: "The complete rMAPS3 SE deliverable for one arm: the authors' released layer, the permutation-calibrated supplement that carries the p and q to report (target-exon cluster permutation, RBP-level min-P), a row-unit sensitivity, foreground-bootstrap rank stability where an audited run exists, the version 4.3.3 region lollipops, a rank workbook and an index page. Invoke for 'full rMAPS run', 'rMAPS with figures', 'calibrated rMAPS', 'rMAPS lollipops for <arm>', 'build the rMAPS deliverable', or 'rMAPS supplement'. For the tool's own answer alone, use rmaps3-quick. Human / GRCh38 (hg38) / GENCODE v49 by default; mouse must be named explicitly."
---

# rMAPS3 full run — both layers, figures, workbooks

Same wrapper as `rmaps3-quick`, `--mode full`. It does everything quick does, then adds the lab
layers, each labelled as ours and kept separate from the authors' layer. Read `rmaps3-quick`
first: every input, default, output and footgun there also applies here. Depth:
`E:\Claude\rMAPS3\LESSONS.md` (sections 2 and 3 for the calibration v2.1 decisions); tool index
`E:\Claude\rMAPS3\tools\README.md`; release tag `lab-v1.0-2026-09-24`.

**Layer order is a hard rule.** Layer 1 is the collaborators' released engine, raw p, no BH, no
calibration: RBP ORDER only. Layer 2 is ours, says so on the figure, in the workbook and in the
sidecar, and is **the p and q to report**. Never decorate a third-party tool's own output with our
statistics.

## When to use

| Use this skill | Use something else |
|---|---|
| A committee, paper or PI deliverable for one arm | Just the tool's own numbers → `rmaps3-quick` |
| A reportable calibrated p/q and the lollipops | RBP-RELI enrichment figures → `rbp-reli-lollipops` |
| Rank stability under foreground resampling | Filtering or merging splicing calls → `splicing-tools` |

## Exact command

```bash
E:/rmaps_venv2/Scripts/python.exe E:/Claude/rMAPS3/tools/rmaps3_skill_run.py \
  --mode full --arm <ARM> --out <OUT_DIR> \
  --up <up.coord.txt> --dn <dn.coord.txt> --bg <bg.coord.txt> --gate-counts <counts.json> --gate-rule <A|B|Beffect> \
  --genome-root E:/references/rmaps_genomes --genome hg38 \
  --engine released --engine-root E:/Claude/rMAPS3_upstream --stat-method mannwhitney \
  --permutations 2000 --refine-perms 100000 --seed 149 \
  --gtf E:/references/gencode_v49/gencode.v49.primary_assembly.annotation.gtf \
  --spliceosome-list "F:/RNA-SEQ-ANALYSIS/RBP-RELI/RBP-RELI/reli/data/spliceosome_census.txt" \
  --broad-binders-list "F:/RNA-SEQ-ANALYSIS/RBP-RELI/RBP-RELI/reli/data/broad_binders.txt" \
  --positive-control QKI --arm-label "<title>"
```

`--arm-label` is the figure and index title and is required whenever figures are drawn; the wrapper
never derives a title from the arm name. For the seven dissertation arms copy the `arm_label` of the
arm's row in `E:\Claude\rMAPS3\data\arm_labels_dissertation.tsv`; for any other arm write a title
that states the event source and rule, e.g. "QKI knockout (rMATS-strict)" for a rule-A arm.
Full mode forwards `--arm-label` and `--positive-control` to the figure builder, which reads no title,
rule, gate record or control from the arm name.

`--mode full` requires `--stat-method mannwhitney` and `--engine released` and refuses before any
work otherwise. `--rmats-se <SE.MATS.JC.txt> --filter <gates.json>` replaces the three set flags
exactly as in `rmaps3-quick`. Pre-split inputs need `--gate-counts`, because the figure footer
prints the gate's expression-unknown counts. Run it with `E:\rmaps_venv2`, the interpreter that
satisfies `requirements-lock.txt`.

The event-set rule comes only from `--gate-rule` or the gate record's `rule` field, never from the arm name; the wrapper prints it in `versions.txt`, the workbook README and the figure footer. Pre-split input without either is refused (exit 2); raw
`--rmats-se` input is rule A.

A run is `complete` (exit 0) only when every artefact its mode, engine and statistic require (`REQUIRED_ARTIFACTS` in the wrapper) exists, is non-empty and has its sha256 in `run_manifest.json`; otherwise it is `INCOMPLETE` (exit 4), or `complete_partial` (exit 0) with `--allow-partial`. A failed positive control gives `complete_with_failed_positive_control` (exit 3), an exception `failed` (exit 1), a refused argument exit 2 before any output exists. In full mode the required set adds both calibration summaries
(`summary/<ARM>/` and `summary_rowunit/<ARM>/`, readouts included), the figure index, the rank workbook,
all 16 supplement and main figures, and the rank-stability workbook when `--rank-stability-run` is given.

## The reportable null, in one paragraph (quote this in methods)

The statistic is the released engine's own one-sided rank-sum on per-event motif hit counts,
reduced to the smallest p over the 50-nt windows of a region. Its null distribution comes from
permuting the changed/background labels over **target exons** (chr, strand, exonStart, exonEnd):
rMATS rows that share a target exon move together and are never deduplicated, and every draw keeps
the observed foreground's cluster-size composition. The minimum over windows is taken inside every
permutation (Westfall–Young min-P), and so is the minimum over an RBP's motifs for the RBP-level p
(min-P: each motif's statistic becomes its own tail count, so motifs weigh equally). Two stages:
2,000 permutations for every test, 100,000 for tests with stage-1 p ≤ 0.005; seed 149. Each test reports its stage-2 p only when its own stage-1 p ≤ `--refine-threshold` (0.005), and its stage-1 p otherwise, so every calibrated p is super-uniform at every α and BH over them controls the FDR under PRDS (`docs/two_stage_validity.md`). BH runs over
unique k-mers × 2 directions × 3 pooled regions at motif level, and over RBPs × 2 × 3 at RBP level.
The null asks whether changed target exons are an exchangeable subset of the tested universe given
their motif counts; exon covariates that track motif content, such as length, enter it.

## Inputs and locked defaults, beyond the quick set

| Flag | Default | Note |
|---|---|---|
| `--calib-unit` | `cluster` | the only reportable unit; the row unit always runs beside it as a sensitivity |
| `--permutations` / `--refine-perms` / `--refine-threshold` | 2000 / 100000 / 0.005 | two-stage, self-triggered: a test takes its stage-2 p only on its own stage-1 p; `calib_stage` / `rbp_calib_stage*` record which |
| `--gate-rule` | none; required for pre-split input unless the gate record states `rule` | the event-set rule (`A`, `B`, `Beffect`), passed to the figure builder |
| `--seed` | 149 | calibration and bootstrap |
| `--bootstraps` / `--rank-stability-run` | 300 / none | needs the audited engine's `positional/*.hits.npz`; otherwise a logged skip |
| `--gtf` | GENCODE v49 primary assembly | gene-name resolution only |
| `--spliceosome-list` / `--broad-binders-list` | RBP-RELI `reli/data` lists | drive `_noSpliceosome_noBroad`; SRSF1 always retained |
| `--underpowered-arms` | none | named arms get the underpowered label on tables and figures, verbatim |
| `--gate-text` / `--gate-record` / `--direction-text` / `--tail-text` | dissertation `reli_v121` text | figure footer and legend text; use them for any other gate |
| `--top-n` / `--method-comparison` | 10 / none | rows per panel; optional method sheet in the rank workbook |

Gates when the wrapper pre-splits: identical to `rmaps3-quick`. Read the live values off
`F:\rMAPS\event_sets\<ARM>\` and the `reli_v121` builder, never from memory.

## Species and assembly

Human / GRCh38 (hg38) / GENCODE v49; FASTA `E:\references\rmaps_genomes\hg38\hg38.fa`. A `mm*`
build or mouse `--species` prints a MOUSE REFERENCE IN PLAY banner. Never merge arms across species.

## Outputs and how to read them

Everything `rmaps3-quick` writes, plus:

| Path | What it is |
|---|---|
| `summary/<ARM>/<ARM>_calibrated_ranksum_v2.xlsx` | the reportable supplement: README (with the null), `condensed`, `per_motif`, `rbp_level`, `positions` |
| `summary/<ARM>/rbp_level.tsv` | RBP-level min-P p/q (primary), max-z and mean-z (sensitivity), row-unit columns |
| `summary/<ARM>/readout.md`, `refinement_report.json` | readout; stage record, BH family sizes, events and target exons per set, both-way exons |
| `summary_rowunit/<ARM>/` | row-unit sensitivity (v1); never reported |
| `figures/<ARM>/*.png` + `.svg` + `_figure_provenance_v43.md` | v4.3.3 lollipops: 2 layers × byRBP/byMotif × main/`_noSpliceosome_noBroad` |
| `figures/index.html`, `figures/<ARM>/<ARM>_rank_comparison.xlsx` | figure index; per-panel ranks under every layer |
| `stability/<ARM>_rank_stability.xlsx` | only with an audited run |

Figure encoding: stem = −log10 of that layer's p (by-RBP supplement: the RBP-level min-P p);
colour = direction hue (included gold `#E69F00`, skipped blue `#0072B2`) deepening with raw p on
the main layer and calibrated BH q on the supplement, grey `#999999` at ≥ 0.05; dot size = the tool's own
motif-score ratio on both layers. The supplement legend ends "these are the p and q to report";
the main legend says "use for RBP ORDER only". Every count reads "n events (rMATS SE rows) over N
target exons".

## Verification

1. Everything in `rmaps3-quick` §Verification, including the blocking archive check.
2. **The calibration re-asserts the native layer.** It aborts if a recomputed sub-region minimum
   disagrees with the released root value, and verifies that motif keys sharing a k-mer carry
   bit-identical archives before testing them once.
3. **Positive control.** With `--positive-control QKI` on a QKI-KO arm, QKI first in INCLUDED × Upstream Intron and SKIPPED ×
   Downstream Intron on both layers: wrapper sheet plus `figures/positive_control_audit.tsv`.
4. The ported calibration reproduced the 2026-09-22 lab run byte for byte on QKI_KO_B (all four
   TSVs); re-check against `refinement_report.json` whenever the engine or archives change.

## Footguns

All of `rmaps3-quick` §Footguns, and:

- **Wall time.** The row-unit sensitivity plus the v2 stage 2 are the long poles; budget from the
  arm's own `refinement_report.json`.
- **Two stages, one contract.** Each test reports its stage-2 p only when its own stage-1 p ≤ `--refine-threshold` (0.005), and its stage-1 p otherwise, so every calibrated p is super-uniform at every α and BH over them controls the FDR under PRDS (`docs/two_stage_validity.md`). A promoted neighbour
  never switches a test.
- **A p at the floor** 1/(B+1) means no permutation reached the observed statistic. Quote the rank.
- **Max-z is a sensitivity, not the RBP statistic.** It penalises an RBP whose second motif has a
  heavier-tailed null; min-P does not.
- **Length.** A length-matched background changed no RBP order in the arms tested; a smaller
  background alone loses calls. Report it as a sensitivity, never as the null.
- **The figure builder refuses a partial or v1 calibration**: it needs `exit=0` in the calibration
  `command.log` and the v2 schema.

## What NOT to claim

- Never quote the released layer's p as a p. It is the RBP ORDER claim.
- Never report the row-unit column: duplicate target exons make the rMATS row an invalid
  permutation unit.
- Never build one figure from two engines, two statistics or two permutation units.
- Underpowered arms are landscape only, not for RBP ranking.
- RBP motif enrichment is not binding and not a mechanism. Re-read every rank, count or q from the
  run's own files.
