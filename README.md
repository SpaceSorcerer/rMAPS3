# rMAPS3

![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue)
![Tests](https://github.com/abdelrahmanvxn/rMAPS-3/workflows/CI/badge.svg)
![Platforms](https://img.shields.io/badge/platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey)

RNA Map Analysis and Plotting Server (rMAPS) generates RNA maps for RNA-binding protein (RBP) analysis around alternative splicing events.

Original rMAPS website: http://rmaps.cecsresearch.org/

## What This Project Provides

- Original rMAPS motif-map and CLIP-map workflows updated for Python 3.
- Unified CLI in `cli.py` for motif maps, CLIP maps, MISO conversion, and exon-set generation.
- Local Web UI launcher in `run_web.py`.
- Selectable p-value methods in both CLI and Web UI: `fisher`, `mannwhitney_greater`, `brunnermunzel_greater`, and `permutation_one_sided`.

## Screenshots

### Web UI

Motif mode (collapsed optional parameters):

![Web UI Motif Map](docs/images/WEBUI_MotifMap.png)

CLIP mode (collapsed optional parameters):

![Web UI CLIP-seq](docs/images/WEBUI_ClipSeq.png)

Optional analysis parameters (expanded):

![Web UI Optional Parameters Expanded](docs/images/WEBUI_OPTIONAL.png)

### CLI

Global help:

![CLI Global Help](docs/images/CLI_HELP.png)

Motif event-specific help:

![CLI Motif Map SE Help](docs/images/CLI_SE.png)

## Requirements

- Python 3.10+
- Perl on `PATH` (for MISO converters)
- Flask (for Web UI)
- Optional: TeX distribution (PyX text rendering)
- Optional: Ghostscript (PNG export)

Install dependencies:

```bash
python -m pip install -r requirements.txt
```

## Documentation

The original rMAPS web help content has been preserved and adapted in the user guide, FAQ, and RBP motif table below.

Original web help:

- User guide and input formats: [`docs/USER_GUIDE.md`](docs/USER_GUIDE.md)
- Curated RBP motif table: [`docs/RBP_MOTIFS.md`](docs/RBP_MOTIFS.md)
- FAQ: [`docs/FAQ.md`](docs/FAQ.md)

rMAPS3 docs:

- Installation and genome setup: [`docs/INSTALL.md`](docs/INSTALL.md)
- Full CLI reference and event-type details: [`docs/CLI_USAGE.md`](docs/CLI_USAGE.md)
- Web UI usage and API: [`webui/README.md`](webui/README.md)
- Testing guide: [`tests/README.md`](tests/README.md)
- Legacy test script notes: [`tests/legacy/README.md`](tests/legacy/README.md)

## Genome Data Layout

Genome FASTA files must be downloaded separately (not included in repo).

Use provided fetch scripts:
- **Windows:** `.\scripts\fetch_genomes.ps1 -Genomes hg19,hg38`
- **Linux/macOS:** `./scripts/fetch_genomes.sh --genomes hg19,hg38`

By default, files download to `genomedata/`. To use a different location, pass `--FastaRoot` (Windows) or `--fasta-root` (Linux/macOS) to the script, or set `$env:RMAPS_FASTA_ROOT` / `RMAPS_FASTA_ROOT` environment variable.

Expected layout after download:
```text
genomedata/
  hg19/hg19.fa
  hg38/hg38.fa
  mm10/mm10.fa
  ...
```

For full setup instructions, see [Installation Guide](docs/INSTALL.md).

## Project Structure

- `cli.py`: CLI entrypoint
- `run_web.py`: local web launcher
- `rmaps_core/`: shared dispatch and utility modules
- `legacy/`: event-specific engines invoked by the CLI
- `bin/`: helper scripts (`miso2rMATS.*.pl`, `RNA.map.noWiggle.py`, `getExonSets.py`)
- `data/`: motif reference inputs
- `data/test/`: bundled sample data
- `tests/`: Python test suites and legacy shell wrappers

## Quick Start

Show CLI help:

```bash
python cli.py --help
```

Run smoke checks:

```bash
python tests/smoke_cli.py
```

Run integration suites:

```bash
python tests/test_clip.py
python tests/test_motif.py --fasta-root genomedata --genome hg19
```

For full command reference (all event types and options), see [`docs/CLI_USAGE.md`](docs/CLI_USAGE.md).

## Web UI (Local)

```bash
python run_web.py
```

Open `http://127.0.0.1:5000`.

## Troubleshooting

See [`docs/FAQ.md`](docs/FAQ.md) for common setup and runtime issues.

## Lab fork: branch `lab/miat-qki`

This branch adds a lab layer on top of the collaborators' engine. `main` mirrors upstream
`b9a9dce` (tag `upstream-base-2026-09-16`). The SE engine repairs are `eaeb303`, `c34776b`
(tag `audited-engine-2026-09-16`) and `3faead9`. Everything else lives in `tools/`, `tests/`,
`data/`, `configs/` and `docs/`. Release tag: `lab-v1.0-2026-09-24`.
[LESSONS.md](LESSONS.md) records every decision below with its date and evidence. Scope is human /
GRCh38 (hg38) / GENCODE v49 unless the inputs say otherwise, and SE events only.

### Layer stack

| Layer | What it is | Report it as | Produced by |
|---|---|---|---|
| Main | The authors' released engine (`b9a9dce`) run with `--stat-method mannwhitney`: a one-sided rank-sum on per-exon motif hit counts, reduced to the smallest p over the 50-nt windows of a region | RBP ORDER only. Its p is scipy's tie-corrected normal approximation on >99 % zeros and is severely anti-conservative | the engine, then `tools/countdist_to_npz.py` and `tools/verify_ranksum_archives.py` |
| Supplement | The same statistic with a Westfall–Young label-permutation p over target-exon clusters (two stages, seed 149), BH q over unique k-mers, and RBP-level min-P over each RBP's motifs | the p and q to report. Every calibrated p is a valid permutation p: stage-1 B = 2,000 (resolution 1/2,001), or stage-2 B = 100,000 (resolution 1/100,001), which a test reports only when its own stage-1 p ≤ 0.005; `calib_perms_used` gives each p's B. BH over them controls the FDR under PRDS-type dependence; PRDS is assumed, not proven ([docs/two_stage_validity.md](docs/two_stage_validity.md)) | `tools/calibrate_ranksum_v2.py` (calibration v2.1) |
| Sensitivities | Row-unit calibration (`*_rowunit` columns), row vs target-exon unit, length-matched background, method comparison, count-aware tests | labelled comparisons, never the p | `calibrate_ranksum.py`, `unit_sensitivity.py`, `length_matched.py`, `compare_stat_methods.py`, `count_aware_stats.py` |
| Stability | Foreground-bootstrap rank stability on an audited Fisher run | how many top RBPs are reproducible | `tools/rank_stability.py` |
| Figures | Region lollipops, version 4.3.3: both layers, by-RBP and by-motif, with and without core-spliceosome and broad binders | — | `tools/build_region_lollipops_v4.py` |

Rules that follow from the stack:

- Lab statistics are never painted onto the tool's own output.
- One figure never mixes two engines, two statistics or two permutation units.
- Every count reads "n events (rMATS SE rows) over N target exons".

The audited engine's `mannwhitney` binarizes before ranking, so it is not an alternative main
layer. The audited Fisher layer of 2026-09-16 is superseded. `tools/summarize_rmaps_regions*.py`
and `tools/rmaps3_lab_run.py` remain for reproduction.

### Tools and skills

[`tools/README.md`](tools/README.md) lists every tool with its purpose, layer and status
(CANONICAL / SENSITIVITY / SUPERSEDED). Two lab skills wrap
[`tools/rmaps3_skill_run.py`](tools/rmaps3_skill_run.py):

- **`rmaps3-quick`** runs `--mode quick`. It gives the authors' layer only: the engine run, the verified count archives, `quick_summary.xlsx`, the four main-layer figures and `index.html`. A run that lacks the verified archives, or the figures when `--no-figures` was not given, exits 4 with status `INCOMPLETE` in `run_manifest.json`; `--allow-partial` records it as `complete_partial` and exits 0.
- **`rmaps3-full`** runs `--mode full`. It adds the row-unit sensitivity, the reportable v2.1 supplement, rank stability when an audited run is given, the v4.3.3 figures of both layers and the rank workbook.

The wrapper carries no site paths. Pass the genome root, the released-engine checkout, the
GENCODE GTF and the two exclusion lists explicitly. The skills hold the lab's values.

### Environment

- **Packages.** `requirements-lock.txt` pins every package the full lab chain imports, for Python 3.11. Install it into a fresh venv with `pip install -r requirements-lock.txt`. `requirements.txt` stays the upstream engine's minimal list.
- **Lab interpreters.** `E:\rmaps_venv2` satisfies the lock: Python 3.11.15 from the conda env `codex_py`, `python -m venv`, then `pip install -r requirements-lock.txt`; the lock is its `pip freeze`. Observed 2026-09-24 at v1.0.1b: full suite 291 passed, 0 skipped. `E:\rmaps_venv` lacks pandas, PyYAML and tzdata; the same day it gave 286 passed, 3 skipped (the pandas-only module and two lock-environment release checks). Details: `_lab/ENV_2026-09-24.md` in the lab checkout.
- **Engine pin.** With `--engine released` the wrapper refuses to start unless `--engine-root` is a git checkout at `--engine-commit` (default `b9a9dce`) with no modified or untracked file. `versions.txt` records the verified SHA.
- **Unversioned inputs.** The FASTA and its index, the GTF, the motif tables, the alias table, the exclusion lists, the gate record and the event sets are not shipped. The wrapper writes the sha256 of each one it was given to `versions.txt` and `run_manifest.json`.

### Reproduce one arm end to end

A checkout of the released engine is needed. For example:

```bash
git worktree add ../rMAPS3_upstream upstream-base-2026-09-16
```

Then, from the repository root, with pre-split inputs and their gate record:

```bash
python tools/rmaps3_skill_run.py --mode full --arm QKI_KO_B --out results/QKI_KO_B \
  --up up.coord.txt --dn dn.coord.txt --bg bg.coord.txt --gate-counts counts.json \
  --genome-root genomedata --genome hg38 \
  --engine released --engine-root ../rMAPS3_upstream --stat-method mannwhitney \
  --permutations 2000 --refine-perms 100000 --seed 149 \
  --gtf gencode.v49.primary_assembly.annotation.gtf \
  --spliceosome-list spliceosome_census.txt --broad-binders-list broad_binders.txt \
  --positive-control QKI
```

`--rmats-se SE.MATS.JC.txt --filter gates.json` replaces the three set flags and
`--gate-counts`. The wrapper then pre-splits with `tools/build_event_sets.py`, which implements
rule A only, so raw input needs an arm named `<NAME>_A`. The wrapper refuses raw input for an arm
whose suffix names rule B or Beffect: rule B sets are frozen concordant files, supplied pre-split.
`--gate-rule` states the rule explicitly and must equal the arm suffix. The output root must be
absent or empty.

The same chain, step by step, reads and writes explicit roots:

```bash
python tools/countdist_to_npz.py --arm QKI_KO_B --released-root runs --out-root counts
python tools/verify_ranksum_archives.py --arm QKI_KO_B --released-root runs --counts-root counts
python tools/calibrate_ranksum.py --arm QKI_KO_B --counts-root counts --released-root runs --out-root summary_rowunit --alias-table data/rbp_alias_hgnc_2026-09-17.tsv --permutation-unit row --seed 149
python tools/calibrate_ranksum_v2.py --arm QKI_KO_B --counts-root counts --released-root runs --out-root summary --alias-table data/rbp_alias_hgnc_2026-09-17.tsv --rowunit-root summary_rowunit --seed 149
python tools/build_region_lollipops_v4.py --arms QKI_KO_B --out-root figures --released-root runs --calibrated-root summary --event-sets-root event_sets --alias-table data/rbp_alias_hgnc_2026-09-17.tsv --gtf gencode.v49.primary_assembly.annotation.gtf --spliceosome-list spliceosome_census.txt --broad-binders-list broad_binders.txt
```

Here `runs/QKI_KO_B/` is the released engine's output directory, run with `--keep-temp`.

Checks to read before trusting a run:

- `counts/<ARM>/VERIFY.md` shows that the archives reproduce every root-table value and every per-position value by exact float equality. Any mismatch, missing table or changed countDist file aborts the run and deletes nothing.
- `summary/<ARM>/refinement_report.json` gives the stage record, the BH family sizes, and the events and target exons per set.
- `figures/positive_control_audit.tsv` checks, on a QKI-KO arm, that QKI ranks first in INCLUDED × Upstream Intron and SKIPPED × Downstream Intron.

### Pre-split input format

Each of `up`, `dn` and `bg` is a nonempty tab-separated file with this mandatory eight-column
header. The spaces shown here must be tabs in the file.

```text
chr strand exonStart exonEnd firstExonStart firstExonEnd secondExonStart secondExonEnd
```

- **Coordinates.** Starts are zero-based and ends exclusive; keep rMATS values unchanged. `firstExon` is genomic-left and `secondExon` genomic-right on both strands. Strand is `+`/`-`. Do not add gene identifiers.
- **Lab gate.** The lab uses the RBP-RELI `reli_v121` ledger gate. Its four named parts are rMATS FDR, |dPSI|, IJC+SJC in every sample, and DESeq2 baseMean.
  - Background takes the same coverage and expression gates plus a high-FDR floor, with the foreground removed.
  - Genes absent from the expression table are `expr_unknown` and retained.
  - Read the values from the gate ledger, not from this README.
- **Chromosome names.** Names are checked against `<root>/<build>/<build>.fa.fai`. A literal `chr` is added or removed only when that resolves an existing key; unresolved names fail. This is name normalisation, not liftover.
- **Engine outputs.** Root `pVal.*.RNAmap.txt` tables hold raw regional minima, not calibrated region p-values. Keep `temp/` until the archives verify. The wrapper deletes exactly the files listed in the verifier's `verified_temporaries.tsv`, only after the verifier exits 0, and only when each still matches its listed md5. It logs each one in `logs/temp_deletion.log`.

Read [the migration note](docs/MIGRATION_2026-09.md) for complete-window geometry, FASTA
crops and the NPZ schema. The repaired engine is not numerically equivalent to the old
web-server calculation merely because both report Fisher p-values.
