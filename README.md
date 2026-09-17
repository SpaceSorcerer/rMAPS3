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

## Running on rMATS-turbo 4.3.0 output (pre-split mode) — lab workflow

Pre-split SE inputs make foreground/background selection explicit and reproducible.
The native rMATS classifier and standalone `exon-sets se` helper do not implement
the lab's coverage/expression gates; the latter also assumes obsolete expression
columns and emits a different coordinate schema. Use
[`tools/build_event_sets.py`](tools/build_event_sets.py) through the wrapper below.
Its portable Rule-A builder adapts the frozen `reli_v121` selection logic, with
source attribution retained in the file. Rule-B rMATS/VAST concordant selection
requires externally prepared pre-split sets; it is not rebuilt by this wrapper.

Each of `up`, `dn` and `bg` is a nonempty tab-separated file with this mandatory
eight-column header (the spaces displayed here must be tabs in the actual file):

```text
chr strand exonStart exonEnd firstExonStart firstExonEnd secondExonStart secondExonEnd
```

- Starts are zero-based and ends exclusive. Keep rMATS starts and ends unchanged.
  `firstExon` is genomic-left and `secondExon` genomic-right, on both strands;
  use `+`/`-`. Do not add gene identifiers to the eight-column file.
- Lab foreground: FDR < 0.05, |dPSI| >= 0.10, IJC+SJC >= 10 in **every sample**,
  and DESeq2 baseMean > 50. Absent expression records are `expr_unknown` and
  retained. Apply coverage/expression gates symmetrically to background,
  require background FDR >= 0.5, and remove foreground events.
- Direction follows the configured treatment group and rMATS dPSI sign; record
  group identity explicitly. Human lab references are GRCh38/hg38, GENCODE v49;
  the motif engine reads the specified FASTA, not a GTF.
- The wrapper checks chromosome names against `<root>/<build>/<build>.fa.fai`.
  It can add/remove `chr` only to match an existing key (including scaffold names);
  unresolved names fail. This is name normalization, not a genome-build conversion.

Edit a copy of [`configs/example_miat_qki.json`](configs/example_miat_qki.json)
or [`configs/example_generic.json`](configs/example_generic.json), then run:

```bash
python tools/rmaps3_lab_run.py --config configs/example_miat_qki.json
```

Paths in config files resolve relative to the config's directory. Set either
`inputs.up/dn/bg` or `inputs.rmats_se` plus gates and optional `inputs.deseq2`;
set `genome.root/build`, `motifs.known/additional`, `output_root`, `engine`,
`blas_threads`, and `summary`. Output roots must be absent or empty. The wrapper
writes prepared `inputs/`, `engine/`, `summary/`, and a `run_manifest.json` with
config hash, revision, versions, input MD5s and step records. The `.yaml` examples
contain the same settings but require an already installed PyYAML; otherwise YAML
is rejected clearly and the JSON companions work without it. The wrapper adds
no package installation step.

The corresponding direct engine invocation from the repository root is:

```bash
python cli.py motif-map se --known-motifs data/knownMotifs.human.mouse.txt --motifs data/ESRP.like.motif.txt --fasta-root genomedata --genome hg38 --rMATS NA --miso NA --up inputs/up.coord.txt --down inputs/dn.coord.txt --background inputs/bg.coord.txt --output results/MIAT_KD/engine --stat-method fisher --fisher-alternative greater --workers 4 --intron 250 --exon 50 --window 50 --step 1 --keep-temp
```

Set `RMAPS_FORCE_MOTIF_FALLBACK=1` when using the tested simplified Pillow plots.
The wrapper configures BLAS thread caps before launching children; the lab example
uses four threads and one summarizer worker. On Windows a venv launcher can spawn
a different base interpreter, so inspect the active child process before diagnosing
an idle job. Use one summarizer worker if spawn fails with `DuplicateHandle`.

Keep `temp/` positional p-values and count distributions, plus
`positional/*.hits.npz`, `exon/` coordinate files and the engine manifest.
The NPZ schema is version 2 and readers must preserve eligibility masks;
`rmaps_core.positional_io` provides the readers. Root `pVal.*.RNAmap.txt` tables
contain raw regional minima, not calibrated region p-values. The wrapper runs
the copied v4 summarizer automatically; a separate invocation is:

```bash
python tools/summarize_rmaps_regions.py --run results/MIAT_KD/engine --out results/MIAT_KD/summary --arm MIAT_KD --perms 2000 --seed 149 --workers 1
```

The summarizer preserves v4's exact log-space hypergeometric-tail arithmetic and
Westfall–Young min-P label permutations. It reports the native layer alongside
calibrated pooled regions (upstream intron, exon body, downstream intron), with BH
across motif × plotted pool × direction. BH on native raw minima is descriptive
and does not correct their within-region selection. Flanking-exon results are
retained outside that calibrated plotting family. Outputs include
`per_motif_regions.tsv`, `condensed_per_rbp.tsv`, `positions_long.tsv` and the
arm's summary workbook. The condensed best-motif rows retain motif-level q-values;
they are not a separate RBP-level hypothesis test. This calibration supports
one-sided greater Fisher runs; the wrapper rejects other methods before execution.
Use `--perms 20` on the wrapper only for a synthetic smoke test.

Read [the migration note](docs/MIGRATION_2026-09.md) for complete-window geometry,
changed FASTA crops, temporary-output retention and schemas; read
[LESSONS.md](LESSONS.md) for source/commit evidence and historical reproduction
limits. The corrected engine is not numerically equivalent to the old web-server
calculation merely because both report Fisher p-values.
