# CLI Usage

This document is the complete CLI reference for rMAPS 3.

## Global Help

```bash
python cli.py --help
python cli.py motif-map --help
python cli.py clip-map --help
python cli.py convert --help
python cli.py exon-sets --help
```

## Version

```bash
python cli.py --version
```

Include the reported version when sharing results, reporting issues, or documenting HPC/module runs.

Use `python cli.py ...` for the most portable invocation. On Linux/macOS, `./cli.py ...` also works when your active environment has the required Python packages.

Path rule: relative input and output paths are resolved from the directory where you run the command. Absolute paths are used as-is, and `NA` still means "no file" for optional inputs. See [`INSTALL.md`](INSTALL.md) for shared-HPC setup notes.

Motif plotting uses PyX when available. On systems without TeX/PyX fonts,
rMAPS3 falls back to simplified Pillow-rendered motif-map PDFs/PNGs. Set
`RMAPS_FORCE_MOTIF_FALLBACK=1` to test that fallback path explicitly.

## Motif Map Commands

### Supported Event Types

- `se`: Skipped Exon
- `a3ss`: Alternative 3' Splice Site
- `a5ss`: Alternative 5' Splice Site
- `ri`: Retained Intron
- `mxe`: Mutually Exclusive Exons

### Event Subcommand Help

```bash
python cli.py motif-map se --help
python cli.py motif-map a3ss --help
python cli.py motif-map a5ss --help
python cli.py motif-map ri --help
python cli.py motif-map mxe --help
```

### Shared Required Options

- `--known-motifs`
- `--fasta-root`
- `--genome`
- `--output`

### Input Modes (choose one)

- rMATS mode: `--rMATS <file>` with `--miso NA --up NA --down NA --background NA`
- MISO mode: `--miso <file>` with `--rMATS NA --up NA --down NA --background NA`
- Coordinate replay mode: `--up <file> --down <file> --background <file>` with `--rMATS NA --miso NA`

### Common Optional Options

- `--motifs` (default `NA`)
- `--label` (default `RBP`)
- `--intron` (default `250`)
- `--exon` (default `50`)
- `--window` (default `50`)
- `--step` (default `1`)
- `--sig-fdr` / `--sigFDR` (default `0.05`)
- `--sig-delta-psi` / `--sigDeltaPSI` (default `0.05`)
- `--stat-method` / `--statMethod` (default `fisher`)
  - Allowed values: `fisher`, `mannwhitney_greater`, `brunnermunzel_greater`, `permutation_one_sided`
- `--stat-permutations` / `--statPermutations` (optional; permutation count for `permutation_one_sided`)
- `--stat-seed` / `--statSeed` (optional; RNG seed for `permutation_one_sided`)
- `--keep-temp` (SE retains positional tables by default; other event types retain the previous opt-in behavior)
- SE only: `--delete-temp` (delete regular files directly inside `output/temp` after success, then remove that directory if empty; preserve `positional/*.hits.tsv`)
- SE only: `--overwrite` (allow an existing non-empty output directory)
- SE only: `--allow-overlap` (allow identical events across up/down/background sets; duplicates within a set remain errors)
- SE only: `--fisher-alternative greater|two-sided` (default `greater`)
- SE only: `--workers N` (positive integer; default `max(1, min(4, cpu_count - 1))`)
- `--separate`

### Statistical Methods

- `fisher` (default): for SE, Fisher exact test on binary events-with-hit versus events-without-hit; see the SE contract below. Other event engines retain their previous counting behavior.
- `mannwhitney_greater`: one-sided Mann-Whitney U test on per-position motif count distributions.
- `brunnermunzel_greater`: one-sided Brunner-Munzel test for stochastic dominance with weaker equal-variance assumptions.
- `permutation_one_sided`: one-sided empirical permutation test on mean differences.

Permutation guidance:
- `permutation_one_sided` is much slower than the other methods.
- Runtime scales with `--stat-permutations`.
- Use `--stat-seed` for reproducibility.
- Practical workflow: develop/debug with `fisher` (or rank-based methods), then run permutation for final robustness checks.

### Event Details and Examples

#### `se` (Skipped Exon)
Use when your input file contains cassette-exon skipping events.

<!-- AUDIT F1, F3, F5, F6, F7, F8, F9, F13, F14, F15, F20: SE numerical and output contract. -->
SE coordinate replay requires a header and exactly eight tab-separated columns:

```text
chr strand exonStart exonEnd firstExonStart firstExonEnd secondExonStart secondExonEnd
```

The displayed spaces must be tabs in the file. Coordinates use zero-based starts and exclusive ends;
`firstExon` is the genomic-left exon and `secondExon` the genomic-right exon on both strands.
Coordinates must be integers, each start must precede its end, and strand must be `+` or `-`.
Each set must be non-empty, chromosomes must exist in the FASTA, and duplicate events are rejected.
Counts and between-set overlaps are reported; overlaps fail unless `--allow-overlap` is explicit.
For human GRCh38/hg38, use `--fasta-root E:\references\rmaps_genomes --genome hg38`, which
loads `E:\references\rmaps_genomes\hg38\hg38.fa`. FASTA keys use the first header token
(`>chr1 1` therefore has key `chr1`). Missing chromosomes or fetch failures are hard errors.

**Windows, eligibility, and matching**

- Every region uses transcript orientation. A window of length `W` starting at region position
  `j` covers `[j, j+W)`; `--step` selects window starts. A motif overlaps when at least one of
  its nucleotides lies inside that interval. Both strands use this same definition.
- Intronic sequence stops at the neighboring exons; genomic sequence stops at chromosome
  boundaries. An event contributes only when it supplies the entire window. A short intron,
  short exon, or chromosome edge can therefore remove it from that window's denominator.
- Each eligible event contributes exactly `1` if any motif hit overlaps the window, otherwise
  `0`. Fisher's table is `[[up_with, up_without], [bg_with, bg_without]]` (likewise down).
  Density is events-with-hit divided by eligible events, hence lies in `[0,1]`.
- Motifs are regular expressions, matched against uppercase DNA after converting `U` to `T`
  in the pattern. Overlapping matches are retained, and every hit retains its own span.
  IUPAC bracket classes such as `[AG]` are already regex-compatible; bare ambiguity letters
  are not automatically expanded. Motif species provenance remains the user's responsibility.
- The engine scans complete true exon/intron features so regex hits can extend beyond the
  requested windows. Very long introns can therefore increase runtime and memory use.
- Statistical errors raise or produce explicit unavailable (`NA`) results with a reason;
  an unusable comparison is never silently replaced by `p=1`.

The eight region labels and junction-relative offsets corresponding to region position `j` are:

| Region | Offset |
|---|---|
| `UpstreamExon_3prime` | `-exon + j` |
| `UpstreamExonIntron` | `j` |
| `UpstreamIntron` | `-intron + j` |
| `TargetExon_5prime` | `j` |
| `TargetExon-3prime` | `-exon + j` |
| `DownstreamIntron` | `j` |
| `DownstreamExonIntron` | `-intron + j` |
| `DownstreamExon_5prime` | `j` |

**Retained output and downstream calibration**

- Per-motif positional p-values and count distributions in `temp/` are preserved by default.
  `--delete-temp` removes regular files directly in `temp/` after success, including older
  files when `--overwrite` is used, and removes the directory if empty. It leaves nested
  directories intact.
- `positional/<RBP>.<motif>.<region>.hits.tsv` contains one row per event, grouped as `up`,
  `dn`, then `bg`, retaining input order within each set. Event identifiers and the `set`
  column identify rows; numeric column headers are region-relative window starts. Values
  are `0`, `1`, or literal `NA` for an ineligible event/window. These matrices are retained
  even with `--delete-temp` and provide the binary observations and missingness required
  by a downstream region-minimum calibrator.
- Root `pVal.{up,dn}.vs.bg.RNAmap.txt` values remain **uncalibrated raw regional minima**.
  The engine does not implement F4 calibration. Applying BH directly to these minima does
  not correct selection across windows; downstream calibration must precede that step.
- Non-empty output directories are refused unless `--overwrite` is supplied. Summaries
  use only results generated in the current run, including when older files remain.
- `run_manifest.json` records relative output paths, SHA256 hashes, byte sizes, effective
  parameters, statistical method, permutations, seed, Python/package versions, and git revision.
  It inventories retained current-run outputs; the manifest excludes its own hash.

```bash
python cli.py motif-map se \
  --known-motifs data/knownMotifs.human.mouse.txt \
  --motifs data/ESRP.like.motif.txt \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/motif_se \
  --rMATS data/test/SE.MATS.ReadsOnTargetAndJunctionCounts.txt \
  --miso NA --up NA --down NA --background NA
```

#### `a3ss` (Alternative 3' Splice Site)
Use for alternative acceptor events.

```bash
python cli.py motif-map a3ss \
  --known-motifs data/testMotifs.txt \
  --motifs NA \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/motif_a3ss \
  --rMATS data/test/clip/A3SS/a3ss.rMATS.txt \
  --miso NA --up NA --down NA --background NA
```

#### `a5ss` (Alternative 5' Splice Site)
Use for alternative donor events.

```bash
python cli.py motif-map a5ss \
  --known-motifs data/testMotifs.txt \
  --motifs NA \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/motif_a5ss \
  --rMATS data/test/clip/A5SS/a5ss.rMATS.txt \
  --miso NA --up NA --down NA --background NA
```

#### `ri` (Retained Intron)
Use for intron-retention events.

```bash
python cli.py motif-map ri \
  --known-motifs data/testMotifs.txt \
  --motifs NA \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/motif_ri \
  --rMATS data/test/clip/RI/ri.rMATS.txt \
  --miso NA --up NA --down NA --background NA
```

#### `mxe` (Mutually Exclusive Exons)
Use for mutually exclusive exon events.

```bash
python cli.py motif-map mxe \
  --known-motifs data/testMotifs.txt \
  --motifs NA \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/motif_mxe \
  --rMATS data/test/clip/MXE/mxe.rMATS.txt \
  --miso NA --up NA --down NA --background NA
```

## CLIP Map Commands

### Supported Event Types

- `se`: Skipped Exon
- `a3ss`: Alternative 3' Splice Site
- `a5ss`: Alternative 5' Splice Site
- `ri`: Retained Intron
- `mxe`: Mutually Exclusive Exons

### Event Subcommand Help

```bash
python cli.py clip-map se --help
python cli.py clip-map a3ss --help
python cli.py clip-map a5ss --help
python cli.py clip-map ri --help
python cli.py clip-map mxe --help
```

### Shared Required Options

- `--peak` / `-p`
- `--output` / `-o`

### Input Modes (choose one)

- rMATS mode: `--rMATS <file>` with `--miso NA --up NA --down NA --background NA`
- MISO mode: `--miso <file>` with `--rMATS NA --up NA --down NA --background NA`
- Coordinate replay mode: `--up <file> --down <file> --background <file>` with `--rMATS NA --miso NA`

### Common Optional Options

- `--label` (default `RBP`)
- `--intron` (default `250`)
- `--exon` (default `50`)
- `--window` (default `10`)
- `--step` (default `1`)
- `--sig-fdr` / `--sigFDR`
- `--sig-delta-psi` / `--sigDeltaPSI`
- `--stat-method` / `--statMethod` (same method set as motif-map; see Statistical Methods above; default `fisher`)
- `--stat-permutations` / `--statPermutations` (optional; permutation count for `permutation_one_sided`)
- `--stat-seed` / `--statSeed` (optional; RNG seed for `permutation_one_sided`)
- `--keep-temp` (keep `output/temp` after success; by default temp is cleaned on success and kept on failures)
- `--separate`

Default thresholds by event:

- `se`, `a3ss`, `ri`, `mxe`: `sigFDR=0.05`, `sigDeltaPSI=0.05`
- `a5ss`: `sigFDR=0.005`, `sigDeltaPSI=0.01`

### Event Details and Examples

#### `se` (Skipped Exon)
Use for cassette-exon CLIP enrichment maps.

```bash
python cli.py clip-map se \
  --peak data/test/clip/PIPE-CLIP.Clusters.bed \
  --output results/clip_se \
  --rMATS data/test/clip/SE/se.rMATS.txt \
  --miso NA --up NA --down NA --background NA \
  --label PTB \
  --sigFDR 0.05 --sigDeltaPSI 0.05
```

#### `a3ss` (Alternative 3' Splice Site)
Use for alternative acceptor CLIP enrichment maps.

```bash
python cli.py clip-map a3ss \
  --peak data/test/clip/PIPE-CLIP.Clusters.bed \
  --output results/clip_a3ss \
  --rMATS data/test/clip/A3SS/a3ss.rMATS.txt \
  --miso NA --up NA --down NA --background NA \
  --label RBFOX2 \
  --sigFDR 0.05 --sigDeltaPSI 0.05
```

#### `a5ss` (Alternative 5' Splice Site)
Use for alternative donor CLIP enrichment maps.

```bash
python cli.py clip-map a5ss \
  --peak data/test/clip/PIPE-CLIP.Clusters.bed \
  --output results/clip_a5ss \
  --rMATS data/test/clip/A5SS/a5ss.rMATS.txt \
  --miso NA --up NA --down NA --background NA \
  --label ESRP \
  --sigFDR 0.005 --sigDeltaPSI 0.01 \
  --separate
```

#### `ri` (Retained Intron)
Use for retained-intron CLIP enrichment maps.

```bash
python cli.py clip-map ri \
  --peak data/test/clip/PIPE-CLIP.Clusters.bed \
  --output results/clip_ri \
  --rMATS data/test/clip/RI/ri.rMATS.txt \
  --miso NA --up NA --down NA --background NA \
  --label MBNL1 \
  --sigFDR 0.05 --sigDeltaPSI 0.05
```

#### `mxe` (Mutually Exclusive Exons)
Use for mutually exclusive exon CLIP enrichment maps.

```bash
python cli.py clip-map mxe \
  --peak data/test/clip/PIPE-CLIP.Clusters.bed \
  --output results/clip_mxe \
  --rMATS data/test/clip/MXE/mxe.rMATS.txt \
  --miso NA --up NA --down NA --background NA \
  --label NOVA1 \
  --sigFDR 0.05 --sigDeltaPSI 0.05
```

## Conversion Commands

```bash
python cli.py convert miso --help
```

Examples:

```bash
python cli.py convert miso --event se --in data/test/ESRP.OE.miso_bf --out temp/se.from_miso.rmats.txt
python cli.py convert miso --event a3ss --in path/to/a3ss.miso_bf --out temp/a3ss.from_miso.rmats.txt
python cli.py convert miso --event a5ss --in path/to/a5ss.miso_bf --out temp/a5ss.from_miso.rmats.txt
python cli.py convert miso --event ri --in path/to/ri.miso_bf --out temp/ri.from_miso.rmats.txt
python cli.py convert miso --event mxe --in path/to/mxe.miso_bf --out temp/mxe.from_miso.rmats.txt
```

## Exon Set Commands

```bash
python cli.py exon-sets se --help
```

Example:

```bash
python cli.py exon-sets se \
  --input path/to/common.txt \
  --sample1 SAMPLE1 \
  --sample2 SAMPLE2 \
  --out temp/exon_sets
```

## Output Notes

Motif-map output typically includes:

- `exon/`
- `fasta/`
- `maps/`
- `temp/`
- `log.motifMap.txt`
- `pVal.up.vs.bg.RNAmap.txt`
- `pVal.dn.vs.bg.RNAmap.txt`

CLIP-map output typically includes:

- `exon/`
- `temp/`
- `*.RNAmap.txt`
- `*.pdf` / `*.eps`
- `log.CLIPSeq*.txt`
