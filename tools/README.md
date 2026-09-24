# tools/ — lab layer of the rMAPS3 fork

Every tool takes its inputs and output roots as arguments. No tool carries a site path.
The lab's own paths live in the `rmaps3-quick` and `rmaps3-full` skills, which call
`rmaps3_skill_run.py`. Scope is human / GRCh38 (hg38) / GENCODE v49 unless the inputs say
otherwise, and SE events only. Background and decisions are in [../LESSONS.md](../LESSONS.md).

Status key:

- **CANONICAL** means the tool produces a reported layer or its inputs.
- **SENSITIVITY** means a labelled comparison that is never reported as the p or q.
- **SUPERSEDED** means the tool is kept for reproduction or as an import only, and its docstring names the replacement.

## Index

| Tool | Purpose | Layer | Status |
|---|---|---|---|
| `rmaps3_skill_run.py` | One command per arm: `--mode quick` runs the engine, verifies the archives and draws the main-layer figures; `--mode full` adds the supplement, the row-unit sensitivity, rank stability and the v4.3.2 figures | wrapper | CANONICAL |
| `build_event_sets.py` | Portable stdlib Rule-A pre-split builder, with eight-column coordinates and a gate ledger | inputs | CANONICAL (portable) |
| `build_event_sets_lab_full.py` | The pandas lab builder behind the lab event sets: rule A, frozen rule B, and `--vast-conf effect-only` for ruleBeffect. Needs `requirements-lab-full.txt` | inputs | CANONICAL (lab) |
| `countdist_to_npz.py` | Packs a released run's `temp/*.countDist.*.txt` into per-motif `*.counts.npz` archives | main layer | CANONICAL |
| `verify_ranksum_archives.py` | Recomputes every released root p and every per-position p from the archives, requires exact float equality, and writes `verified_temporaries.tsv`, the only list the wrapper deletes from | main layer | CANONICAL |
| `rmaps_countdist_io.py` | Shared readers and the released rank-sum kernel | library | CANONICAL |
| `calibrate_ranksum_v2.py` | Calibration v2.1, which gives the p and q to report: target-exon cluster permutation, two stages, RBP-level min-P, and the unique-k-mer family | supplement | CANONICAL |
| `rmaps_calib_v2_lib.py` | Cluster drawer, RBP min-P / max-z / mean-z and k-mer grouping for v2 | library | CANONICAL |
| `build_region_lollipops_v4.py` | Region lollipops, figure version 4.3.2: both layers, by-RBP and by-motif, main and `_noSpliceosome_noBroad` | figures | CANONICAL |
| `rank_stability.py` | Foreground-bootstrap stability of the RBP order on an audited Fisher run (`positional/*.hits.npz`) | stability | CANONICAL |
| `calibrate_ranksum.py` | v1 row-unit calibration. It writes the `*_rowunit` sensitivity columns, and v2 imports its `MotifModel` | supplement | SUPERSEDED by `calibrate_ranksum_v2.py`; SENSITIVITY only |
| `unit_sensitivity.py` | Row versus target-exon unit comparison. A is the row unit, B deduplicates to one row per target exon, C keeps rows and permutes clusters | sensitivity | SENSITIVITY |
| `length_matched.py` | Recalibration against a length-matched background, with a size-matched random control. Writes `lengthmatched_report.json` for the figure builder's `--length-root` | sensitivity | SENSITIVITY |
| `compare_stat_methods.py` | Rank agreement across every method whose root tables exist | comparison | SENSITIVITY |
| `count_aware_stats.py` | Count-aware rank and Poisson-rate window tests from the audited engine's archives | comparison | SENSITIVITY |
| `motif_scores_ranksum.py` | Motif scores, score verification and map links added to v1 summaries on 2026-09-21 | supplement | SUPERSEDED by the `calibrate_ranksum_v2.py` `count_ratio` columns |
| `summarize_rmaps_regions.py` | Portable Fisher regional-minima summarizer for the audited engine, from the 2026-09-16 layer | Fisher layer | SUPERSEDED by `calibrate_ranksum_v2.py` |
| `summarize_rmaps_regions_lab_v5.py` | Lab v5 Fisher summarizer with two-stage refinement. `rank_stability.py` imports its readers | Fisher layer | SUPERSEDED by `calibrate_ranksum_v2.py` |
| `rmaps3_lab_run.py` | Config-driven audited-engine Fisher wrapper. `rmaps3_skill_run.py` imports its coordinate normalisation | wrapper | SUPERSEDED by `rmaps3_skill_run.py` |

## Order for one arm (what `rmaps3_skill_run.py --mode full` runs)

1. Pre-split the inputs with `build_event_sets.py`, or supply `up` / `dn` / `bg` from the lab builder.
2. Run the released engine with `--stat-method mannwhitney`.
3. Run `countdist_to_npz.py`, then `verify_ranksum_archives.py`, which blocks the run on a mismatch.
4. Run `calibrate_ranksum.py --permutation-unit row` for the sensitivity, then `calibrate_ranksum_v2.py --rowunit-root …` for the supplement.
5. Run `rank_stability.py`, only when an audited Fisher run exists.
6. Run `build_region_lollipops_v4.py`.

`unit_sensitivity.py`, `length_matched.py`, `compare_stat_methods.py` and
`count_aware_stats.py` are run by hand when a sensitivity is wanted. Their outputs are labelled
and never replace the supplement.
