# Lessons from the September 2026 SE review

- Scope: `motif-map se`, compared with upstream `b9a9dce`. Other event engines have not received the same numerical repairs.
- Branch repairs: `eaeb303` (initial audit), `c34776b` (first review), `3faead9` (second review). These commits are cumulative; the initial repair alone is insufficient.
- Paths and line numbers in the cause column refer to **upstream `b9a9dce`**, retrievable with `git show b9a9dce:<path>`. `SE` abbreviates `legacy/motifMapSE_MP.py` below.
- `tools/` additions accompanying this document are **uncommitted proposals**, not a fourth engine-fix commit. Environment workarounds have no upstream code fix unless explicitly stated.
- Verification entries are reproducible checks, not claims that production runs were repeated for this documentation change.

## Statistical lessons

| Symptom | Cause at upstream b9a9dce | Fix and revision | How to verify |
|---|---|---|---|
| Repeated motif hits inflate Fisher counts; clipped density curves hide the problem. | `SE:529-550,1187-1195` sums hits then subtracts them from event counts, clipping negative cells; `SE:1092-1104` caps density at one. | `eaeb303`: each eligible event/window contributes 0 or 1; Fisher uses hit/no-hit counts with eligible denominators. | `tests/test_engine_synthetic.py::test_planted_enrichment_matches_hand_computed_fisher_table` and `test_cli_every_pvalue_roots_alternative_and_na`; independently reconstruct each 2x2 table from positional NPZ. |
| Long exons lose junction-proximal matches. | `SE:535-538` combines the two end exclusions with OR, requiring a match to satisfy both ends. | `eaeb303`: independent feature-end windows. | `test_known_hit_positions_and_half_open_overlap` plants known hits and checks literal expected windows. |
| Positive/negative windows differ; terminal positions have artificial zeros. | `SE:554-556` and analogous loops at `569-571,583-585,599-601` omit range endpoints. | `eaeb303` establishes half-open overlap; `c34776b` restricts starts to `range(0,L-W+1,step)`. Default L=250/W=50/step=1 gives 201 intron windows; L=50 gives one exon window. | `test_default_region_complete_window_counts`, `test_cli_sparse_schema_stream_and_parallel_agreement`, and `test_cli_exon_window_forwarding`. The singleton exon has a test value but no line segment. |
| Small regional minima look like valid region-level p-values. | `SE:1274-1275,1328-1345` selects the smallest positional p-value without calibrating the selection. | **Not fixed in engine commits.** Uncommitted `tools/summarize_rmaps_regions.py` preserves v4 Westfall–Young min-P label-permutation calibration, including pooled regions. | Compare fixed-seed calibration against enumerated assignments on synthetic matrices; inspect `native_p`, `calib_p`, and `calib_p_pooled` separately. `--perms 20` is a smoke test, not a production precision setting. |
| An apparently significant raw minimum remains selection-biased after BH. | `SE:1335-1345` writes raw minima; upstream implements neither selection calibration nor a multiple-testing correction there. | **No engine correction.** Uncommitted summarizer reports descriptive `native_q` across motif × eight regions × direction and calibrated `calib_q` across motif × three plotted pools × direction. BH on raw minima is not a remedy for within-region selection. | Inspect summarizer BH families and workbook columns; flanking-exon pools are outside the calibrated plotting family. Condensed best-motif rows retain motif-family q-values and are not an additional RBP-level test. |
| Statistical failures are indistinguishable from null results. | `SE:1196-1200` replaces exceptions/nonfinite results with p=1. | `eaeb303`: raise errors or emit explicit NA with a reason. | `tests/test_audit_shared.py` and `test_cli_every_pvalue_roots_alternative_and_na`; preserve NA and its reason downstream. |
| Native rMATS selection does not implement the lab universe. | `SE:351-367` omits flanking exons from identity and uses its own PSI/background gate; `SE:370-417` can substitute a fallback background. `bin/getExonSets.py:60-76` emits ten columns and applies obsolete expression-column assumptions. | **Native classifiers remain unchanged.** Uncommitted `tools/build_event_sets.py` supplies eight-column Rule-A inputs with explicit gates and frozen-logic attribution. Rule B requires already-built pre-split sets. | Compare accepted event identities against the gate ledger; include equality-boundary, missing-expression and duplicate-context cases. Never infer the native classifier is repaired from wrapper success. |

## Coordinates and sequence

| Symptom | Cause at upstream b9a9dce | Fix and revision | How to verify |
|---|---|---|---|
| Short introns include exon sequence; chromosome-edge events remain in every denominator. | `SE:452-467` extends flanks without neighboring-exon clipping; `rmaps_core/genome_access.py:35-47` hides fetch failures behind Ns. | `eaeb303`: clip at true feature/chromosome bounds, require complete-window eligibility, fail missing chromosomes/fetch errors. | `test_short_intron_and_exon_exclude_ineligible_windows`, `test_fetch_error_names_event`, `test_chromosome_edge_and_minus_padding`. |
| Minus-strand padding/labels can be displaced. | `rmaps_core/genome_access.py:38-42` reverse-complements before padding; `SE:460-474` swaps genomic flanks into transcript orientation. The swap itself is required, not a defect. | `eaeb303`: pad before orientation conversion; `3faead9`: asymmetric integrated strand checks. | `test_cli_asymmetric_strands_preserve_all_region_assignments`; retain genomic-left/right coordinates in the input on both strands. |
| Scaffold events fetch no sequence although the assembly contains them. | `rmaps_core/genome_access.py:36,46-47` uses exact FASTA keys and previously returned Ns on mismatch. | `eaeb303` makes absence fatal; uncommitted wrapper maps an unmatched name only when adding/removing `chr` gives a key in the supplied `.fai`. No assembly conversion is performed. | Check mapped-name counts, input hashes and unresolved-name refusal. Lab prep evidence below records `chrGL…`/`GL…` normalization with no dropped events. |
| Headerless or malformed coordinates silently alter the event set. | `SE:446-458` discards the first line and interprets positional columns without validation. | `eaeb303`: mandatory eight-column header, valid half-open coordinates and strand, nonempty groups, duplicate checks and explicit cross-set overlap handling. | `test_invalid_coordinate_sets_raise` and `test_duplicate_and_between_set_overlap_validation`. rMATS starts are zero-based; keep ends unchanged. |
| Overlapping/lowercase/RNA-pattern hits disappear or variable-length hits have wrong spans. | `SE:481-483,529-531` uses nonoverlapping case-sensitive matching and the first hit's length. | `eaeb303`: overlapping uppercase-DNA matching, U→T and individual spans. `c34776b`: finite-width, context-independent motif scans restricted to region plus halo. | `test_overlapping_case_rna_and_variable_span_regex`, `test_long_introns_fetch_only_region_plus_motif_padding`, and unsupported-pattern checks. Bare IUPAC ambiguity letters are not expanded. |
| Retained FASTA is mistaken for a full biological feature. | `SE:452-475` formerly exported its full extraction intervals. | `c34776b` crops to requested region plus maximum motif width minus one; `3faead9` documents the changed contract. | Read `docs/MIGRATION_2026-09.md`; compare crop bounds against feature/chromosome bounds. Do not reuse these crops as full intron/exon sequences. |

## I/O and reproducibility

| Symptom | Cause at upstream b9a9dce | Fix and revision | How to verify |
|---|---|---|---|
| Successful runs lose positional evidence. | `rmaps_core/motif_map_core.py:191-193` recursively deletes `temp/` unless retention is requested. | `eaeb303`: SE retains temp by default; `c34776b`: explicit deletion affects registered current-run files only, after hashing; sparse positional matrices remain. | `test_cli_delete_temp_preserves_unowned_and_hashes_summary_inputs` and `test_cli_xlsx_delete_temp_preserves_converted_input_hash`; use `--keep-temp` explicitly for the summarizer workflow. |
| Reruns mix stale motif tables or overwrite unowned files. | `SE:192-205` reuses directories; `SE:1280-1289` globs previous tables; `SE:1335` truncates outputs. | `eaeb303`: fresh-output guard and current-run aggregation; `c34776b`: registered archival; `3faead9`: no-manifest/collision refusal. Uncommitted wrapper always refuses nonempty output roots. | Run `test_cli_delete_temp_preserves_unowned_and_hashes_summary_inputs`; additionally probe no-manifest reuse and unregistered target collisions and verify sentinel hashes. |
| Logs cannot reconstruct effective statistical settings or deleted summary inputs. | `SE:173-181` and `rmaps_core/motif_map_core.py:183-187` pass settings via environment; `SE:210-214` only configures a log. | `eaeb303`: parameter/version/revision/output manifest; `c34776b`: pre-deletion hashes and retained/deleted status; `3faead9`: actual worker PIDs. Uncommitted wrapper adds config SHA256, input MD5 and step provenance. | `test_cli_manifest_matches_explicit_inventory` and `test_cli_optional_motifs_manifest_inventory`; validate every listed retained artifact and expected inventory, not only manifest self-consistency. |
| Old consumers silently misread new intermediates. | `SE:1133-1143` writes four-column countDist with dense lists; `SE:1205-1223` emits positional p-values without the later NA-reason contract. | `c34776b`: schema 2, sparse NPZ and six-column countDist; `3faead9`: reject unsupported sparse schemas in both readers. | `test_sparse_readers_reject_schema_one`; import `rmaps_core.positional_io.load_hits` or `iter_hits`, retaining eligibility. See the migration note before adapting downstream code. |
| Exit zero or plot presence is treated as scientific validation. | `tests/test_motif.py:99-103,116-131,179-190` lacks expected numeric assertions; `SE:1400-1405` can log rendering failures without aborting. | `eaeb303`, `c34776b`, `3faead9`: synthetic numeric, strand, inventory and worker-PID regressions. Plot-failure exit-status defect remains outside these repairs. | Recompute Fisher values and minima, inspect NA reasons and manifests, and verify required PDF/PNG artifacts explicitly. |

## Windows operation

| Symptom | Cause / upstream relation | Fix and revision | How to verify |
|---|---|---|---|
| A venv process appears idle while the job consumes CPU elsewhere. | Environment launcher-stub behavior; no causal upstream source line exists. The running interpreter may be the venv's base executable. | Operational workaround only; no branch commit: inspect the process tree, command line and executable path. | Match the actual interpreter child to the job and inspect its CPU/artifacts; do not conclude a hang from launcher counters. |
| Spawn workers fail with `PermissionError: [WinError 5] DuplicateHandle`. | Environment-dependent handle permissions; `SE:1423-1426` is the upstream spawn site, not proof of an engine arithmetic defect. | `eaeb303` exposes bounded engine workers; summarizer workaround is `--workers 1`, optionally separate processes per arm. `3faead9` records actual engine worker PIDs; it does not repair OS permissions. | Start a small job and inspect exceptions and completed artifacts. Multiworker engine tests establish tested-environment behavior only. |
| Several calibration jobs oversubscribe CPU through BLAS. | Native-library thread pools, outside upstream source; `SE:1421` also formerly requested CPU-count-minus-one workers. | `eaeb303` caps engine workers; uncommitted wrapper sets `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` from `blas_threads` before launching children. Lab setting is four, with one summarizer worker. | Inspect the configured environment and actual process/thread workload. Set separate budgets for engine workers and library threads. |
| PyX/TeX/Ghostscript rendering fails or differs across hosts. | `SE:1407-1416` invokes a Pillow fallback on unavailable fonts/raised draw errors; export failures elsewhere may only be logged. | Fallback predates the audit baseline; no new renderer fix is claimed. Set `RMAPS_FORCE_MOTIF_FALLBACK=1` for the tested simplified plots. | Compare numeric artifacts independently of graphics and require expected PDF/PNG files. Native PyX rendering remains unverified in the cited reviews. |
| Windows paths acquire control characters in generated scripts. | Shell/string escaping, outside upstream source. A non-raw `\r` in a drive path is a carriage return. | Operational workaround only: use `pathlib`, raw strings or forward slashes and config files; no branch commit. | Print escaped representations during a small path validation and require input existence before running. |

## Historical reproduction boundary

| Observation | Cause / upstream location | Fix or boundary | How to verify |
|---|---|---|---|
| March tables are labelled `wilcoxon.ranksum.pVal`, but Mann–Whitney does not reproduce them. | Upstream `SE:1167-1168` implements actual Mann–Whitney; `SE:1187-1195` implements the Fisher calculation that reproduces March. The reproduction report traces the earlier misleading label to historical `computeWilcoxonP` before `96d16c6`. | No relabelling or arithmetic change here. March reproduction uses unfixed `b9a9dce`, exact March pre-split inputs/motifs/genome and Fisher greater. | The cited reproduction report records equality of 2,016 root p-values, 302,400 positional p-values and 378 count-distribution files. These are read report results, not new checks performed for this document. |
| Corrected maps differ from the collaborator output despite the same method name. | Upstream hit-count table (`SE:1187-1195`), windows (`SE:554-556`) and extraction (`SE:452-467`) are part of the historical calculation. | `eaeb303` and `c34776b` intentionally change those numerical contracts. `3faead9` preserves the round-1 numerical outputs in the cited smoke comparison. | Do not claim fixed-engine historical equivalence. Root filenames/order and raw-minimum presentation remain compatible; numbers need not agree. |
| A web-server comparison is assumed to validate every dataset or assembly. | No universal web-server equivalence follows from upstream source or one replay. | No branch fix or blanket equivalence claim. The verified replay is the March MIAT-KD artifact only, with its original pre-strand-fix inputs; the lab's current gates differ. | Verify each reference's event set, strand convention, assembly, motifs and test separately. The lab-local MIAT-OE hg19 web result is illustrative, not an hg38 validation. |

## Evidence read for this document

All drive paths below are **lab-local evidence**, unavailable in a standalone clone:

- `F:\rMAPS\audit\rMAPS3_audit_astra_2026-09-16.md`: initial findings F1–F25.
- `F:\rMAPS\audit\review_fix_sol_2026-09-16.md` and `F:\rMAPS\audit\review2_fix_sol_2026-09-16.md`: first and second review findings.
- `F:\rMAPS\audit\CHANGELOG_fixes_2026-09-16.md`, `F:\rMAPS\audit\CHANGELOG_fix_round1_2026-09-16.md`, `F:\rMAPS\audit\CHANGELOG_fix_round2_2026-09-16.md`: repair evidence and limitations, including the final report's 102 matching numerical files against round 1.
- `F:\rMAPS\runs\repro_march_MIAT_KD_2026-09-16\REPRO_REPORT.md`: exact historical Fisher reproduction and failed Mann–Whitney equivalence.
- `F:\rMAPS\runs\launch_fisher_c34776b.sh`: production launch at `c34776b`; Fisher greater, four workers, forced plotting fallback.
- `F:\rMAPS\runs\inputs_engine_ready\PREP_MANIFEST.json`: scaffold-name mapping, hashes and no unresolved names.
- `F:\rMAPS\scripts\build_event_sets.py`: source attribution and gate-builder CLI; frozen `reli_v121` source MD5 `32be07fd0e6e2edff5f077384f45b352`.
- `F:\rMAPS\CLAUDE.md`, locked conventions and operational footgun bullets: launcher, spawn permissions, thread budgets and path escaping.
- `E:\rmaps_summ\scripts\summarize_rmaps_regions_v4.py`: downstream calibration implementation; exact hypergeometric upper tails are evaluated in log space, preserving the v4 arithmetic.
- Portable contracts: [CLI usage](docs/CLI_USAGE.md), [migration notes](docs/MIGRATION_2026-09.md), and branch history `git log --oneline upstream-base-2026-09-16..HEAD`.

## Addendum (September 17, 2026): ordering versus significance, and the authors' layer

- Show the authors' method first, untouched. The released engine at the upstream base commit is kept as a
  separate layer (raw per-region minimum Fisher p, no adjustment, constant dot size because the tool reports no
  enrichment ratio). Lab additions (exon-level counting, BH, permutation calibration) are separate, labelled layers.
- Multiple-testing adjustment never changes the ORDER of RBPs (BH is monotone) and the Westfall-Young calibration
  barely does (rank correlation 0.98 to 0.99 with the raw-p order). What changes the order is the counting fix:
  motif hits per exon (released) versus exons carrying the motif (audited) agree at rank correlation 0.75 to 0.86 and
  share only 4 to 6.5 of each panel's top 10. The QKI positive control is first in its two control panels under every layer.
- Permutation floor: with 2,000 label permutations the calibrated p cannot go below 1/2001. A two-stage scheme
  (100,000 fresh permutations for tests with stage-1 p <= 0.005) lifts strong signals to q ~ 0.002 at a few
  dozen refined tests per arm; each refined p is a valid Monte Carlo p on its own.
- Rank stability by resampling the changed exons (300 bootstraps of the foreground, background fixed): at about
  100 to 600 events per direction only the top 1 to 3 RBPs per panel are reproducible (top-10 frequency >= 0.8);
  positions 4 to 10 are not. Report order with this stability measure, not with p-values alone.
- Legacy RBP names in the shipped motif table (9G8, BRUNOL4/5/6, HNRPLL, HuR, PTB, SF2-ASF, SRp20, SRp40, SRp55,
  Tra2-beta) resolve through HGNC alias/previous-symbol records (tools/ users: see the alias table format in
  the lab workflow); SF2-ASF, 9G8 and PTB merge into SRSF1, SRSF7 and PTBP1. 9G8 and SRp40 are HGNC-ambiguous.
- VAST-tools human database Hs2 is hg38 on Ensembl v88 (README); gene labels differ from GENCODE v49 for a few
  percent of exons, so a coordinate-based join with a gene requirement loses ~5 % of matches. This is annotation
  vintage, not a coordinate error (no +/-1 bp offset spike; direction agreement 92 to 96 % on reciprocal-unique matches).
