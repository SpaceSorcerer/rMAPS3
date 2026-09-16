# September 2026 SE migration

<!-- AUDIT R9: intermediate-format changes versus upstream b9a9dce. -->
Scope: `motif-map se`. Other event engines retain their earlier contracts.
`run_manifest.json` now declares `schema_version: 2`.

| Change from upstream b9a9dce | Consumer action |
|---|---|
| Temporary output retained by default | Use `--delete-temp` for explicit cleanup; only current registered temp files are deleted after hashing. |
| Non-empty output refused | Use a fresh directory or `--overwrite`; previous registered retained outputs plus manifest move to `_previous_<UTC timestamp>/`. Unowned files are preserved. |
| Default worker cap of four | Set `--workers N` explicitly when another allocation is intended; default is `max(1,min(4,cpu_count-1))`. |
| Strict coordinate replay inputs | Supply a header and exactly eight tab-separated columns, zero-based half-open coordinates, valid strands and chromosomes; duplicates are invalid and cross-set overlaps require `--allow-overlap`. |
| Complete windows only | Starts are `range(0,L-W+1,step)`: defaults give 201 intron windows and one exon window. `--exon-window` defaults to `--window` and permits an explicit smaller exon window. A single exon window is tested numerically but produces no curve segment in the current line renderer, which requires at least two points. |
| Eligibility-aware binary event observations | Each eligible event contributes 0 or 1 for any overlapping motif; missing windows are excluded from denominators. |
| Six-column countDist schema | Columns are `Region, position, sum, eligible, density, values`; values is the compact histogram `0:n,1:n,NA:n`, replacing event-length lists. |
| Per-motif pVal reason column | Parse headers and retain the `reason` for NA comparisons; do not replace NA with p=1. |
| Sparse positional NPZ | Dense `.hits.tsv` is replaced by `.hits.npz` with int32 event indexes, int16 half-open hit spans, eligibility bounds, labels, and window metadata. Use `load_hits` or `iter_hits`; field details are in CLI_USAGE.md. |
| Bounded motif scans | Regexes need positive minimum match width, finite maximum lengths, and no context-dependent assertions; requested regions include a maximum-hit-length-minus-one halo. Unsupported patterns fail explicitly. |
| Shared worker sequence inputs | `temp/sequence_region_0.npy` through `sequence_region_7.npy` and `temp/sequence_metadata.npz` provide read-only worker inputs and follow registered-temp cleanup. |
| Manifest provenance and cleanup status | All registered outputs, including deleted summary inputs, retain SHA256 and size; `status` records retained/deleted. Do not assume every manifest path still exists. |

Root `pVal.up.vs.bg.RNAmap.txt` and `pVal.dn.vs.bg.RNAmap.txt` filenames and
nine-column order are unchanged. Regional minima are raw and uncalibrated; no
region-minimum calibration is introduced by these repairs.
