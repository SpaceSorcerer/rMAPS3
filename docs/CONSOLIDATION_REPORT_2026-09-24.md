# rMAPS3 consolidation — 2026-09-24

**Result:** there is one lab branch (`lab/miat-qki` @ `fe7d32b`, pushed) plus `main` (upstream mirror `b9a9dce`). The worktrees are the main checkout and `rMAPS3_upstream`. Tag `lab-v1.0-2026-09-24` = `920c36c` was pushed and was not moved for the later tie fix. The suite has 223 passed and 1 skipped (a pandas-only test); the baseline was 183. `git status` is clean, and `tools/` has no drive paths.

| Commit | What |
|---|---|
| `fde272a` | Adds the S1 overwrite tests from `rMAPS3_fix`, byte-identical (md5 0eb3b307). 17 passed at HEAD. |
| `87b8230` | Adds `tools/unit_sensitivity.py`, `length_matched.py` and `motif_scores_ranksum.py`, plus 22 synthetic tests. On lab data, 9 deterministic TSVs are byte-identical to the published files. Treatment B, treatment C and length-matched are byte-identical to the originals at 200/1000 permutations. Motif scores are byte-identical on QKI_KO_B. |
| `425d193` | Removes site paths from `tools/`. `rank_stability.py` and the v5 summarizer read this checkout's engine (unchanged since `3faead9`). `rmaps3_skill_run.py` requires `--genome-root`, `--engine-root`, and `--gtf` plus the lists when figures are drawn; its figure label follows `FIG_VERSION` 4.3.2. `build_event_sets_lab_full.py` takes its roots as arguments. A QKI_KO_B `rank_stability.py` rerun exited 0 with byte-identical TSVs and npz arrays. |
| `59a20ff` | Adds the `tools/README.md` index and SUPERSEDED headers to `calibrate_ranksum`, both Fisher summarizers and `rmaps3_lab_run`. |
| `920c36c` | Rewrites `LESSONS.md` as one document with 5 topic sections plus provenance, and rewrites the README lab section. |
| `fe7d32b` | Applies the 4e5d717 tie fix to `build_event_sets_lab_full.py` (orientation `agreement` excludes dPSI == 0 ties) and adds `tests/test_orientation_ties_lab_full.py`. It gives 3 passed under `codex_py` with pandas; against the pre-fix file the tie test fails. It skips under `rmaps_venv`, which has no pandas. The portable builder is the only other `tools/` builder with this guard, and it was already fixed. **The frozen dissertation sets would not change:** selection never used this ratio, and `informative_agreement` was 1.0 in all 3 arms. A rebuild would only move the recorded `agreement` (0.961 / 0.981 / 0.980, `E:\rmaps_work\logs\<ARM>\orientation.json`) to 1.0. Not rebuilt. |

**Git steps.** Both `merge-base --is-ancestor` checks exited 0.
- `tag -a audited-engine-2026-09-16 c34776b` succeeded. `rMAPS3_run` was clean, including ignored files, so `worktree remove` succeeded.
- `branch -d fix/audit-2026-09-16` deleted the branch (was `3faead9`). `push origin lab/miat-qki --tags` updated `e7e623d..920c36c` and pushed the two new tags. `upstream-base-2026-09-16` already existed at `b9a9dce`.

**Deviation: `rMAPS3_fix` was not `worktree remove`d.** It held the untracked S1 test, identical to `fde272a`. It also held ignored `_lab/` (396 MB of remediation records), `.agent-handoff/` and `.venv/`, which removal would have deleted. The directory was instead renamed on the same volume to `E:\Claude\_archive\rMAPS3_fix_worktree_2026-09-24\` (832 MB), then `git worktree prune` removed the stale registration. Nothing was deleted.

**Files moved on F:.** Three files went to `F:\rMAPS\scripts\_superseded_by_fork_2026-09-24\` with md5 unchanged: `build_rmaps_region_lollipops_v4.py` (FIG 4.1), `_v42.py` and `_v43.py`. The other 11 entries stay. `README_POINTERS.md` names each file's fork counterpart and says why it stays; `build_event_sets.py` is the lab copy that built the event sets. `E:\rmaps_*\scripts\` was untouched.

**Skills.** Both SKILL.md files now pass `--engine-root` and name 4.3.2 and the tool index; quick also passes `--gtf` and the lists. Every concrete path resolves, and both commands pass `validate` and `require_site_paths`.

**Skipped or open.**
- `origin/fix/audit-2026-09-16` still exists on the remote; deleting it is Brian's call. The spaceflight `RMAPS3_jc10` callers may have relied on the removed wrapper defaults (unverified).
- `build_event_sets_lab_full.py` was not re-run, because it writes into `E:\rmaps_work`. Only its root binding and the synthetic tie tests were checked. The frozen RBP-RELI `build_reli_foregrounds_v2.py:510` has the same tie formula; that shared engine was not touched.

**file_map entries.**
- `F:\rMAPS\scripts\README_POINTERS.md` maps each out-of-repo rMAPS script to its fork tool and commit. `F:\rMAPS\scripts\_superseded_by_fork_2026-09-24\` holds the three superseded figure builders, byte-unchanged.
- `E:\Claude\_archive\rMAPS3_fix_worktree_2026-09-24\` is the former fix worktree at `3faead9` with its `_lab/` records. `E:\Claude\rMAPS3\_lab\CONSOLIDATION_REPORT_2026-09-24.md` is this report. Update `F:\rMAPS\file_map.md` lines 12 and 71 to `fe7d32b`, 223 passed + 1 pandas-only skip.
