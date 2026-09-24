# Two-stage calibration: why every calibrated p is a valid permutation p (v1.0.1b, 2026-09-24)

**Rule (self-triggered promotion).** Each reported test reports its stage-2 p only when its own stage-1 p <= t (`--refine-threshold`, default 0.005); otherwise it reports its stage-1 p. The tests are the motif-level test per sub-region and per pooled region, and RBP min-P, max-z and mean-z per pooled region. Stage-2 draws computed because another test of the same pair or RBP was promoted never switch a test. `calib_stage`, `calib_perms_used` and `rbp_calib_stage[_maxz|_meanz]`, `rbp_calib_perms_used[_maxz|_meanz]` say which stage each p comes from. Code: `tools/rmaps_calib_v2_lib.py` `two_stage_p`.

**Each stage p is super-uniform.** p_k = (b_k+1)/(B_k+1) = (b+1)/(B+1) with b_k = #{draws with statistic >= observed}. Under the null the observed labelling and the B_k random cluster draws are exchangeable. The observed statistic's rank among the B_k + 1 values is then uniform, with ties counted against it, so P(p_k <= a) <= a for every a.

**The switched p.** Let p = p2 if p1 <= t, else p1.
- a <= t: {p <= a} = {p1 <= t, p2 <= a}, so P(p <= a) <= P(p2 <= a) <= a.
- a > t: P(p <= a) = P(p1 <= t, p2 <= a) + P(t < p1 <= a) <= P(p1 <= t) + P(t < p1 <= a) = P(p1 <= a) <= a.
- Neither line uses the joint law of p1 and p2. The stage streams are independent (`SeedSequence([seed, stage, direction])`), but p1 and p2 share the observed data, so they are independent only given the data. Under full independence with p1 exactly uniform, the lines reduce to P(p1 <= t) P(p2 <= a) <= t a for a <= t, and t a + (a - t) <= a for a > t.

**Why group triggering was dropped (v1.0.1 rule).** If a test can take its stage-2 p because another test was promoted, then for a > t the event {p <= a} is no longer inside {p1 <= a}. A test with p1 > a can be switched and land below a. `tests/test_two_stage_fdr.py` builds valid partner p-values that promote a test on slices of its own p1 above 0.2. Under the group rule P(p <= 0.2) then exceeds 0.2; under the self-triggered rule it does not.

**FDR, and nothing stronger.** BH over p-values that are super-uniform under their nulls controls the FDR at the nominal level when those p-values are PRDS. The rule above gives super-uniformity. PRDS is assumed, not proven, for these p: shared permutations within an arm and direction, overlapping k-mers and the shared background all induce dependence. The simulation checks FDR on global-null arms with correlated motifs (`tests/test_two_stage_fdr.py`).
