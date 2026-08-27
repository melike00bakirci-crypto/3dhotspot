# DECISION-STAGE-D-FIXED-RFP-0001

```
DECISION_ID    = DECISION-STAGE-D-FIXED-RFP-0001
DATE           = 2026-08-18
AUTHORIZED_BY  = user
STATUS         = IMPLEMENTED
SCOPE          = Stage D (geometric footprint robustness) only
```

## Decision

Every Stage D robustness iteration reconstructs the perturbed footprint (`S_iter`) directly
at the FROZEN Stage C `r_fp` via `hotspot3d.footprint.api.reconstruct_at`. The per-iteration
domain re-derivation for selection purposes, multi-scale sweep, QC-driven admissibility, and
Pareto / distance-to-ideal selection are removed from the iteration path. What was previously
called the "secondary fixed-radius reconstruction" is now the sole per-iteration
reconstruction; the previous "primary" reconstruction (which re-derived `r_fp^iter` via the
full frozen methodology on `S_iter`) no longer exists.

This is a **user-authorized scientific design decision**, not a bug fix and not a
QA-initiated repair. `pipeline-qa-debugger` has no authority to make this class of change
(scientific methodology firewall) and did not make it; the user made it directly.

## Rationale (as given by the user)

> The robustness question is "does the selected footprint survive loss of its supporting
> centers?" — not "would Stage C have selected a different radius after each perturbation?"
> Radius selection is Stage C's job, performed once on the full data with the multi-scale
> sweep, ARI neighboring-radius stability and Pareto criteria.

Cost was **not** the driver of this decision — Stage D was already within its 40.0 h budget
(measured 18.8 h, see below) at the time this decision was authorized. The change was made
for scientific correctness: re-deriving a fresh footprint radius inside a robustness
iteration answers a different, unintended question (would Stage C's own selection procedure
be unstable under center loss) rather than the intended one (is the footprint Stage C
actually selected robust to center loss).

## What was measured, in order, this session

1. **Cost analysis of the pre-existing design** (single real KCNA2 iteration,
   `20260818T114125Z_f5c0984d_bb8ccb85`, k=1): the recorded Stage D estimate of 42.1 h was
   built from `_baseline_reproduction`'s wall time alone (one `build_footprint` sweep), which
   undercounted a real iteration's actual cost (primary sweep + fixed-radius reconstruction +
   two comparisons). A direct single-iteration measurement gave 16.696 s/iteration
   (46.4 h projected at N_CAP=10,000) — already over budget, more than the recorded estimate,
   not less.
2. **`bottleneck_radius` performance defect** (`QA-09-BOTTLENECK-PERF-0001`, filed and fixed
   separately, unrelated to this decision): a pure-Python union-find dominated 36% of one
   iteration's cost. Fixed with a provably bit-identical binary-search + `scipy.ndimage.label`
   replacement. Re-measured single-iteration cost: 6.7715 s/iteration (18.8 h projected) — the
   first point at which Stage D was comfortably within its 40.0 h budget.
3. **This decision's premise**: with cost no longer the constraint, the pre-existing design's
   *scientific* shape was reconsidered on its own terms, independent of the performance work,
   and found to be answering the wrong question inside each iteration (see Rationale). See
   step 4 below for the post-decision re-profile.

## What changed

- `src/hotspot3d/robustness/perturb.py` — `run_iteration()` calls only `reconstruct_at`
  (never `build_footprint`) for the perturbed center set. `IterationResult.robustness_r_fp`
  replaces `r_fp_iter`; it is always the frozen Stage C `r_fp` for a completed iteration,
  asserted at its own origin.
- `src/hotspot3d/robustness/stage.py` — `assert_no_radius_optimization_in_perturb()` (new,
  AST-based) is a mandatory structural assertion on every Stage D run, alongside a
  stage-level check that every completed iteration's `robustness_r_fp` equals `phase_c.r_fp`.
  The per-iteration cost estimate that feeds the `max_wall_seconds` budget check is now probed
  by timing one real `reconstruct_at` call (matching what an iteration actually does) instead
  of reusing the (now much more expensive and unrepresentative) baseline-reproduction sweep
  time. `rfp_diagnostics.tsv` and its `RFP_DIAGNOSTIC_COLUMNS` schema are removed (nothing left
  to diagnose: the radius no longer varies). The `_fixed`-suffixed duplicate geometry columns
  in `perturbation_results.tsv` are removed; the (formerly "fixed", now sole) geometry block is
  reported unsuffixed.
- `src/hotspot3d/robustness/sensitivity.py` — `center_recurrence_rate_fixed` (a duplicate of
  `center_recurrence_rate` once there is only one reconstruction) is removed.
- `src/hotspot3d/robustness/report.py` — the `r_fp^iter` diagnostics section of
  `stage_d_report.md` is replaced with a statement of the frozen-radius invariant.
- `src/hotspot3d/orchestration/contracts.py`, `src/hotspot3d/footprint/stage.py` — the
  `handoff_04` required key `rfp_diagnostics` is renamed `robustness_r_fp_invariant`.
- `docs/3D_Hotspot_Footprint_Workflow_v2_2.md` §8 — a `[CLARIFIED]` paragraph states that "the
  entire hotspot analysis should be repeated from the beginning" does not include re-selecting
  `r_fp`, which is frozen after Stage C.
- `.claude/agents/footprint-robustness.md` — §3 items 9–10, §4, §9 item 19, the file-listing
  section, the reporting-protocol paragraph, and §16's prohibited-actions list updated to
  describe fixed-radius-only reconstruction and the new invariant; the moot "`r_fp^iter ≠ r_fp`
  is a diagnostic, never an error" language is replaced by "re-selecting `r_fp` inside a
  robustness iteration" as a prohibited action.
- Tests: `tests/robustness/test_robustness_phase_d.py`,
  `test_robustness_v2_underpowered_investigation.py`, `test_robustness_determinism.py`, and
  `tests/fixtures/synthetic_annotation.py` updated for the new schema/status vocabulary
  (`STATUS_EMPTY_DOMAIN` / `STATUS_NO_ADMISSIBLE` are removed — no longer reachable, since no
  selection is attempted inside an iteration). Two new pins added: a static AST-based check
  (`assert_no_radius_optimization_in_perturb`, invoked on every real run) and a dynamic
  monkeypatch spy proving `build_footprint` is called exactly twice in a real Phase C + Phase D
  run (Stage C's own selection, Phase D's mandatory baseline reproduction) regardless of how
  many perturbation iterations execute.

## What was preserved, unchanged

`N_CAP` (10,000), `K_MAX` rule (`min(n_S-1, floor(n_S/2))`), `B`, the erosion fraction (0.15),
every QC threshold including the 0.50 coverage rule, the voxel grid (`h=0.5 Å`, fallback
`1.0 Å`), the subset design and sampling policy (deterministic, recorded before iteration 1,
seeded combinatorial unranking), and every comparison/classification metric (Jaccard, Dice,
volume, surface area, components, Hausdorff, centroid shift, TP/FP/TN/FN, MCC, sensitivity,
specificity, PPV, NPV, F1, balanced accuracy, preservation frequencies) — now computed once per
iteration, on the sole fixed-`r_fp` reconstruction, instead of twice (once per reconstruction).
`config/pipeline.yaml` is byte-unchanged; no `FROZEN_ASSERTIONS` value was touched.

## Post-decision re-profile

See the session's final report for the single-iteration measurement taken after this change
and the resulting Stage D total estimate against the 40.0 h budget.
