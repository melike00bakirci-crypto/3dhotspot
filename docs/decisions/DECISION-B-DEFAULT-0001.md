# DECISION-B-DEFAULT-0001

```
DECISION_ID    = DECISION-B-DEFAULT-0001
DATE           = 2026-08-18
AUTHORIZED_BY  = user
STATUS         = IMPLEMENTED
SCOPE          = permutation.B_default (Stage B primary + secondary-null permutation count) —
                  pipeline-wide, all genes
```

## Decision

`permutation.B_default` is raised from **10,000 to 100,000**, effective for every gene run under
this pipeline from this point forward.

This is a **user-authorized scientific design decision**, not a bug fix and not a QA-initiated
repair. `pipeline-qa-debugger` has no authority to make this class of change (scientific
methodology firewall) and did not make it. `hotspot-statistics` itself is explicitly forbidden
from changing `B` ("Never change `B` or the FDR method yourself... state the ambiguity, options,
consequences, recommendation — then stop" — `.claude/agents/hotspot-statistics.md` §15). The user
made this call directly, in response to the Lead's report on the first live CACNA1A run
(`RUN_ID 20260818T215927Z_58546cc4_1e029297`), which terminated `UNDERPOWERED` because the
frozen `B = 10,000` could not resolve a small enough p-value for that run's realized test family
(`m = 792`).

## Rationale (as given by the user)

> B is a resolution parameter, not a threshold, null model or cohort definition — raising it
> changes no statistical criterion, only the precision of the empirical p-value. CACNA1A requires
> B ≥ 15,839; 100,000 clears that with margin and matches KCNA2's own `b_recommended`.

This tracks the pipeline's own classification of `B`: the permutation-resolution diagnostic
(II.4) exists precisely because `B` sets the *achievable resolution* of the empirical p-value
(`p_res = 1/(B+1)`) rather than any scientific criterion — `q`, the null model, `fdr.method`, the
cohort definition, and every threshold in the methodology are untouched by this change. KCNA2's
own `06_FINAL_HOTSPOTS` run (`20260818T163102Z_58546cc4_68f90db3`, B=10,000, certificate PASSED)
already reported `b_recommended = 100,000` in its own diagnostic — i.e. the pipeline was already
recommending this exact value for full resolution before CACNA1A ever ran; CACNA1A's failure just
made the recommendation binding rather than advisory.

## Cross-gene comparability

**All genes run under this pipeline must use the same `B` for results to be comparable.**
KCNA2's completed run (`20260818T163102Z_58546cc4_68f90db3`, B=10,000) is superseded for
comparability purposes by this decision and will be re-run at B=100,000 under a new RUN_ID in
this same session, alongside CACNA1A. The original KCNA2 run directory is left in place
(nothing is overwritten — RUN_ID is content-addressed via `config_sha256`, so the new B produces
a distinct RUN_ID automatically) but should not be cited as the current production result for
KCNA2 once its B=100,000 counterpart exists.

## What changed

- `config/pipeline.yaml` — `permutation.B_default: 10000` → `100000` (still FROZEN, new pinned
  value; comment updated to cite this decision).
- `src/hotspot3d/utils/config.py` — the `FROZEN_ASSERTIONS` entry for `permutation.B_default`
  updated from `10000` to `100000`, so `assert_frozen_methodology` pins the new value and a run
  configured with any other `B` (outside a synthetic context) is still rejected before any stage
  executes — the enforcement mechanism itself is unchanged, only the pinned constant is.
- `src/hotspot3d/hotspot/stage.py` — Stage B's own explicit guard (independent of
  `assert_frozen_methodology`, checked again at Stage B entry) updated from `!= 10000` /
  `"frozen default 10000"` to `!= 100000` / `"frozen default 100000 (DECISION-B-DEFAULT-0001)"`,
  in both the non-synthetic `BlockedError` path and the synthetic `SYNTHETIC_INPUT_MODE` warning
  path. `config_sha256` changes as a direct, expected consequence (verified:
  `4ec9d6a82a9dcf42…`, was `58546cc41a21585a…`), which is exactly what gives every subsequent
  RUN_ID its own identity — no run_id collision with any B=10,000 run is possible.

## What was preserved, unchanged

`primary_null` (`structure_aware_positional`), `local_null`, `p_value_estimator`
(`phipson_smyth_add_one`), `exclude_uninformative_centers`, `fdr.method` (`BH`), `fdr.q` (0.05),
the power certificate mechanics (§5.4) and its formulas, the permutation-resolution diagnostic's
own formulas/ladder/safety factor (`s=10`, `k_target=1`, `auto_change_B: false`),
`robustness.N_CAP`, `robustness.max_wall_seconds`, every radius-selection and radius-domain
constant, `r_hot`/`r_fp` methodology, and every other `FROZEN_ASSERTIONS` entry not named above.
Only the single scalar `B_default` moved, in the three places above.

## Known follow-up NOT done as part of this decision

The test suite was not updated. In particular:

- `tests/hotspot_stats/test_stage.py::test_production_default_is_ten_thousand_permutations` now
  asserts a stale literal (`== 10_000`) against the shipped config and will fail.
- `tests/integration/test_config_provenance.py` (`test_d_production_config_is_unchanged_and_
  unwritten`) asserts `shipped.get("permutation.B_default") == 10000` and will fail for the same
  reason.
- `tests/integration/test_underpowered_terminal_state.py` is calibrated to the OLD frozen value
  in a way that is not a mechanical constant swap: its synthetic `permutation_resolution` fixture
  was deliberately tuned (via a swept `run_id`) to realize a test family `m ≈ 529`, chosen to
  clear the old threshold `m > q·(B+1) = 500.05` by a narrow, checked margin. At B=100,000 that
  threshold becomes `m > 5000.05`, which this fixture's realized family size does not clear — the
  module's entire premise (a genuinely underpowered case at the frozen `B`) needs new
  domain-level recalibration (a new or resized fixture, a new swept `run_id`, or an argument that
  no such case is reachable at this scale), which is `hotspot-statistics`' owned judgment call
  (`tests/hotspot_stats/**`-adjacent), not a mechanical edit. Routing this to `hotspot-statistics`
  is a follow-up, not done here, so as not to delay the re-runs this decision exists to unblock.

This does not affect the validity of any real-gene run — `config/pipeline.yaml`,
`utils/config.py`, and `hotspot/stage.py` (the three places that actually govern what a live
`hotspot3d run` does) are fully updated and internally consistent, confirmed by
`hotspot3d validate-config` (`config_sha256 4ec9d6a82a9dcf42…`).

## Post-decision re-runs (this session)

Both genes re-run at B=100,000 under new RUN_IDs, same `config_sha256 4ec9d6a82a9dcf42…`:

- **CACNA1A** `20260818T233025Z_4ec9d6a8_ad8f32d6` — `COMPLETED_WITH_WARNINGS`, exit 0. The
  power certificate that failed at B=10,000 (`m=792`, `p_floor/c_1=1.584`) now passes with margin
  (`p_res=1/100,001=9.9999e-6 << c_1=6.313e-5`); 72 significant centers, `PERMUTATION_RESOLUTION_
  LIMITED=False`. Full report: `results/CACNA1A/CACNA1A_20260818T233025Z_4ec9d6a8_ad8f32d6/
  LEAD_REPORT.md`.
- **KCNA2** `20260819T010956Z_4ec9d6a8_ad8f32d6` — `COMPLETED_WITH_WARNINGS`, exit 0. `r_hot`
  moved 10.0 Å → 9.5 Å and significant centers 50 → 46 versus the B=10,000 run — an expected
  consequence of `B` changing the precision of every permutation-derived selection objective, not
  only the final significance test. Full report: `results/KCNA2/
  KCNA2_20260819T010956Z_4ec9d6a8_ad8f32d6/LEAD_REPORT.md`.

Both genes now share `config_sha256 4ec9d6a82a9dcf42…`, satisfying the cross-gene comparability
requirement this decision states. The B=10,000 run directories for both genes remain on disk for
provenance but are superseded and should not be cited as current results.
