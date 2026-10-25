# DECISION-UNDERPOWERED-FIXTURE-0001

**Status:** IN PROGRESS
**Scope:** one synthetic fixture, and lifting the skip it was blocking. No statistic,
threshold, null or cohort policy changes.
**Closes:** the `DECISION-B-DEFAULT-0001` follow-up recorded in
`tests/integration/test_underpowered_terminal_state.py`.

## The hole

UNDERPOWERED means the test **could not have rejected anything regardless of the data**
— the attainable p-value floor sits above the rank-1 bar `c_1 = q/m`. It is not "we
tested and found nothing", and Workflow v2 Appendix A forbids reporting it as a negative
just as firmly as it forbids reporting it as a fault.

Both misreports have happened. The one that motivated the orchestrator hook: Stage B
returns a complete handoff rather than raising, so before the hook existed the pipeline
called Stage CD next, which refused on `qc_status = FAIL` with a `BlockedError` — and the
run was typed `BLOCKED / TECHNICAL_FAILURE`. An uninformative run reported as a pipeline
fault.

The hook was added and a module was written to guard it. **That module has been skipped
in its entirety**, so the protection has not been running:

* the fixture it used, `permutation_resolution` (M = 620), reached the state through the
  **permutation-resolution floor** `p_res = 1/(B+1)`, which held at `B = 10,000` because
  the family size cleared `q*(B+1) = 500.05`;
* `DECISION-B-DEFAULT-0001` raised `B` to 100,000. The line moved to `m > 5000.05`, and
  raising `B` *improves* resolution — the floor gets harder to hit, not easier. No
  synthetic fixture reaches it;
* an earlier attempt to force the state by overlaying `B` down to 200 was correctly
  rejected by `assert_frozen_methodology`;
* so the module was skipped rather than weakened, with a note assigning the missing
  fixture to `hotspot-statistics`.

## The fixture

`single_pathogenic_residue` — `length = 180, n_plp = 1, n_blb = 19`, uniform spread.

It reaches TEST_CANNOT_REJECT through the **combinatorial** floor, which does not move
with `B`. The arithmetic: with `N_P = 1` the most extreme observation any sphere can
produce is its single pathogenic residue, so the smallest attainable tail is `n_U/M` —
of order 1e-2 on a protein this size — while `c_1 = q/m` is of order 1e-3 or below once
the family has tens of members. `p_floor > c_1` therefore holds from the cohort alone,
at any `B`.

Measured on a real run at the frozen `B = 100,000`:

```
binding_floor            combinatorial
p_floor                  0.01666666666666784
c_1_rank1_critical_value 0.0006849315068493151
m_test_family_size       73
increasing_B_is_futile   True
PASSES                   False
TEST_CANNOT_REJECT       True
```

Run status `UNDERPOWERED / UNINFORMATIVE`, `negative_result` absent,
`05_HOTSPOT_RADIUS = UNDERPOWERED/UNINFORMATIVE`, `06_FINAL_HOTSPOTS = NOT_RUN/NOT_APPLICABLE`.

**`increasing_B_is_futile = True` is the property that matters.** It is what makes this
fixture stable against any future change to `B`, which is precisely how the previous one
decayed.

**This is not a parameter chosen to force an outcome.** It is the degenerate cohort, and
the ACMG SF v3.2 batch contains it for real: **DSC2 (N_P = 1, N_B = 19)** and
**TMEM43 (N_P = 1, N_B = 15)** both terminate UNDERPOWERED. The fixture mirrors DSC2.
`PLN (N_P = 2, N_B = 1)` and `SDHAF2 (N_P = 2, N_B = 6)` are the two other real cases.

## What changed in the test module

* The module-level `pytest.mark.skip` is **lifted**; all 6 tests run and pass.
* The fixture call points at `single_pathogenic_residue`.
* `_certificate()` reads the certificate from `05_HOTSPOT_RADIUS` or
  `06_FINAL_HOTSPOTS`, because with the combinatorial floor the PREFLIGHT certificate
  fires and Stage B stops at 05 — so 06 is NOT_RUN and holds none. The old route filled
  06.
* Two tests were rewritten for the route the fixture now takes. They had asserted
  `binding_floor == "permutation_resolution"`, `increasing_B_is_futile is False` and a
  `B_to_clear_c_1` ladder value — all correct for the resolution branch and all the
  opposite here. They now assert the combinatorial branch's own contract, including that
  `p_floor == p_comb_best` and that no `B` clears it.
* A `_binding_block()` helper returns whichever certificate half carries a floor. The
  contract is that the certificate names which floor binds, not which half is filled.
* The obsolete skip rationale is kept as history rather than deleted: it records why the
  fixture was hard to build, which is why the protection lapsed.

## Not claimed

* Nothing about whether the four real UNDERPOWERED genes in `results_v2/` are classified
  correctly — they already are. This restores a guard against a future regression; it
  does not change any current result.
* The resolution-floor branch of the certificate is now exercised by no fixture. It is
  still covered by unit tests in `tests/hotspot_stats/`, but no end-to-end run reaches
  it at the frozen `B`, and by the arithmetic above none can.
