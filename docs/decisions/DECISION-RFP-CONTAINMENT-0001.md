# DECISION-RFP-CONTAINMENT-0001

```
DECISION_ID    = DECISION-RFP-CONTAINMENT-0001
DATE           = 2026-10-24
AUTHORIZED_BY  = user
STATUS         = IN PROGRESS
SCOPE          = Stage C (II.9/II.10) — how r_fp is determined. Pipeline-wide.
REVOKES        = F10 (in full), DECISION-FOOTPRINT-DOMAIN-0001 (operative effect)
```

## Decision

`r_fp` is no longer selected by an independent multi-objective search. It is
**determined** by a containment constraint:

```
evidence set  E  = { i in P/LP : min_{c in S} d(i,c) <= r_hot }
r_fp             = max_{i in E} min_{c in S} d(i,c)
footprint   FP   = { i in U_struct : min_{c in S} d(i,c) <= r_fp }
```

`S` = SIGNIFICANT_HOTSPOT_CENTERS, `d` = Ca-Ca Euclidean, `P/LP` = the pathogenic
half of the classified cohort.

Two non-selecting rules accompany it:

- `coverage > 0.50` -> ESCALATION and `COMPLETED_WITH_WARNINGS`. The radius is
  **never reduced** to satisfy it; reducing it would break the containment
  guarantee, which is the whole content of the decision.
- `r_fp < 7.0 A` -> flag "detection scale below 3D resolution". The radius is
  **never raised**.

`B/LB` enters no step. The geometric descriptors (compactness, convexity,
connectivity, bridging index, expansion rate, ARI stability, component count) are
computed and reported but decide nothing.

## Why

The statistical claim this pipeline makes is about `r_hot` neighbourhoods: a centre
is significant because *its `r_hot` sphere* is enriched in P/LP beyond the positional
null. Delimiting the result with a radius smaller than that draws a boundary sharper
than the evidence supports.

**TP53RK is the measured instance of the defect.** Archived run
`TP53RK_20260819T084210Z_4ec9d6a8_902e1b03`: `r_hot = 9.5 A`, one significant centre
(residue 62), three P/LP residues (42, 65, 81) at 9.473 / 9.298 / 7.813 A from it —
all three inside the detection neighbourhood, all three in the test family, none
significant on its own. The independent search selected `r_fp = 4.5 A`, giving a
4-residue footprint (45, 61, 62, 63) that contained **zero of the three pathogenic
residues**. Under this decision `r_fp = 9.473 A`, the footprint is 20 residues and
contains **3/3**.

## Properties

**L1 — `E` is never empty when `S` is non-empty.** A centre is significant because its
sphere is P/LP-enriched; a sphere with `n_P = 0` cannot reject. `r_fp` is therefore
always defined.

**L2 — `r_fp <= r_hot`, structurally.** Members of `E` satisfy `d <= r_hot` by
definition, so their maximum cannot exceed it. The footprint is never wider than the
detection neighbourhood.

**L3 — containment.** `E` is a subset of `FP`. Every supported pathogenic residue is
inside. This is the defining guarantee.

**L4 — minimality.** `r_fp` is the *smallest* radius satisfying L3: for any
`r < r_fp` some member of `E` falls outside. The radius is the unique solution of a
constraint, not the argmax of an objective. No weight, tolerance or normalisation
appears anywhere in it.

**L5 — `B/LB` independence.** `r_fp` is computed without reading any B/LB residue.
This is what makes "the footprint contains N benign residues" a *finding* rather than
the output of an objective function.

## What is revoked

### F10 — the config flag and the declarative guard, NOT the type-level interlock

F10 stated that `r_fp` "is neither set equal to nor bounded below by `r_hot`". It was
enforced in two places and the two are treated differently.

**Revoked:** the declarative check in `footprint/params.py`, which only read
`footprint_domain.bound_below_by_r_hot`. A guard that inspects a config flag cannot
see what the code actually does — implementing this rule without touching the flag
would have left it silent. It is replaced by a **behavioural** guard (AST + signature
+ field inspection, modelled on the Stage A star audit in `data/audit.py`).

**Retained and narrowed, NOT removed:** `ReadOnlyHotspotRadius` in
`footprint/inputs.py`. Its numeric dunders — `__float__`, `__int__`, `__lt__`,
`__le__`, `__gt__`, `__ge__`, `__add__`, `__sub__`, `__mul__`, `__truediv__` — all
keep raising. The wrapper exists to make accidental leakage of `r_hot` into Stage C
impossible, and that purpose survives this decision intact.

Exactly one new authorised method is added:

```
ReadOnlyHotspotRadius.derive_evidence_set(
    *, centers, universe, cohort_plp) -> EvidenceSet
```

It reads `self._value` only to produce a per-residue `in_evidence` boolean. **The
numeric value never leaves the method.** `EvidenceSet` carries, per mapped P/LP
residue, `(residue_id, nearest_center_id, distance_A, in_evidence)` and nothing else.
`r_contain` is computed outside the wrapper, from that output.

A consequence worth stating: because `derive_containment_domain` takes the
`EvidenceSet` rather than a radius, **`footprint/domain.py`'s original guarantee —
"`r_hot` is not an argument to any function in this module" — remains literally
true.** This decision opens a narrow, named door in `inputs.py`; it does not weaken
the module the guarantee was written about.

`value_for_coverage_query` is also retained. Audited: its five production call sites
all write the value into handoff or provenance JSON and none of them feeds a
decision.

### DECISION-FOOTPRINT-DOMAIN-0001

Its operative content — the search domain `[5.0, min(0.25*D_max, 25.0)]` — has no
subject once the search is gone. The record is retained for provenance: the batch in
`results_v2/` was produced under it, and its measured findings (6 of 8 BLOCKED genes
failed on the MST-derived bound; 18 of 50 genes had `r_fp` pinned at the floor) are
the evidence that an independent search was unreliable in the first place.

### Partially: DECISION-NEAR-TIE-0001

The `footprint_selection.near_tie_*` half becomes inert — there is no runner-up to be
near. The `radius_selection.near_tie_*` half, which governs `r_hot`, is untouched and
remains in force.

## Separate item: the voxel grid stopped being padded by a retired bound

This is not part of the selection rule. It is a defect the rule *exposed*, fixed in the
same change because the retired keys could not honestly be labelled inert while three
live code paths were still reading one.

**What was wrong.** The grid that the distance transform is sampled on was padded by
`domain.rho_max`. Under the retired search that was `min(0.25*D_max, 25.0)` — a bound
unrelated to the radius being measured. Under containment the domain is a single point,
so `rho_max == r_fp` and the sweep's own grid was already correct; but two other callers
rebuilt the retired domain to get their padding:

| Call site | What it is | Consequence |
|---|---|---|
| `footprint/api.py` `reconstruct_at` | Stage D's reconstruction, the sensitivity band | box inflated 3–10x; measured geometry depended on the domain the radius sat in |
| `footprint/regenerate_grid.py` | the Output Contract P5 regeneration recipe | regenerated lattice offset from the run's — the grid it promises to reproduce byte-for-byte, it did not |

`reconstruct_at` is production code, not a reproducer: `robustness/perturb.py` calls it
once per Stage D iteration and `robustness/stage.py` calls it for the baseline
reproduction check.

**The fix.** `sweep.build_sweep_grid` now takes the radius as a number instead of reading
it off a domain object, and all three callers pass the radius they are about to evaluate
(`max(domain.grid)` for the sweep, `rho` for the reconstruction, `radius` for the CLI).
The grid is a function of `(centers, rho, h)` alone. No second grid builder was added —
that is what the regeneration module's docstring forbids, and the reuse is now what makes
the byte-identity claim true rather than merely asserted.

Two things went with it. `FixedReconstruction.inside_iteration_domain` reported membership
of a domain that no longer exists, and nothing read it. `regenerate_grid`'s guard refused
any radius outside `[5.0, min(0.25*D_max, 25.0)]`; under containment `r_fp` may legitimately
fall below `5.0`, and the band's radii are the resolution floor up to `r_hot`, so the guard
could only reject requests that were in fact valid. It is replaced by a positivity check.
`--universe` is still accepted and validated but no longer reaches the geometry.

**The gate, and what the measurement found.** The fix was made conditional on residue
membership and component counts being untouched — a grid-origin shift must not move the
footprint's content. Old padding vs new padding, each gene at its own `r_fp`:

| gene | `r_fp` (A) | residue set | components | volume | surface | compactness |
|---|---|---|---|---|---|---|
| CACNA1A | 8.0 | identical | 5 -> 5 | 56656 -> 56656 (+0.00%) | +0.00% | +0.00% |
| FGFR3 | 10.5 | identical | 2 -> 2 | 47487 -> 47487 (+0.00%) | +0.00% | +0.00% |
| KCNA2 (`r_hot` 9.5) | 9.0 | identical | 2 -> 2 | 34392 -> 34392 (+0.00%) | +0.00% | +0.00% |
| KCNA2 (`r_hot` 10.0) | 10.0 | identical | 1 -> 1 | 43573 -> 43573 (+0.00%) | +0.00% | +0.00% |
| PTPN11 | 6.0 | identical | 1 -> 1 | 7915 -> 7947 (+0.41%) | +0.51% | -0.24% |
| TP53RK | 9.5 | identical | 1 -> 1 | 3584 -> 3584 (+0.00%) | +0.00% | +0.00% |

Stage D's perturbed path was measured separately, since that is where `reconstruct_at`
actually runs: 12 center-removal iterations each for CACNA1A and PTPN11, old padding vs
new. **24 of 24 reconstructions returned an identical residue set and an identical
component count.** Volume moved only in PTPN11, by at most 0.49%.

So the gate holds: only the voxel-derived descriptors move, in one gene, by well under a
percent. The numbers above are the record of what this change did, separate from what the
selection rule did, so that the re-run's differences can be attributed.

### Why PTPN11 is the only one that moves — traced, not inferred

Two plausible-sounding explanations were checked and both are **wrong**, so they are
recorded here to stop them being repeated. It is not that `r_fp = 6.0 A` sits furthest
below the retired ceiling and so had the largest padding share: CACNA1A's `r_fp` (8.0 A)
is also well below it and moves by exactly nothing. It is not that PTPN11 has the
smallest footprint or the highest surface-to-volume ratio: TP53RK's footprint is smaller
(3584 vs 7915 A^3) and does not move at all.

The actual cause is a floating-point tie in voxel assignment:

| | |
|---|---|
| PTPN11 center 507, y | `-0.104` |
| old origin y (pad 27.0) | `-27.104` -> `(p-o)/h` = `54.0` exactly -> `floor` = **54** |
| new origin y (pad 8.0) | `-8.104` -> `(p-o)/h` = `15.999999999999998` -> `floor` = **15** |

The center lies **exactly on a voxel boundary**. `GridSpec.index_of` assigns a voxel with
`floor((p - origin) / h)`, which is discontinuous there, so in exact arithmetic the two
origins give 54 and 16 — a difference of exactly the 38-voxel origin shift, i.e. the same
cell — but in float64 the second division lands one ULP below 16 and floors to 15. One
seed voxel therefore moves by one voxel. The EDT around it reclassifies a shell of voxels
at the `r_fp` isosurface: 434 voxels differ, net +260 occupied, and 260 * 0.5^3 = 32.5 A^3
against an observed volume change of +32 A^3. The differing voxels are all interior (3 to
23 voxels from the box edge, none adjacent to it), which is what rules out clipping or any
boundary artefact of the smaller box.

So the discriminator is a coordinate coincidence, not footprint geometry. CACNA1A: 0 of 72
centers change cell. PTPN11: 3 of 23 centers sit exactly on a lattice boundary and one of
them flips. The surface-to-volume ratio governs *how many* voxels flip once a seed moves;
it does not govern *whether* one moves.

**This is a latent fragility, independent of this change.** Any edit that moves the grid
origin — this padding fix, a different `bbox_padding_A`, a coordinate revision upstream —
can flip such a tie, and `index_of` has no tie-breaking rule. It is recorded rather than
fixed: adding one (`floor(x + eps)`, or snapping the origin to a multiple of `h`) changes
the voxel assignment of every tied center in every gene, which is a geometry change
needing its own decision. Flagged for that decision; see weakness 7.

**Not claimed.** Six archived runs and two genes' perturbation samples are not the 68-gene
batch. Nothing here establishes that no gene anywhere shifts a residue; it establishes
that none of the six did, that the two measured Stage D paths did not, and that the
mechanism — a sub-voxel origin shift — cannot move a residue whose distance to its nearest
center is not within half a voxel edge of `r_fp`.

## Separate item: three run artefacts were still describing the retired search

Found while wiring the exports, and fixed in the same change. These are reporting
defects, not methodology — but the first one made a run's own records state the
opposite of the rule that produced it, which is worse than a missing field.

**`footprint_decision.json` asserted the inverse of the rule.** The file named
"decision" carried `r_hot_used_in_selection: false`, `r_fp_bounded_by_r_hot: false`
and an `f10_statement` reading "r_fp was selected from geometry alone. r_hot is not an
argument to the footprint API and r_fp is neither set equal to nor bounded below by
it." Under containment `r_hot` defines the evidence set and L2 makes `r_fp <= r_hot`
structural, so an auditor reading that file would have been told the dependence does
not exist. It is rewritten around the determination: `r_contain` exact and rounded,
the grid slack, the binding residue and whether it is a distance outlier, the evidence
counts, and an explicit `r_hot_dependence_statement` that says what the dependence is
and what still bounds it (one boolean, no numeric value, `ReadOnlyHotspotRadius`
intact). The retired keys are removed rather than left to decay, and a test asserts
their absence so a later edit cannot quietly restore a false claim.

**`footprint_domain.json` did not contain the record it is declared to contain.**
`EXPECTED_07` calls it "containment record: r_contain, evidence counts, binding
residue, outlier fence, MST merge scales"; the containment record was in fact only in
`handoff_03`. It now carries `containment` in full. Its `note` field also still cited
F10 and claimed "Merge scales inform rho_max" — both revoked.

**The regeneration recipe written into every run reproduced the wrong lattice.** The
recipe string said `pad=rho_max + 2.0`. That is the same defect as the padding item
above, recorded into the artefact a reader is told to follow. It now goes through
`build_sweep_grid`, so recipe and run share one implementation. `regenerate_grid`
itself also refused any radius outside `[5.0, min(0.25*D_max, 25.0)]` — a bound that
no longer exists and that `r_fp` may legitimately fall below — and padded from its
ceiling, so the grid it returned was not the run's grid. Verified after the fix on
PTPN11 at `r_fp = 6.0`: same shape `(70, 65, 66)`, same origin, same voxel edge,
occupancy identical element for element.

## Separate item: Stage D now reports the radius each iteration would have determined

`r_fp` is a function of the center set, so "how much does the footprint depend on any
one center" became a question with a direct numeric answer — and Stage D is where it
is already removing centers. Each iteration now records `robustness_r_contain`,
`robustness_r_fp_would_be` and the residue that would bind it, aggregated in
`robustness_profile.json` under `containment_sensitivity`.

Three properties make it safe to report:

* It is computed over the **frozen** Stage C evidence set, re-measured against the
  surviving centers. Membership is not re-derived, which would need `r_hot` inside
  Stage D and would mix "a residue left the evidence set" into a number meant to
  isolate "the surviving residues got further from a center".
* It never reaches a reconstruction. `DECISION-STAGE-D-FIXED-RFP-0001` holds the
  radius at Stage C's `r_fp` for every iteration, `robustness_r_fp_invariant` still
  asserts it, and the AST guard still forbids the determination inside `perturb.py`.
* Removing centers can only increase a nearest-center distance, so
  `would_be >= r_fp_original` by construction. A value below it is a fault, and
  `monotonicity_holds` reports that rather than assuming it.

On the four-center synthetic case every iteration would have moved the radius, which is
what a four-center set should do and is not a statement about real genes; the archived
leave-one-out measurement on the six real runs moved `r_contain` by 0.030-0.836 A, four
of six below the 0.5 A grid step. The point of the field is that this stops being a
one-off measurement and becomes a per-gene number in the record.

## Verified before implementation

Reconstructed from archived runs with `footprint.api.reconstruct_at`, Stage B not
re-run. Six runs with at least one significant centre; MEFV excluded (0 significant
centres in its archived ClinVar run).

| gene | r_hot | old r_fp | old residues | old P/LP | **r_fp** | residues | P/LP | B/LB | components |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| CACNA1A | 8.0 | 10.0 | 413 | 85/115 | **8.0** | 274 | 78/115 | 1/98 | 5 |
| FGFR3 | 12.0 | 12.5 | 133 | 24/34 | **10.5** | 118 | 24/34 | 3/26 | 2 |
| KCNA2 (rh 9.5) | 9.5 | 10.0 | 130 | 34/38 | **9.0** | 117 | 34/38 | 0/12 | 2 |
| KCNA2 (rh 10) | 10.0 | 10.0 | 143 | 36/38 | **10.0** | 143 | 36/38 | 0/12 | 1 |
| PTPN11 | 6.0 | 3.0 | 23 | 21/56 | **6.0** | 51 | 28/56 | 0/11 | 1 |
| TP53RK | 9.5 | 4.5 | 4 | **0/3** | **9.5** | 20 | **3/3** | 0/8 | 1 |

`coverage <= 0.50` held in all six (highest 0.287, KCNA2). The `r_fp < 7.0` flag fires
for PTPN11 (5.850) — one of six, not a rare event.

Stability: removing the single residue that determines `r_fp` moves it by 0.030 /
0.120 / 0.175 / 0.211 / 0.720 / 0.836 A. Four of six are below the 0.5 A grid step.
The determining residue is a distributional outlier (beyond `Q3 + 1.5*IQR`) in **one**
case, PTPN11 (5.850 against a fence of 2.386, median 0.000).

## Known weaknesses, recorded rather than fixed

**1. `r_fp` inherits `r_hot`'s instability with no buffer.** `r_fp` is now a function of
`r_hot`, and `r_hot` is still chosen by the multi-objective machinery this decision
rejects for `r_fp`. Measured on KCNA2, whose `r_hot` has flipped between 9.5 and 10.0
across runs at a 0.18-1.58 % margin:

```
r_hot 9.5  -> r_fp 9.0  -> 117 residues, 2 components, 34/38 P/LP
r_hot 10.0 -> r_fp 10.0 -> 143 residues, 1 component,  36/38 P/LP
```

26 residues and the "single contiguous region" claim turn on that coin flip. The
independent search partially buffered this; this decision removes the buffer.
**Mitigation adopted: `r_contain` is also computed at the runner-up `r_hot` and both
are reported.** It does not change the selection.

**2. Density is traded for coverage.** Density falls monotonically with `r`; in PTPN11
it is 0.913 at the centres and 0.549 at `r_fp`. Containment was privileged over
concentration. That is a choice, not a result.

**3. TP53RK sits on a knife edge.** At `r_fp = 9.5` it captures 3/3; at 9.0, 1/3. The
strongest piece of evidence for this decision is also its most radius-critical case.

**4. Component merging is radius-critical.** KCNA2 (rh 10) drops from 2 components to 1
at exactly 10.0. "One contiguous region" is therefore a statement about the exact
radius, and must be reported as such.

**5. Up to half the pathogenic residues can lie outside the evidence set** — CACNA1A
32.2 %, FGFR3 29.4 %, PTPN11 50.0 %. These are real pathogenic variants that no
centre-anchored footprint can contain. Reported, never hidden.

### Weakness 6 — `r_contain` can be 0, and the rule as written does not forbid it

If every mapped P/LP residue is itself a significant center, every nearest-center
distance is 0, so `r_contain = 0` and `r_fp = 0`. The rule produces it without
complaint; the geometry then cannot be represented — a zero-radius union of balls has
no volume, the mesh fails QC-F3, and the run terminates without a footprint.

Reached while recalibrating the tests, on a synthetic case where all 64 universe
residues were centers. Whether it is reachable on a real gene is not established here:
it needs a gene whose entire mapped P/LP set is significant, and none of the six
archived runs is close. `DETECTION_SCALE_BELOW_3D_RESOLUTION` already flags a radius
below 7.0 A, so the condition is *reported*; what is not decided is whether a
determined radius below the voxel edge should be a positive lower bound, an immediate
`NegativeResult` naming the degeneracy, or left exactly as it is. **Left open for a
decision rather than patched**, because every option changes the rule and the rule is
the user's to set. The behaviour today is: reported, then terminated by QC with no
footprint.

### Weakness 7 — voxel assignment has no tie-breaking rule

`GridSpec.index_of` uses `floor((p - origin) / h)`. A center whose coordinate sits exactly
on a voxel boundary is assigned by whichever side the float64 division happens to land on,
so the same physical point can occupy different cells under different grid origins. Found
while tracing PTPN11's volume change (above): 3 of its 23 centers sit on exact lattice
boundaries, and one of them flipped under a 19 A origin shift, moving the reported volume
by 32 A^3.

Consequences, all bounded: residue membership was identical in every case measured here,
because a residue flips only if its distance to the nearest center is within half a voxel
edge of `r_fp`; `r_fp` itself is immune, since containment is computed on C-alpha distances
and never touches the grid. What is affected is the voxel-derived descriptors — volume,
surface area, compactness — and therefore their comparability across any run whose grid
origin differs.

Left open deliberately. The candidate fixes (an epsilon in `floor`, or snapping the origin
to a multiple of `h` so ties resolve identically everywhere) each reassign every tied
center in every gene, which changes reported geometry pipeline-wide and is a decision in
its own right, not a wiring fix.

## Mandatory reporting

Produced per run, from `evidence_set.tsv` and the sensitivity band:

1. `n_plp_outside_evidence` and its fraction
2. **`07/sensitivity_band.tsv`** — `r` from 7.0 A to `r_hot` in 0.5 A steps (or from
   `r_contain` to `r_hot` when `r_hot < 7.0`, possibly a single row), columns
   `r | n_residues | n_components | n_plp | n_blb | coverage`, with the selected
   `r_fp` row marked. 6-11 extra `reconstruct_at` calls per gene. **This, not the
   runner-up radius alone, is the real mitigation for weakness 1.**
3. binding residue id, the evidence-distance distribution (median, Q1, Q3, IQR), and
   `binding_residue_is_outlier` = distance > Q3 + 1.5*IQR. The verdict is recorded in
   `handoff_03`, not left for the reader to derive. It fires for PTPN11 (5.850
   against a fence of 2.386).
4. B/LB content — an observation, since B/LB entered no step
5. component count and the residue ranges of each component
6. `r_contain` recomputed at the runner-up `r_hot`, together with
   `r_hot_runner_up_margin_rel`. The radius alone says "it could have been 10.0"; the
   margin says "it was a coin flip" (KCNA2: 0.178 %). Stage B now writes both into
   `handoff_02` unconditionally — reporting only, its selection is untouched.
7. **`containment_sensitivity`** in `robustness_profile.json`, with
   `robustness_r_contain` / `robustness_r_fp_would_be` / the binding residue on every
   row of `perturbation_results.tsv`. `r_fp` is a function of the center set, so "how
   much of the footprint rests on any one center" has a numeric answer and Stage D is
   already removing centers. Report `fraction_would_differ` and `max_increase_A`
   alongside `r_fp`; a run where most iterations would have moved the radius is a run
   whose footprint rests on a handful of centers, whatever the Jaccard profile says.

## Out of scope for this decision

**The candidate list of unclassified footprint residues.** `classified_cohort.tsv`
carries only `PLP` and `BLB`, so a VUS-bearing residue cannot be distinguished from a
residue with no variant at all. Reporting them merged would present two different
things as one. Deferred to a separate change; the cohort schema has to carry a VUS
class first.

## Required follow-up, NOT done here

**`.claude/agents/footprint-robustness.md` still instructs the retired search.** 637
lines, with roughly 34 referencing concepts this decision removed — including its
`description:` frontmatter, which is the text the Agent tool surfaces when the agent is
dispatched. As it stands, a future Stage C/D agent invocation would be told to own "the
independent r_fp search domain and optimization ... Pareto plus distance-to-ideal
selection of r_fp", to treat F10 as in force, and to write `r_fp_scan.tsv`. §4
("Selection"), §6 (the output list), the F10 clauses at lines 47/59/125/256 and the
Pareto re-verification step all describe the retired procedure.

Not touched here deliberately: it is a substantial rewrite of a charter, not wiring, and
charters are the user's to approve. It is the single largest remaining inconsistency and
should be done before the next Stage C/D agent run. The same pass should pick up the
near-tie and domain descriptions already flagged as stale in
`hotspot-statistics.md:123` and `footprint-robustness.md:55-57`.

## Scope limits

- `r_hot` selection (Stage B, II.6) is untouched
- `permutation.B_default`, the null models, `fdr.*` untouched
- `QC_F1_max_coverage = 0.50` keeps its value; only its role changes, from a
  candidate filter to an escalation threshold
- Stage E annotation untouched
- Seed derivation untouched
