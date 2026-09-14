# DECISION-FOOTPRINT-DOMAIN-0001

```
DECISION_ID    = DECISION-FOOTPRINT-DOMAIN-0001
DATE           = 2026-09-12
AUTHORIZED_BY  = user
STATUS         = IMPLEMENTED
SCOPE          = footprint_domain (II.9) — the candidate r_fp search domain,
                 pipeline-wide, all genes
```

## Decision

The candidate `r_fp` search domain changes from

```
[3.0,  min(1.25 * rho_all, 0.25 * D_max, 20.0)]
```

to

```
[5.0,  min(0.25 * D_max, 25.0)]
```

Three separate edits, one of which is the substantive one:

| Change | From | To | Class |
|---|---|---|---|
| `rho_min_A` | 3.0 | **5.0** | floor raised |
| `hard_ceiling_A` | 20.0 | **25.0** | ceiling raised |
| **`1.25 * rho_all` term** | caps `rho_max` | **no longer caps `rho_max`** | **substantive** |

The MST-derived term is *removed from the upper bound only*. Merge scales continue to
be computed and reported as merge events — that reporting is untouched.

`dmax_fraction_cap = 0.25` is **unchanged**: a scale-relative cap ("the footprint may not
span more than a quarter of the protein") is well-motivated and was never the defect.

The effective ceiling is now `min(0.25 * D_max, 25.0)`, backed by the **unchanged**
`QC_F1_max_coverage = 0.50` hard gate — candidates swallowing more than half of
`U_struct` are still rejected outright.

## Rationale

### The defect

`rho_all` is half the longest edge of the Euclidean MST over the significant centres. It
measures **how tightly the centres are packed** — which is not what the footprint radius
is for. The footprint's purpose is to delimit the volume that best represents the region
around the selected hotspots; tying its upper bound to centre packing is a category error.

The failure mode is mechanical. When significant centres are sequence-adjacent — the
normal case for a contiguous hotspot, e.g. 681, 682, 683, 684 — the longest MST edge is
one Cα–Cα step (~3.8 Å), so `rho_all ≈ 1.9 Å` and `1.25 * rho_all ≈ 2.4 Å`, which is
**below the old 3.0 Å floor**. The domain is empty and the run terminates
`NO_ADMISSIBLE_FOOTPRINT_DOMAIN`.

### Measured incidence (ACMG SF v3.2, 68 genes, 2026-09-11)

Six genes were BLOCKED by exactly this, with near-identical numbers:

| Gene | `rho_max` | Significant centres |
|---|---:|---:|
| SDHD | 2.417 | 2 |
| MAX | 2.416 | 8 |
| TNNT2 | 2.414 | 6 |
| DSG2 | 2.393 | 2 |
| PCSK9 | 2.389 | 2 |
| PMS2 | — | 4 |

In all six the recorded binding constraint was `1.25*rho_all`. That is **6 of the 8
BLOCKED genes in the batch**; the remaining two are unrelated and legitimate
(`WT1` upstream `qc_status=FAIL`; `CYP27A1` `MULTI_FRAGMENT_AFDB_ENTRY`).

The non-fatal form of the same constraint affected the majority of completed genes:

```
r_fp <= 3.5 A  :  18/50 genes (36%)   <- pinned to the old 3.0 floor
r_fp <  5.0 A  :  27/50 genes (54%)
```

The earlier PTPN11 case is the same pathology: its domain came out `[3.0, 3.405]` — a grid
of one point — so `r_fp` was forced to the floor, the footprint covered exactly the 23
significant centres, and the footprint stage produced no information.

### Why 5.0 and 25.0

These are not new constants. `radius_domain` (II.2) already uses exactly this envelope for
`r_hot`:

```yaml
R_FLOOR: 5.0                              # FROZEN (II.2)
R_CEIL: 25.0                              # FROZEN (II.2)
r_cap_rule: "min(R_CEIL, 0.25 * D_max)"   # FROZEN (II.2)
```

The change therefore makes the `r_fp` domain **structurally identical in form** to the
`r_hot` domain rather than introducing a novel rule.

**This does not weaken F10.** `bound_below_by_r_hot` remains `false` and the guard in
`params.py` is untouched. The 5.0 Å floor is a fixed constant that happens to equal the
II.2 constant; it is *not* derived from the selected `r_hot` of any run. `r_hot` is still
not an argument to any function in `footprint/domain.py`.

## Admissibility verification (pre-implementation)

The risk of raising the floor is emptying the domain from the other direction. Verified
over all 68 batch genes that `min(0.25 * D_max, 25.0) >= 5.0` holds in every case:

```
smallest effective ceiling:  TP53RK  17.83 A   (D_max = 71.3 A)
genes failing the 5.0 floor: 0
```

No gene loses an admissible domain under this change.

## Expected effect

- **6 genes unblocked** — SDHD, MAX, TNNT2, DSG2, PCSK9, PMS2
- **18 genes released** from a floor-pinned `r_fp`
- Every gene's `config_sha256` changes, therefore every `run_id` changes. This is correct
  and intended: results produced under the old domain are **not comparable** to results
  produced under the new one and must not be pooled.

## Cross-gene comparability

**All genes must be re-run under the new domain for results to be comparable.** Every
completed run in the tree predating this decision is superseded for cross-gene comparison.
Individual pre-decision runs remain valid as provenance records of what was computed under
the old methodology and are not deleted.

Re-running is required to confirm two things this decision does *not* assume:

### 1. KCNA2 — MEASURED 2026-09-12: `r_fp` MOVES, 10.0 -> 10.5 A

The pre-implementation expectation recorded here was **that `r_fp = 10 A` would survive.
It does not.** That expectation rested on a replay of KCNA2's own `r_fp_scan.tsv` with the
sub-5.0 A candidates dropped, which kept 10.0 A winning by ~12% (0.706 vs 0.790). The
replay could not simulate the *upper* half of the change — the 12.5–25 A candidates were
never scored in that run — and those are exactly what moves the winner.

**Isolation experiment.** KCNA2 was re-run with `--run-id 20260901T230055Z_4ec9d6a8_eac88f86`,
pinning it to the pre-decision run's identity so that `seed(ctx) =
sha256(MASTER_SEED|run_id|ctx)` reproduces every draw. Verified: **73 of 73 derived seeds
identical**, and Stage B came out bit-identical.

> **CORRECTED 2026-09-13** — three factual errors in the first version of this section,
> found by re-deriving every number from the archived artefacts.
>
> **(a) The old domain was `[3.0, 12.087]`, not `[3, ~4.8]`.** 19 grid points, binding
> constraint `1.25*rho_all = 12.087`. The `~4.8` figure was the synthetic `clustered`
> fixture's old ceiling and does not belong here.
>
> **(b) The inputs were not byte-identical.** `handoff_01` differs in `clinvar_release`
> (Mon 31 Aug 2026 -> Sun 06 Sep 2026) and `uniprot_entry_version` (226 -> 227). The
> derived cohort sizes and the whole of Stage B are nonetheless identical, so the
> conclusion stands — but "the domain is therefore the *only* variable" was too strong
> and is withdrawn.
>
> **(c) Widening the domain perturbs the voxel geometry of every candidate.** At the same
> radius, with identical seeds and an identical residue set, the mesh measurements moved:
>
> | r_fp | old volume (A^3) | new volume | covered residues |
> |---:|---:|---:|---:|
> | 5.0 | 12660.2 | 12584.2 | 75 = 75 |
> | 10.0 | 43816.4 | 43573.1 | 143 = 143 |
> | 10.5 | 48246.1 | 48015.0 | 150 = 150 |
> | 12.0 | 61873.5 | 61640.2 | 176 = 176 |
>
> **0 of 15 shared radii are bit-identical.** No voxel or padding constant was edited:
> diffing the archived `config_base.yaml` against the current config shows the only
> differences are this decision's own keys plus `near_tie_major_threshold`. The occupancy
> grid is sized at run time from the domain's extent, so raising `hard_ceiling_A` to 25
> grows the bounding box, moves the grid origin, and shifts every distance-transform
> measurement. This is a **defect in its own right**, recorded below; this decision
> exposed it rather than created it.

| | old domain `[3.0, 12.087]` | new domain `[5.0, 25.0]` |
|---|---:|---:|
| grid points | 19 | 41 |
| `r_hot` | 10.0 | 10.0 |
| significant centres | 50 | 50 |
| hotspot regions | 1 | 1 |
| covered residues | 143 | 143 |
| **`r_fp`** | **10.0** | **10.5** |
| footprint residues | 143 | 150 |

**Cause — min-max re-scaling, and the magnitudes support it.** In the old domain
`r_fp = 10.0` beat `10.5` by **6.449%** (0.750971 vs 0.799400) — a decisive win, not a
coin flip. In the new domain `10.5` wins by 1.431% and `10.0` is no longer in the top
four. The voxel shift in (c) is ~0.55% in volume and cannot reverse a 6.4% margin; the
added 12.5-25 A candidates change the min and max of every objective, and min-max
normalisation runs over the admissible set, so every normalised value and therefore every
distance moves.

**Under DECISION-NEAR-TIE-0001 the new selection is flagged ADVISORY (1.431%)** — the new
`r_fp = 10.5 A` is itself not a robust optimum, and that is now reported.

**Robustness of the new answer.** The unpinned re-run (new `run_id`, hence different seeds,
hence `r_hot = 9.5` via KCNA2's known 9.5-vs-10.0 coin flip and 46 centres in 2 regions)
*also* selected `r_fp = 10.5 A`. So 10.5 is the new domain's answer under two independent
seed sets and two different centre sets.

**Consequence for the manuscript.** `CODE_AVAILABILITY.md` and any text citing
`r_fp = 10 A` for KCNA2 describe the pre-decision methodology. Under this decision the
value is **10.5 A** with 150 footprint residues. Those documents must be updated or must
state explicitly which domain they were computed under. This is not a defect — it is the
re-run requirement in this section doing its job.

### 2. The six unblocked genes must reach a terminal state, not merely get past II.9

**Verified 2026-09-12** on `SDHD` and `MAX` (both previously BLOCKED by this exact
constraint); both now run to `exit=0` with a selected footprint:

| Gene | old `rho_max` | new domain | binding | `r_fp` | footprint residues |
|---|---:|---|---|---:|---:|
| SDHD | 2.417 (empty) | `[5.0, 22.11]` | `0.25*D_max` | 9.5 | 18 |
| MAX | 2.416 (empty) | `[5.0, 25.0]` | `hard_ceiling` | 10.5 | 20 |

In both the candidate set is `{0.25*D_max, hard_ceiling}` — the `1.25*rho_all` term is gone
from the bound, as intended. `TNNT2`, `DSG2`, `PCSK9` and `PMS2` were not re-run
individually; they share the identical failure signature and are expected to behave the
same, to be confirmed by the full batch re-run.

## Defect exposed, not caused, by this decision

**A candidate's voxel geometry depends on the extent of its domain** — correction (c)
above. The occupancy grid's bounding box is derived at run time from the domain, so the
volume, surface area and hence `compactness` of a radius are not intrinsic to that radius:
they shift when unrelated candidates enter or leave. Two consequences.

1. `compactness` and `connectivity` are not comparable across runs with different domains,
   even for the same radius on the same protein.
2. Any future change to `rho_min_A`, `hard_ceiling_A` or `dmax_fraction_cap` silently
   perturbs every candidate's geometry, not merely the candidate set.

The magnitude here is small (~0.55% in volume; the residue-level footprint was identical at
every shared radius). **Fixing it is not authorized by this decision.** The fix would be to
derive the grid from a domain-independent reference — the protein bounding box plus a fixed
pad — and needs its own decision record plus a re-run to quantify.


## Scope limits — what this decision does NOT authorize

- No change to `QC_F1_max_coverage` (0.50), `QC_F2_max_bridging_index`, or any QC rule
- No change to the four r_fp objectives (`compactness`, `connectivity`, `stability`,
  `parsimony`) or their equal weighting — the II.10 guard in `params.py` stands
- No change to `step_fp_A` (0.5), `dmax_fraction_cap` (0.25), or `bound_below_by_r_hot`
- No change to merge-scale computation or merge-event reporting
- No label-aware objective of any kind (see `hotspot3d-secenek-C-tasarim.md`)

## Test suite after implementation

`tests/footprint`, `tests/data_structure`, `tests/annotation`, `tests/robustness`,
`tests/regression`: pass. `tests/integration`: **1 failed, 134 passed, 6 errors**.

All seven remaining failures are **pre-existing and unrelated to this decision** — they are
the stale `B = 10_000` assertions left behind by DECISION-B-DEFAULT-0001:

- `test_config_provenance.py::test_d_production_config_is_unchanged_and_unwritten` asserts
  `permutation.B_default == 10000`; `HEAD:config/pipeline.yaml` already ships `100000`
- `test_fixture_contract.py::test_permutation_resolution_case_makes_the_floor_reachable`
  asserts `m > q*(B+1)`, i.e. `620 > 5000` — unreachable at `B = 100000`
- `test_underpowered_terminal_state.py` (6 errors) — its module fixture cannot construct a
  genuinely underpowered run at `B = 100000`

This diff touches no `permutation.*` or `robustness.N_CAP` key. Fixing these three is a
separate item and is **not authorized by this decision**.

An earlier run of the suite reported 52 failures. That was an **incomplete application of
this decision**, not a defect in it: `config/pipeline.yaml` had been updated but
`FROZEN_ASSERTIONS` in `utils/config.py` had not, so `assert_frozen_methodology` blocked
every pipeline invocation at run start. The gate behaved exactly as designed.

## Files changed

- `config/pipeline.yaml` — `footprint_domain` block
- `src/hotspot3d/utils/config.py` — `FROZEN_ASSERTIONS`: `rho_min_A` 3.0 -> 5.0 and
  `hard_ceiling_A` 20.0 -> 25.0 repinned; `post_merge_factor_caps_rho_max` newly pinned
  at `False` so the defect cannot be reintroduced by a silent config edit. This registry
  is the run-start gate (`assert_frozen_methodology`) and is enforced before any stage
  executes — it refused the first verification run, which is the gate working correctly.
- `src/hotspot3d/footprint/params.py` — new `post_merge_factor_caps_rho_max` key
- `src/hotspot3d/footprint/domain.py` — module docstring; MST candidate gated;
  `hard_ceiling_20A` candidate key renamed to `hard_ceiling` (the old name hard-coded a
  value that this decision changes)
- `tests/footprint/test_footprint_domain.py`, `tests/data_structure/test_structure.py` —
  updated to the new envelope
