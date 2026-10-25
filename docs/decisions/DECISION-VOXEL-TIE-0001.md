# DECISION-VOXEL-TIE-0001

**Status:** IN PROGRESS
**Scope:** `GridSpec.index_of` only. No statistic, threshold, null, cohort or radius.
**Supersedes:** weakness 7 of `DECISION-RFP-CONTAINMENT-0001`, which recorded this and
deliberately left it open.

## The defect

`index_of` assigns a point to a cell with `floor((p - origin) / h)`. That function is
discontinuous exactly at a cell boundary, and in float64 a quotient that is an integer
in exact arithmetic can land a few ULP below it. The same physical point then occupies
different cells under different grid origins.

Found while tracing why PTPN11 was the only gene whose volume moved under the padding
repair (`DECISION-RFP-CONTAINMENT-0001`):

| | |
|---|---|
| PTPN11 center 507, y | `-0.104` |
| origin y = `-27.104` | `(p-o)/h` = `54.0` -> `floor` = **54** |
| origin y = `-8.104` | `(p-o)/h` = `15.999999999999998` -> `floor` = **15** |

The origins differ by exactly 19.0 A = 38 cells, so the two indices should differ by
exactly 38; 54 - 38 = 16, not 15. One seed moved one cell, the EDT reclassified a shell
of 434 voxels at the `r_fp` isosurface, and the reported volume changed by 32 A^3
(260 net voxels x 0.5^3 = 32.5) for a purely arithmetic reason.

## Prevalence, measured on the batch

Across 50 genes and 3894 significant centres, **170 axis-level boundary ties**. They are
not one population:

| | count | behaviour |
|---|---|---|
| bounding-box extremes | 150 (88 %) | **stable.** `pad` is always a multiple of `h`, so for the extremal residue `(p - origin)/h = pad/h` is an exact integer under any pad. These never flip. |
| coincidental | 20 (12 %) | **the risk.** Spread over 13 of 50 genes; at most 5 in one gene (HNF1A, out of 187 centres). |

## Decision

Resolve the tie deterministically: a point lying exactly on a cell boundary belongs to
the **upper** cell. Implemented as a tolerance of `1e-9` in index units added before the
floor.

**Why that tolerance is safe, with a number.** Coordinates arrive at 1e-3 A precision
and the cell edge is 0.5 A, so the smallest distinction that can be *real* in index
units is 2e-3 — five hundred thousand times the tolerance. It can only move a point the
arithmetic placed a few ULP below an integer it sits exactly on; it cannot reassign a
point with a genuine sub-cell offset. Both halves of that are asserted by test.

**Rejected alternative: snapping the origin to a multiple of `h`.** It would also make
assignment origin-independent, but it moves the origin itself and therefore every cell
position in the grid — a far larger perturbation to every reported descriptor than the
defect it repairs. The tolerance touches only the points that are actually tied.

## Measured impact

Simulated both ways on 14 archived runs across 7 genes, at each run's own `r_fp`:

| | result |
|---|---|
| residue set | **identical in 14 of 14** |
| component count | **identical in 14 of 14** |
| volume / surface | move in 2 genes only: FGFR3 -0.16 % / -0.12 %, HNF1A -0.18 % / -0.21 % and -0.20 % / -0.22 % |
| compactness | at most +0.09 % |
| every other run | 0.00 % |

**This retracts an estimate made in the same session.** The recommendation against doing
this change argued it would "reassign every tied centre in every gene, so every gene's
reported geometry moves". That is wrong: the 150 structural ties are already stable, so
the repair only touches the 20 coincidental ones, and most of those do not straddle the
`rho` isosurface. The measured cost is two genes moving by under a quarter of a percent,
and no residue or component anywhere.

`r_fp` and footprint membership were never exposed: containment is computed on C-alpha
distances and never touches the grid.

## Tests

`tests/footprint/test_footprint_geometry_core.py`:

* `test_a_point_on_a_cell_boundary_lands_in_the_same_cell_under_any_origin` — shifts the
  origin by 0, 8, 19, 27 and 100 A (all exact multiples of `h`) and requires every index
  to shift by exactly that many cells. Verified to FAIL without the tolerance, with the
  off-by-one-cell error the defect produces.
* `test_the_tie_tolerance_cannot_move_a_genuine_sub_cell_offset` — a point 1e-3 A below a
  boundary stays in the lower cell; a point exactly on it takes the upper.

## What this does not fix

Comparability across configurations is still governed by `config_sha256` and
`code_version`, as before. This repair removes one arithmetic source of drift in the
voxel-derived descriptors; it makes no claim that runs under different digests are
comparable.
