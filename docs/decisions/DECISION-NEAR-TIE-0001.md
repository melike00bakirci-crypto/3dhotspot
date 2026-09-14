# DECISION-NEAR-TIE-0001

```
DECISION_ID    = DECISION-NEAR-TIE-0001
DATE           = 2026-09-13
AUTHORIZED_BY  = user
STATUS         = IMPLEMENTED
SCOPE          = near-tie reporting for BOTH r_hot (II.6) and r_fp (II.10).
                 Reporting only — no selection changes.
```

## Decision

Two changes to the near-tie flag. **Neither changes which radius is selected.**

### 1. The radius-separation clause is removed

Before:

```python
if rel < near_tie_rel and abs(ranked[0] - ranked[1]) > near_tie_sep_steps * step:
```

The second clause suppressed the flag whenever the top two candidates were **adjacent**
on the 0.5 A grid, on the premise that neighbouring radii give similar answers.

**KCNA2 refutes the premise.** Its top two candidates are 10.0 and 9.5 A — adjacent, so
the flag never fired — and the two are 0.178% apart in distance-to-ideal:

| | r_hot = 9.5 A | r_hot = 10.0 A |
|---|---:|---:|
| significant centres | 46 | 50 |
| hotspot regions | **2** | **1** |

"One contiguous hotspot region" versus "two regions" was decided by 0.178%. Adjacency is
precisely where a tie is *least* visible in the scan table and therefore *most* worth
flagging: a tie between distant radii is obvious to anyone reading the scan, a tie between
neighbours is not.

`near_tie_radius_separation_steps` is retained in config and in the `select()` signature
for provenance and API stability. It no longer gates anything, and is marked INERT.

### 2. Severity is graded by the margin

| Margin | Severity | Meaning |
|---|---|---|
| `< 1%` | **MAJOR** | the seed, not the data, decides the winner |
| `1% - 5%` | **ADVISORY** | worth knowing, not a coin flip |
| `> 5%` | silent | a real optimum |

The MAJOR tier is not a guess. It was demonstrated in-session: changing
`config_sha256` re-derived every seed (`seed(ctx) = sha256(MASTER_SEED|run_id|ctx)`),
and KCNA2's `r_hot` moved from 10.0 to 9.5 A **on the same gene with the same cohort**,
taking the region count from 1 to 2 with it. At a sub-1% margin the winner is a property
of the permutation draw, not of the protein.

New config key, both selection blocks: `near_tie_major_threshold: 0.01`.

## Measured incidence (69 genes with an r_hot scan)

| Severity | Genes | |
|---|---:|---|
| **MAJOR (<1%)** | **6** | CACNA1S 0.77% · COL3A1 0.57% · HNF1A 0.26% · KCNA2 0.18% · MSH2 0.44% · **TP53 0.53%** |
| ADVISORY (1-5%) | 17 | ABCD1, BRCA1, CACNA1A, ENG, FGFR3, HRAS, KCNH2, LDLR, MYBPC3, MYH7, PRKAG2, RBM20, SDHB, TGFBR1, TPM1, TSC2, VHL |
| silent (>5%) | 44 | |
| not applicable | 2 | PKP2, TMEM127 |

The old rule fired on 8 of 62 genes and on **none** of the six sub-1% cases — it was silent
exactly where the selection was least robust.

**`TP53` is in the MAJOR tier at 0.53%.** It is one of the genes used to validate the
pipeline against known biology, and its `r_hot = 15.5 A` beat the runner-up by half a
percent. That was not reported anywhere before this decision.

## Verified

Replayed through the real `select()` on archived scan tables:

```
KCNA2   selected 10.0   near_tie True   MAJOR      0.178%
TP53    selected 15.5   near_tie True   MAJOR      0.527%
HNF1A   selected 23.5   near_tie True   MAJOR      0.265%
MSH2    selected  8.0   near_tie True   MAJOR      0.436%
BRCA1   selected 12.0   near_tie True   ADVISORY   1.934%
MYBPC3  selected  8.5   near_tie True   ADVISORY   1.695%
PTEN / MEN1 / SCN5A     near_tie False  silent
```

**Every selected radius is unchanged.** The flag is reporting, not selection.

## What this decision does NOT do

- It does **not** change any selected radius, and no archived result is invalidated
- It does **not** change the tie chain, `w_k`, the metric, normalisation, or the utopia point
- It does **not** make the selection depend on the flag. Overriding a selection because it
  was close would be outcome-dependent and is prohibited.
- It does **not** touch `near_tie_relative_threshold` (0.05), which keeps its meaning as
  the ADVISORY boundary

## Files changed

- `src/hotspot3d/utils/multiobjective.py` — separation clause removed; `severity`,
  `major_threshold`, `advisory_threshold` added to `near_tie_detail`; new
  `near_tie_major_rel` parameter
- `config/pipeline.yaml` — `near_tie_major_threshold: 0.01` in `radius_selection` and
  `footprint_selection`; `near_tie_radius_separation_steps` marked INERT
- `src/hotspot3d/hotspot/stage.py` — graded severity passed to the warning; the message
  no longer claims "more than 2 steps apart" (which is now false by construction)
- `src/hotspot3d/footprint/params.py`, `src/hotspot3d/footprint/selection.py` — threshold
  plumbed through
