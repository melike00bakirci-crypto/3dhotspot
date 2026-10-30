# DECISION-HIGH-COVERAGE-GENES-0001

**Status:** DECIDED (user, 2026-10-25)
**Scope:** how a class of results is PRESENTED. No code, no statistic, no threshold.
**Decision:** the 14 high-coverage genes are presented as findings.

## The observation that raised the question

In the `results_v2/` batch, 14 of the 51 genes that produced significant centres did so
from spheres covering more than half the structure, at fold enrichment ~1.0:

| gene | `r_hot` (A) | coverage | fold enrichment | centres |
|---|---|---|---|---|
| LDLR | 22.5 | 0.660 | 1.043 | 346 |
| GAA | 22.0 | 0.872 | 1.016 | 334 |
| MLH1 | 20.0 | 0.729 | 1.081 | 248 |
| MEN1 | 20.0 | 0.849 | 1.032 | 230 |
| KCNQ1 | 23.0 | 0.536 | — | 218 |
| HNF1A | 23.5 | 0.509 | — | 187 |
| PTEN | 16.0 | 0.754 | 1.032 | 166 |
| TP53 | 15.5 | 0.555 | — | 157 |
| ACVRL1 | 12.0 | 0.644 | 1.029 | 137 |
| ABCD1 | 16.5 | 0.644 | 1.149 | 132 |
| STK11 | 16.5 | 0.642 | 1.311 | 100 |
| RPE65 | 17.5 | 0.872 | 0.999 | 88 |
| SDHB | 13.5 | 0.621 | 1.112 | 49 |
| OTC | 11.5 | 0.576 | — | 43 |

The batch's largest centre counts are concentrated in this group. The obvious reviewer
question follows: at coverage 0.87, what does "hotspot" mean?

## Why this is not a contradiction

The primary test and the enrichment metric answer different questions, and the code says
so rather than conflating them (`hotspot/scan.py`: fold enrichment is "deliberately kept
apart from significance"):

* the **structure-aware positional null** asks whether P/LP residues are placed more
  tightly than uniform placement over `U_struct`, conditioned on `N_P`;
* **fold enrichment** asks whether a sphere holds a higher P/LP share than the labelled
  background does.

A gene can satisfy the first while scoring ~1.0 on the second, and both facts are true
simultaneously. Significance here is a claim about positional clustering, not about
pathogenic-versus-benign enrichment, and the manuscript's wording must say that — a
sentence implying enrichment would not be supported by this test.

Nor is it a detection failure: **all 14 carry a MAJOR `EXCESSIVE_STRUCTURAL_COVERAGE`**
and 11 also carry `HIGH_COHORT_ABSORPTION`. The pipeline identifies the condition
without exception; `DECISION-REPORTING-CONTEXT-0001` makes the coverage and enrichment
travel with the count so the context cannot be separated from the number again.

## What supports presenting them as findings

The clustering signal is **not created** by the saturating radius. At the largest
admissible radius whose coverage stays below 0.30, every one of the 14 still has
significant centres:

| gene | at `r_hot` | at coverage < 0.30 |
|---|---|---|
| LDLR | 346 @ 0.660 | **56** @ 9.5 A, 0.277 |
| GAA | 334 @ 0.872 | **65** @ 9.0 A, 0.270 |
| HNF1A | 187 @ 0.509 | **100** @ 9.5 A, 0.284 |
| KCNQ1 | 218 @ 0.536 | **75** @ 8.0 A, 0.294 |
| TP53 | 157 @ 0.555 | **56** @ 7.0 A, 0.288 |
| MLH1 | 248 @ 0.729 | 42 @ 9.5 A, 0.233 |
| ABCD1 | 132 @ 0.644 | 37 @ 10.0 A, 0.239 |
| ACVRL1 | 137 @ 0.644 | 34 @ 7.5 A, 0.284 |
| PTEN | 166 @ 0.754 | 29 @ 7.0 A, 0.213 |
| STK11 | 100 @ 0.642 | 28 @ 9.0 A, 0.261 |
| MEN1 | 230 @ 0.849 | 26 @ 8.5 A, 0.221 |
| RPE65 | 88 @ 0.872 | 16 @ 10.5 A, 0.263 |
| SDHB | 49 @ 0.621 | **2** @ 9.0 A, 0.118 |
| OTC | 43 @ 0.576 | **5** @ 8.5 A, 0.127 |

So for 12 of the 14 the finding survives a conservative coverage bound with a
substantial centre count; the large radius grows the result rather than producing it.

**SDHB and OTC are the exceptions and must be qualified separately.** They collapse to
2 and 5 centres once coverage is held below 0.30, so in those two the reported result
does rest on saturation. They are not to be presented in the same sentence as the other
twelve.

## Mandatory reporting for this class

1. Coverage at `r_hot` and fold enrichment at `r_hot`, beside the centre count. Already
   implemented (`DECISION-REPORTING-CONTEXT-0001`); they ride on `handoff_02` and in the
   run-summary metrics.
2. The compact-radius sensitivity above. **No code is required**: `r_hot_scan.tsv`
   already carries `coverage_fraction` and `n_significant_centers` at every scanned
   radius, and that file is deposited — a reviewer can reproduce the table, and it is
   the `r_hot` counterpart of what `sensitivity_band.tsv` does for `r_fp`.
3. SDHB and OTC qualified as resting on saturation.
4. The MAJOR `EXCESSIVE_STRUCTURAL_COVERAGE` for each of the 14, not suppressed.

**Precision of the compact-radius counts.** They are the scan's PROVISIONAL `S(r)`,
which comes from a different permutation draw than the final detection (II.12 / v2
§5.5). Measured across these 14 at `r_hot`, provisional versus final differs by at most
**1 centre** (mean absolute difference 0.3), so the table is reliable to +/-1 and must
be described that way rather than as the final count.

## What is NOT claimed

* Not that the spheres are enriched for pathogenic variants relative to benign ones —
  fold enrichment is ~1.0 and is reported as such.
* Not that `r_hot` is the biologically correct scale for these genes. The radius was
  selected by the frozen multi-objective procedure; for RPE65 the pair-correlation
  interval was `[8.0, 11.0]` and the fallback envelope `[5.0, 17.86]` was used instead,
  with the selection landing at 17.5. That is recorded, not defended.
* Nothing about the 6 genes in the overlap group that produced no centres at all.

---

## Addendum, 2026-10-30 — the regime is renamed, and a seventh defect is recorded

The user objected to the label, and the objection was right. Calling these genes
"REFUSED — coverage" read as a failure. It is not one: **the hotspot centres stand and
are the finding.** What QC-F1 withholds is only the footprint *geometry*. Verified across
all 11: centres retained 11/11, footprint produced 0/11.

The regime is therefore renamed `refused_excessive_coverage` →
**`diffuse_no_localised_region`**, and the per-gene caveat now reads "pathogenicity is
diffuse, not localised … Hotspot centres STAND and are a finding; only the footprint
geometry is withheld". No code under `src/` changed; this is the reporting layer.

### Why withholding the geometry is nonetheless correct here

The user's concern was that in some genes pathogenic variants genuinely are spread across
the whole protein while benign ones are not, and that discarding such a gene would throw
away a real result. That case would show as fold enrichment above 1. It does not occur in
any of the 11:

| gene | centres | coverage @ r_hot | fold |
|---|---:|---:|---:|
| `GAA` | 340 | 0.882 | 1.02 |
| `OTC` | 137 | 0.915 | 1.01 |
| `RPE65` | 85 | 0.872 | 1.00 |
| `MEN1` | 233 | 0.849 | 1.03 |
| `LDLR` | 346 | 0.664 | 1.04 |

The full range across the 11 is **1.00–1.16**. At fold ≈ 1.0 the pathogenic density inside
the sphere equals the density outside it: the sphere is not distinguishing anything, it is
merely large enough to contain everything. Reporting a region that covers 88 % of the
protein adds nothing to the statement "the variants are spread", which the centres already
carry. So nothing was lost in these 11 — but that is an empirical fact about this batch,
not a property of the rule.

### Defect 7, now recorded

**QC-F1 decides on coverage alone and never consults fold enrichment.** A gene at
coverage > 0.50 **and** fold > 1.5 — spread *and* enriched — would have its footprint
withheld wrongly. That has not happened in 70 genes, so the rule has not yet cost a
finding, but it cannot distinguish the two cases and that is a real gap.

The principled repair is to make the withholding conditional: withhold only when
`coverage > 0.50 AND fold < 1.5`. It changes `config`, every number in the batch, and
needs a re-run, so it joins the other six open methodology items rather than being
patched. It is listed as item 7 in the batch report's defect table.
