# FINDING-B1-DISORDER-CONFOUND-0001 — the positional null is confounded in six genes

**Status:** MEASURED. No change applied. Acting on this is a methodology decision
requiring a re-run, and is not authorized.
**Measured:** 2026-10-30, read-only, against `results_v3/`.
**Closes:** B1 of the 2026-09-13 audit, which had been open and untested through every
subsequent session.

## The mechanism

AlphaFold places intrinsically disordered regions as extended, low-pLDDT sprawl; the
folded core is compact. Pathogenic missense variants concentrate in the folded core,
because that is where they disrupt structure and function. The structure-aware positional
null, however, redraws `N_P` residues **uniformly over all of `U_struct`** — core and
sprawl alike.

So the observed cohort and its own null can be drawn from different structural regimes.
When they are, observed P/LP pairs are closer together than null pairs for a reason that
has nothing to do with clustering: they are all inside the compact part and the null is
not. That inflates `g(r)` across the whole width of the folded core rather than at one
radius, which is exactly what a single wide elevated block looks like.

The premise holds in this batch. A median of **93.3 %** of P/LP residues sit at
pLDDT ≥ 70. In `BRCA1`, 80.4 % of `U_struct` is below pLDDT 50 and mean `U_struct` pLDDT
is 41.6, while P/LP mean 89.6 and 93.8 % of them are above 70.

## The test

Null against null. For each gene, draw `N_P` residues **uniformly from the folded subset
only** (pLDDT ≥ 70) — a set containing **no clustering whatever** — and run the same pair
correlation against the same positional null over `U_struct`. Five independent draws per
gene, diagnostic `B = 600` (the quantity being located is the elevated band, not a
p-value).

Two measurements per gene:

- `b1_explained` — the confound's maximum excess over the null, as a fraction of the
  observed maximum excess
- `b1_overlap` — the fraction of the observed elevated radii that the confound also
  elevates

**The ratio alone is not a finding.** In genes whose observed excess is tiny it explodes:
`PALB2` reaches 40.96 with `b1_overlap = 0.00` and `n_E = 1`. A confound that elevates
*different* radii is not an explanation of the observed elevation. The ratio is therefore
gated on overlap throughout.

## Result

`CONFOUNDED` — `b1_explained ≥ 0.50` **and** `b1_overlap ≥ 0.90`. The disorder contrast
alone reproduces both where the elevation sits and most of how large it is:

| gene | regime | centres | explained | overlap |
|---|---|---:|---:|---:|
| `HNF1A` | WORKS — CAVEAT REQUIRED | 187 | 1.16 | 0.97 |
| `LDLR` | REFUSED — coverage | 346 | 0.76 | 1.00 |
| `BRCA1` | WORKS, tier B | 226 | 0.75 | 0.98 |
| `LMNA` | WORKS | 6 | 0.60 | 0.92 |
| `MYH7` | WORKS | 25 | 0.58 | 0.92 |
| `KCNQ1` | REFUSED — coverage | 217 | 0.54 | 0.95 |

`HNF1A`'s confound produces **more** elevation than was observed, and the 2026-09-13
audit's separate observation — that a pLDDT ≥ 70 restriction collapses `HNF1A` from 187
centres to 0 — is the same phenomenon seen from the other side.

`BAND_CONFOUNDED` — overlap ≥ 0.90 but under half the magnitude: `KCNH2` (1.00 / 0.32)
and `TP53` (1.00 / 0.37). The confound fixes *where* the elevated band falls, and
therefore the derived radius, even where it does not account for the signal's size.

`CLEAN` — overlap < 0.50 and explained < 0.30: **25 genes**, including `KCNA2`, the case
gene, at 0.29 / 0.13; `TSC1` at 0.09 / 0.13; `MSH2`, the cleanest in the batch, at
0.05 / 0.05; and `SMAD4`, `STK11`, `PMS2`, `TGFBR1`, `TGFBR2`.

`PARTIAL` 22 genes. `NOT_TESTED` 15 — no centres, or fewer than three P/LP residues.

`COL3A1` could not be tested: its folded subset is **smaller than its `N_P` of 306**. A
collagen triple helix carries low pLDDT along its whole length, so there is no folded core
to draw from. That is also why it found zero centres, and it means the gene's negative
result is not vulnerable to this confound.

## What this corrects

Two statements in earlier reporting were wrong and are retracted here:

- `BRCA1` was presented as needing no caveat — it meets every one of the four quality
  criteria, 226 centres at coverage 0.106, tier B. It is nonetheless **not quotable as
  clustering evidence**: 75 % of its excess is reproduced by the confound at 0.98 overlap.
  The quality criteria cannot see this failure mode, which is the point of recording it as
  a separate column rather than folding it into the tier.
- `MYH7` was in `works` with no qualification. It is compromised twice over: its near-tie
  decides 25 against 73 centres on a 0.82 % distance margin, and 58 % of its elevation is
  confound at 0.92 overlap.

## What was NOT done, and why

This is a **diagnostic, not a filter.** No residue is excluded for low pLDDT, no cohort is
altered, and no result in `results_v3/` is rewritten. The user cancelled pLDDT-based
exclusion in an earlier session and nothing here reverses that: the folded subset is used
only as a *source of draws for a synthetic negative control*, never as an admission rule.

The repair, if authorized, is at the null and not at the cohort: stratify the positional
null so it draws from the same structural regime as the observed cohort — for example
within pLDDT bands, or within annotated domains — so that compactness contrast cancels
instead of registering as signal. That changes `config`, every number in the batch, and
requires a re-run. It should be decided together with the other open methodology items
(block-selection key, the `>` vs `>=` E criterion, near-tie severity,
`BH_BY_DISAGREEMENT`, `CASQ2`'s `r_contain = 0`) so that one re-run settles all of them.

## Where the numbers live

- `results_v3/BATCH_DATA_ACMG_SF_v3.2_v3_classified.tsv` — columns `b1_explained`,
  `b1_overlap`, `b1_status`
- `results_v3/BATCH_REPORT_ACMG_SF_v3.2_v3.md` — section "The measurement itself"
