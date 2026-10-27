# DECISION-MEFV-RESULT-0001

**Status:** RECORDED
**Scope:** how MEFV's result is to be read and cited. No code, no methodology change.
**Run:** `results_v3/MEFV/MEFV_20261027T211825Z_ff721543_ac7cce0f`

MEFV was run on request after the 68-gene ACMG SF v3.2 batch, under the same
`config_sha256 ff721543` / `code_version ac7cce0f`, so it is directly comparable with
those 68 and with KCNA2. Exit 0, 14 minutes.

## The result

**`COMPLETED_NEGATIVE / SCIENTIFIC_NEGATIVE`** — no significant hotspot centre.

| | |
|---|---|
| UniProt / MANE | O15553 / `ENST00000219596.6` ↔ `NM_000243.3` |
| `M` (U_struct) | 781 |
| `N` / `N_P` / `N_B` | 19 / 8 / 11, no conflicts |
| P/LP residues | 111, 242, 244, 653, 680, 694, 726, 761 |
| `r_hot` | 16.0 A (runner-up 19.5, margin 19 % — not a near-tie) |
| Significant centres | **0** |
| Global Ripley's K, P/LP | `T_obs = 2.517`, `p_global = 0.1842` — not significant |
| Global Ripley's K, B/LB | `T_obs = 1.252`, `p_global = 0.7739` — not significant |
| Smallest observed per-centre p | 0.00727 against a rank-1 bar of 0.000139 |
| `B` / `q` | 100,000 / 0.05, BH |

## This is a legitimate negative, not a failure to compute

The distinction matters and MEFV sits on the informative side of it. A run that *cannot*
compute terminates `UNDERPOWERED` (the attainable p-floor exceeds `q/m` regardless of the
data — DSC2, TMEM43, SDHAF2 and PLN in the batch) or `BLOCKED` (inputs unusable —
CYP27A1, TSC1). MEFV is neither. Its power certificate passes at both instances:

```
preflight   PASSES=True   TEST_CANNOT_REJECT=False   m=41    p_floor 1.0e-05 < c_1 1.22e-03
posthoc     PASSES=True   TEST_CANNOT_REJECT=False   m=359   p_floor 1.0e-05 < c_1 1.39e-04
```

The test had the resolution to reject and did not. The smallest per-centre p it found,
0.00727, is **52 times** the rank-1 critical value — not a marginal miss.

## The geometry, for the record

C-alpha pairwise distances among the eight P/LP residues: median 45.4 A, min 5.9, max
99.8. Nearest same-class neighbour per residue:

| residue | nearest P/LP |
|---|---|
| 242, 244 | 5.9 A (an adjacent pair) |
| 694, 761 | 12.4 A |
| 653, 726 | 12.5 A |
| 680 | 15.9 A |
| **111** | **60.8 A** |

So MEFV carries one tight pair, a loose group of five in the C-terminal region, and one
residue 60 A from any other. No single sphere concentrates enough of the eight to survive
BH across 359 candidate centres, which is what the per-centre result says.

**A correction belongs here.** When this run was first reported in conversation, a quick
median-pairwise-distance resampling was quoted as putting the eight at the 6th percentile
of random draws — "tighter than chance, just outside 5 %". That statistic is cruder than
the pipeline's and overstated the signal. The pipeline's own global test, `max_r |Z_K(r)|`
against the structure-aware positional null, gives **p = 0.184**. That is the number to
cite; the percentile figure should not be repeated.

## Reconciliation with the 22-residue INFEVERS analysis

An earlier side analysis of MEFV, reviewed in session on 2026-09-12, worked from a
22-residue P/LP cohort and found global clustering at `p = 0.0013`. That does not
contradict this run. The two use different evidence bases:

* **this run** admits ClinVar only, under the frozen binary significance policy
  (`compound_significance_policy: whole_field`, no review-star filtering). Of MEFV's
  1410 classified ClinVar records, **23 (1.6 %)** are P/LP; 553 are B/LB, 531 VUS and 270
  conflicting. The 23 records collapse to 8 residues.
* **the side analysis** drew on INFEVERS curation. The user ruled on 2026-09-12 that
  INFEVERS will never enter the pipeline, so that cohort is out of scope by decision, not
  by oversight.

Same gene, two evidence bases, two answers, both correct for what they measure. The
pipeline reports what ClinVar supports. **Neither supersedes the other, and the two
numbers must never be presented as one result.**

The limiting factor is therefore the cohort, not the method: 98 % of MEFV's ClinVar
content falls outside the binary policy, most of it VUS and conflicting calls, which is
the expected shape for a gene whose variants are common in the population. What would
change this result is more classified pathogenic variants, not a different statistic.

## Also recorded: two claims made and withdrawn in session

Both concerned MEFV and both were wrong; they are written down so they do not resurface.

* **A "disorder confound halving MEFV's clustering evidence" was asserted and
  retracted.** Decomposition showed the cohort change 22 -> 16 moved the global p from
  0.0015 to 0.0016 — nothing. Only shrinking the universe moved it, and a stratified null
  preserving 16/6 gave p = 0.0013. The confound claim had no support.
* **"7 of 11 out-of-domain P/LP residues sit near a centre" was wrong; the pipeline's
  answer is 2 of 11.** The error came from using a 25-centre list that included residue
  245 at pLDDT 34. The user's own table confirmed 2 of 11.

Separately, and on the user's instruction: pLDDT sensitivity was cancelled as a filter —
no residue is ever excluded for low pLDDT (F2). pLDDT is recorded, never filtering.

## Not claimed

* Nothing about whether MEFV *has* a pathogenic hotspot. A negative under this cohort is
  a statement about the evidence available, not about the protein.
* No comparison with the 68-gene batch's positive genes is implied by the shared
  `config_sha256`; comparability means the runs are reproducible against one another, not
  that their outcomes are commensurate.
* One reporting gap noticed and not fixed: `handoff_02` carries
  `global_clustering_plp_p: null` while the value, 0.1842, is present in
  `04_GLOBAL_CLUSTERING/ripleys_k_summary.json`. The handoff drops a number it has.
