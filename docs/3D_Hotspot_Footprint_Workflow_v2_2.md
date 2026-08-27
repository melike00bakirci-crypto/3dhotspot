# 3D Hotspot / Footprint Workflow with Functional Mechanism — v2

> **Revision note.** This is a revision of the original workflow specification following the KCNA2 run
> `20260817T155124Z_df8ec59f_f30b0d52`, which terminated with zero significant hotspot centers and was
> reported as a scientific negative. Post-hoc analysis showed the negative was not scientific: under the
> configuration the agent derived from v1, the minimum p-value attainable by the per-center test
> (7.1e-04) exceeded the rank-1 Benjamini-Hochberg critical value (1.1e-04), so no center could have
> been declared significant regardless of the data. The revisions below close the specification gaps
> that permitted this. Text carried over from v1 is unchanged except where marked; all new normative
> requirements are marked **[NEW]** or **[REVISED]**.

---

## 0. Scope and automation contract

The analysis aims to identify regions within the target protein where pathogenic missense variants
exhibit statistically significant three-dimensional spatial clustering and to determine the optimal
sphere radius for hotspot analysis using an objective, reproducible, and fully automated decision
framework. Every stage of the analysis should follow the same standardized workflow for all genes
without requiring user intervention.

**[NEW] Automation boundary.** "Without user intervention" governs the *execution* of the pre-registered
analysis, not the *amendment* of it. The agent must run the entire workflow end to end and must never
pause for confirmation of a step that is specified here. However, the agent must never silently widen a
pre-registered search domain, increase the permutation count, loosen a multiple-testing correction, or
drop a residue in order to obtain a positive result. Where the specification requires such an amendment,
the agent halts, emits the diagnostic that motivates it, and escalates. Escalation is a defined terminal
state of an automated run, not a failure of automation.

**[NEW] Terminal states.** Every run must terminate in exactly one of:

| state | meaning |
|---|---|
| `COMPLETED` | Significant hotspots found; full downstream chain executed. |
| `COMPLETED_NEGATIVE` | No hotspots found **and** the design was demonstrably capable of finding them (see §5.4). |
| `UNDERPOWERED` | No hotspots found **and** the design could not have found them. This is **not** a scientific result about the gene. |
| `BLOCKED` | A structural, data, or QC precondition was not met. |

A run may be reported as `COMPLETED_NEGATIVE` only after the power certificate of §5.4 has been
computed and has passed. Reporting an underpowered run as a negative result is a specification violation.

---

## 1. Variant cohort construction

All missense variants associated with the target gene should be retrieved from the ClinVar database.
Only variants classified as **Pathogenic/Likely Pathogenic (P/LP)** or **Benign/Likely Benign (B/LB)**
should be included in the analysis. Variants classified as **Variants of Uncertain Significance (VUS)**,
**conflicting interpretations**, or records lacking reliable clinical classification should be excluded.
Duplicate records representing the same amino acid position should be collapsed so that the analysis is
performed at the residue level rather than the variant-record level. The resulting P/LP and B/LB
reference cohorts should be comprehensively summarized and reported.

**[NEW] Exclusion accounting must be faithful.** Every excluded record must be reported under the reason
that actually made it ineligible. Records excluded because they are synonymous, nonsense, frameshift, or
indel are *not* parse failures and must not be aggregated under a parse-failure category. The exclusion
ledger must reconcile exactly: `n_retrieved = n_eligible + n_excluded`, with a per-reason breakdown at
the granularity of the true cause.

**[NEW] Compound significance labels.** The policy for compound ClinVar labels
(e.g. `Pathogenic, low penetrance`, `Likely pathogenic|drug response`) must be declared explicitly and
its cost reported: the run must state how many records each policy excluded. A whole-field exact-match
policy is permitted but must be accompanied by the count of records it discarded, so the reader can
judge whether the benign cohort was materially reduced.

**[NEW] Benign cohort adequacy.** The size of the B/LB cohort (`N_B`) determines the resolution of any
test that contrasts pathogenic against benign residues. `N_B` must be reported prominently in the cohort
summary, and §5.4 must be evaluated against it before any negative result is issued.

---

## 2. Structure retrieval and quality control

The AlphaFold structural model of the target protein should be obtained. Before any downstream analysis,
the structure must undergo rigorous quality control to verify model integrity, chain composition, missing
coordinates, consistency of residue numbering with the UniProt sequence, and exact correspondence between
AlphaFold residue positions and ClinVar variant positions.

**[REVISED] Oligomeric state.** Before selecting a structural model, the biological assembly of the target
protein must be determined from UniProt/ComplexPortal annotation. If the protein functions as a homo- or
hetero-oligomer, the monomeric AlphaFold model does not represent the spatial neighbourhoods that variants
actually occupy, because residues contributed by neighbouring subunits are absent. In that case the run
must either use a multimer model of the biological assembly, or, if none is available, halt with
`BLOCKED: OLIGOMERIC_ASSEMBLY_UNAVAILABLE`. Proceeding on a monomer for an obligate oligomer is permitted
only if it is declared as an explicit, reported limitation and the resulting hotspots are labelled
intra-subunit only.

**[REVISED] Low-confidence regions — explicit and binding criteria.** Regions with low structural confidence
must be identified, and the specification now fixes their disposition rather than leaving it to the
implementing agent:

- Residues with `pLDDT < 50` are **excluded from the candidate center universe** `U_center`. Coordinates in
  these regions are not meaningful and their inclusion inflates the multiple-testing family without
  contributing interpretable geometry.
- Residues with `pLDDT < 50` are **retained** in the cohort `L` and in the positional reference universe
  `U_struct` if they carry a classified variant, so that no variant is silently discarded.
- The two universes are therefore distinct and must be reported separately: `|U_struct|` (positional
  reference) and `|U_center|` (test family before further exclusions).
- The number of residues removed from `U_center` by this rule, and the contiguous regions they form, must
  be reported.

After all variants have been accurately mapped onto their corresponding three-dimensional coordinates, the
spatial analyses may proceed.

---

## 3. Global spatial assessment

Before initiating hotspot detection, the overall spatial distribution of variants across the protein should
be evaluated to determine whether significant clustering exists at the global level. To accomplish this,
**Global Ripley's K analysis** should be performed separately for the P/LP and B/LB variant sets. The
objective of this analysis is not to determine the optimal sphere radius but rather to assess whether the
observed spatial distribution deviates from complete spatial randomness and whether significant clustering
exists across the protein.

Subsequently, **Pair Correlation Analysis** should be performed to identify the spatial scales at which
variants tend to cluster. However, peaks identified by the pair correlation function should **not** be
interpreted as the optimal analysis radius. Instead, these results should be used solely to define a
biologically plausible candidate radius range that will be evaluated during the subsequent optimization
process.

**[NEW] The null model used here is normative for the whole workflow.** The null model adopted for the
global analysis — uniform without-replacement placement of the P/LP point set over the positional
universe `U_struct`, conditioned on the observed cohort size ("structure-aware positional null") — is the
reference null for this specification. §5.2 requires the per-center test to use the same null. Any
deviation must be justified in the report and accompanied by the power certificate of §5.4 computed
under the deviating null.

---

## 4. Radius scan and multi-objective selection

A systematic radius scan should then be performed across the candidate radius range. For each candidate
radius, the hotspot analysis should be repeated independently. Each run should calculate **Leave-One-Out
Cross Validation (LOO)** performance, **Matthews Correlation Coefficient (MCC)**, permutation-based
statistical significance, fold enrichment, the number of detected hotspots, hotspot size, spatial
distribution of hotspots, and stability across neighboring radii.

**[REVISED] Objectives versus admissibility constraints.** The quantities listed above are *objectives*
entering the Pareto optimization. They are not admissibility filters. In particular:

- **The number of significant centers at a candidate radius must never be used to eliminate that radius
  from the Pareto set.** A radius that yields zero significant centers is a legitimate Pareto candidate
  that simply scores zero on that objective. Eliminating it converts the radius scan into a search
  conditioned on significance, which is circular: significance then determines the radius, and the radius
  determines significance.
- The only permitted admissibility constraints are those that make a radius *computationally or
  geometrically undefined*, and each must be declared in the configuration with its threshold. A radius may
  be ruled inadmissible only if a stated quantity is undefined (e.g. division by an empty set), never
  because a result was not significant.
- If more than half the scanned radii are ruled inadmissible, the run must emit
  `ADVISORY: ADMISSIBILITY_DOMINATES_SELECTION` and report which constraint did the eliminating.
- **Any objective that is degenerate across the admissible set (identical value at every admissible radius)
  must be reported as non-informative and excluded from the distance-to-utopia computation.** A Pareto set
  of size one, or a selection in which all objectives are degenerate, must emit
  `MAJOR: VACUOUS_PARETO_SELECTION` — the radius was not selected, it was the only survivor.

**[NEW] Biological plausibility as a computed objective.** The specification's requirement that each radius
be evaluated for biological plausibility must be operationalized numerically rather than left as prose.
The following are computed at every candidate radius and enter selection:

| quantity | definition | requirement |
|---|---|---|
| `structural_coverage` | fraction of `U_struct` lying within `r` of any significant center | soft penalty above 0.25; `MAJOR` warning above 0.50 |
| `cohort_absorption` | fraction of the classified cohort `L` contained in the single largest sphere | soft penalty above 0.33 |
| `r_to_domain_ratio` | `r` divided by the radius of gyration of the modelled structure | soft penalty above 0.5 |

A radius whose sphere absorbs most of the cohort does not localize anything and must not be selected merely
because no smaller radius reached significance.

Radius selection should never rely on a single statistical metric. In particular, selecting the radius
solely because it produces the highest MCC value or solely because it corresponds to the principal peak of
the pair correlation function should be considered inappropriate. Instead, radius selection should be
formulated as a **multi-objective optimization problem**.

Accordingly, every candidate radius should be evaluated simultaneously with respect to classification
performance, statistical significance, hotspot enrichment, stability, biological plausibility, and any
additional relevant performance metrics. A **Pareto optimization** procedure should then be applied to
identify all **non-dominated** candidate radii, namely those for which no alternative radius performs
better across all evaluation criteria simultaneously. This approach ensures that the selected candidates
represent optimal trade-offs among multiple competing objectives rather than maximizing a single metric.

In the final selection stage, the Pareto-optimal solution set should be evaluated using predefined,
objective decision rules to determine the final analysis radius. This decision framework should jointly
consider stability across neighboring radii, biological relevance, statistical robustness, and hotspot
reliability. The objective is not simply to identify the radius with the highest numerical performance but
rather to select the most balanced, reliable, and biologically meaningful solution.

**[NEW] Boundary optima.** If the selected radius lies in the outermost band of the candidate domain, the
run emits `MAJOR: BOUNDARY_OPTIMUM`. The agent must not widen the domain automatically. It must, however,
report what the diagnostic implies: a boundary optimum at the upper end combined with a vacuous Pareto set
is evidence that no interior radius was ever competitive, which is a symptom of the power failure in §5.4
rather than a property of the radius domain.

---

## 5. Primary hotspot detection

Using the selected optimal radius, a protein-wide **sliding-sphere hotspot analysis** should be performed.
All detected hotspot regions should be validated using permutation testing and appropriate multiple-testing
correction procedures. Statistically significant hotspots should then be comprehensively characterized and
reported in a format suitable for downstream biological interpretation.

### 5.1 [NEW] Test family

Candidate centers are the residues of `U_center` as defined in §2, further restricted to centers whose
sphere contains at least one classified residue. The family size `m` used for multiple-testing correction
is the number of centers actually tested. `m` must be reported alongside every corrected p-value, because
it sets the significance bar.

### 5.2 [REVISED] Null model for the per-center test

The per-center permutation test must use the **same structure-aware positional null as §3**: the P/LP point
set is redrawn uniformly without replacement over `U_struct`, conditioned on `N_P`, and the observed count
of P/LP residues within the sphere is compared with its permutation distribution.

A label-permutation null — shuffling P/LP and B/LB labels among the residues of the cohort `L` while holding
their positions fixed — tests a different hypothesis ("are pathogenic residues more clustered than benign
residues, given the set of variant-bearing positions") and is permitted **only** as a secondary, reported
analysis, never as the primary test, and never without the §5.4 certificate. The two nulls answer different
questions and must not be mixed across stages of a single run: a workflow that declares global clustering
under one null and then fails to detect it under another has not produced a negative result, it has produced
an inconsistency.

If both nulls are computed, both must be reported side by side with an explicit statement of which
hypothesis each addresses.

### 5.3 Multiple-testing correction

Benjamini-Hochberg at the declared level `q`, with the Benjamini-Yekutieli result reported alongside.
Boundary ties are all rejected. The correction is fixed before the data are seen and is never loosened.

### 5.4 [NEW] Power certificate — mandatory precondition for any negative result

Before the per-center test is run, and again after it, the agent must compute the **minimum p-value the
test could possibly return**, and compare it with the hardest bar the correction imposes. This is the
single check whose absence permitted the KCNA2 run to be misreported.

Two distinct floors exist and **both** must be evaluated:

1. **Permutation resolution floor.** `p_res = 1 / (B + 1)`. Governed by the number of permutations `B`.
   Remediable by increasing `B`.
2. **Combinatorial floor.** The smallest p-value attainable under the chosen null given the cohort
   composition, achieved when the observed configuration is maximally extreme. **This floor is not
   remediable by increasing `B`.** For a label-permutation null with `N_P` pathogenic and `N_B` benign
   residues and a sphere containing `n_L` classified residues, it is the hypergeometric tail
   `C(N_P, n_L) / C(N_P + N_B, n_L)`. For a positional null over `M` reference residues with `n_U`
   residues inside the sphere, it is the corresponding hypergeometric tail over `M`.

Define `p_floor = max(p_res, p_comb_best)`, where `p_comb_best` is the combinatorial floor evaluated at the
most favourable sphere occupancy attainable at the selected radius. Define the rank-1 critical value
`c_1 = q / m`.

- If `p_floor > c_1`, **no center can be declared significant under any realization of the data.** The run
  must terminate as `UNDERPOWERED` and emit `BLOCKING: TEST_CANNOT_REJECT`, reporting `p_floor`, `c_1`,
  their ratio, `N_P`, `N_B`, `m`, and which of the two floors is binding. It must **not** be reported as a
  negative result, and no claim about the absence of hotspots in the gene may be made.
- If the binding floor is `p_res`, the report must state the value of `B` that would clear `c_1`, and the
  run must be re-executed at that `B` under the same RUN_ID policy.
- If the binding floor is `p_comb_best`, increasing `B` is futile and must not be recommended. The report
  must instead state which of the three remedies applies: a larger benign cohort, a different null (§5.2),
  or a smaller test family (§2, §5.1).
- If `p_floor <= c_1` and no center is nonetheless rejected, the result **is** a valid
  `COMPLETED_NEGATIVE`, and this must be stated with the certificate attached as evidence.

The pre-flight instance of this check uses the largest attainable `n_L` over the candidate domain; the
post-hoc instance uses the realized values. A failure at either instance is binding.

### 5.5 [NEW] Scan/detection consistency

The radius scan and the final detection must draw permutations from independent seed contexts. If their
center sets disagree, the final detection is authoritative, and the disagreement must be reported together
with the §5.4 certificate — because a disagreement in the presence of a failed certificate is not evidence
that the result sits near the significance boundary, it is evidence that the test is operating at its floor.

---

## 6. Documentation of the primary analysis

The complete analysis should be fully documented. The report should include all analysis parameters,
intermediate calculations, quality control results, statistical tests, optimization procedures, decision
criteria, candidate radii that were not selected, and a detailed justification for the final radius
selection. This documentation should ensure that the entire workflow is fully transparent, reproducible,
and readily applicable to any protein or gene using the same standardized analytical framework.

**[NEW] Every report must open with a one-page decision summary** stating, in this order: terminal state;
`N_P`, `N_B`, `|U_struct|`, `|U_center|`, `m`; the null used for the primary test; `p_floor`, `c_1` and
their ratio; the selected radius and whether the Pareto set was vacuous; and the number of significant
centers. A reader must be able to determine from this page alone whether a negative result is
interpretable.

---

## 7. Footprint radius as a multi-scale process

Footprint radius determination should not be treated as the independent evaluation of a single radius.
Instead, the evolution of the footprint geometry as the radius increases should be analyzed as a
**multi-scale process**. In this framework, each radius does not represent an isolated solution but rather a
continuation of the previous one, allowing the geometric evolution of the footprint across increasing
spatial scales to be systematically characterized.

The analysis should be performed across a predefined minimum-to-maximum radius range. At each radius,
spheres should be constructed around statistically significant center residues, and overlapping spheres
should be merged to generate the corresponding footprint. The radius should then be incrementally
increased, repeating this procedure while recording how the footprint geometry changes as a function of
radius.

At each successive radius, the analysis should evaluate not only the resulting footprint itself but also how
it differs from the footprint obtained at the previous radius. In particular, the following geometric events
should be monitored:

- The radius at which previously separate geometric regions begin to merge.
- Whether newly merged regions are connected naturally or through narrow geometric connections
  (artificial bridges).
- How the footprint volume changes as the radius increases.
- Whether the footprint geometry remains stable across specific radius intervals.
- Whether small changes in radius produce disproportionately large structural changes in the footprint.

This approach allows the footprint to be interpreted not as a static structure observed at a single radius,
but as a continuously evolving geometric object whose behavior across multiple spatial scales can be
analyzed.

Geometrically stable regions should be defined as those in which the footprint maintains essentially the
same structural organization over a range of neighboring radii. Conversely, situations in which numerous new
mergers occur over a small radius increment, or where previously distinct regions abruptly collapse into a
single footprint, should be interpreted as evidence of geometric instability.

Accordingly, the optimal footprint radius should not be selected solely on the basis of geometric quality
metrics computed at a single radius. Instead, footprint continuity and geometric stability across multiple
radii should also be incorporated into the decision process. This ensures that the selected footprint is not
merely a structure that appears optimal at one specific radius, but rather represents a robust and
reproducible three-dimensional region whose geometric characteristics are preserved across multiple spatial
scales.

Candidate radii identified through this multi-scale evaluation should subsequently be compared using
geometric quality metrics including **compactness**, **convexity**, **connectivity**, **artificial
bridging**, **expansion rate**, and other relevant geometric descriptors. Finally, a **Pareto-based
multi-objective optimization** should be performed to identify the footprint that provides the best overall
geometric balance across all evaluation criteria.

**[NEW]** The objective/constraint distinction of §4 applies verbatim to footprint selection: no candidate
footprint radius may be eliminated on the grounds that a downstream quantity failed to reach significance,
and a Pareto set of size one must emit `VACUOUS_PARETO_SELECTION`.

---

## 8. Footprint robustness

After the final footprint has been selected, its robustness should be evaluated to determine whether the
observed geometric region reflects a reproducible biological signal rather than a result driven by the
specific set of variants included in the analysis. This should be accomplished using a **Leave-One-Out
(LOO)-based multi-validation framework**. However, unlike the classical LOO approach, removing only a single
residue should not be considered sufficient. Instead, all statistically significant center residues
contributing to hotspot formation should first be identified, and systematic subsets containing different
combinations of these residues should be generated. Whenever computationally feasible, all possible
combinations should be evaluated; otherwise, representative subset sampling strategies should be employed.

For each validation iteration, the selected residue combination should be removed, after which the entire
hotspot analysis should be repeated from the beginning. A new footprint should then be reconstructed and
quantitatively compared with the original solution. This procedure allows assessment of whether the
footprint depends excessively on a small number of influential residues or instead represents a robust
spatial structure supported by multiple independent pathogenic variants.

**[CLARIFIED, DECISION-STAGE-D-FIXED-RFP-0001, AUTHORIZED_BY=user]** "The entire hotspot analysis should be
repeated from the beginning" describes the footprint *construction* procedure applied to the reduced center
set — it does not include re-selecting `r_fp`. The footprint radius is chosen once, by Stage C, using the
full multi-scale sweep, ARI neighboring-radius stability and Pareto/distance-to-ideal criteria over the
complete, unperturbed center set; that choice is frozen before any validation iteration runs. Each iteration
reconstructs the footprint for `S_iter` directly at the frozen `r_fp` — no per-iteration domain re-derivation
for selection purposes, multi-scale sweep, QC-driven admissibility, or Pareto/distance-to-ideal re-selection.
The scientific question this section answers is whether the *selected* footprint survives loss of its
supporting centers, not whether a fresh sweep over `S_iter` would have chosen a different radius. See
`docs/decisions/DECISION-STAGE-D-FIXED-RFP-0001.md` for the measurement that motivated this and what was
preserved.

During each validation iteration, footprint volume, surface area, the number of connected components,
overall geometric integrity, hotspot center locations, and spatial relationships among hotspots should all be
recalculated. In addition, quantitative similarity between footprints should be evaluated using measures
such as **Jaccard similarity**, **Dice coefficient**, or other volume-overlap metrics. Changes in hotspot
center positions, geometric deviation, and structural continuity should also be assessed. Reliable
footprints are expected to remain largely preserved across validation iterations, whereas dramatic geometric
changes following minor modifications of the variant set should be interpreted as evidence of limited
robustness.

In addition to evaluating footprint robustness, hotspot classification performance should also be
recalculated during every validation iteration. The distribution of pathogenic and benign variants located
inside and outside the footprint should be used to determine **true positives (TP), false positives (FP),
true negatives (TN), and false negatives (FN)**. These values should then be used to calculate **Matthews
Correlation Coefficient (MCC)**, sensitivity, specificity, accuracy, positive predictive value, negative
predictive value, F1 score, and balanced accuracy. This enables assessment of whether the selected footprint
remains stable not only geometrically but also in terms of classification performance.

After all validation iterations have been completed, the mean, standard deviation, confidence interval, and
coefficient of variation should be calculated for every performance metric. In addition, the proportion of
iterations in which the footprint is preserved, the stability of hotspot center locations, and the frequency
with which individual hotspots reappear should all be reported. These analyses provide a quantitative
assessment of footprint reproducibility.

**[NEW] Re-running the analysis inside a validation iteration must re-run §5.4.** If removing a subset of
centers causes the power certificate to fail, that iteration is recorded as `UNDERPOWERED`, not as a
footprint that failed to reproduce. Underpowered iterations must be counted and reported separately from
genuine non-reproductions; pooling them understates robustness.

---

## 9. Biological annotation

Subsequently, the final hotspot regions should undergo comprehensive biological annotation. For each
hotspot, the report should include the covered amino acid interval, the numbers of pathogenic and benign
variants, fold enrichment, permutation-based *p*-value, multiple-testing-adjusted significance level,
AlphaFold confidence (pLDDT), secondary structural features, functional regions, ligand- or
protein-interaction sites, known protein domains, evolutionary conservation, and previously reported disease
associations from the literature. This comprehensive annotation ensures that statistically identified
hotspots are also supported by biological evidence.

---

## 10. Functional mechanism overlay

As an additional layer of biological interpretation, all pathogenic variants for which functional effects
have been experimentally characterized should be systematically retrieved from the published literature and
annotated according to their reported functional mechanism, including **Gain-of-Function (GOF)**,
**Loss-of-Function (LOF)**, **Dominant-Negative (DN)**, **Mixed**, **Unclear**, or **Not Experimentally
Characterized**. For each annotated variant, the supporting reference, experimental system, assay type,
principal functional finding, and strength of evidence should be recorded to distinguish directly
demonstrated mechanisms from inferred or weakly supported classifications. These functional labels should
**not** be used during the primary hotspot or footprint discovery process, but should instead be overlaid
onto the independently identified spatial regions after hotspot detection and validation, thereby avoiding
circular reasoning. The distribution of functionally characterized variants across individual hotspots and
footprints should then be summarized to determine whether particular three-dimensional regions are
preferentially associated with specific functional mechanisms. When a sufficient number of experimentally
characterized variants is available, an additional mechanism-specific spatial analysis may be performed to
evaluate whether GOF, LOF, DN, or mixed-effect variants exhibit significant spatial segregation or
enrichment within distinct structural regions of the protein.

---

## 11. Output structure and reproducibility

At the conclusion of the analysis, all intermediate and final outputs should be saved within a standardized,
reproducible directory structure that can be directly applied to different genes without modification.
Numerical results should be exported in machine-readable tabular formats (CSV or TSV), while analysis
parameters, software versions, execution settings, and other metadata should be stored in structured JSON
files. Detailed descriptions of the analytical methods, intermediate results, and final findings should be
generated in Markdown and, where appropriate, PDF format. All visualizations should be saved in both
high-resolution raster (PNG) and vector (SVG or PDF) formats. Protein structures, hotspot annotations, and
footprint models should be exported in appropriate structural biology formats (PDB and/or mmCIF). In
addition, all intermediate outputs — including quality control results, radius optimization analyses, Pareto
optimization results, footprint evaluations, Leave-One-Out (LOO) validation results, hotspot analyses, and
all calculations underlying the final decision process — should be stored as separate outputs, with every
analysis parameter and execution setting comprehensively documented. This ensures that the workflow provides
not only the final hotspot results but also a complete record of the entire analytical decision-making
process, yielding a fully transparent, reproducible, and readily transferable framework applicable to any
protein or gene.

**[NEW] The power certificate of §5.4 is a first-class output** and must be written to disk as a structured
JSON artefact in the primary hotspot stage directory, regardless of the run's terminal state.

---

## Appendix A — [NEW] Prohibited inferences

The following statements must never appear in a report unless the §5.4 certificate passed:

- "No 3D hotspot is detectable in this gene."
- "This is a valid, complete scientific negative."
- "This result licenses no modification of the method."

When the certificate fails, the correct statement is: *"Under the configuration executed, the test could
not have rejected any center regardless of the data; the analysis is uninformative about the presence or
absence of hotspots in this gene."*

## Appendix B — [NEW] Worked example (KCNA2, run 20260817T155124Z)

Retained as a regression case for any implementation of this specification.

| quantity | value |
|---|---|
| `N_P` / `N_B` | 38 / 12 |
| `\|U_struct\|` | 499 |
| `m` (family size, as executed) | 452 |
| null used for per-center test | label permutation |
| best observed configuration | center 335, r = 21 Å, 20 of 20 classified residues P/LP, zero benign |
| observed smallest p | 9.999e-04 |
| combinatorial floor `C(38,20)/C(50,20)` | 7.125e-04 |
| rank-1 critical value `q/m` | 1.106e-04 |
| verdict under §5.4 | `p_floor > c_1` by a factor of 6.4 → `UNDERPOWERED` |
| same data, positional null, r = 10 Å, center 404 | p = 1.2e-07 (clears `c_1` by three orders of magnitude) |

A conforming implementation must return `UNDERPOWERED` on the executed configuration and must detect
hotspots on the §5.2-conforming configuration.
