# Code and Data Availability

Draft paragraph for the manuscript. Every number below is taken verbatim from the run
artefacts named at the end of this file; none of it is re-derived here.

> **Provenance warning — read before citing.** The Zenodo archive
> [10.5281/zenodo.22237976](https://doi.org/10.5281/zenodo.22237976) is pinned to release
> `v1.0.0` (commit `c734256e`). **That archive does not contain the code that produced the
> numbers below.** Between `v1.0.0` and the code version used here there are 34 changed
> files under `src/`, +1502 / −485 lines, including the replacement of the footprint
> radius selection procedure (`DECISION-RFP-CONTAINMENT-0001`). The code version that
> produced these numbers, `ac7cce0f`, is reachable as commit `8cf973a` — now merged to
> `main` — but it is not tagged and not archived, and the tree no longer reproduces that
> digest: `DECISION-RENAME-SPHARE3D-0001` renamed the pipeline from `hotspot3d` to
> `SPHARE3D`, which changed 602 lines of file content and therefore the source digest
> (`ac7cce0f` → `4eb89577`). See **Code version mapping** below. **A new tagged release
> and a new Zenodo version DOI are required before this paragraph can be submitted.**
> Until then the availability statement is not satisfiable by a reader.

---

## Draft paragraph

> **Code and data availability.** The analysis pipeline (`sphare3d`) is openly available at
> https://github.com/melike00bakirci-crypto/SPHARE3D, released under the MIT licence, and
> archived at Zenodo under DOI [10.5281/zenodo.22237976](https://doi.org/10.5281/zenodo.22237976).
> *[PENDING: cite the version DOI of the release that archives commit `8cf973a` (code
> version `ac7cce0f`), or its renamed equivalent — see the provenance warning above and
> the code version mapping below. Do not submit with the `v1.0.0` DOI.]* The pipeline is
> gene-independent: no step branches on a specific gene or protein, and every scientific
> constant is fixed in a frozen, pre-registered configuration file whose SHA-256 digest is
> recorded in each run.
>
> Variant data were obtained from the ClinVar `variant_summary` weekly release of
> 26 October 2026, with no review-status (star) filtering applied at any point. Missense
> variants were mapped onto the MANE Select transcript ENST00000316361.10 / NM_004974.4
> and the corresponding UniProt canonical sequence P16389 (KCNA2, entry version 227). The
> structural model was the AlphaFold DB prediction AF-P16389-F1 model v6 (database version
> 2025-08-01); residues were represented by their Cα atoms and all distances are Cα–Cα
> Euclidean. pLDDT was recorded throughout and was never used as a filter on the variant
> cohort; it restricts only the set of residues admissible as sphere *centres*.
>
> Spatial clustering was assessed under a structure-aware positional null in which the
> observed number of pathogenic residues is redrawn uniformly without replacement over
> the positional universe of residues with a usable Cα coordinate. Per-residue tests used
> B = 100,000 permutations, and significance was controlled at a Benjamini–Hochberg false
> discovery rate of q < 0.05.
>
> The hotspot radius (r_hot = 10.0 Å) was selected by a pre-registered multi-objective
> procedure over a deterministically derived candidate domain: Pareto filtering on four
> objectives (leave-one-out MCC, permutation evidence, fold enrichment, neighbouring-radius
> stability), min–max normalisation over the admissible set, and L2 distance to the utopia
> point. The selection is reported as a near-tie: the runner-up (9.5 Å) is 0.0139 % behind
> in distance-to-utopia, which the pipeline records at MAJOR severity — below the 1 %
> threshold at which the winner is determined by permutation noise rather than by the data.
> The selected radius is therefore not a sharp optimum and the derived region count is
> sensitive to it; at the runner-up radius the containment radius is 9.0 Å. Both the full
> decision trace and the near-tie diagnostic are deposited.
>
> The footprint radius was **not** selected by an optimisation. It is determined by a
> containment rule (`DECISION-RFP-CONTAINMENT-0001`): the evidence set E is the set of
> P/LP residues lying within r_hot of a significant centre, and r_fp is the maximum
> nearest-significant-centre distance over E, rounded up to a 0.5 Å grid. For KCNA2 this
> gives an exact containment radius of 9.9562 Å and r_fp = 10.0 Å, with every residue of E
> contained in the footprint by construction. Benign/likely-benign residues take no part in
> the determination.
>
> All run artefacts — the complete decision trace for every candidate radius, per-residue
> test statistics for every tested centre including non-significant ones, the perturbation
> robustness profile, and the full provenance record (configuration digest, code version,
> derived random seeds, and software environment) — are deposited with the archived code
> and are sufficient to reproduce the reported analysis exactly.

---

## Parameters, as executed

| Item | Value |
| --- | --- |
| Gene | KCNA2 |
| UniProt accession | P16389 (entry version 227) |
| MANE Select transcript | ENST00000316361.10 / NM_004974.4 |
| ClinVar release | 26 October 2026 |
| Review-status (star) filter | none, at any stage |
| Structural model | AlphaFold DB, `AF-P16389-F1` model v6 (DB 2025-08-01) |
| Residue representation | Cα; Cα–Cα Euclidean distance |
| pLDDT | recorded; never a cohort filter |
| Positional universe \|U_struct\| | 499 residues |
| Candidate-centre universe \|U_center\| | 378 residues (pLDDT ≥ 50; cohort never filtered) |
| Residues excluded as centres | 121, in 4 contiguous regions |
| Protein diameter, D_max | 129.26 Å |
| Classified cohort \|L\| | 50 residues (38 P/LP, 12 B/LB, 0 conflict) |
| Unmapped cohort residues | 0 |
| Primary null | structure-aware positional, without replacement |
| Permutations, B | 100,000 |
| Multiple-testing correction | Benjamini–Hochberg, q < 0.05 |
| Hotspot radius, r_hot | **10.0 Å** |
| r_hot runner-up | 9.5 Å, 0.0139 % behind — near-tie at **MAJOR** severity |
| Significant hotspot centres | **50, in 1 region** |
| Coverage of U_struct at r_hot | 0.2866 |
| Footprint rule | containment (`r_fp = r_contain`), **not** an optimisation |
| Evidence set \|E\| | 36 of 38 P/LP residues (2 outside, fraction 0.0526) |
| Exact containment radius | 9.9562 Å |
| Footprint radius, r_fp | **10.0 Å** (9.9562 rounded up to the 0.5 Å grid; slack 0.0438 Å) |
| Footprint residues | **143** |
| Footprint components | 1 |
| Footprint volume | 43 573.1 Å³ |
| Footprint surface area | 8 285.7 Å² |
| Voxel edge, h | 0.5 Å (no resolution fallback applied) |
| Binding residue | 362 — not an outlier |
| Stage QC status | PASS_WITH_WARNINGS |

## Provenance of these numbers

Every value above is copied from the artefacts of a single run:

```
RUN_ID         KCNA2_20261027T084336Z_ff721543_ac7cce0f
config_sha256  ff7215430c7a23aec095135c63248c111269976f14fa8549ce1e7b930414f3fd
code_version   ac7cce0f8f38efbe09e60b7e0af512a8f9695c6f970bf5c01f45f925693d525f
results root   results_v3/
git            commit 8cf973a on main (NOT tagged, NOT archived)
```

This is the run under the **current** frozen configuration, which incorporates
`DECISION-B-DEFAULT-0001`, `DECISION-NEAR-TIE-0001`, `DECISION-STAGE-E-SCOPE-0001`,
`DECISION-MANE-MATCH-0001`, `DECISION-RFP-CONTAINMENT-0001`, `DECISION-STAGE-D-FIXED-RFP-0001`,
`DECISION-VOXEL-TIE-0001` and `DECISION-UNDERPOWERED-FIXTURE-0001`. The same
`config_sha256` produced the 70-run batch in `results_v3/` (see
`results_v3/BATCH_REPORT_ACMG_SF_v3.2_v3.md`), so the case gene and the reference batch
are directly comparable.

Source files within that run:

- `FULL_RESULTS/03_STRUCTURE_QC/handoff_01.json` — accession, transcript, ClinVar release, cohort sizes, D_max
- `FULL_RESULTS/05_HOTSPOT_RADIUS/center_universe.json` — U_center and the centre exclusions
- `FULL_RESULTS/06_FINAL_HOTSPOTS/handoff_02.json` — r_hot, B, q, FDR method, centre count, near-tie
- `FULL_RESULTS/08_FINAL_FOOTPRINT/handoff_03.json` — r_contain, r_fp, components, coverage, volume, surface
- `FULL_RESULTS/12_REPRODUCIBILITY/` — configuration snapshot, derived seeds, software versions

The `RUN_ID` encodes its own provenance: the timestamp, the first eight characters of the
configuration digest, and the first eight of the code version. Two runs are comparable —
or knowably not — from the directory name alone.

## Code version mapping

`code_version` is a content digest of the source tree. `DECISION-RENAME-SPHARE3D-0001`
renamed the pipeline, which changed file contents and therefore the digest. Every run in
`results_v3/` was produced by the pre-rename tree, so the digest its `RUN_ID` carries is
not reproducible from the current tree. The mapping is recorded here so a reader can still
obtain the exact code that produced the reported numbers:

| `code_version` | tree | commit | what it produced |
|---|---|---|---|
| `ac7cce0f` | last `hotspot3d`-named tree | `8cf973a` | every run in `results_v3/`, including the KCNA2 run cited above |
| `4eb89577` | first `sphare3d`-named tree | the rename commit | nothing yet |

`config_sha256` is unaffected and remains `ff721543`: the digest covers the effective
parameter tree, and the only occurrence of the old name in `config/pipeline.yaml` is a
comment. This was verified by recomputation, not assumed.

To reproduce the reported run exactly, check out `8cf973a` — the package is `hotspot3d`
there, not `sphare3d`.

## Gene set

The reported batch is the 70 runs enumerable as directories under `results_v3/`: 68 genes
from the ACMG SF v3.2 secondary-findings list, plus KCNA2 (the case gene) and MEFV, all
under one `config_sha256` and one `code_version`.

**Gap:** the criterion by which these 68 genes were selected from the full ACMG SF v3.2
list is not recorded anywhere in this repository — not in `config/`, not in `docs/`. The
run set is reconstructible from the result directories, but the *selection rule* is not,
and it is therefore not currently citable. It must be written down before submission.

## Superseded runs

Earlier KCNA2 runs remain in the tree for provenance and **must not be cited**:

| What changed | Decision | Affected runs |
|---|---|---|
| `permutation.B_default` 10,000 → 100,000 | `DECISION-B-DEFAULT-0001` | every run before 2026-08-18 |
| `r_fp` domain → `[5.0, min(0.25·D_max, 25.0)]` | `DECISION-FOOTPRINT-DOMAIN-0001` | all of `results/` (`config_sha256 = 4ec9d6a8`) |
| `r_fp` Pareto search **retired**, replaced by the containment rule | `DECISION-RFP-CONTAINMENT-0001` | all of `results_v2/` (`config_sha256 = 42e59721`) |
| near-tie severity graded; radius-separation suppression removed | `DECISION-NEAR-TIE-0001` | all of `results_v2/` |

The third is not a small correction. Under the retired multi-objective search KCNA2's
footprint radius was reported as 10.5 Å over 150 residues; under the containment rule it is
**10.0 Å over 143 residues**. The 10.5 Å / 150-residue figures belong to the superseded
methodology and must not be cited.

A run's `RUN_ID` carries its own `config_sha256`, so comparability is decidable from the
directory name alone: **cite only runs whose digest begins `ff721543`.**
