# Code and Data Availability

Draft paragraph for the manuscript. Every number below is taken verbatim from the run
artefacts named at the end of this file; none of it is re-derived here.

Archived at Zenodo under DOI [10.5281/zenodo.22237976](https://doi.org/10.5281/zenodo.22237976)
— the version DOI, fixed to release v1.0.0.

---

## Draft paragraph

> **Code and data availability.** The analysis pipeline (`hotspot3d`) is openly available at
> https://github.com/melike00bakirci-crypto/3dhotspot and archived at Zenodo under DOI
> [10.5281/zenodo.22237976](https://doi.org/10.5281/zenodo.22237976) (version 1.0.0), released under the MIT
> licence. The pipeline is gene-independent: no step branches on a specific gene or
> protein, and every scientific constant is fixed in a frozen, pre-registered
> configuration file whose SHA-256 digest is recorded in each run.
>
> Variant data were obtained from the ClinVar `variant_summary` weekly release of
> 17 August 2026, with no review-status (star) filtering applied at any point. Missense
> variants were mapped onto the MANE Select transcript ENST00000316361.10 and the
> corresponding UniProt canonical sequence P16389 (KCNA2). The structural model was the
> AlphaFold DB prediction AF-P16389-F1-model_v6; residues were represented by their Cα
> atoms and all distances are Cα–Cα Euclidean. pLDDT was recorded throughout and was
> never used as a primary filter.
>
> Spatial clustering was assessed under a structure-aware positional null in which the
> observed number of pathogenic residues is redrawn uniformly without replacement over
> the positional universe of residues with a usable Cα coordinate. Per-residue tests used
> B = 100,000 permutations, and significance was controlled at a Benjamini–Hochberg false
> discovery rate of q < 0.05. The hotspot radius (r_hot = 10.0 Å) and the footprint radius
> (r_fp = 10.5 Å) were each selected by an independent, pre-registered multi-objective
> procedure (Pareto filtering, min–max normalisation over the admissible set, and L2
> distance to the utopia point), never by a single metric and never post hoc. The r_hot
> selection is reported as a near-tie: the runner-up (9.5 Å) is 1.58 % behind in
> distance-to-utopia, and the pipeline records this as an ADVISORY-level diagnostic. The
> selected radius is therefore not a sharp optimum, and the derived region count is
> sensitive to it; both radii and their objective vectors are in the deposited decision
> trace.
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
| MANE Select transcript | ENST00000316361.10 |
| ClinVar release | 14 September 2026 |
| Review-status (star) filter | none, at any stage |
| Structural model | AlphaFold DB, `AF-P16389-F1` model v6 (DB 2025-08-01) |
| Residue representation | Cα; Cα–Cα Euclidean distance |
| pLDDT | recorded; never a primary filter |
| Positional universe \|U_struct\| | 499 residues |
| Candidate-centre universe \|U_center\| | 378 residues (pLDDT ≥ 50; cohort never filtered) |
| Protein diameter, D_max | 129.26 Å |
| Classified cohort \|L\| | 50 residues (38 P/LP, 12 B/LB, 0 conflict) |
| Primary null | structure-aware positional, without replacement |
| Permutations, B | 100,000 |
| Multiple-testing correction | Benjamini–Hochberg, q < 0.05 |
| Hotspot radius, r_hot | **10.0 Å** (near-tie, ADVISORY: 9.5 Å is 1.58 % behind) |
| Significant hotspot centres | **50, in 1 region** |
| Residues covered by hotspot spheres | 143 |
| Footprint search domain | [5.0, 25.0] Å; binding constraint `hard_ceiling` |
| Footprint radius, r_fp | **10.5 Å** |
| Footprint residues | **150** |
| Footprint components | 1 |

## Provenance of these numbers

Every value above is copied from the artefacts of a single run:

```
RUN_ID         KCNA2_20260916T214717Z_42e59721_f459f3a2
config_sha256  42e5972130bccda13f1d644be0e7325eb359d255c8564019fe529ba3d193e9c1
code_version   f459f3a22f00cea20db97e4b5adae95ffbf5e32e8f256f39d071aedf5e6bbd74
results root   results_v2/
```

This is the run under the **current** frozen configuration, which incorporates
`DECISION-FOOTPRINT-DOMAIN-0001`, `DECISION-NEAR-TIE-0001` and
`DECISION-STAGE-E-SCOPE-0001`. The same `config_sha256` produced the 68-gene ACMG SF v3.2
batch in `results_v2/` (see `results_v2/BATCH_REPORT_ACMG_SF_v3.2_v2.md`), so the case gene
and the reference batch are directly comparable.

Source files within that run:

- `FULL_RESULTS/03_STRUCTURE_QC/handoff_01.json` — accession, transcript, ClinVar release, cohort sizes
- `FULL_RESULTS/06_FINAL_HOTSPOTS/handoff_02.json` — r_hot, B, q, FDR method, centre count
- `FULL_RESULTS/08_FINAL_FOOTPRINT/handoff_03.json` — r_fp, components, coverage
- `FULL_RESULTS/12_REPRODUCIBILITY/` — configuration snapshot, derived seeds, software versions

The `RUN_ID` encodes its own provenance: the timestamp, the first eight characters of the
configuration digest, and the first eight of the code version. Two runs are comparable —
or knowably not — from the directory name alone.

## Superseded runs

Earlier KCNA2 runs remain in the tree for provenance and **must not be cited**:

| What changed | Decision | Affected runs |
|---|---|---|
| `permutation.B_default` 10,000 -> 100,000 | `DECISION-B-DEFAULT-0001` | every run before 2026-08-18 |
| `r_fp` domain -> `[5.0, min(0.25*D_max, 25.0)]` | `DECISION-FOOTPRINT-DOMAIN-0001` | every run under `config_sha256 = 4ec9d6a8`, i.e. all of `results/` |

The second is not a small correction: across the 68-gene ACMG SF batch, **the footprint
radius changed in every gene that produced one under both configurations** (43 of 43), and
six genes that had terminated `BLOCKED` on the old domain now complete. The old
`r_fp = 10 Å` / 143-residue figures for KCNA2 belong to the superseded methodology.

A run's `RUN_ID` carries its own `config_sha256`, so comparability is decidable from the
directory name alone: cite only runs whose digest begins `42e59721`.
