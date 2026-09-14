# Code and Data Availability

Draft paragraph for the manuscript. Every number below is taken verbatim from the run
artefacts named at the end of this file; none of it is re-derived here.

Archived at Zenodo under DOI [10.5281/zenodo.22237976](https://doi.org/10.5281/zenodo.22237976)
— the version DOI, fixed to release v1.0.0.

> ## ⚠ NUMBERS BELOW ARE PRE-DECISION — DO NOT CITE YET
>
> Every value in this file comes from `KCNA2_20260819T010956Z_4ec9d6a8_ad8f32d6`, run under
> `config_sha256 = 4ec9d6a8`. Three methodology decisions have since been authorized and
> implemented, so that configuration is no longer what the pipeline ships:
>
> | Decision | Effect on this file |
> |---|---|
> | `DECISION-FOOTPRINT-DOMAIN-0001` | the `r_fp` domain is now `[5.0, min(0.25*D_max, 25.0)]`. Measured on KCNA2 with seeds pinned: **`r_fp` moves 10.0 -> 10.5 A, footprint residues 143 -> 150.** |
> | `DECISION-NEAR-TIE-0001` | KCNA2's `r_hot` is a **MAJOR near-tie: 0.178%** from the runner-up. The region count (1 vs 2) and the centre count (50 vs 46) flip with it. This must be stated in the manuscript, not omitted. |
> | `DECISION-STAGE-E-SCOPE-0001` | annotation reporting only; no number here changes. |
>
> `config_sha256` has changed, so **results produced under the old configuration are not
> comparable with results produced under the new one and must not be pooled.** This file
> will be regenerated in full from a single re-run of the frozen pipeline under the current
> configuration; until that run completes, treat every figure below as superseded.
>
> The prose paragraph and the parameter table are kept verbatim meanwhile so the diff after
> the re-run shows exactly which quantities moved.

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
> discovery rate of q < 0.05. The hotspot radius (r_hot = 9.5 Å) and the footprint radius
> (r_fp = 10 Å) were each selected by an independent, pre-registered multi-objective
> procedure (Pareto filtering, min–max normalisation over the admissible set, and L2
> distance to the utopia point), never by a single metric and never post hoc.
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
| UniProt accession | P16389 |
| MANE Select transcript | ENST00000316361.10 |
| ClinVar release | 17 August 2026 |
| Review-status (star) filter | none, at any stage |
| Structural model | AlphaFold DB, `AF-P16389-F1-model_v6` |
| Residue representation | Cα; Cα–Cα Euclidean distance |
| pLDDT | recorded; never a primary filter |
| Positional universe \|U_struct\| | 499 residues |
| Classified cohort \|L\| | 50 residues (38 P/LP, 12 B/LB) |
| Primary null | structure-aware positional, without replacement |
| Permutations, B | 100,000 |
| Multiple-testing correction | Benjamini–Hochberg, q < 0.05 |
| Hotspot radius, r_hot | 9.5 Å |
| Significant hotspot centres | 46, in 2 regions |
| Footprint radius, r_fp | 10 Å |
| Footprint components | 1 |

## Provenance of these numbers

Every value above is copied from the artefacts of a single run:

```
RUN_ID         KCNA2_20260819T010956Z_4ec9d6a8_ad8f32d6
config_sha256  4ec9d6a82a9dcf423364944806f99a7070e215804de2b9ea719c2bb7b9689184
code_version   ad8f32d6...
```

Source files within that run:

- `FULL_RESULTS/03_STRUCTURE_QC/handoff_01.json` — accession, transcript, ClinVar release, cohort sizes
- `FULL_RESULTS/06_FINAL_HOTSPOTS/handoff_02.json` — r_hot, B, q, FDR method, centre count
- `FULL_RESULTS/08_FINAL_FOOTPRINT/handoff_03.json` — r_fp, components, coverage
- `FULL_RESULTS/12_REPRODUCIBILITY/` — configuration snapshot, derived seeds, software versions

The `RUN_ID` encodes its own provenance: the timestamp, the first eight characters of the
configuration digest, and the first eight of the code version. Two runs are comparable —
or knowably not — from the directory name alone.

## Note on B

`permutation.B_default` was raised from 10,000 to 100,000 before this run, under the
user-authorised decision recorded in `docs/decisions/DECISION-B-DEFAULT-0001.md`. Earlier
KCNA2 runs at B = 10,000 exist in the run tree for provenance but are superseded and must
not be cited as the production result.
