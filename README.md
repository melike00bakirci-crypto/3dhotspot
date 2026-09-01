# hotspot3d

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22237976.svg)](https://doi.org/10.5281/zenodo.22237976)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Gene-independent 3D spatial hotspot and footprint pipeline for pathogenic missense variants.**

`hotspot3d` takes an HGNC gene symbol, pulls the ClinVar missense record and the AlphaFold
structure for its MANE Select transcript, and asks whether pathogenic variants cluster in
three-dimensional space more than the structure itself would produce by chance. Where they do,
it delimits the hotspot, builds a geometric *footprint* around it, stress-tests that footprint
against perturbation, and annotates the surviving regions with curated functional-mechanism
evidence.

Nothing in the pipeline branches on a specific gene or protein. The same code path runs TP53,
KCNA2 and a synthetic fixture identically.

---

## Table of contents

- [Design commitments](#design-commitments)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Command-line interface](#command-line-interface)
- [Pipeline stages](#pipeline-stages)
- [Output layout](#output-layout)
- [Terminal states](#terminal-states)
- [Configuration and frozen methodology](#configuration-and-frozen-methodology)
- [Reproducibility](#reproducibility)
- [Testing](#testing)
- [Repository layout](#repository-layout)
- [Research data](#research-data)
- [Decision records](#decision-records)

---

## Design commitments

These are structural, not stylistic. Most of the code exists to enforce them.

**Gene-independence.** No module branches on gene identity. A gene-specific result must come
from gene-specific *data*, never from gene-specific *code*.

**Frozen, pre-registered methodology.** Every scientific constant that could be tuned to
influence an outcome is listed in `FROZEN_ASSERTIONS` (`src/hotspot3d/utils/config.py`) and
re-checked against `config/pipeline.yaml` before any stage executes. Drift aborts the run with
`BLOCKED` rather than silently producing a different analysis. Changing one of these requires an
authorised decision record in `docs/decisions/`.

**Separation of discovery from interpretation.** Biological annotation runs strictly downstream
of, and cannot feed back into, hotspot discovery — `annotation.may_influence_discovery` is
frozen to `False`. Stage E cannot move, resize, re-threshold or validate a hotspot.

**Structure-aware null.** Significance is never assessed against a uniform-in-space null. Both
the global clustering test and the final residue-centred test permute within the observed Cα
positional universe, so the protein's own shape is not mistaken for signal.

**No filtering that quietly reshapes the cohort.** ClinVar review stars are recorded but never
used for inclusion (`review_star_filter: None`). pLDDT is recorded but never used for primary
filtering — it enters only as a *sensitivity* check at ≥ 70.

**Honest terminal states.** A run that lacks power reports `UNDERPOWERED`; a run with no signal
reports `COMPLETED_NEGATIVE`. Neither is repaired into something more publishable.

---

## Installation

Requires Python ≥ 3.10. This project uses [`uv`](https://docs.astral.sh/uv/).

```bash
git clone git@github-melike:melike00bakirci-crypto/3dhotspot.git
cd 3dhotspot

uv venv
source .venv/bin/activate
uv pip install -e ".[dev]"
```

Runtime dependencies: `numpy`, `scipy`, `pandas`, `pyyaml`, `scikit-image`, `scikit-learn`,
`matplotlib`, `networkx`, `requests`, `biopython`. Exact resolved versions are pinned in
`uv.lock`.

Verify the install and the methodology guard:

```bash
hotspot3d version
hotspot3d validate-config     # → CONFIG VALID — frozen methodology constants intact
```

`validate-config` also prints the `config_sha256` recorded in every run's provenance.

---

## Quick start

```bash
# Full run against live ClinVar + AlphaFold
hotspot3d run --gene TP53

# Offline run against a bundled synthetic fixture — no network
hotspot3d run --gene SYNTH --synthetic-case clustered
```

Results land in `results/<GENE>/<GENE>_<RUN_ID>/`. Start with `REVIEW_PACK/key_metrics.tsv`
and `00_RUN_SUMMARY/`.

---

## Command-line interface

```
hotspot3d run --gene GENE [options]
hotspot3d validate-config
hotspot3d show-config [KEY]
hotspot3d version
```

### `run` options

| Flag | Default | Purpose |
|---|---|---|
| `--gene` | *required* | HGNC gene symbol, e.g. `TP53` |
| `--config` | `config/pipeline.yaml` | Frozen configuration |
| `--gene-config` | auto | Per-gene overlay; auto-detected at `config/genes/<GENE>.yaml` |
| `--results-root` | `results` | Output root |
| `--run-id` | generated | Pin the `RUN_ID` for deterministic reproduction checks |
| `--synthetic-case` | — | Run a named fixture instead of live data (`clustered`, `no_cluster`, `sparse`, `conflict`, `low_plddt`, …) |
| `--quiet` | off | Suppress the terminal summary |

`show-config` with no key dumps the fully resolved configuration; with a dotted key
(`hotspot3d show-config fdr.q`) it prints one value.

Exit codes: `0` success, `2` blocked (invalid config, or a stage that could not proceed).

---

## Pipeline stages

Five scientific stages write into thirteen numbered output directories. Each directory has
exactly one owning agent, and write scope is enforced at runtime by `assert_may_write`
(`src/hotspot3d/utils/runctx.py`) — a stage physically cannot write outside its own territory.

### Stage A — input preparation and QC · `data-structure`

`01_INPUT_RAW` · `02_CLINVAR` · `03_STRUCTURE_QC`

Acquires ClinVar missense records on the MANE Select transcript with **no review-star
filtering**, classifies residues binary P/LP vs B/LB (handling `RESIDUE_CLASS_CONFLICT` where a
residue carries both), fetches and QCs the AlphaFold model, maps ClinVar residues onto the
structure, and establishes the Cα positional universe `U_struct` and classified cohort `L`.
pLDDT is *recorded, never filtered on*.

### Stage B — spatial hotspot discovery · `hotspot-statistics`

`04_GLOBAL_CLUSTERING` · `05_HOTSPOT_RADIUS` · `06_FINAL_HOTSPOTS` · `11_SENSITIVITY`

Global Ripley's K and pair correlation under the structure-aware positional null, then a
deterministic radius scan over `[R_FLOOR, R_CEIL] = [5.0, 25.0] Å` in 0.5 Å steps, scoring
LOO-MCC, permutation evidence, fold enrichment and neighbouring-radius stability. `r_hot` is
selected by Pareto filtering plus min–max normalisation and L2 distance-to-ideal. Final
residue-centred testing draws from the *same* null, then BH-FDR at `q = 0.05` yields
`SIGNIFICANT_HOTSPOT_CENTERS`. Emits a power certificate and a permutation-resolution
diagnostic.

### Stage C — footprint construction · `footprint-robustness`

`07_FOOTPRINT_RADIUS` · `08_FINAL_FOOTPRINT`

An **independent** `r_fp` search — not bounded below by `r_hot`
(`bound_below_by_r_hot: False`) — over MST merge scales and a multi-scale occupancy-grid sweep,
scoring components, merge events, artificial bridging, compactness, convexity, connectivity,
expansion rate and ARI stability. Candidates must satisfy the coverage ≤ 0.50 QC rule.
Selection is again Pareto plus distance-to-ideal, producing `FP_original`.

### Stage D — geometric robustness · `footprint-robustness`

`09_ROBUSTNESS`

Perturbs the `SIGNIFICANT_HOTSPOT_CENTERS` set and rebuilds `FP_iteration` using the **frozen**
footprint methodology only. It never re-runs Ripley's K, the radius scan, permutation testing
or FDR (`rerun_hotspot_pipeline: False`), and never modifies ClinVar records
(`modify_clinvar_records: False`). Reports a continuous preservation profile against a 0.50
primary threshold rather than a categorical verdict.

### Stage E — biological annotation · `biological-annotation`

`10_ANNOTATION`

Strictly downstream and non-feedback. Per hotspot: amino-acid interval, pathogenic/benign
content copied verbatim, pLDDT, secondary structure, functional regions, ligand and
protein-interaction sites, domains, conservation, disease associations — plus systematically
curated experimental functional mechanism (GOF / LOF / DN / Mixed / Unclear / Not
Experimentally Characterized) with references, systems, assays and graded evidence.

### Lead

`00_RUN_SUMMARY` · `12_REPRODUCIBILITY`

Run assembly, manifest, provenance, escalations and the determinism record.

---

## Output layout

```
results/<GENE>/<GENE>_<RUN_ID>/
├── README.txt
├── run_summary.json
├── gene_summary.tsv
├── manifest.tsv              # one row per file — produced AND expected-but-absent
├── warnings.tsv
├── FULL_RESULTS/
│   ├── 00_RUN_SUMMARY/  01_INPUT_RAW/     02_CLINVAR/       03_STRUCTURE_QC/
│   ├── 04_GLOBAL_CLUSTERING/  05_HOTSPOT_RADIUS/  06_FINAL_HOTSPOTS/
│   ├── 07_FOOTPRINT_RADIUS/   08_FINAL_FOOTPRINT/ 09_ROBUSTNESS/
│   └── 10_ANNOTATION/   11_SENSITIVITY/   12_REPRODUCIBILITY/
├── REVIEW_PACK/              # the human entry point
│   ├── key_metrics.tsv / .json
│   ├── selected_radii.tsv
│   ├── decision_trace/
│   ├── key_figures/
│   └── structures/
└── ARCHIVES/
```

`manifest.tsv` deliberately records files that were *expected but not produced*, with a reason —
absence is evidence, so it is logged rather than left implicit.

---

## Terminal states

Reported as-is. The pipeline never adjusts methodology, thresholds or cohort definitions to
move a run into a more favourable state.

| State | Meaning |
|---|---|
| `COMPLETED` | Ran to completion with significant findings |
| `COMPLETED_WITH_WARNINGS` | Completed; see `warnings.tsv` |
| `COMPLETED_NEGATIVE` | Ran correctly, no significant spatial clustering |
| `UNDERPOWERED` | Statistically inconclusive given cohort size |
| `UNDERPOWERED_COHORT` | Cohort too small to test meaningfully |
| `ROBUSTNESS_NOT_EVALUABLE` | Too few centres for Stage D perturbation |
| `BLOCKED` | Could not proceed — config drift, missing input, unrecoverable stage failure |

---

## Configuration and frozen methodology

All parameters live in `config/pipeline.yaml`, optionally overlaid per gene from
`config/genes/<GENE>.yaml` (not present by default — create it only if a gene genuinely needs
an overlay). A subset is **FROZEN** and asserted before every run:

| Parameter | Frozen value | Why it matters |
|---|---|---|
| `representation.coordinate_atom` | `CA` | One atom per residue; no side-chain ambiguity |
| `representation.distance_metric` | `euclidean` | — |
| `clinvar.review_star_filter` | `None` | Review stars recorded, never used for inclusion |
| `clinvar.use_review_stars_for_inclusion` | `False` | — |
| `plddt.primary_filtering_enabled` | `False` | Confidence never silently reshapes the cohort |
| `plddt.sensitivity_threshold` | `70.0` | Enters only as a sensitivity check |
| `permutation.B_default` | `100000` | Raised from 10 000 — see `DECISION-B-DEFAULT-0001` |
| `fdr.method` / `fdr.q` | `BH` / `0.05` | BY reported for comparison only, never as a substitute |
| `loo_mcc.role` | `radius_selection_only` | MCC selects a radius; it never classifies |
| `radius_domain.R_FLOOR` / `R_CEIL` / `step_hot_A` | `5.0` / `25.0` / `0.5` Å | Deterministic search domain |
| `radius_domain.automatic_domain_expansion` | `False` | The domain cannot silently widen to find a hit |
| `boundary_diagnostic.auto_widen_domain` | `False` | — |
| `radius_selection` / `footprint_selection` | `min_max`, `L2`, `w_k = 1.0` | Identical selection rule both times |
| `footprint_domain.bound_below_by_r_hot` | `False` | Stage C is genuinely independent of Stage B |
| `footprint_qc.QC_F1_max_coverage` | `0.50` | A footprint cannot swallow the protein |
| `robustness.rerun_hotspot_pipeline` | `False` | Stage D tests geometry, not significance |
| `robustness.modify_clinvar_records` | `False` | — |
| `robustness.emit_categorical_verdict` | `False` | Continuous profile, not pass/fail |
| `robustness.preservation_threshold_primary` | `0.50` | — |
| `annotation.may_influence_discovery` | `False` | Interpretation cannot feed back into discovery |
| `seeding.MASTER_SEED` | `20250101` | Fixed seed for reproducibility |

Editing any of these without a decision record causes the next run to abort:

```
FROZEN methodology violation: 'permutation.B_default' is 10000, expected 100000
```

---

## Reproducibility

- **Fixed master seed** (`20250101`), with per-stage seeds derived deterministically.
- **`config_sha256`** and a **`code_version`** hash recorded in every run's provenance.
- **`uv_pip_freeze`** captured per run, so the resolved environment is recoverable.
- **`12_REPRODUCIBILITY`** holds the determinism record; `--run-id` pins the `RUN_ID` so a run
  can be repeated byte-for-byte and diffed.

The run directory name carries this provenance in the clear:

```
FGFR3_20260819T094511Z_4ec9d6a8_902e1b03
 |     |                |        └── code_version prefix
 |     |                └── config_sha256 prefix
 |     └── UTC timestamp
 └── gene
```

So two runs are directly comparable — or knowably not — from their directory names alone,
without opening a single file. Current values are printed by `hotspot3d validate-config`
(`config_sha256`) and `hotspot3d version` (`code_version`).

---

## Testing

705 test functions across unit, integration and determinism suites.

```bash
uv run pytest -q                       # everything
uv run pytest -q -m unit               # scientific-core units
uv run pytest -q -m integration        # end-to-end synthetic
uv run pytest -q -m "not slow"         # skip long-running
```

Markers: `unit`, `integration`, `slow`. Note that `filterwarnings = ["error::RuntimeWarning"]`
is set — a numerical `RuntimeWarning` (an overflow, an invalid divide) fails the suite rather
than scrolling past.

`scripts/run_definitive_suite.sh` runs the full definitive suite.

---

## Repository layout

```
src/hotspot3d/
├── data/           Stage A — ClinVar, classification, AlphaFold, QC, mapping
├── spatial/        Ripley's K, pair correlation, the structure-aware null
├── hotspot/        Stage B — radius scan, selection, permutation, FDR
├── footprint/      Stages C & D — r_fp search, FP_original, robustness
├── robustness/     Perturbation and geometric profiling
├── annotation/     Stage E — annotation and functional-mechanism curation
├── reporting/      Manifest, summaries, review pack, figures
├── orchestration/  CLI, pipeline driver, stage contracts
└── utils/          Config guard, seeding, hashing, run context, QA taxonomy
config/pipeline.yaml    Frozen configuration
docs/                   Workflow spec v2.2 and decision records
tests/                  Unit, integration, determinism
scripts/                Definitive suite runner
.claude/agents/         Per-stage agent definitions
```

---

## Research data

Pipeline outputs are **not tracked in git** — `results/` alone reaches ~9.3 GB, with individual
run archives around 445 MB, well past GitHub's 100 MB per-file hard limit. `logs/` and the
per-gene `*_analiz.tar.gz` bundles are excluded for the same reason. This is a size constraint,
not a confidentiality one: the data is ordinary research output.

Genes analysed to date: **CACNA1A, FGFR3, KCNA2, MEFV, PTPN11, TP53RK**.

To share these outputs, use a data repository built for the size — Zenodo, Figshare or an
institutional store — and cite the DOI here. Git LFS is possible but would need paid data packs
at this volume.

Any run is regenerable from a clean checkout with `hotspot3d run --gene <GENE>`, since the
configuration, seed and code version are all pinned.

---

## Decision records

Changes to frozen methodology are documented in `docs/decisions/`:

- `DECISION-B-DEFAULT-0001` — permutation `B_default` raised 10 000 → 100 000.
- `DECISION-STAGE-D-FIXED-RFP-0001` — Stage D uses the frozen footprint methodology.

The full specification is `docs/3D_Hotspot_Footprint_Workflow_v2_2.md`; the mechanism-annotation
extension is `3D_Hotspot_Footprint_Workflow_with_Functional_Mechanism.md`.
