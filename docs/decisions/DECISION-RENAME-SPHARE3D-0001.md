# DECISION-RENAME-SPHARE3D-0001 — the pipeline is renamed to SPHARE3D

**Status:** APPLIED
**Authorized by:** user, 2026-10-29
**Scope:** identity only. No scientific constant, threshold, rule or result is changed.

## Why the name had to change

The pipeline was called `hotspot3d`. That name is taken by published software in the
same field: **HotSpot3D** (Ding Lab, Washington University), which identifies mutation
hotspots in protein 3D structures and was published in *Nature Genetics* in 2016
(`10.1038/ng.3586`, PMID 27294619). It is distributed as `ding-lab/hotspot3d` on GitHub
and as the package `hotspot3d` on Bioconda, so `conda install hotspot3d` already resolves
to their tool.

This is not a cosmetic clash. Same name, same GitHub repository name, same package name,
same problem domain. It would have been raised in review and would have made the two
tools indistinguishable in citation and installation.

## The new name

**SPHARE3D** — **S**patial **P**athogenic **H**otspot **A**nalysis and **R**egion
**E**xtraction, in **3D**.

Convention: the package, the CLI command and the repository are lowercase `sphare3d`;
prose, figures and the manuscript use `SPHARE3D`.

Availability was checked on 2026-10-29 before the name was adopted, because a name
embedded in a git tag and a Zenodo DOI is expensive to change afterwards:

| channel | result |
|---|---|
| GitHub, exact repository name `sphare3d` | 0 repositories |
| PyPI `sphare3d` | 404 — unregistered |
| Literature (PubMed, bioRxiv, Nature, OUP) | no occurrence |
| Bioconda / anaconda.org | no occurrence |

Near neighbours, recorded so they are not mistaken for collisions later: `sphAIre`
(`sphaire3d.com`, AI 3D modelling for industrial design — extra `i`, unrelated field);
`sphae` (a phage-therapy toolkit, *Bioinformatics Advances* 2025, on PyPI and Bioconda —
five letters against our eight); `SPHARM-PDM` (spherical-harmonics shape analysis in
medical imaging).

One property of the name is worth stating plainly rather than discovering later:
`sphare` is a common misspelling of `sphere` — 208 GitHub repositories carry it, almost
all of them small web projects that meant `sphere`. The exact name `sphare3d` is clear,
but readers will mistype it. The mitigation is to spell the acronym out on first use, so
the name reads as deliberate rather than as an error; no defensive typo-squat package is
registered and none is planned.

## What changed mechanically

- `src/hotspot3d/` → `src/sphare3d/` (`git mv`, history preserved)
- 602 occurrences of `hotspot3d` across 119 tracked files
- 5 occurrences of `3dhotspot` / `3Dhotspot` (repository URLs, the clone example in
  `README.md`, the `CLAUDE.md` heading)
- `pyproject.toml`: `name`, the `project.scripts` entry point, and `version`
  `0.1.0` → `1.0.0` (it had drifted from the `v1.0.0` git tag)
- display strings raised to `SPHARE3D` where a human reads them: the CLI version banner,
  the run index `<title>`, the two footprint-view headers, the `CITATION.cff` title and
  abstract
- `git remote set-url origin git@github-melike:melike00bakirci-crypto/sphare3d.git`

Nothing under `config/` changed except the word in a comment on line 2.

## The one consequence, and how it is contained

`code_version` is the SHA-256 of the source tree's contents and intra-package relative
paths. Renaming the package **directory** does not affect it — paths are hashed relative
to the package root — but renaming 602 lines of file **content** does:

```
before the rename   ac7cce0f8f38efbe09e60b7e0af512a8f9695c6f970bf5c01f45f925693d525f
after the rename    4eb895777f569630987025638d86499c0eda34d1442d1a3009f906db1b17b225
```

Every one of the 71 runs in `results_v3/` records `code_version ac7cce0f…` and carries it
in its `RUN_ID`. After this commit no state of the tree reproduces that digest, so the
`RUN_ID` promise — that comparability is decidable from the directory name alone — would
have become unresolvable. It is preserved by recording the mapping here and in
`CODE_AVAILABILITY.md`:

| digest | tree | commit |
|---|---|---|
| `ac7cce0f` | the last `hotspot3d`-named tree; produced every run in `results_v3/` | `8cf973a` |
| `4eb89577` | the first `sphare3d`-named tree | this commit |

`config_sha256` is **unchanged** at `ff721543`. This was verified rather than assumed: the
digest is taken over the effective parameter tree, the only `hotspot3d` string in
`config/pipeline.yaml` is a comment on line 2, and recomputing the digest with that
comment altered returns the same value, including the eight characters that enter the
`RUN_ID`.

## What is NOT changed

No result. `results_v3/` is not rewritten, and its numbers, reports, regimes and caveats
stand exactly as they were. The rename is an identity change above the methodology, not
within it.

Seed derivation is also untouched, but a related property deserves recording because it
surfaced while checking this: seeds are `sha256(master_seed|run_id|context)` and `run_id`
contains a timestamp, so **two runs of identical code at different times already draw
different permutations.** The pipeline's reproducibility guarantee is "same `RUN_ID` →
same result", not "same code → same result". The rename does not weaken that; the
timestamp already determined it.

## Follow-through still owed

- the GitHub repository must be renamed `3dhotspot` → `sphare3d` by the account owner;
  GitHub then redirects the old URL
- a new tag, and a new Zenodo version under the existing record — see the provenance
  warning in `CODE_AVAILABILITY.md`, which is a separate and larger problem: the archived
  `v1.0.0` does not contain the code that produced the reported results
- the local working directory is still `~/Sandbox/Melike/3Dhotspot`; renaming it was not
  authorized in this session and no file depends on the path
