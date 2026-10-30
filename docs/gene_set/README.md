# Gene-set provenance

How the 68 ACMG SF v3.2 genes in `results_v3/` were chosen, kept here so the selection is
reproducible rather than reconstructible only from the result directories.

| file | what it is |
|---|---|
| `acmg_sf_v3.2_official_84.txt` | the ACMG SF v3.2 secondary-findings list as published, 84 gene symbols |
| `prescreen.py` | the pre-screen, exactly as run. Reads the gene list and a ClinVar `variant_summary` TSV on stdin |
| `prescreen_runnable_68.json` | its output: the 68 that passed, each as `[gene, n_PLP, n_BLB, L, n_PLP/L]` |

`L` is the largest missense residue position observed in ClinVar for that gene, used as a
proxy for protein length.

## The three exclusion tests

1. `L > 2700` — AlphaFold DB serves no model at that length
2. `n_PLP < 1` or `n_BLB < 1` — one binary class empty; the run would terminate BLOCKED
3. `L = 0` — no parseable missense record

Missense means `Assembly = GRCh38`, `Type = single nucleotide variant`, a
`(p.Xxx###Yyy)` three-letter change in the ClinVar `Name` field, and `ref != alt`. The
significance sets match the pipeline's binary policy. **No review-status (star) filter** at
this or any later stage.

## Reproducing it

```bash
zcat variant_summary.txt.gz \
  | python docs/gene_set/prescreen.py docs/gene_set/acmg_sf_v3.2_official_84.txt OUTDIR
```

The printed summary names which gene failed which test — that output was **not persisted**
for the batch actually run, which is why `CODE_AVAILABILITY.md` marks the per-gene
attribution of the 16 exclusions as inference rather than measurement. Re-running the
pre-screen against the same ClinVar release settles it.

## The 16 excluded

```
ACTA2  APC  APOB  BRCA2  CALM1  CALM2  CALM3  DSP
FBN1   FLNC RYR1  RYR2   MYL2   MYL3   TNNC1  TTN
```

They are **not negative results.** Nine are a limitation of the structural source, seven
are cohorts that cannot support a binary contrast. Any "worked in N of M genes" statement
takes **68** as its denominator, not 84.

KCNA2 (the case gene) and MEFV were run alongside the 68 under the same `config_sha256`
and `code_version`, bringing `results_v3/` to 70 runs. Neither came from this list.
