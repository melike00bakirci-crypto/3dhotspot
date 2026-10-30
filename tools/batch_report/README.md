# Batch report generators

The scripts that turn a `results_v3/`-style run tree into the batch data table, the regime
classification and the consolidated report. They are kept here because the report is
generated output and is never committed — without them it could not be rebuilt.

Run from the repository root, in this order:

```bash
python tools/batch_report/collect.py          # -> BATCH_DATA_*.tsv            (35 cols)
python tools/batch_report/b1_confound.py      # -> b1_confound.json            (slow)
python tools/batch_report/classify.py         # -> BATCH_DATA_*_classified.tsv (41 cols)
python tools/batch_report/consolidated.py     # -> BATCH_REPORT_*.md
python tools/batch_report/shell_tie_audit.py  # -> SHELL_TIE_AUDIT_*.tsv / .md
```

`classify.py` reads `results_v3/b1_confound.json` if it is there and leaves the three `b1_*`
columns empty if it is not, so the chain still completes without the slow step.

| script | what it does |
|---|---|
| `collect.py` | one row per gene from the run artefacts. Picks the newest `RUN_ID` per gene and skips `.partial` trees, so a re-run (e.g. TSC1) does not produce two rows |
| `b1_confound.py` | the B1 diagnostic: draws `N_P` residues uniformly from the folded subset — a set with no clustering — and runs the same pair correlation against the same null. Five draws per gene, `B = 600`. Imports `sphare3d.spatial.ripley`; **read-only**, writes nothing into any run tree |
| `classify.py` | assigns `regime`, `tier`, `caveat` and the three `b1_*` columns |
| `consolidated.py` | the report: regime key, every gene with its caveat, the B1 section, the pipeline assessment, QA list, terminal states |
| `shell_tie_audit.py` | shell sparsity and the `>` vs `>=` sensitivity of the E criterion |

## Two things in here that are easy to get wrong

**The B1 ratio is meaningless without the overlap.** `b1_explained` is the confound's
excess over the observed excess, and in genes whose observed excess is tiny it explodes —
`PALB2` reaches 40.96 at `b1_overlap = 0.00` with one elevated shell. A confound that
elevates *different* radii explains nothing. Everything is gated on overlap; don't report
the ratio alone.

**`diffuse_no_localised_region` is not a refusal.** The centres stand and are the finding;
only the footprint geometry is withheld. See `DECISION-HIGH-COVERAGE-GENES-0001` and its
2026-10-30 addendum.
