# DECISION-STAGE-E-SCOPE-0001

```
DECISION_ID    = DECISION-STAGE-E-SCOPE-0001
DATE           = 2026-09-12
AUTHORIZED_BY  = user (in principle 2026-09-08; implementation 2026-09-12)
STATUS         = IMPLEMENTED
SCOPE          = Stage E (II.13) — how the out-of-scope literature-mechanism half
                 is reported. Pipeline-wide, all genes.
```

## Decision

The literature/mechanism half of Stage E is **formally scoped out** rather than left as
live-but-empty machinery. Nothing is deleted: the curation code, its config rubric (A23),
the post hoc gate (A24) and the schema columns all stay, because they are implemented,
tested and correct — what is missing is a *curated source*, not the logic.

Two reporting defects are fixed.

### 1. `LITERATURE_MECHANISMS_NOT_RUN`: MAJOR -> ADVISORY

`Severity.MAJOR` means *"may materially change interpretation; must be read"* and is meant
for run-specific anomalies. `LiveAnnotationSource` deliberately does not implement
`literature_mechanisms` (`live_sources.py:429-440`), so this warning fires on **every run
of every gene** — verified in KCNA2 and TP53 run artefacts.

A warning with 100% incidence carries no discriminating information and actively harms the
reader: it desensitises them to the real MAJOR warnings sitting beside it.

`Severity.ADVISORY` — *"worth knowing; does not affect validity"* — is the accurate level.
The caveat genuinely matters (a reader must not read the mechanism labels as negative
findings) but it does not invalidate the structural annotation Stage E does produce:
pLDDT, secondary structure, functional regions, conservation, disease associations, star
composition, robustness summary. **The warning text is unchanged.**

### 2. `mechanism_counts`: zeros -> null

Before this decision, an unrun search reported:

```json
"mechanism_counts": {"GOF": 0, "LOF": 0, "DN": 0, "Mixed": 0,
                     "Unclear": 0, "Not_Experimentally_Characterized": 0}
```

This is the same defect class as `fold_enrichment = 0.0` for an empty selection: **a number
standing in for a quantity that was never measured.** All-zero counts read as *"we curated
the literature and found no mechanism evidence"* — the opposite of the truth. They also
**contradict this stage's own warning**, which states that every label *defaults to*
`Not_Experimentally_Characterized`; that category's count was nevertheless `0`.

When the search was not run, every category is now `None` — the pipeline's existing marker
for a quantity that was never measured (cf. `degenerate_objectives: None` on the
UNDERPOWERED path). The **key** is untouched: `HANDOFF_05_KEYS` is frozen and
`mechanism_counts` remains present; only the value shape changes, which is not a contract
change.

### 3. An explicit caveat in the handoff

`interpretation_caveats` (already a frozen `HANDOFF_05_KEYS` member) now carries a sentence
stating that no curated source was queried, that `mechanism_counts` is null to mark the
quantity NOT MEASURED rather than measured-as-zero, and that the post hoc gate therefore
cannot pass. No new top-level key was introduced.

## What this decision does NOT do

- It does **not** delete the mechanism curation code, `config` rubric (A23) or gate (A24)
- It does **not** change `min_evidence_for_label`, `conflicting_strong_resolution`,
  `weak_only_resolution`, `posthoc_gate.*`, or any mechanism category
- It does **not** implement `literature_mechanisms`; the source remains out of scope
- It does **not** touch `annotation.may_influence_discovery` (FROZEN) or the
  mechanism-blindness of Stages B–D. Discovery remains mechanism-blind.
- It does **not** relax `POSTHOC_GATE_FAILED`, which is already INFO and correct: a gate
  that cannot pass on null input is reporting exactly that.

## Verification

`tests/annotation`: **172 passed** after the change.

## Files changed

- `src/hotspot3d/annotation/stage.py` — severity; `mechanism_counts` null-on-not-run;
  scope sentence appended to `interpretation_caveats`

## Follow-up, not authorized here

If a curated mechanism source is ever supplied, `literature_mechanisms_not_run` goes false
on its own and all three behaviours above revert automatically — the severity, the counts
and the caveat are all conditioned on that one flag. No further edit is required to
re-enable the half.
