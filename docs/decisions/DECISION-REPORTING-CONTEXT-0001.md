# DECISION-REPORTING-CONTEXT-0001

**Status:** IN PROGRESS
**Scope:** reporting only. No statistic, threshold, null, cohort or radius changes.
**Motivation:** a read-only Stage A->E assessment of the `results_v2/` batch (69 runs).

Two changes, both about what a reader is shown. Neither alters a number, a decision or
a terminal state of any gene that produced a result.

## 1. The centre count no longer travels alone

`n_significant_centers` reached `handoff_02`, `run_summary.json` and the review-pack
metrics with the radius beside it and nothing else. Measured across the batch, that is
not enough context to read it:

| gene | coverage at `r_hot` | fold enrichment | centres |
|---|---|---|---|
| LDLR | 0.660 | 1.043 | 346 |
| GAA | 0.872 | 1.016 | 334 |
| MLH1 | 0.729 | 1.081 | 248 |
| MEN1 | 0.849 | 1.032 | 230 |
| PTEN | 0.754 | 1.032 | 166 |
| RPE65 | 0.872 | 0.999 | 88 |

**14 of the 51 genes that produced centres did so from spheres covering more than half
the structure, at fold enrichment ~1.0** — i.e. no enrichment for P/LP relative to the
labelled background. The batch's largest centre counts are concentrated there.

This is not a detection failure. **All 14 carry a MAJOR `EXCESSIVE_STRUCTURAL_COVERAGE`,
and 11 also carry `HIGH_COHORT_ABSORPTION`** — the pipeline identifies the condition
without exception. The defect is that the warning lived in a separate file while the
number travelled on, so "334 significant hotspot centres in GAA" could be read, cited
and carried into a manuscript without the fact that the spheres cover 87 % of the
protein.

Nor is it a methodological inconsistency: the primary null asks whether P/LP residues
are placed more tightly than uniform placement over `U_struct`, which is a different
question from whether a sphere is enriched for P/LP against B/LB. `hotspot/scan.py`
says so explicitly — fold enrichment is "deliberately kept apart from significance".
Both facts can hold at once, and both are now reported together.

**Implemented:** `coverage_fraction_at_r_hot`, `fold_enrichment_at_r_hot` and
`fold_enrichment_defined_at_r_hot` on `handoff_02`; `hotspot_coverage_at_r_hot` and
`hotspot_fold_enrichment_at_r_hot` in the run-summary metrics, positioned immediately
after the count, and registered in the review-pack metric list.

Not added to `HANDOFF_REQUIRED_KEYS`: the UNDERPOWERED and NEGATIVE terminal paths have
no selected radius, so the two quantities are undefined there, and the count they would
qualify is 0 — no claim is being made for them to qualify.

## 2. An unanalysable cohort is a negative, not a crash

Stage A, on a cohort with an empty binary class, does not raise. It records
`BLOCKED / SCIENTIFIC_NEGATIVE` in its own `stage_status.json`, puts a
`negative_result` on `handoff_01`, sets `qc_status = FAIL` and returns — its own
interpretation text states "hotspot-statistics must NOT run".

The orchestrator already has a hook for that route, and its comment names the defect
class verbatim: *"A stage may DECLARE a settled scientific negative in its handoff and
return normally... Only the raised route was handled, so a stage that completed
COMPLETED_NEGATIVE still had its successor invoked; that successor correctly refused to
build on an empty result, but it refused with a BlockedError, which types as
TECHNICAL_FAILURE."*

The hook keys on `negative_result["settled_negative"]`, and **only Stage B was setting
that flag.** So Stage A's case fell through: the run reached Stage B, `Handoff.verify()`
raised the generic "upstream reports `qc_status=FAIL`", and because `error` outranks
`negative` in `_finalize`, the run was typed `BLOCKED / TECHNICAL_FAILURE` with
`negative_result: null`.

**WT1 is the case:** classified cohort `N = 1, N_P = 0, N_B = 1`. A statement about
ClinVar's content for that gene, filed as a pipeline fault.

**Implemented:** `"settled_negative": True` on Stage A's negative dict — one key, using
the mechanism that already exists. A first attempt added a parallel reinterpretation in
`_finalize` and was reverted: a second mechanism for the same seam is how the two drift
apart.

**Guarded:** `tests/integration/test_end_to_end.py::test_unanalysable_cohort_is_a_negative_not_a_technical_failure`
drives a real `run_pipeline` on a cohort with the benign class dropped. Verified to FAIL
without the flag (`RunStatus.BLOCKED`, error "upstream handoff_01 reports
qc_status=FAIL") and pass with it. The module-level invariant already in that file — "a
scientific negative is never reported as a technical failure" — is guarded by
`if result.negative_result is not None`, which is precisely what this defect nulls out,
so it could not have caught this.

Not placed in `test_underpowered_terminal_state.py`, where the sibling defect is
documented: that whole module is currently `pytest.mark.skip` pending an UNDERPOWERED
fixture reachable at the frozen `B = 100,000`, so a test added there would never run.

## Addendum — four cleanups, same class

None changes a number, a decision or a terminal state of any gene that produced a
result. All four were found in the Stage A->E assessment of the `results_v2/` batch.

**1. `src/clinspot3d/` deleted.** 1.3 MB of orphaned bytecode from a previous package
name: 100 `.pyc` files, zero `.py`, untracked and invisible to `git status` because
`*.pyc` is ignored, last written 2026-09-01. It imported as a namespace package while
every submodule raised `ModuleNotFoundError`, so it could not shadow `sphare3d` and
was outside `code_version`'s root — harmless, but a trap for any future
`import clinspot3d`. `import clinspot3d` now fails cleanly. The stale
`_editable_impl_clinspot3d.pth` remains in the venv; it points at `src`, which still
exists for `sphare3d`, so it is inert, and `.venv` is machine-local rather than repo
state.

**2. Stage E reported one fact as two warnings.** `LITERATURE_MECHANISMS_NOT_RUN` and
`POSTHOC_GATE_FAILED` each fired on **49 of 49** completed runs, and the second is only
the first's consequence: with the curated source out of scope, every mechanism label
defaults to `Not_Experimentally_Characterized`, so no variant reaches any evidence level
and the gate has nothing to weigh.

Worse than redundant, the record was wrong about its own cause. `gate_failure_record`
wrote `outcome_type: SCIENTIFIC_NEGATIVE` with the detail "it says the literature does
not yet support the test" — but the literature was never searched. That is exactly the
misreport `DECISION-STAGE-E-SCOPE-0001` fixed one function over, when all-zero
`mechanism_counts` read as "we curated and found nothing". That record's own text
anticipates this case: it states the post hoc gate "therefore cannot pass". **A gate that
cannot pass has not failed; it was not evaluable.**

So `gate_failure_record` now takes `mechanisms_not_run` and, when set, writes
`status: NOT_EVALUABLE`, `outcome_type: NOT_APPLICABLE`, condition
`MECHANISM_SOURCE_NOT_QUERIED` and `upstream_cause: LITERATURE_MECHANISMS_NOT_RUN`. The
second warning is suppressed on that path only — nothing is lost, because the record
names the cause and `LITERATURE_MECHANISMS_NOT_RUN` still fires. When a curated source
is eventually supplied, the genuine-negative branch is untouched and reports as before.

**3. An escalation was stamped `TECHNICAL_FAILURE` by the stage that raised it.** Stage A
caught `BlockedError` and `EscalationRequired` in one `except` and wrote
`BLOCKED / TECHNICAL_FAILURE` for both. The orchestrator types an escalated run
`PARTIAL / NOT_APPLICABLE`, so the stage contradicted its own run: **BTD** (AlphaFold and
UniProt numbering disagree by -20 positions, alignment identity 1.0000 over 523 residues)
and **GLA** (gene resolves to two reviewed UniProt accessions, Q9UHE5 and P06280) both
read `BLOCKED/TECHNICAL_FAILURE` at the stage and `PARTIAL/NOT_APPLICABLE` at the run.
`_mark_not_run` cannot correct it, because it skips a status file the stage has already
written.

The branches are now split: `EscalationRequired` writes `BLOCKED / NOT_APPLICABLE` and a
`recommended_action` saying a Lead DECISION is required rather than "resolve it and
re-run"; `BlockedError` is unchanged. An escalation is an ambiguity in the inputs, not a
fault in the code, and the two must not be read the same way. This also corrects a claim
made earlier in the assessment — that CYP27A1 and BTD/GLA showed "one condition typed two
ways". They are different conditions; what was wrong was the label on one of them.

**4. `n_binary_records_matched_component_wise` was zero by configuration, not by data.**
It counts records admitted as binary whose significance matched component-wise. Under the
frozen `compound_significance_policy: whole_field`, `component_matching` is False, so a
compound label such as "Pathogenic, risk factor" classifies as `OTHER` and is never
eligible — the path that could make the counter non-zero cannot execute. It read 0 in all
66 genes, which reads as "we checked and found none". It is now `None` when the policy
makes it unmeasurable, with a note naming the policy; the same correction, for the same
reason, as zeros -> null in `DECISION-STAGE-E-SCOPE-0001`.

## Deliberately NOT done

**`BH_BY_DISAGREEMENT` severity.** It fires as MAJOR whenever BH and BY reject different
counts. BY is strictly more conservative (batch totals: BH 4123 centres, BY 2661, ratio
0.65), so it fires in **every gene that finds anything** — 50 of 50. Of the 13 genes
where it does not fire, 12 found nothing and the 13th is MAX (2 = 2, the knife-edge
gene). A MAJOR that always fires desensitises readers to the ones that matter, and this
repository has already made exactly this downgrade once, for
`LITERATURE_MECHANISMS_NOT_RUN`, for exactly this reason.

It is nonetheless **left for the user to decide**, because of what it would do: three
genes have `BH_BY_DISAGREEMENT` as their ONLY MAJOR — ATP7B, SCN5A and **KCNA2, the
manuscript's case gene**. Downgrading it moves those three from
`COMPLETED_WITH_WARNINGS` to `COMPLETED`. The change is defensible on its merits and its
effect is to make the case gene's headline status cleaner; that is not a call the
implementer should make unilaterally.

**Carrying a "lost" failure reason into the run summary — withdrawn, it was wrong.**
BTD and GLA have an empty `error` field in `run_summary.json`, which was read as the
reason having been lost. It is not: both raised `EscalationRequired`, and the full
cause is in the summary's `escalations` list — BTD "the AlphaFold model and the UniProt
canonical sequence disagree about residue numbering by -20 positions (alignment
identity 1.0000 over 523 residues)", GLA "gene 'GLA' resolves to 2 reviewed UniProt
accessions: ['Q9UHE5', 'P06280']". An escalation is not an error, and it is already
reported in the field that is for it. The related claim that CYP27A1 and BTD/GLA show
"the same stage condition typed two different ways" was also wrong: CYP27A1 raised a
`BlockedError` (technical) and BTD/GLA raised escalations (a Lead decision). Those are
different conditions, correctly distinguished.

One small untidiness does survive that retraction, recorded and not fixed: on the
escalation route the stage writes `TECHNICAL_FAILURE` into its own `stage_status.json`
while the orchestrator types the run `NOT_APPLICABLE`. The reason is not lost, but the
two labels disagree.
