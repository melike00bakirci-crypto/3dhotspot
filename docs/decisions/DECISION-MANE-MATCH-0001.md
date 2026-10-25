# DECISION-MANE-MATCH-0001

**Status:** IN PROGRESS
**Supersedes:** nothing. This is a defect repair, not a methodology change.
**Scope:** Stage A transcript eligibility only.

## The defect

`on_mane_transcript` is `FALSE` for every ClinVar record in every gene, in every run
ever produced. Measured on `results_v2/`:

| gene | records | `on_mane_transcript = TRUE` |
|---|---|---|
| ABCD1 | 1780 | 0 |
| BRCA1 | 15331 | 0 |
| SCN5A | 4999 | 0 |

The cause is a type mismatch across two sources that was never compared in a test:

* UniProt's `MANE-Select` cross-reference puts an **Ensembl** accession in its `id`
  field — ABCD1: `ENST00000218104.6`. `sources.parse_entry` assigns that to
  `GeneResolution.mane_transcript`.
* ClinVar's `Name` field names a **RefSeq** accession — ABCD1: `NM_000033.4`.
* `clinvar.transcript_matches` strips versions and compares the two strings, so it
  compares `ENST00000218104` against `NM_000033`. It cannot ever be true.

**What the pipeline has actually been filtering on.** With `on_mane` always false, the
only surviving eligibility test is step 3 of `evaluate_records` — the reference amino
acid must match the UniProt canonical sequence at that position. Where the MANE RefSeq
transcript and the UniProt canonical agree on numbering (the common case) this costs
nothing: ABCD1 and BRCA1 record **zero** `non_canonical_transcript` exclusions. Where
they disagree it is a weaker filter than intended: SCN5A excludes 1178 records that
way, and any record from a non-MANE transcript whose reference AA happens to match at
that position is silently admitted.

**Why no test caught it.** The fixtures encode the misconception. Both places that
build a synthetic UniProt entry put a RefSeq accession in the MANE-Select `id`:

```
tests/data_structure/test_live_source.py:63
    {"database": "MANE-Select", "id": "NM_000000.1"}
tests/data_structure/test_stage_a.py:470
    {"database": "MANE-Select", "id": "NM_004333.6"}
```

and `test_stage_a.py:475` then *asserts* `resolution.mane_transcript == "NM_004333.6"`.
Real UniProt never returns that shape. So the suite compared `NM_` against `NM_`,
matched, and reported green while production compared `ENST` against `NM_`.

## Decision

Compare the record's transcript against **either identifier of the MANE pair**, taking
the RefSeq side from the same cross-reference that supplies the Ensembl side.

UniProt's MANE-Select xref carries the pairing directly (verified live against
`rest.uniprot.org` for P33897):

```
db = MANE-Select   id = ENST00000218104.6
    ProteinId           = ENSP00000218104.3
    RefSeqNucleotideId  = NM_000033.4
    RefSeqProteinId     = NP_000024.2
```

So `GeneResolution` gains `mane_refseq_nucleotide`, read from that property, and a
record is on-MANE when its accession matches the Ensembl id **or** the RefSeq
nucleotide id, versions ignored as before.

**Why this source and not another.** MANE *is* the ENST/NM_ pairing, so the xref is the
definition rather than a lookup that could disagree with it. It arrives in a payload the
pipeline already fetches, so there is no new network dependency, no new failure mode,
and no third-party mapping table to version. `mane_transcript` keeps reporting the
Ensembl accession, so nothing downstream changes shape.

Rejected: resolving `NM_` -> `ENST` through an external MANE summary table (a second
network source and a second thing to pin), and abandoning the MANE claim to document the
reference-AA rule as the real methodology (defensible, but it would require rewriting the
manuscript's transcript sentence and would leave the weaker filter in place).

## Second defect in the same function, fixed with it

The fallback path — taken when a UniProt entry names no MANE-Select xref — reads

```python
transcript = sorted(x.get("id") or "" for x in refseq)[-1]
```

A RefSeq cross-reference's `id` is the **protein** accession (`NP_000024.2`); the
nucleotide accession is in its `NucleotideSequenceId` property. So the fallback assigns
an `NP_` to `transcript`, which can no more match ClinVar's `NM_` than the Ensembl id
could. Its comment also claims "the RefSeq transcript matching the canonical sequence
with the longest CDS", but the code sorts accession strings and takes the last — CDS
length is not available in this payload and is never consulted.

**Dormant, not active:** all 66 genes in `results_v2/` resolved through `MANE_SELECT`
with `mane_fallback_used = false`, so this path has never run in a recorded batch. It is
corrected because the repair touches the same six lines, and the source label is changed
to describe what the code does instead of what it does not.

## What this changes — measured, and smaller than it first looks

**The cohort does not change at all.** `on_mane` is read in exactly two places in
`evaluate_records`, and neither can turn an admitted record into an excluded one:

* on the `not ref_match` branch it selects the *reason label* —
  `ref_aa_mismatch` when on-MANE, `non_canonical_transcript` otherwise. Both call
  `ev.exclude(...)`, so the record is excluded either way.
* on the branch where position and reference AA both agree, `not on_mane` only sets
  `remapped_to_canonical = True`. The record is admitted either way.

So eligibility is decided by the reference-AA rule alone, before and after this repair.
Measured over all 200,074 missense records in `results_v2/`:

| | records |
|---|---|
| admission changes | **0** |
| exclusion reason may flip `non_canonical_transcript` -> `ref_aa_mismatch` | 4,767 |
| false `remapped_to_canonical = TRUE` replaced by `on_mane_transcript = TRUE` | 85,313 |

`N_P`, `N_B`, `|L|`, every centre, every radius and every footprint are therefore
untouched, and no terminal state can move. **This repair does not gate the batch
re-run.** That corrects two earlier statements: the 2026-09-13 audit's warning that
re-running before A1 is fixed "would reproduce the errors", and this record's own first
draft, which claimed the filter would start biting and could push a thin-benign gene
such as LMNA (`N_B = 2`) into UNDERPOWERED. Neither holds — the ref-AA rule admits
exactly the same records with or without a working MANE comparison.

What the repair does buy, which is still worth having:

1. `CODE_AVAILABILITY.md`'s claim that variants "were mapped onto the MANE Select
   transcript" becomes true instead of aspirational.
2. `remapped_to_canonical` stops being a false positive on 85,313 records. Today it
   says "this record was remapped from a non-canonical transcript" about records that
   sit on the MANE transcript and were never remapped — the provenance trail asserts an
   operation that did not occur.
3. The two exclusion reasons become distinguishable, so "this record is on the right
   transcript but its reference AA disagrees" stops being reported as "this record is on
   the wrong transcript".
4. The eligibility gate becomes a real gate, so a future tightening (e.g. refusing
   remapping) has something correct to tighten.

A re-run is needed only to regenerate these labels in the artefacts; it is not needed
for any number.

## Out of scope

* `p_emp` vs the exact hypergeometric p as BH's input — a separate, undecided matter.
* The reference-AA rule itself, the significance policy, the star policy.
* `CODE_AVAILABILITY.md`'s transcript sentence becomes true once this lands; the text
  still needs re-checking against the gene it cites.
