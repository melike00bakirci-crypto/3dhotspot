"""Record evaluation, residue-level collapse and the four ClinVar strata (F12/F3).

This is the accounting core of Stage A. Every retrieved record ends up in exactly
one place:

  stratum 1  ``variants_missense_all.tsv``            every retrieved record
  stratum 2  ``variants_residue_level.tsv``           residues (the collapse, ¶3)
  stratum 3  ``variants_excluded_from_primary.tsv``   every exclusion + reason
  stratum 4  ``review_star_distribution.tsv``         stars, described never used

and the identity ``stratum 1 = eligible(stratum 2) union stratum 3`` is proven in
code, not asserted in prose. Review stars are carried through all four strata and
key no decision anywhere (F12).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Sequence

from ..utils.errors import NegativeResult
from .classify import (
    CLASS_BLB,
    CLASS_CONFLICT,
    CLASS_PLP,
    REASON_NO_CA,
    REASON_NON_CANONICAL_TX,
    REASON_NON_MISSENSE_VARIANT_TYPE,
    REASON_REF_AA_MISMATCH,
    REASON_RESIDUE_CONFLICT,
    REASON_UNPARSEABLE_HGVS,
    SIG_BLB,
    SIG_PLP,
    SignificanceCall,
    SignificancePolicy,
    non_missense_exclusion_reason,
    residue_class,
)
from .clinvar import ClinVarRecord, transcript_matches
from .review_status import star_or_na, star_sort_key

# Reason precedence: a record is reported under the FIRST fact that made it
# unusable, reading from "cannot be placed on a residue" to "placed but not
# binary". Fixed here so the strata are deterministic and reproducible; the
# untouched significance_class column keeps every other fact recoverable.
REASON_PRECEDENCE: tuple[str, ...] = (
    REASON_UNPARSEABLE_HGVS,
    REASON_NON_MISSENSE_VARIANT_TYPE,
    REASON_NON_CANONICAL_TX,
    REASON_REF_AA_MISMATCH,
    "vus",
    "conflicting",
    "non_binary_significance",
    REASON_RESIDUE_CONFLICT,
    REASON_NO_CA,
)

STRATUM1_COLUMNS: tuple[str, ...] = (
    "record_id", "variation_id", "allele_id", "gene_symbol", "name", "variant_type",
    "clinical_significance", "significance_class", "significance_exact_match",
    "significance_modifiers", "review_status", "star_level",
    "transcript", "on_mane_transcript", "remapped_to_canonical",
    "hgvs_c", "hgvs_p", "aa_ref_hgvs", "aa_alt_hgvs", "residue_index",
    "aa_ref_uniprot", "ref_aa_match", "condition_text", "last_evaluated",
    "n_submitters", "eligible_primary", "exclusion_reason", "exclusion_detail",
)

STRATUM2_COLUMNS: tuple[str, ...] = (
    "residue_index", "aa_ref", "class", "n_records_plp", "n_records_blb",
    "clinvar_ids", "star_levels", "max_star", "eligible_primary",
)

STRATUM3_COLUMNS: tuple[str, ...] = (
    "record_id", "variation_id", "residue_index", "clinical_significance",
    "significance_class", "review_status", "star_level", "transcript", "hgvs_p",
    "exclusion_reason", "exclusion_detail",
)

STRATUM4_COLUMNS: tuple[str, ...] = (
    "review_status", "star_level", "class", "inclusion_stratum", "n_records",
)

CONFLICT_COLUMNS: tuple[str, ...] = (
    "residue_index", "aa_ref", "n_records_plp", "n_records_blb",
    "clinvar_ids_plp", "clinvar_ids_blb", "star_levels_plp", "star_levels_blb",
    "max_star", "exclusion_reason", "policy",
)

ELIGIBLE_STRATUM = "eligible_primary"
EXCLUDED_STRATUM = "excluded_from_primary"
ALL_GROUP = "ALL"


# --- per-record evaluation --------------------------------------------------

@dataclass
class RecordEval:
    """A normalized record plus every Stage A judgement made about it."""

    record: ClinVarRecord
    significance: SignificanceCall
    on_mane_transcript: bool
    remapped_to_canonical: bool
    residue_index: int | None
    aa_ref_uniprot: str | None
    ref_aa_match: bool | None
    eligible_primary: bool
    exclusion_reason: str | None = None
    exclusion_detail: str = "NA"

    # -- exclusion is always additive and always reasoned --------------------
    def exclude(self, reason: str, detail: str = "NA") -> None:
        """Record an exclusion. Never overwrites an earlier, higher-precedence one."""
        incoming = REASON_PRECEDENCE.index(reason)
        current = (REASON_PRECEDENCE.index(self.exclusion_reason)
                   if self.exclusion_reason in REASON_PRECEDENCE else len(REASON_PRECEDENCE))
        self.eligible_primary = False
        if incoming < current:
            self.exclusion_reason = reason
            self.exclusion_detail = detail

    @property
    def record_id(self) -> str:
        return self.record.record_id

    @property
    def binary_class(self) -> str | None:
        return {SIG_PLP: CLASS_PLP, SIG_BLB: CLASS_BLB}.get(
            self.significance.significance_class)

    def stratum1_row(self) -> dict:
        rec = self.record
        return {
            "record_id": rec.record_id, "variation_id": rec.variation_id,
            "allele_id": rec.allele_id, "gene_symbol": rec.gene_symbol,
            "name": rec.name, "variant_type": rec.variant_type,
            "clinical_significance": rec.clinical_significance,
            "significance_class": self.significance.significance_class,
            "significance_exact_match": self.significance.exact_match,
            "significance_modifiers": list(self.significance.modifiers) or None,
            "review_status": rec.review_status, "star_level": rec.star_level,
            "transcript": rec.transcript, "on_mane_transcript": self.on_mane_transcript,
            "remapped_to_canonical": self.remapped_to_canonical,
            "hgvs_c": rec.hgvs_c, "hgvs_p": rec.hgvs_p,
            "aa_ref_hgvs": rec.aa_ref_hgvs, "aa_alt_hgvs": rec.aa_alt_hgvs,
            "residue_index": self.residue_index,
            "aa_ref_uniprot": self.aa_ref_uniprot, "ref_aa_match": self.ref_aa_match,
            "condition_text": rec.condition_text, "last_evaluated": rec.last_evaluated,
            "n_submitters": rec.n_submitters,
            "eligible_primary": self.eligible_primary,
            "exclusion_reason": self.exclusion_reason,
            "exclusion_detail": self.exclusion_detail,
        }

    def stratum3_row(self) -> dict:
        rec = self.record
        return {
            "record_id": rec.record_id, "variation_id": rec.variation_id,
            "residue_index": self.residue_index,
            "clinical_significance": rec.clinical_significance,
            "significance_class": self.significance.significance_class,
            "review_status": rec.review_status, "star_level": rec.star_level,
            "transcript": rec.transcript, "hgvs_p": rec.hgvs_p,
            "exclusion_reason": self.exclusion_reason,
            "exclusion_detail": self.exclusion_detail,
        }


def evaluate_records(records: Sequence[ClinVarRecord], *, policy: SignificancePolicy,
                     canonical_sequence: str, mane_transcript: str,
                     require_ref_aa_match: bool = True) -> list[RecordEval]:
    """Evaluate every retrieved record against the frozen inclusion policy.

    Order of judgement (METHOD_SPEC II.13):
      1. protein HGVS must parse to a missense substitution;
      2. the record must sit on the MANE Select transcript, or be remappable —
         remapping is permitted **only** when protein position and reference AA
         are both consistent with the UniProt canonical sequence;
      3. the reference AA must match the canonical sequence;
      4. the significance must be P/LP or B/LB.
    No step looks at review status or star level.
    """
    length = len(canonical_sequence)
    out: list[RecordEval] = []

    for rec in records:
        call = policy.classify(rec.clinical_significance)
        on_mane = transcript_matches(rec.transcript, mane_transcript)
        position = rec.residue_index
        aa_ref_uniprot = (canonical_sequence[position - 1]
                          if position is not None and 1 <= position <= length else None)
        ref_match = (None if (aa_ref_uniprot is None or rec.aa_ref_hgvs is None)
                     else aa_ref_uniprot == rec.aa_ref_hgvs)

        ev = RecordEval(
            record=rec, significance=call, on_mane_transcript=on_mane,
            remapped_to_canonical=False,
            residue_index=position if rec.is_missense else None,
            aa_ref_uniprot=aa_ref_uniprot, ref_aa_match=ref_match,
            eligible_primary=True,
        )

        # 1 — usable missense protein change. A change that parsed CLEANLY into
        # a non-missense consequence (synonymous, nonsense, frameshift, in-frame
        # indel/dup/ext) is reported under REASON_NON_MISSENSE_VARIANT_TYPE, never
        # aggregated with a genuine parse failure (Workflow v2 §1).
        if not rec.is_missense:
            detail = rec.change.reason or "no_protein_hgvs"
            ev.exclude(non_missense_exclusion_reason(detail), detail)

        # 2 / 3 — transcript and reference-AA consistency
        elif aa_ref_uniprot is None:
            ev.exclude(REASON_REF_AA_MISMATCH,
                       f"position_out_of_range:{position}>{length}")
        elif not ref_match:
            reason = REASON_REF_AA_MISMATCH if on_mane else REASON_NON_CANONICAL_TX
            ev.exclude(reason,
                       f"hgvs_ref={rec.aa_ref_hgvs};uniprot_ref={aa_ref_uniprot};"
                       f"transcript={rec.transcript or 'NA'}")
        elif not on_mane:
            # Remap permitted: position AND reference AA both agree with canonical.
            ev.remapped_to_canonical = True

        # 4 — binary significance
        if not call.eligible_primary:
            ev.exclude(call.exclusion_reason,
                       f"significance_class={call.significance_class}")

        if require_ref_aa_match is False:      # pragma: no cover - frozen TRUE
            raise ValueError(
                "transcript.require_ref_aa_match is FROZEN at TRUE; disabling the "
                "reference-AA check is a methodological revision, not a run option."
            )
        out.append(ev)

    return out


# --- residue-level collapse -------------------------------------------------

@dataclass
class ResidueRow:
    """One residue after the ¶3 collapse. Class is assigned only by F3."""

    residue_index: int
    aa_ref: str
    residue_class: str
    n_records_plp: int
    n_records_blb: int
    evals_plp: list[RecordEval] = field(default_factory=list)
    evals_blb: list[RecordEval] = field(default_factory=list)
    eligible_primary: bool = True
    exclusion_reason: str | None = None

    @property
    def evals(self) -> list[RecordEval]:
        return sorted(self.evals_plp + self.evals_blb, key=lambda e: e.record_id)

    @property
    def clinvar_ids(self) -> list[str]:
        return [e.record.variation_id for e in self.evals]

    @property
    def star_levels(self) -> list[int | None]:
        return [e.record.star_level for e in self.evals]

    @property
    def max_star(self) -> int | None:
        known = [s for s in self.star_levels if s is not None]
        return max(known) if known else None

    def stratum2_row(self) -> dict:
        return {
            "residue_index": self.residue_index, "aa_ref": self.aa_ref,
            "class": self.residue_class,
            "n_records_plp": self.n_records_plp, "n_records_blb": self.n_records_blb,
            "clinvar_ids": self.clinvar_ids,
            "star_levels": [star_or_na(s) for s in self.star_levels],
            "max_star": self.max_star,
            "eligible_primary": self.eligible_primary,
        }

    def conflict_row(self) -> dict:
        return {
            "residue_index": self.residue_index, "aa_ref": self.aa_ref,
            "n_records_plp": self.n_records_plp, "n_records_blb": self.n_records_blb,
            "clinvar_ids_plp": [e.record.variation_id for e in self.evals_plp],
            "clinvar_ids_blb": [e.record.variation_id for e in self.evals_blb],
            "star_levels_plp": [star_or_na(e.record.star_level) for e in self.evals_plp],
            "star_levels_blb": [star_or_na(e.record.star_level) for e in self.evals_blb],
            "max_star": self.max_star,
            "exclusion_reason": REASON_RESIDUE_CONFLICT,
            "policy": "exclude_from_primary;preserve_all_source_records",
        }


def collapse_to_residues(evals: Sequence[RecordEval], *,
                         canonical_sequence: str) -> list[ResidueRow]:
    """¶3 — collapse duplicate records onto residues and apply F3.

    A residue carrying both P/LP and B/LB evidence becomes RESIDUE_CLASS_CONFLICT:
    it keeps all of its source records, is excluded from the primary comparison,
    and is never arbitrarily assigned to either class.
    """
    buckets: dict[int, tuple[list[RecordEval], list[RecordEval]]] = {}
    for ev in evals:
        if not ev.eligible_primary or ev.residue_index is None:
            continue
        plp, blb = buckets.setdefault(ev.residue_index, ([], []))
        {CLASS_PLP: plp, CLASS_BLB: blb}[ev.binary_class].append(ev)

    rows: list[ResidueRow] = []
    for index in sorted(buckets):
        plp, blb = buckets[index]
        row = ResidueRow(
            residue_index=index, aa_ref=canonical_sequence[index - 1],
            residue_class=residue_class(len(plp), len(blb)),
            n_records_plp=len(plp), n_records_blb=len(blb),
            evals_plp=sorted(plp, key=lambda e: e.record_id),
            evals_blb=sorted(blb, key=lambda e: e.record_id),
        )
        if row.residue_class == CLASS_CONFLICT:
            row.eligible_primary = False
            row.exclusion_reason = REASON_RESIDUE_CONFLICT
            for ev in row.evals:
                ev.exclude(REASON_RESIDUE_CONFLICT, f"residue_index={index}")
        rows.append(row)
    return rows


def apply_structural_exclusions(rows: Sequence[ResidueRow],
                                unusable_residues: Iterable[int]) -> list[int]:
    """Exclude cohort residues that have no usable CA (F1). Never silent.

    Returns the residue indices actually excluded, which the mapping report and
    stratum 3 both carry.
    """
    unusable = set(unusable_residues)
    excluded: list[int] = []
    for row in rows:
        if not row.eligible_primary or row.residue_index not in unusable:
            continue
        row.eligible_primary = False
        row.exclusion_reason = REASON_NO_CA
        for ev in row.evals:
            ev.exclude(REASON_NO_CA, f"residue_index={row.residue_index}")
        excluded.append(row.residue_index)
    return excluded


# --- stratum 4: review-star description (never a filter) --------------------

def review_star_distribution(evals: Sequence[RecordEval]) -> list[dict]:
    """Stratum 4 — stars overall, per class, per inclusion stratum, and crossed.

    Emitted so a reader can see exactly what evidence quality the cohort rests on.
    It is a *description*. No row of this table has ever removed a record.
    """
    counts: dict[tuple[str, object, str, str], int] = {}

    def bump(status: str, level, klass: str, stratum: str) -> None:
        key = (status, level, klass, stratum)
        counts[key] = counts.get(key, 0) + 1

    for ev in evals:
        status = ev.record.review_status
        level = ev.record.star_level
        klass = ev.significance.significance_class
        stratum = ELIGIBLE_STRATUM if ev.eligible_primary else EXCLUDED_STRATUM
        bump(status, level, ALL_GROUP, ALL_GROUP)
        bump(status, level, klass, ALL_GROUP)
        bump(status, level, ALL_GROUP, stratum)
        bump(status, level, klass, stratum)

    def sort_key(item):
        (status, level, klass, stratum), _ = item
        return (star_sort_key(level), status, klass, stratum)

    return [
        {"review_status": status, "star_level": star_or_na(level), "class": klass,
         "inclusion_stratum": stratum, "n_records": n}
        for (status, level, klass, stratum), n in sorted(counts.items(), key=sort_key)
    ]


def star_summary(evals: Sequence[RecordEval]) -> dict:
    """Compact star tally for ``cohort_summary.json`` and the Lead's report."""
    tally: dict[str, int] = {}
    for ev in evals:
        key = str(star_or_na(ev.record.star_level))
        tally[key] = tally.get(key, 0) + 1
    return dict(sorted(tally.items()))


# --- conservation proof (QC rule 1) -----------------------------------------

@dataclass(frozen=True)
class ConservationProof:
    n_stratum1: int
    n_eligible: int
    n_excluded: int
    conserved: bool
    overlap: tuple[str, ...]
    unaccounted: tuple[str, ...]
    reasonless: tuple[str, ...]

    def as_dict(self) -> dict:
        return {
            "identity": "stratum1 == eligible(stratum2) union stratum3",
            "n_stratum1": self.n_stratum1, "n_eligible": self.n_eligible,
            "n_excluded": self.n_excluded, "conserved": self.conserved,
            "records_in_both": list(self.overlap),
            "records_unaccounted": list(self.unaccounted),
            "excluded_rows_without_reason": list(self.reasonless),
        }


def prove_conservation(evals: Sequence[RecordEval]) -> ConservationProof:
    """Prove record-count conservation across the strata (QC rule 1)."""
    all_ids = {e.record_id for e in evals}
    eligible = {e.record_id for e in evals if e.eligible_primary}
    excluded = {e.record_id for e in evals if not e.eligible_primary}
    reasonless = tuple(sorted(e.record_id for e in evals
                              if not e.eligible_primary and not e.exclusion_reason))
    overlap = tuple(sorted(eligible & excluded))
    unaccounted = tuple(sorted(all_ids - (eligible | excluded)))
    conserved = not overlap and not unaccounted and not reasonless
    return ConservationProof(
        n_stratum1=len(all_ids), n_eligible=len(eligible), n_excluded=len(excluded),
        conserved=conserved, overlap=overlap, unaccounted=unaccounted,
        reasonless=reasonless,
    )


# --- cohort assembly --------------------------------------------------------

@dataclass(frozen=True)
class CohortCounts:
    N: int
    N_P: int
    N_B: int
    n_conflict: int
    n_residues_total: int


def cohort_counts(rows: Sequence[ResidueRow]) -> CohortCounts:
    eligible = [r for r in rows if r.eligible_primary]
    return CohortCounts(
        N=len(eligible),
        N_P=sum(1 for r in eligible if r.residue_class == CLASS_PLP),
        N_B=sum(1 for r in eligible if r.residue_class == CLASS_BLB),
        n_conflict=sum(1 for r in rows if r.residue_class == CLASS_CONFLICT),
        n_residues_total=len(rows),
    )


def assert_cohort_viable(counts: CohortCounts) -> None:
    """QC rule 10 — ``|L| >= 2`` with BOTH classes non-empty.

    An empty class is a fact about ClinVar's content, not a bug: it raises
    :class:`NegativeResult`, which Stage A converts into a complete
    COMPLETED_NEGATIVE handoff. It is never repaired by admitting VUS or
    conflicting records, by extending to a non-canonical transcript, by
    reassigning a conflict residue or by introducing a star filter.
    """
    if counts.N_P == 0 or counts.N_B == 0 or counts.N < 2:
        raise NegativeResult(
            "INSUFFICIENT_CLASSIFIED_RESIDUES",
            f"classified cohort L is not analysable: N={counts.N}, N_P={counts.N_P}, "
            f"N_B={counts.N_B}, n_conflict={counts.n_conflict}. Both binary classes "
            f"must be non-empty and |L| >= 2.",
            N=counts.N, N_P=counts.N_P, N_B=counts.N_B, n_conflict=counts.n_conflict,
        )
