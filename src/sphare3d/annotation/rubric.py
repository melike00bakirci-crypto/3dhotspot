"""The FROZEN evidence rubric — METHOD_SPEC II.13 (A23).

One deterministic function decides every mechanism label. It is isolated here, with
no I/O and no knowledge of hotspots, for two reasons:

  * it is the only place a mechanism may be assigned, so it can be unit-tested
    exhaustively over the evidence truth table; and
  * it *cannot* see whether a variant lies inside a hotspot, which makes the
    circularity guardrail (agent §8) structural rather than a matter of discipline.

Frozen rules (config ``annotation.*``)::

    label requires          >= MODERATE
    conflicting STRONG      -> Mixed          (both references retained)
    WEAK-only / conflicting-WEAK -> Unclear
    no evidence (NONE)      -> Not_Experimentally_Characterized

``Not_Experimentally_Characterized`` is a real category and is never filled by
inference. Evidence is never upgraded to populate a pattern.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from ..utils.errors import BlockedError

# --- frozen vocabularies ----------------------------------------------------


class Evidence(str, Enum):
    """Graded evidence strength. Order matters; see ``RANK``."""

    STRONG = "STRONG"
    MODERATE = "MODERATE"
    WEAK = "WEAK"
    NONE = "NONE"


class Mechanism(str, Enum):
    GOF = "GOF"
    LOF = "LOF"
    DN = "DN"
    MIXED = "Mixed"
    UNCLEAR = "Unclear"
    NOT_CHARACTERIZED = "Not_Experimentally_Characterized"


RANK: dict[Evidence, int] = {
    Evidence.NONE: 0,
    Evidence.WEAK: 1,
    Evidence.MODERATE: 2,
    Evidence.STRONG: 3,
}

#: Mechanisms that actually *assert* a functional effect. ``Unclear`` and
#: ``Not_Experimentally_Characterized`` are outcomes of the rubric, not claims a
#: source can contribute, so a source reporting them adds evidence but no assertion.
ASSERTABLE: frozenset[Mechanism] = frozenset(
    {Mechanism.GOF, Mechanism.LOF, Mechanism.DN, Mechanism.MIXED}
)

#: Rubric definitions, reproduced verbatim into the report so the grade a reader
#: sees is the grade that was applied.
RUBRIC_DEFINITIONS: dict[Evidence, str] = {
    Evidence.STRONG: (
        "Direct assay of the specific variant, relevant system, controlled "
        "comparison, peer-reviewed primary source"
    ),
    Evidence.MODERATE: (
        "Direct assay of the specific variant but limited — single assay, "
        "non-physiological system, or uncontrolled"
    ),
    Evidence.WEAK: (
        "Indirect — a different variant at the same residue, a computational "
        "prediction reported in the source, or clinical inference"
    ),
    Evidence.NONE: "No experimental evidence -> Not_Experimentally_Characterized",
}


def parse_evidence(value: str | Evidence) -> Evidence:
    if isinstance(value, Evidence):
        return value
    try:
        return Evidence(str(value).strip().upper())
    except ValueError as exc:
        raise BlockedError(
            f"BLOCKED — unknown evidence grade {value!r}. The rubric is FROZEN at "
            f"{[e.value for e in Evidence]}; an unrecognised grade is never coerced."
        ) from exc


def parse_mechanism(value: str | Mechanism) -> Mechanism:
    if isinstance(value, Mechanism):
        return value
    text = str(value).strip()
    for mech in Mechanism:
        if mech.value.lower() == text.lower():
            return mech
    raise BlockedError(
        f"BLOCKED — unknown mechanism {value!r}. Categories are FROZEN at "
        f"{[m.value for m in Mechanism]}."
    )


# --- rubric configuration (read from the frozen config, never hard-coded) ----


@dataclass(frozen=True)
class Rubric:
    """The rubric as configured. Absent config is BLOCKING (agent §6.5)."""

    min_evidence_for_label: Evidence
    conflicting_strong_resolution: Mechanism
    weak_only_resolution: Mechanism
    categories: tuple[Mechanism, ...]

    #: Not in the frozen text. The methodology fixes conflicting-STRONG -> Mixed and
    #: conflicting-WEAK -> Unclear but is silent on conflicting-MODERATE. We take the
    #: conservative reading — ``Mixed`` is an affirmative claim that a variant genuinely
    #: has two effects, and two *limited* assays do not support it, so the conflict
    #: resolves to ``Unclear``. Every row decided this way is flagged in
    #: ``curator_note`` and reported to the Lead as an open rubric point.
    conflicting_moderate_resolution: Mechanism = Mechanism.UNCLEAR

    @classmethod
    def from_config(cls, cfg) -> "Rubric":
        categories = tuple(
            parse_mechanism(c) for c in cfg.get("annotation.mechanism_categories")
        )
        levels = [parse_evidence(v) for v in cfg.get("annotation.evidence_levels")]
        if set(levels) != set(Evidence):
            raise BlockedError(
                f"BLOCKED — config annotation.evidence_levels is {levels}, but the "
                f"FROZEN rubric is {[e.value for e in Evidence]}."
            )
        return cls(
            min_evidence_for_label=parse_evidence(
                cfg.get("annotation.min_evidence_for_label")
            ),
            conflicting_strong_resolution=parse_mechanism(
                cfg.get("annotation.conflicting_strong_resolution")
            ),
            weak_only_resolution=parse_mechanism(
                cfg.get("annotation.weak_only_resolution")
            ),
            categories=categories,
        )


# --- evidence records -------------------------------------------------------


@dataclass
class MechanismRecord:
    """One curated literature record about one variant.

    This is the raw evidence unit. ``mechanism`` is what *the source reports*; the
    label the pipeline assigns is decided only by :func:`resolve_mechanism`.
    """

    residue_index: int
    aa_change: str
    mechanism: Mechanism
    reference: str
    experimental_system: str
    assay_type: str
    principal_finding: str
    evidence_strength: Evidence
    # retrieval provenance — makes literature_search_log.tsv reproducible by construction
    query: str = "NA"
    database: str = "NA"
    search_date: str = "NA"
    hit_count: int = 0
    included: bool = True
    exclusion_reason: str = "NA"

    def __post_init__(self) -> None:
        self.mechanism = parse_mechanism(self.mechanism)
        self.evidence_strength = parse_evidence(self.evidence_strength)

    @property
    def rank(self) -> int:
        return RANK[self.evidence_strength]


@dataclass
class Resolution:
    """The rubric's verdict for one variant, with the reasoning kept attached."""

    mechanism: Mechanism
    evidence_strength: Evidence
    references: list[str] = field(default_factory=list)
    curator_note: str = "NA"
    conflict: bool = False
    rubric_gap: bool = False

    @property
    def is_labelled(self) -> bool:
        """True when an affirmative mechanism was assigned."""
        return self.mechanism in ASSERTABLE


# --- the rubric -------------------------------------------------------------


def resolve_mechanism(records: list[MechanismRecord], rubric: Rubric) -> Resolution:
    """Apply the FROZEN rubric to all evidence for a single variant.

    Deterministic and order-independent: records are considered by evidence tier,
    and ties are broken by the frozen category order, never by input ordering.
    Contradicting sources are *retained*, never dropped (agent §16).
    """
    usable = [r for r in records if r.included]
    if not usable:
        return Resolution(
            mechanism=Mechanism.NOT_CHARACTERIZED,
            evidence_strength=Evidence.NONE,
            curator_note="No included literature record.",
        )

    with_evidence = [r for r in usable if r.evidence_strength is not Evidence.NONE]
    if not with_evidence:
        return Resolution(
            mechanism=Mechanism.NOT_CHARACTERIZED,
            evidence_strength=Evidence.NONE,
            references=_refs(usable),
            curator_note=(
                "Records retrieved but none carries experimental evidence "
                "(all NONE); category assigned by rubric, never by inference."
            ),
        )

    asserting = [r for r in with_evidence if r.mechanism in ASSERTABLE]
    if not asserting:
        return Resolution(
            mechanism=Mechanism.UNCLEAR,
            evidence_strength=max(
                (r.evidence_strength for r in with_evidence), key=lambda e: RANK[e]
            ),
            references=_refs(with_evidence),
            curator_note=(
                "Evidence exists but no source asserts a mechanism category."
            ),
        )

    top_rank = max(r.rank for r in asserting)
    top_tier = next(e for e in Evidence if RANK[e] == top_rank)
    at_top = [r for r in asserting if r.rank == top_rank]
    distinct = _ordered_unique(r.mechanism for r in at_top)

    # Below the labelling bar: WEAK-only, or conflicting WEAK -> Unclear (frozen).
    if top_rank < RANK[rubric.min_evidence_for_label]:
        return Resolution(
            mechanism=rubric.weak_only_resolution,
            evidence_strength=top_tier,
            references=_refs(asserting),
            curator_note=(
                f"Highest available evidence is {top_tier.value}, below the frozen "
                f"{rubric.min_evidence_for_label.value} bar for a mechanism label."
            ),
            conflict=len(distinct) > 1,
        )

    dissent = _dissent_note(asserting, at_top, distinct)

    if len(distinct) == 1:
        return Resolution(
            mechanism=distinct[0],
            evidence_strength=top_tier,
            references=_refs(asserting),
            curator_note=dissent,
            conflict=bool(dissent != "NA"),
        )

    # Conflict at or above the labelling bar.
    if top_tier is Evidence.STRONG:
        return Resolution(
            mechanism=rubric.conflicting_strong_resolution,
            evidence_strength=Evidence.STRONG,
            references=_refs(asserting),
            curator_note=(
                "Conflicting STRONG sources ("
                + ", ".join(m.value for m in distinct)
                + "); frozen rubric resolves to Mixed with all references retained."
            ),
            conflict=True,
        )

    return Resolution(
        mechanism=rubric.conflicting_moderate_resolution,
        evidence_strength=top_tier,
        references=_refs(asserting),
        curator_note=(
            "Conflicting "
            + top_tier.value
            + " sources ("
            + ", ".join(m.value for m in distinct)
            + "); the frozen rubric specifies conflicting-STRONG and conflicting-WEAK "
            "only. Resolved conservatively to "
            + rubric.conflicting_moderate_resolution.value
            + " because Mixed is an affirmative dual-mechanism claim that limited "
            "assays do not support. OPEN RUBRIC POINT — escalated to the Lead."
        ),
        conflict=True,
        rubric_gap=True,
    )


def _refs(records: list[MechanismRecord]) -> list[str]:
    return _ordered_unique(r.reference for r in records)


def _dissent_note(
    asserting: list[MechanismRecord],
    at_top: list[MechanismRecord],
    distinct: list[Mechanism],
) -> str:
    """Record lower-tier disagreement instead of discarding it."""
    lower = [r for r in asserting if r not in at_top and r.mechanism not in distinct]
    if not lower:
        return "NA"
    return (
        "Higher-tier evidence decides; contradicting lower-tier source(s) retained: "
        + "; ".join(
            f"{r.mechanism.value}@{r.evidence_strength.value} ({r.reference})"
            for r in lower
        )
    )


def _ordered_unique(values) -> list:
    """Deterministic de-duplication preserving first appearance."""
    seen: list = []
    for value in values:
        if value not in seen:
            seen.append(value)
    return seen
