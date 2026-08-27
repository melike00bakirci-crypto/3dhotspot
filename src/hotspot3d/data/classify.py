"""Clinical-significance normalization — the ONLY place a binary class is assigned.

FROZEN policy (¶2, F12, ``config/clinvar``):
  * ``P`` + ``LP``  -> ``PLP``
  * ``B`` + ``LB``  -> ``BLB``
  * everything else -> ``eligible_primary = FALSE`` with a reason from
    ``clinvar.exclusion_reason_enum``, preserved in strata 1 and 3 and never
    relabelled into a binary class to enlarge a cohort.

The label vocabulary is not hard-coded here: it is read from the frozen config
(``clinvar.pathogenic_labels`` / ``clinvar.benign_labels``) and passed in.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Sequence

# --- residue classes (F3) ---------------------------------------------------
CLASS_PLP = "PLP"
CLASS_BLB = "BLB"
CLASS_CONFLICT = "CONFLICT"

# --- record-level significance classes --------------------------------------
SIG_PLP = "PLP"
SIG_BLB = "BLB"
SIG_VUS = "VUS"
SIG_CONFLICTING = "CONFLICTING"
SIG_OTHER = "OTHER"
SIG_RECORD_CONFLICT = "RECORD_SIGNIFICANCE_CONFLICT"

BINARY_SIG_CLASSES = frozenset({SIG_PLP, SIG_BLB})

# --- exclusion reasons (must stay inside clinvar.exclusion_reason_enum) ------
REASON_VUS = "vus"
REASON_CONFLICTING = "conflicting"
REASON_NON_BINARY = "non_binary_significance"
REASON_REF_AA_MISMATCH = "ref_aa_mismatch"
REASON_NON_CANONICAL_TX = "non_canonical_transcript"
REASON_UNPARSEABLE_HGVS = "unparseable_hgvs"
#: [NEW — Workflow v2 §1, Lead ruling 2026-08-18] a protein change that PARSED
#: CLEANLY into a non-missense consequence: synonymous, nonsense, frameshift, or
#: an in-frame indel/duplication/extension/stop-loss. Distinct from
#: REASON_UNPARSEABLE_HGVS, which is reserved for changes that could NOT be
#: interpreted at all (unmatched pattern, empty/non-protein-level HGVS, a
#: non-standard reference/alternate amino acid, or a malformed position). Never
#: aggregated together: v2 §1 requires these to be reported under the reason
#: that actually made them ineligible, not lumped under a parse-failure label.
REASON_NON_MISSENSE_VARIANT_TYPE = "non_missense_variant_type"
REASON_NO_CA = "no_ca_coordinate"
REASON_RESIDUE_CONFLICT = "residue_class_conflict"

SIG_CLASS_TO_REASON: dict[str, str] = {
    SIG_VUS: REASON_VUS,
    SIG_CONFLICTING: REASON_CONFLICTING,
    SIG_OTHER: REASON_NON_BINARY,
    SIG_RECORD_CONFLICT: REASON_CONFLICTING,
}

#: hgvs.parse_protein_hgvs() reason values that represent a CLEAN interpretation
#: of a non-missense protein consequence (Lead ruling, verbatim categories):
#: synonymous, nonsense, frameshift, in-frame indel/duplication/extension. Also
#: includes ``reference_is_termination`` (stop-loss/read-through, p.Ter600Glu) —
#: the same "interpreted, and it's not missense" character as the others, so it
#: is grouped with them rather than with the true parse failures below.
_NON_MISSENSE_VARIANT_TYPE_REASONS = frozenset({
    "synonymous", "nonsense", "reference_is_termination",
    "not_a_substitution:fs", "not_a_substitution:del", "not_a_substitution:ins",
    "not_a_substitution:dup", "not_a_substitution:ext", "not_a_substitution:delins",
})


def non_missense_exclusion_reason(change_reason: str | None) -> str:
    """Workflow v2 §1 — the top-level exclusion reason for a non-missense change.

    ``change_reason`` is ``hotspot3d.data.hgvs.ProteinChange.reason`` (or the
    ``"no_protein_hgvs"`` fallback for a record with no protein HGVS at all).
    Returns :data:`REASON_NON_MISSENSE_VARIANT_TYPE` for a change that parsed
    cleanly into a non-missense consequence, and :data:`REASON_UNPARSEABLE_HGVS`
    for everything else — including a non-standard reference/alternate amino
    acid, an unmatched HGVS pattern, an empty or non-protein-level HGVS string,
    a malformed (non-positive) position, or an uncertain/bracketed notation
    (``not_a_substitution:=/`` / ``:?`` / ``:[`` / ``:]``) that never resolved to
    a definite consequence. This function makes no other decision: it never
    changes a record's eligibility, only which of the two reasons it is filed
    under.
    """
    text = str(change_reason or "")
    if text in _NON_MISSENSE_VARIANT_TYPE_REASONS:
        return REASON_NON_MISSENSE_VARIANT_TYPE
    return REASON_UNPARSEABLE_HGVS

# ClinVar wordings that denote uncertainty or aggregate conflict. Both the
# pre-2024 "interpretations" and the current "classifications" wordings appear.
_VUS_PATTERNS = (
    "uncertain significance",
    "uncertain risk allele",
    "uncertain significance/uncertain risk allele",
)
_CONFLICTING_PATTERNS = (
    "conflicting interpretations of pathogenicity",
    "conflicting classifications of pathogenicity",
    "conflicting interpretations",
    "conflicting classifications",
    "conflicting data from submitters",
)

# ClinicalSignificance is a compound field: a primary assertion plus optional
# modifiers, separated by ',' or ';'. '/' is NOT a separator — it is part of the
# aggregate labels "Pathogenic/Likely pathogenic" and "Benign/Likely benign".
_COMPONENT_SPLIT_RE = re.compile(r"[;,]")


def _norm(text: str | None) -> str:
    return " ".join(str(text or "").strip().lower().split())


def split_significance(text: str | None) -> list[str]:
    """Split a ClinicalSignificance string into its asserted components."""
    return [part.strip() for part in _COMPONENT_SPLIT_RE.split(str(text or "")) if part.strip()]


@dataclass(frozen=True)
class SignificanceCall:
    """The normalized reading of one ClinicalSignificance string."""

    raw: str
    significance_class: str
    exact_match: bool                       # whole field equalled a frozen label
    matched_label: str | None = None
    modifiers: tuple[str, ...] = ()         # components outside the frozen labels
    exclusion_reason: str | None = None

    @property
    def eligible_primary(self) -> bool:
        return self.significance_class in BINARY_SIG_CLASSES


# --- compound-label reading policy (config: clinvar.compound_significance_policy)
# Real ClinVar aggregate labels carry orthogonal modifiers ("Pathogenic, risk
# factor"). The frozen rubric lists only the bare labels, so the compound case is
# resolved by an explicit config ruling rather than here.
POLICY_COMPONENT_WISE = "component_wise"
POLICY_WHOLE_FIELD = "whole_field"
COMPOUND_POLICIES = (POLICY_COMPONENT_WISE, POLICY_WHOLE_FIELD)

_POLICY_NOTES = {
    POLICY_COMPONENT_WISE: (
        "A ClinicalSignificance field is read component-wise: it is split on ',' and "
        "';' and a component equal to a frozen label assigns the class, so "
        "'Pathogenic, risk factor' enters P/LP. '/' is never a separator — it is part "
        "of the aggregate labels. Records matched through a component rather than the "
        "whole field carry significance_exact_match=FALSE and keep their extra "
        "components in significance_modifiers, so the whole_field reading is fully "
        "recoverable from variants_missense_all.tsv without re-running."
    ),
    POLICY_WHOLE_FIELD: (
        "Only an exact bare-label match on the whole ClinicalSignificance field enters "
        "the binary classes. Every compound label becomes non_binary_significance and "
        "is preserved in strata 1 and 3 with its components in significance_modifiers."
    ),
}


@dataclass
class SignificancePolicy:
    """Frozen label vocabulary and compound-label reading, from ``pipeline.yaml``."""

    pathogenic_labels: tuple[str, ...]
    benign_labels: tuple[str, ...]
    compound_policy: str = POLICY_COMPONENT_WISE
    _plp: dict[str, str] = field(default_factory=dict, init=False, repr=False)
    _blb: dict[str, str] = field(default_factory=dict, init=False, repr=False)

    def __post_init__(self) -> None:
        if self.compound_policy not in COMPOUND_POLICIES:
            raise ValueError(
                f"clinvar.compound_significance_policy must be one of "
                f"{list(COMPOUND_POLICIES)}, got {self.compound_policy!r}. A cohort "
                f"rule is never guessed."
            )
        self._plp = {_norm(label): label for label in self.pathogenic_labels}
        self._blb = {_norm(label): label for label in self.benign_labels}

    @classmethod
    def from_config(cls, cfg) -> "SignificancePolicy":
        return cls(
            pathogenic_labels=tuple(cfg.get("clinvar.pathogenic_labels")),
            benign_labels=tuple(cfg.get("clinvar.benign_labels")),
            compound_policy=cfg.get("clinvar.compound_significance_policy"),
        )

    @property
    def component_matching(self) -> bool:
        return self.compound_policy == POLICY_COMPONENT_WISE

    def as_dict(self) -> dict:
        return {
            "pathogenic_labels": list(self.pathogenic_labels),
            "benign_labels": list(self.benign_labels),
            "compound_significance_policy": self.compound_policy,
            "component_matching": self.component_matching,
            "component_separators": [",", ";"],
            "note": _POLICY_NOTES[self.compound_policy],
        }

    # -- the single classification entry point -------------------------------
    def classify(self, significance: str | None) -> SignificanceCall:
        raw = str(significance or "").strip() or "NA"
        whole = _norm(raw)
        components = split_significance(raw)
        normalized = [_norm(c) for c in components] if self.component_matching else [whole]

        plp_hits = [self._plp[n] for n in normalized if n in self._plp]
        blb_hits = [self._blb[n] for n in normalized if n in self._blb]
        modifiers = tuple(c for c, n in zip(components, [_norm(c) for c in components])
                          if n not in self._plp and n not in self._blb)

        if plp_hits and blb_hits:
            return SignificanceCall(
                raw=raw, significance_class=SIG_RECORD_CONFLICT, exact_match=False,
                modifiers=tuple(components),
                exclusion_reason=SIG_CLASS_TO_REASON[SIG_RECORD_CONFLICT],
            )
        if plp_hits:
            return SignificanceCall(raw=raw, significance_class=SIG_PLP,
                                    exact_match=whole in self._plp,
                                    matched_label=plp_hits[0], modifiers=modifiers)
        if blb_hits:
            return SignificanceCall(raw=raw, significance_class=SIG_BLB,
                                    exact_match=whole in self._blb,
                                    matched_label=blb_hits[0], modifiers=modifiers)

        if any(p in whole for p in _CONFLICTING_PATTERNS):
            sig = SIG_CONFLICTING
        elif any(n in _VUS_PATTERNS for n in normalized) or "uncertain" in whole:
            sig = SIG_VUS
        else:
            sig = SIG_OTHER
        return SignificanceCall(raw=raw, significance_class=sig, exact_match=False,
                                modifiers=tuple(components),
                                exclusion_reason=SIG_CLASS_TO_REASON[sig])


def compound_policy_cost(records: Sequence, *, pathogenic_labels: Sequence[str],
                         benign_labels: Sequence[str],
                         configured_policy: str) -> dict:
    """Workflow v2 §1 — the counterfactual cost of the frozen compound-significance
    policy, for THIS run's actual ClinVar payload.

    ``configured_policy`` (``clinvar.compound_significance_policy``, currently
    ``whole_field`` by Lead ruling) is never touched here and remains the only
    reading that assigns a class anywhere in Stage A. This function classifies
    every retrieved record's raw ``ClinicalSignificance`` text under BOTH policies
    and reports only aggregate counts of where they disagree, so the reader can
    judge whether the frozen ruling materially reduced either cohort — never a
    second inclusion path.
    """
    other_policy = (POLICY_COMPONENT_WISE if configured_policy == POLICY_WHOLE_FIELD
                    else POLICY_WHOLE_FIELD)
    configured = SignificancePolicy(pathogenic_labels=tuple(pathogenic_labels),
                                    benign_labels=tuple(benign_labels),
                                    compound_policy=configured_policy)
    counterfactual = SignificancePolicy(pathogenic_labels=tuple(pathogenic_labels),
                                        benign_labels=tuple(benign_labels),
                                        compound_policy=other_policy)

    n_excluded_by_configured_but_kept_by_counterfactual = 0
    n_kept_by_configured_but_excluded_by_counterfactual = 0
    examples: list[str] = []
    for rec in records:
        text = rec.clinical_significance
        a = configured.classify(text)
        b = counterfactual.classify(text)
        if (not a.eligible_primary) and b.eligible_primary:
            n_excluded_by_configured_but_kept_by_counterfactual += 1
            if text not in examples and len(examples) < 20:
                examples.append(text)
        elif a.eligible_primary and not b.eligible_primary:
            n_kept_by_configured_but_excluded_by_counterfactual += 1

    return {
        "configured_policy": configured_policy,
        "counterfactual_policy": other_policy,
        "n_records_evaluated": len(records),
        "n_excluded_by_configured_but_kept_by_counterfactual":
            n_excluded_by_configured_but_kept_by_counterfactual,
        "n_kept_by_configured_but_excluded_by_counterfactual":
            n_kept_by_configured_but_excluded_by_counterfactual,
        "note": (
            f"Counterfactual comparator only — {other_policy!r} is never applied "
            f"anywhere else in Stage A and assigns no class. "
            f"{n_excluded_by_configured_but_kept_by_counterfactual} record(s) that "
            f"{configured_policy!r} excludes from the binary cohort would have been "
            f"read into a P/LP or B/LB class under {other_policy!r}: this is the "
            f"reportable cost of the frozen ruling "
            f"(clinvar.compound_significance_policy)."
        ),
        "example_clinical_significance_values_affected": sorted(set(examples)),
    }


def residue_class(n_plp: int, n_blb: int) -> str:
    """F3 — a residue carrying both kinds of evidence is CONFLICT, never assigned.

    The conflict branch is deliberately unresolvable: there is no majority rule,
    no star-weighted vote and no "strongest evidence wins". Such a residue leaves
    the primary comparison with all of its source records preserved.
    """
    table = {
        (True, True): CLASS_CONFLICT,
        (True, False): CLASS_PLP,
        (False, True): CLASS_BLB,
    }
    key = (int(n_plp) > 0, int(n_blb) > 0)
    if key not in table:
        raise ValueError("residue_class requires at least one binary record")
    return table[key]
