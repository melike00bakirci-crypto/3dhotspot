"""The frozen evidence rubric, exhaustively.

The rubric decides every mechanism label, so it is tested over the whole evidence
truth table rather than at a few representative points.
"""
from __future__ import annotations

import itertools

import pytest

from hotspot3d.annotation.rubric import (
    ASSERTABLE,
    Evidence,
    Mechanism,
    MechanismRecord,
    RANK,
    resolve_mechanism,
)
from hotspot3d.utils.errors import BlockedError

pytestmark = pytest.mark.unit


def rec(mech: Mechanism, ev: Evidence, ref: str = "PMID:1", residue: int = 10,
        **kw) -> MechanismRecord:
    return MechanismRecord(
        residue_index=residue, aa_change="p.Ala10Val", mechanism=mech,
        reference=ref, experimental_system="sys", assay_type="assay",
        principal_finding="finding", evidence_strength=ev, **kw,
    )


# --- the four frozen rules --------------------------------------------------

def test_no_records_is_not_experimentally_characterized(rubric):
    res = resolve_mechanism([], rubric)
    assert res.mechanism is Mechanism.NOT_CHARACTERIZED
    assert res.evidence_strength is Evidence.NONE


def test_evidence_none_is_not_experimentally_characterized(rubric):
    """A real category — never filled by inference, even though a source named LOF."""
    res = resolve_mechanism([rec(Mechanism.LOF, Evidence.NONE)], rubric)
    assert res.mechanism is Mechanism.NOT_CHARACTERIZED
    assert res.evidence_strength is Evidence.NONE
    assert not res.is_labelled


def test_conflicting_strong_is_mixed_and_retains_both_references(rubric):
    res = resolve_mechanism(
        [rec(Mechanism.GOF, Evidence.STRONG, "PMID:A"),
         rec(Mechanism.LOF, Evidence.STRONG, "PMID:B")],
        rubric,
    )
    assert res.mechanism is Mechanism.MIXED
    assert res.evidence_strength is Evidence.STRONG
    assert set(res.references) == {"PMID:A", "PMID:B"}, "both sources must be retained"
    assert res.conflict


def test_weak_only_is_unclear(rubric):
    res = resolve_mechanism([rec(Mechanism.GOF, Evidence.WEAK)], rubric)
    assert res.mechanism is Mechanism.UNCLEAR


def test_conflicting_weak_is_unclear(rubric):
    res = resolve_mechanism(
        [rec(Mechanism.GOF, Evidence.WEAK, "PMID:A"),
         rec(Mechanism.LOF, Evidence.WEAK, "PMID:B")],
        rubric,
    )
    assert res.mechanism is Mechanism.UNCLEAR
    assert res.conflict


@pytest.mark.parametrize("mech", sorted(ASSERTABLE, key=lambda m: m.value))
@pytest.mark.parametrize("ev", [Evidence.MODERATE, Evidence.STRONG])
def test_single_source_at_or_above_moderate_is_labelled(rubric, mech, ev):
    res = resolve_mechanism([rec(mech, ev)], rubric)
    assert res.mechanism is mech
    assert res.evidence_strength is ev
    assert res.is_labelled


@pytest.mark.parametrize("mech", sorted(ASSERTABLE, key=lambda m: m.value))
def test_weak_never_produces_a_label(rubric, mech):
    """A label requires >= MODERATE. WEAK evidence can never assert a mechanism."""
    res = resolve_mechanism([rec(mech, Evidence.WEAK)], rubric)
    assert res.mechanism is Mechanism.UNCLEAR


# --- the full truth table ---------------------------------------------------

@pytest.mark.parametrize(
    "ev_a,ev_b", list(itertools.product(list(Evidence), list(Evidence)))
)
def test_two_source_truth_table_is_total_and_never_upgrades(rubric, ev_a, ev_b):
    """Every pair of grades resolves deterministically, and never above its evidence."""
    res = resolve_mechanism(
        [rec(Mechanism.GOF, ev_a, "PMID:A"), rec(Mechanism.LOF, ev_b, "PMID:B")],
        rubric,
    )
    assert res.mechanism in set(Mechanism)

    top = max(RANK[ev_a], RANK[ev_b])
    if top == RANK[Evidence.NONE]:
        assert res.mechanism is Mechanism.NOT_CHARACTERIZED
    elif top < RANK[Evidence.MODERATE]:
        assert res.mechanism is Mechanism.UNCLEAR, "WEAK evidence may not be upgraded"
    else:
        # At or above the bar the outcome is a real label or the conservative Unclear,
        # and the recorded grade never exceeds the evidence actually available.
        assert RANK[res.evidence_strength] <= top


@pytest.mark.parametrize("ev_a,ev_b", list(itertools.product(list(Evidence), repeat=2)))
def test_rubric_is_order_independent(rubric, ev_a, ev_b):
    a = rec(Mechanism.GOF, ev_a, "PMID:A")
    b = rec(Mechanism.LOF, ev_b, "PMID:B")
    forward = resolve_mechanism([a, b], rubric)
    backward = resolve_mechanism([b, a], rubric)
    assert forward.mechanism is backward.mechanism
    assert forward.evidence_strength is backward.evidence_strength
    assert set(forward.references) == set(backward.references)


def test_agreeing_sources_keep_the_label(rubric):
    res = resolve_mechanism(
        [rec(Mechanism.LOF, Evidence.STRONG, "PMID:A"),
         rec(Mechanism.LOF, Evidence.MODERATE, "PMID:B")],
        rubric,
    )
    assert res.mechanism is Mechanism.LOF
    assert res.evidence_strength is Evidence.STRONG


def test_higher_tier_decides_but_dissent_is_retained(rubric):
    """Contradicting sources are never dropped (agent §16)."""
    res = resolve_mechanism(
        [rec(Mechanism.GOF, Evidence.STRONG, "PMID:A"),
         rec(Mechanism.LOF, Evidence.MODERATE, "PMID:B")],
        rubric,
    )
    assert res.mechanism is Mechanism.GOF
    assert "PMID:B" in res.references
    assert "LOF" in res.curator_note


def test_conflicting_moderate_is_conservative_and_flagged(rubric):
    """Not covered by the frozen text; resolved conservatively and marked as a gap."""
    res = resolve_mechanism(
        [rec(Mechanism.GOF, Evidence.MODERATE, "PMID:A"),
         rec(Mechanism.LOF, Evidence.MODERATE, "PMID:B")],
        rubric,
    )
    assert res.mechanism is Mechanism.UNCLEAR
    assert res.rubric_gap, "the open rubric point must be surfaced, not hidden"
    assert "OPEN RUBRIC POINT" in res.curator_note
    assert set(res.references) == {"PMID:A", "PMID:B"}


def test_excluded_records_do_not_contribute(rubric):
    res = resolve_mechanism(
        [rec(Mechanism.GOF, Evidence.STRONG, "PMID:A", included=False)], rubric
    )
    assert res.mechanism is Mechanism.NOT_CHARACTERIZED


# --- vocabulary is frozen ---------------------------------------------------

def test_unknown_grade_is_blocking():
    with pytest.raises(BlockedError, match="unknown evidence grade"):
        rec(Mechanism.GOF, "PROBABLY_STRONG")


def test_unknown_mechanism_is_blocking():
    with pytest.raises(BlockedError, match="unknown mechanism"):
        rec("PARTIAL_LOF", Evidence.STRONG)


def test_config_matches_the_frozen_rubric(rubric):
    assert rubric.min_evidence_for_label is Evidence.MODERATE
    assert rubric.conflicting_strong_resolution is Mechanism.MIXED
    assert rubric.weak_only_resolution is Mechanism.UNCLEAR
