"""Residue collapse, F3 conflicts and the four-stratum conservation proof."""
from __future__ import annotations

import pytest

from hotspot3d.data.classify import (
    CLASS_BLB,
    CLASS_CONFLICT,
    CLASS_PLP,
    REASON_NO_CA,
    REASON_RESIDUE_CONFLICT,
    SignificancePolicy,
)
from hotspot3d.data.clinvar import normalize_records
from hotspot3d.data.cohort import (
    ELIGIBLE_STRATUM,
    EXCLUDED_STRATUM,
    apply_structural_exclusions,
    assert_cohort_viable,
    cohort_counts,
    collapse_to_residues,
    evaluate_records,
    prove_conservation,
    review_star_distribution,
)
from hotspot3d.utils.errors import NegativeResult

pytestmark = pytest.mark.unit


@pytest.fixture
def built(config):
    policy = SignificancePolicy.from_config(config)

    def _build(case):
        evals = evaluate_records(
            normalize_records(case.records, gene=case.gene), policy=policy,
            canonical_sequence=case.sequence, mane_transcript=case.mane_transcript)
        rows = collapse_to_residues(evals, canonical_sequence=case.sequence)
        return evals, rows

    return _build


# --- collapse ---------------------------------------------------------------

def test_duplicate_records_collapse_onto_one_residue(cases, built):
    case = cases["clustered"]
    evals, rows = built(case)
    by_index = {r.residue_index: r for r in rows}

    assert len(by_index) == len(rows), "residue indices must be unique after collapse"
    multi = [r for r in rows if r.n_records_plp + r.n_records_blb > 1]
    assert multi, "the fixture must exercise the duplicate-record collapse"
    for row in multi:
        assert len(row.clinvar_ids) == row.n_records_plp + row.n_records_blb


def test_residue_classes_match_the_fixture_ground_truth(cases, built):
    for name, case in cases.items():
        _, rows = built(case)
        plp = {r.residue_index for r in rows if r.residue_class == CLASS_PLP}
        blb = {r.residue_index for r in rows if r.residue_class == CLASS_BLB}
        assert plp == set(case.expected["plp_residues"]), name
        assert blb == set(case.expected["blb_residues"]), name


def test_every_residue_row_has_a_class_and_at_least_one_record(cases, built):
    for case in cases.values():
        _, rows = built(case)
        for row in rows:
            assert row.residue_class in (CLASS_PLP, CLASS_BLB, CLASS_CONFLICT)
            assert len(row.clinvar_ids) >= 1


# --- F3 ---------------------------------------------------------------------

def test_conflict_residues_keep_every_record_and_join_no_class(cases, built):
    case = cases["conflict"]
    _, rows = built(case)
    conflicts = [r for r in rows if r.residue_class == CLASS_CONFLICT]

    assert {r.residue_index for r in conflicts} == set(case.expected["conflict_residues"])
    for row in conflicts:
        assert row.eligible_primary is False
        assert row.exclusion_reason == REASON_RESIDUE_CONFLICT
        assert row.n_records_plp > 0 and row.n_records_blb > 0
        assert len(row.clinvar_ids) == row.n_records_plp + row.n_records_blb
        for ev in row.evals:
            assert ev.eligible_primary is False
            assert ev.exclusion_reason == REASON_RESIDUE_CONFLICT


def test_conflict_residues_are_excluded_from_the_counts(cases, built):
    case = cases["conflict"]
    _, rows = built(case)
    counts = cohort_counts(rows)
    assert counts.n_conflict == case.expected["n_conflict"]
    assert counts.N == case.expected["N"]
    assert counts.N_P + counts.N_B == counts.N


# --- conservation (QC rule 1) ----------------------------------------------

def test_stratum_conservation_holds_for_every_case(cases, built):
    for name, case in cases.items():
        evals, _ = built(case)
        proof = prove_conservation(evals)
        assert proof.conserved, f"{name}: {proof.as_dict()}"
        assert proof.n_stratum1 == len(case.records)
        assert proof.n_eligible + proof.n_excluded == proof.n_stratum1
        assert proof.overlap == () and proof.unaccounted == ()


def test_every_excluded_record_carries_a_reason_from_the_frozen_enum(cases, built, config):
    allowed = set(config.get("clinvar.exclusion_reason_enum"))
    for case in cases.values():
        evals, _ = built(case)
        excluded = [e for e in evals if not e.eligible_primary]
        assert excluded
        for ev in excluded:
            assert ev.exclusion_reason in allowed


def test_exclusion_reason_counts_match_the_fixture(cases, built):
    for name, case in cases.items():
        evals, _ = built(case)
        tally: dict[str, int] = {}
        for ev in evals:
            if not ev.eligible_primary:
                tally[ev.exclusion_reason] = tally.get(ev.exclusion_reason, 0) + 1
        assert tally == case.expected["exclusion_reason_counts"], name


def test_record_ids_stay_unique_even_when_variation_ids_repeat():
    rows = [{"VariationID": "42", "GeneSymbol": "G",
             "Name": "NM_1.1(G):c.1A>G (p.Met1Val)",
             "ClinicalSignificance": "Pathogenic", "ReviewStatus": "practice guideline"}
            for _ in range(3)]
    records = normalize_records(rows, gene="G")
    assert len({r.record_id for r in records}) == 3


# --- structural exclusion ---------------------------------------------------

def test_residues_without_a_usable_ca_are_excluded_with_a_reason(cases, built):
    case = cases["clustered"]
    evals, rows = built(case)
    victim = case.expected["plp_residues"][0]

    excluded = apply_structural_exclusions(rows, [victim])
    assert excluded == [victim]

    row = next(r for r in rows if r.residue_index == victim)
    assert row.eligible_primary is False
    assert row.exclusion_reason == REASON_NO_CA
    assert all(e.exclusion_reason == REASON_NO_CA for e in row.evals)
    assert prove_conservation(evals).conserved


# --- stratum 4 --------------------------------------------------------------

def test_star_distribution_covers_overall_class_and_stratum(cases, built):
    evals, _ = built(cases["clustered"])
    rows = review_star_distribution(evals)
    assert rows

    overall = sum(r["n_records"] for r in rows
                  if r["class"] == "ALL" and r["inclusion_stratum"] == "ALL")
    assert overall == len(evals)

    by_stratum = sum(r["n_records"] for r in rows
                     if r["class"] == "ALL"
                     and r["inclusion_stratum"] in (ELIGIBLE_STRATUM, EXCLUDED_STRATUM))
    assert by_stratum == len(evals)

    assert {r["class"] for r in rows} - {"ALL"}, "per-class rows must be present"
    assert "NA" in {str(r["star_level"]) for r in rows}, "unknown wordings appear as NA"


# --- viability --------------------------------------------------------------

def test_viable_cohort_passes(cases, built):
    _, rows = built(cases["clustered"])
    assert_cohort_viable(cohort_counts(rows)) is None


@pytest.mark.parametrize("dropped", ["PLP", "BLB"])
def test_an_empty_class_is_a_negative_result_not_a_crash(cases, built, dropped):
    case = cases["sparse"].without_class(dropped)
    _, rows = built(case)
    counts = cohort_counts(rows)
    with pytest.raises(NegativeResult) as excinfo:
        assert_cohort_viable(counts)
    assert excinfo.value.condition == "INSUFFICIENT_CLASSIFIED_RESIDUES"
    assert (counts.N_P == 0) or (counts.N_B == 0)
