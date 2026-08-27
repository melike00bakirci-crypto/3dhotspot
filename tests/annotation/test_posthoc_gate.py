"""The post hoc gate: boundaries, near-misses, and recorded numbers on failure."""
from __future__ import annotations

import numpy as np
import pytest

from hotspot3d.annotation.mechanisms import curate
from hotspot3d.annotation.posthoc import (
    benjamini_hochberg,
    evaluate_gate,
    gate_failure_record,
    run_posthoc_spatial,
)
from hotspot3d.annotation.rubric import Evidence

from tests.fixtures.synthetic_annotation import coords_for, make_mechanism_records

pytestmark = pytest.mark.unit


def curated_with(rubric, *, n_gof: int, n_lof: int,
                 evidence: Evidence = Evidence.MODERATE):
    records = make_mechanism_records(
        n_gof=n_gof, n_lof=n_lof, evidence=evidence, include_special=False
    )
    return curate(records, rubric)


# --- the variant-count boundary ---------------------------------------------

def test_gate_passes_at_exactly_ten_in_two_categories(rubric, gate):
    ev = evaluate_gate(curated_with(rubric, n_gof=10, n_lof=10), gate)
    assert ev.passed
    assert ev.qualifying_categories == ["GOF", "LOF"]
    assert ev.counts["GOF"] == 10 and ev.counts["LOF"] == 10


def test_nine_is_a_near_miss_and_a_near_miss_is_not_a_pass(rubric, gate):
    ev = evaluate_gate(curated_with(rubric, n_gof=9, n_lof=10), gate)
    assert not ev.passed
    assert ev.qualifying_categories == ["LOF"]
    assert ev.counts["GOF"] == 9, "the near-miss count is recorded, not rounded up"


def test_both_categories_at_nine_fails(rubric, gate):
    ev = evaluate_gate(curated_with(rubric, n_gof=9, n_lof=9), gate)
    assert not ev.passed
    assert ev.n_qualifying_categories == 0


# --- the category-count boundary --------------------------------------------

def test_one_qualifying_category_fails(rubric, gate):
    ev = evaluate_gate(curated_with(rubric, n_gof=10, n_lof=0), gate)
    assert not ev.passed
    assert ev.n_qualifying_categories == 1
    assert ev.qualifying_categories == ["GOF"]


def test_two_qualifying_categories_pass(rubric, gate):
    ev = evaluate_gate(curated_with(rubric, n_gof=10, n_lof=10), gate)
    assert ev.n_qualifying_categories == 2
    assert ev.passed


# --- the evidence bar -------------------------------------------------------

def test_weak_evidence_does_not_count_toward_the_gate(rubric, gate):
    """WEAK-only variants resolve to Unclear, which is not a mechanism category."""
    ev = evaluate_gate(
        curated_with(rubric, n_gof=20, n_lof=20, evidence=Evidence.WEAK), gate
    )
    assert not ev.passed
    assert ev.counts.get("GOF", 0) == 0 and ev.counts.get("LOF", 0) == 0


def test_no_evidence_does_not_count_toward_the_gate(rubric, gate):
    ev = evaluate_gate(
        curated_with(rubric, n_gof=20, n_lof=20, evidence=Evidence.NONE), gate
    )
    assert not ev.passed
    assert ev.n_eligible_variants == 0


def test_gate_parameters_match_the_frozen_config(gate):
    assert gate.min_variants_per_category == 10
    assert gate.min_evidence is Evidence.MODERATE
    assert gate.min_categories == 2
    assert gate.separate_correction is True
    assert gate.feeds_upstream is False


# --- failure records numbers, not just a verdict ----------------------------

def test_gate_failure_records_the_numbers(rubric, gate):
    ev = evaluate_gate(curated_with(rubric, n_gof=9, n_lof=3), gate)
    record = gate_failure_record(ev, gate)

    numbers = record["gate"]
    assert numbers["counts_by_category_at_or_above_min_evidence"]["GOF"] == 9
    assert numbers["counts_by_category_at_or_above_min_evidence"]["LOF"] == 3
    assert numbers["threshold_min_variants_per_category"] == 10
    assert numbers["threshold_min_categories"] == 2
    assert numbers["threshold_min_evidence"] == "MODERATE"
    assert numbers["passed"] is False
    assert record["status"] == "NOT_RUN"
    assert record["outcome_type"] == "SCIENTIFIC_NEGATIVE"
    assert record["near_miss_is_not_a_pass"] is True
    assert record["feeds_upstream"] is False


def test_gate_failure_is_a_scientific_negative_not_an_omission(rubric, gate):
    ev = evaluate_gate(curated_with(rubric, n_gof=2, n_lof=1), gate)
    record = gate_failure_record(ev, gate)
    assert "INSUFFICIENT" in record["gate"]["reason"].upper()
    assert "valid, informative result" in record["negative_result"]["detail"]


def test_running_the_analysis_on_a_failed_gate_is_refused(rubric, gate):
    curated = curated_with(rubric, n_gof=9, n_lof=9)
    ev = evaluate_gate(curated, gate)
    coords = {cv.residue_index: coords_for(cv.residue_index) for cv in curated}
    with pytest.raises(ValueError, match="not advisory"):
        run_posthoc_spatial(
            curated, coords, gate, ev, rng=np.random.default_rng(0), seed=0,
            n_permutations=10, q=0.05,
        )


# --- the analysis itself ----------------------------------------------------

def test_posthoc_analysis_is_labelled_and_terminal(rubric, gate):
    curated = curated_with(rubric, n_gof=10, n_lof=10)
    ev = evaluate_gate(curated, gate)
    coords = {cv.residue_index: coords_for(cv.residue_index) for cv in curated}

    result = run_posthoc_spatial(
        curated, coords, gate, ev, rng=np.random.default_rng(20250101), seed=20250101,
        n_permutations=200, q=0.05,
    )
    assert result["analysis_type"] == "post_hoc"
    assert result["hypothesis_generating"] is True
    assert result["feeds_upstream"] is False
    assert result["correction"]["separate_from_stage_b"] is True
    assert {c["mechanism"] for c in result["categories"]} == {"GOF", "LOF"}
    for cat in result["categories"]:
        assert 0.0 < cat["p_emp"] <= 1.0
        assert 0.0 < cat["q_bh"] <= 1.0


def test_posthoc_analysis_is_deterministic_for_a_fixed_seed(rubric, gate):
    curated = curated_with(rubric, n_gof=10, n_lof=10)
    ev = evaluate_gate(curated, gate)
    coords = {cv.residue_index: coords_for(cv.residue_index) for cv in curated}

    def run():
        return run_posthoc_spatial(
            curated, coords, gate, ev, rng=np.random.default_rng(7), seed=7,
            n_permutations=100, q=0.05,
        )

    assert run()["categories"] == run()["categories"]


# --- the separate correction ------------------------------------------------

def test_bh_is_monotone_and_bounded():
    adjusted, significant = benjamini_hochberg([0.001, 0.01, 0.04, 0.9], 0.05)
    assert all(0.0 <= a <= 1.0 for a in adjusted)
    assert adjusted == sorted(adjusted), "step-up must not invert the p-value order"
    assert significant[0] and not significant[-1]


def test_bh_on_empty_family():
    assert benjamini_hochberg([], 0.05) == ([], [])


def test_bh_single_pvalue_is_unchanged():
    adjusted, significant = benjamini_hochberg([0.03], 0.05)
    assert adjusted == [0.03] and significant == [True]
