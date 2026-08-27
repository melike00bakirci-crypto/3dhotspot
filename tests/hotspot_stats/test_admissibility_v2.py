"""Workflow v2 §4 — significance is an objective, never an admissibility filter.

v1 struck a radius from the Pareto set when it produced no significant center
(``QC-H1``), and computed ``QC-H2``/``QC-H4`` over ``S(r)`` as well. Significance
then chose the radius and the radius chose significance. On KCNA2 that left 21 A as
effectively the only survivor. These tests pin the fix:

  * a zero-significance radius stays admissible and Pareto-eligible;
  * ``QC_H1``/``QC_H2``/``QC_H4`` are recorded as diagnostics and change nothing;
  * only genuinely undefined quantities may rule a radius inadmissible;
  * degenerate objectives are excluded from distance-to-utopia and reported;
  * ``VACUOUS_PARETO_SELECTION`` and ``ADMISSIBILITY_DOMINATES_SELECTION`` fire.
"""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from hotspot3d.hotspot.scan import (
    ADMISSIBILITY_RULES,
    OBJECTIVE_NAMES,
    SIGNIFICANCE_CONDITIONED_DIAGNOSTICS,
    RadiusResult,
    scan_one_radius,
)
from hotspot3d.hotspot.selection import run_selection

pytestmark = pytest.mark.unit

SELECT = dict(tie_chain=["distance", "neighbor_stability", "perm_evidence_zg",
                         "smaller_radius"],
              w_k=1.0, metric="L2", near_tie_rel=0.05, near_tie_sep_steps=2, step=0.5,
              band_fraction=0.10, bw2_fraction=0.50, bw3_top_n=3, redundancy_rho=0.98)


def _result(radius, mcc, zg, fe, stab, admissible=True, reason="NA",
            diagnostics=()) -> RadiusResult:
    return RadiusResult(
        radius_A=radius, centers=[], center_rows=[],
        loo=SimpleNamespace(mcc=mcc), perm=SimpleNamespace(zg=zg),
        positional=SimpleNamespace(null_model="structure_aware_positional"),
        fold_enrichment=fe, fold_enrichment_defined=True,
        n_plp_in=0, n_labeled_in=0, coverage_fraction=0.1, n_singleton_centers=0,
        median_n_labeled=3.0, qc_failures=([] if admissible else [reason]),
        diagnostics_fired=list(diagnostics), diagnostic_values={},
        boundary_neutral_verified=True, bh_boundary_p=None, n_in_family=10,
        neighbor_stability=stab)


# --- the constraint / diagnostic split ---------------------------------------

def test_the_two_rule_sets_are_disjoint_and_named():
    """The declared constraint list must not contain a significance-conditioned rule."""
    assert set(ADMISSIBILITY_RULES) == {"FE_UNDEF", "QC_H5_DOMAIN",
                                        "QC_H3_UNDEFINED_MEDIAN"}
    assert set(SIGNIFICANCE_CONDITIONED_DIAGNOSTICS) == {"QC_H1", "QC_H2", "QC_H4"}
    assert not set(ADMISSIBILITY_RULES) & set(SIGNIFICANCE_CONDITIONED_DIAGNOSTICS)


def test_config_declares_exactly_the_rules_the_code_applies():
    from hotspot3d.utils.config import load_config

    cfg = load_config("config/pipeline.yaml")
    assert cfg.get("radius_qc.significance_may_determine_admissibility") is False
    assert set(cfg.get("radius_qc.admissibility_constraints")) == set(ADMISSIBILITY_RULES)
    assert set(cfg.get("radius_qc.significance_conditioned_diagnostics")) == \
        set(SIGNIFICANCE_CONDITIONED_DIAGNOSTICS)


@pytest.mark.integration
def test_a_declared_rule_the_code_does_not_apply_blocks_the_stage(tmp_path):
    """Declaration and enforcement must not drift apart silently (v2 §4)."""
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from synthetic_hotspot import full_synthetic_run

    from hotspot3d.hotspot.stage import run_stage_b
    from hotspot3d.utils.errors import BlockedError

    ctx, handoff, _ = full_synthetic_run(
        tmp_path, B=200,
        overrides={"radius_qc": {"admissibility_constraints": ["FE_UNDEF", "QC_H1"]}})
    with pytest.raises(BlockedError, match="admissibility_constraints"):
        run_stage_b(ctx, upstream=handoff)


# --- a zero-significance radius stays in the optimization ---------------------

def _dispersed_inputs(n_side=5, spacing=6.0, seed=3):
    grid = np.arange(n_side) * spacing
    coords = np.array([[x, y, z] for x in grid for y in grid for z in grid],
                      dtype=np.float64)
    rng = np.random.default_rng(seed)
    positions = np.sort(rng.choice(len(coords), size=16, replace=False))
    y = np.array([1.0] * 8 + [0.0] * 8)
    return coords, np.arange(1, len(coords) + 1), positions, coords[positions], y


def test_a_radius_with_no_significant_center_remains_admissible():
    """v2 §4 — it simply scores zero on the enrichment objective."""
    coords, index, positions, labeled, y = _dispersed_inputs()
    out = scan_one_radius(
        4.0, coords, index, positions, labeled, y, np.random.default_rng(1),
        label_rng=np.random.default_rng(2), B=200, q=0.05, kappa=2.0,
        tie_tolerance=1e-12, sparse_cutoff=3, undefined_mcc=0.0,
        qc_h2_max_coverage=0.50, qc_h3_min_median=2, qc_h4_max_isolated=0.50,
        r_floor=3.0, r_cap=30.0, seed_context="scan|r=4")

    assert out.centers == []
    assert out.admissible is True                      # the v1 QC-H1 failure is gone
    assert out.qc_failure_reason == "NA"
    assert "QC_H1" in out.diagnostics_fired            # recorded, not binding
    assert out.diagnostic_values["QC_H1_zero_significant_centers"] is True
    assert out.objectives()["fold_enrichment"] == 0.0
    assert out.fold_enrichment_defined is False


def test_the_significance_conditioned_diagnostics_are_still_computed():
    """Demoting a rule must not delete its evidence."""
    coords, index, positions, labeled, y = _dispersed_inputs()
    out = scan_one_radius(
        4.0, coords, index, positions, labeled, y, np.random.default_rng(1),
        label_rng=np.random.default_rng(2), B=200, q=0.05, kappa=2.0,
        tie_tolerance=1e-12, sparse_cutoff=3, undefined_mcc=0.0,
        qc_h2_max_coverage=0.50, qc_h3_min_median=2, qc_h4_max_isolated=0.50,
        r_floor=3.0, r_cap=30.0, seed_context="scan|r=4")
    for key in ("QC_H1_zero_significant_centers", "QC_H2_coverage_fraction",
                "QC_H3_median_n_labeled", "QC_H4_isolated_center_fraction"):
        assert key in out.diagnostic_values


def test_an_off_grid_radius_is_still_inadmissible():
    """QC_H5 is a genuine domain assertion and survives the demotion of the others."""
    coords, index, positions, labeled, y = _dispersed_inputs()
    out = scan_one_radius(
        4.0, coords, index, positions, labeled, y, np.random.default_rng(1),
        label_rng=np.random.default_rng(2), B=100, q=0.05, kappa=2.0,
        tie_tolerance=1e-12, sparse_cutoff=3, undefined_mcc=0.0,
        qc_h2_max_coverage=0.50, qc_h3_min_median=2, qc_h4_max_isolated=0.50,
        r_floor=5.0, r_cap=30.0, seed_context="scan|r=4")     # 4.0 < R_FLOOR
    assert out.admissible is False
    assert out.qc_failure_reason == "QC_H5_DOMAIN"
    assert "outside" in out.admissibility_reason


def test_a_zero_significance_radius_can_win_the_selection():
    """The circularity is only broken if such a radius can actually be selected."""
    results = [
        _result(5.0, 0.20, 0.5, 1.0, 0.2, diagnostics=["QC_H2"]),
        _result(5.5, 0.90, 3.0, 2.5, 0.9, diagnostics=["QC_H1"]),   # found nothing
        _result(6.0, 0.30, 1.0, 1.2, 0.3),
    ]
    sel = run_selection(results, [5.0, 5.5, 6.0], **SELECT)
    assert sel.selected == 5.5
    assert sel.n_admissible == 3
    assert sel.admissibility_dominates is False


# --- degenerate objectives ----------------------------------------------------

def test_degenerate_objectives_are_excluded_from_the_distance_and_reported():
    results = [
        _result(5.0, 0.2, 1.0, 2.0, 0.1),
        _result(5.5, 0.6, 1.0, 2.0, 0.5),
        _result(6.0, 0.9, 1.0, 2.0, 0.9),
    ]
    sel = run_selection(results, [5.0, 5.5, 6.0], **SELECT)
    assert set(sel.degenerate_objectives) == {"perm_evidence_zg", "fold_enrichment"}
    assert sel.effective_objectives == ["loo_mcc", "neighbor_stability"]
    assert sel.degenerate_exclusion_is_numerically_neutral is True
    # Exclusion changes no distance: a degenerate axis normalizes to 1 and adds 0.
    assert sel.result.distances[6.0] == pytest.approx(0.0)
    assert sel.selected == 6.0


def test_excluding_degenerate_objectives_can_be_switched_off_for_audit():
    results = [_result(r, m, 1.0, 2.0, m) for r, m in ((5.0, 0.2), (5.5, 0.9))]
    kept = run_selection(results, [5.0, 5.5],
                         **{**SELECT, "exclude_degenerate_objectives": False})
    dropped = run_selection(results, [5.0, 5.5], **SELECT)
    assert set(kept.effective_objectives) == set(OBJECTIVE_NAMES)
    assert kept.selected == dropped.selected           # the numbers are identical
    assert kept.result.distances == dropped.result.distances


# --- the two selection guards -------------------------------------------------

def test_vacuous_pareto_fires_when_one_candidate_dominates_everything():
    results = [
        _result(5.0, 0.9, 3.0, 3.0, 0.9),              # dominates the rest outright
        _result(5.5, 0.2, 1.0, 1.0, 0.2),
        _result(6.0, 0.1, 0.5, 0.5, 0.1),
    ]
    sel = run_selection(results, [5.0, 5.5, 6.0], **SELECT)
    assert sel.result.pareto_members == [5.0]
    assert sel.vacuous_pareto is True
    assert "only non-dominated survivor" in sel.vacuous_reason


def test_vacuous_pareto_fires_when_every_objective_is_degenerate():
    results = [_result(r, 0.5, 1.0, 1.5, 0.4) for r in (5.0, 5.5, 6.0)]
    sel = run_selection(results, [5.0, 5.5, 6.0], **SELECT)
    assert set(sel.degenerate_objectives) == set(OBJECTIVE_NAMES)
    assert sel.effective_objectives == []
    assert sel.vacuous_pareto is True
    assert "every objective" in sel.vacuous_reason
    assert sel.degenerate_exclusion_is_numerically_neutral is True


def test_vacuous_pareto_is_silent_on_a_genuine_trade_off():
    results = [
        _result(5.0, 0.9, 0.1, 1.0, 0.1),
        _result(5.5, 0.1, 0.9, 1.0, 0.5),
        _result(6.0, 0.5, 0.5, 2.0, 0.9),
    ]
    sel = run_selection(results, [5.0, 5.5, 6.0], **SELECT)
    assert len(sel.result.pareto_members) > 1
    assert sel.vacuous_pareto is False


def test_admissibility_dominates_fires_above_the_configured_fraction():
    results = [
        _result(5.0, 0.5, 1.0, 1.0, 0.5, admissible=False, reason="QC_H5_DOMAIN"),
        _result(5.5, 0.5, 1.0, 1.0, 0.5, admissible=False, reason="QC_H5_DOMAIN"),
        _result(6.0, 0.5, 1.0, 1.0, 0.5, admissible=False, reason="FE_UNDEF"),
        _result(6.5, 0.6, 1.2, 1.4, 0.6),
    ]
    sel = run_selection(results, [5.0, 5.5, 6.0, 6.5], **SELECT)
    assert sel.inadmissible_fraction == pytest.approx(0.75)
    assert sel.admissibility_dominates is True
    assert sel.inadmissible_by_rule == {"FE_UNDEF": 1, "QC_H5_DOMAIN": 2}


def test_admissibility_dominates_is_silent_at_or_below_the_threshold():
    results = [
        _result(5.0, 0.5, 1.0, 1.0, 0.5, admissible=False, reason="QC_H5_DOMAIN"),
        _result(5.5, 0.6, 1.2, 1.4, 0.6),
    ]
    sel = run_selection(results, [5.0, 5.5], **SELECT)
    assert sel.inadmissible_fraction == pytest.approx(0.50)
    assert sel.admissibility_dominates is False        # strictly MORE than half
