"""Workflow v2 §4/§7 — objective/constraint separation, mirrored from Stage B.

Footprint admissibility (QC_F1/F2/F3) was already purely geometric before v2 —
r_fp operates on the already-frozen SIGNIFICANT_HOTSPOT_CENTERS set and never
re-tests significance, so the "significance is not admissibility" defect v2 §4
fixed in Stage B never existed here. What these tests pin is the DIAGNOSTIC layer
v2 §7 adds on top: degenerate objectives excluded from distance-to-utopia and
reported, ``VACUOUS_PARETO_SELECTION`` and ``ADMISSIBILITY_DOMINATES_SELECTION``.
"""
from __future__ import annotations

import json

import pytest

from synthetic_footprint import lattice

from hotspot3d.footprint.api import CenterSet, Universe, build_footprint
from hotspot3d.footprint.qc import QC_F1, QC_F2, QC_F3, QCVerdict
from hotspot3d.footprint.selection import SelectionDiagnostics, diagnose_selection
from hotspot3d.utils.multiobjective import Candidate, select

pytestmark = pytest.mark.unit

OBJECTIVES = ("compactness", "connectivity", "stability", "parsimony")
TIE_CHAIN = ("distance", "stability", "compactness", "smaller_radius")


def _candidate(rho, compactness, connectivity, stability, parsimony,
              admissible=True, reason="NA") -> Candidate:
    return Candidate(
        key=rho,
        objectives={"compactness": compactness, "connectivity": connectivity,
                   "stability": stability, "parsimony": parsimony},
        admissible=admissible, qc_failure_reason=reason)


def _select(candidates):
    return select(candidates, OBJECTIVES, TIE_CHAIN, w_k=1.0, metric="L2",
                 near_tie_rel=0.05, near_tie_sep_steps=2, step=0.5)


def _verdict(rho, admissible=True, reason="NA") -> QCVerdict:
    reasons = () if admissible else (reason,)
    return QCVerdict(rho=rho, admissible=admissible, reasons=reasons,
                     excessive_coverage=False)


# --- admissibility stays purely geometric, never significance-conditioned ----

def _sublattice_run(n_side, spacing, every=2):
    """Centers on a coarser sublattice, which drives coverage upward on purpose."""
    ids, coords = lattice(n_side, spacing)
    selected = [i for i in range(len(coords))
               if all((c / spacing) % every < 0.5 for c in coords[i])]
    centers = CenterSet(
        ids=tuple(ids[i] for i in selected), coords=coords[selected],
        universe_rows=tuple(selected))
    return centers, Universe(ids=tuple(ids), coords=coords)


def test_geometrically_inadmissible_candidates_never_strike_the_footprint(params):
    """v2 §7: several QC_F1/F2/F3-inadmissible radii coexist with a found footprint.

    SIGNIFICANT_HOTSPOT_CENTERS is untouched by admissibility (it is not a
    significance concept at all for r_fp) and a footprint is still selected —
    partial QC_F failure is never treated as 'no footprint found'.
    """
    centers, universe = _sublattice_run(5, 4.0)
    n_centers_before = len(centers)

    solution = build_footprint(centers, universe, params)

    inadmissible = [v for v in solution.verdicts.values() if not v.admissible]
    admissible = [v for v in solution.verdicts.values() if v.admissible]
    assert inadmissible and admissible, "fixture must produce a genuine mix"
    assert solution.r_fp is not None, "a footprint must still be selected"
    assert len(centers) == n_centers_before, "centers are never struck by admissibility"


# --- degenerate objectives ----------------------------------------------------

def test_degenerate_objective_is_excluded_and_neutral(params):
    candidates = [
        _candidate(5.0, 0.2, 1.0, 0.1, 0.9),
        _candidate(5.5, 0.6, 1.0, 0.5, 0.5),
        _candidate(6.0, 0.9, 1.0, 0.9, 0.1),
    ]
    result = _select(candidates)
    verdicts = {c.key: _verdict(c.key) for c in candidates}

    diag = diagnose_selection(result, verdicts, params)
    assert diag.degenerate_objectives == ["connectivity"]
    assert set(diag.effective_objectives) == {"compactness", "stability", "parsimony"}
    assert diag.degenerate_exclusion_is_numerically_neutral is True
    assert diag.vacuous_pareto is False


def test_degenerate_exclusion_selection_result_is_unchanged_by_inclusion(params):
    """Numerical-neutrality: excluding a degenerate axis changes no distance."""
    candidates = [
        _candidate(5.0, 0.2, 1.0, 0.9, 0.1),
        _candidate(5.5, 0.9, 1.0, 0.2, 0.9),
    ]
    result = _select(candidates)
    verdicts = {c.key: _verdict(c.key) for c in candidates}
    diag = diagnose_selection(result, verdicts, params)

    assert diag.degenerate_objectives == ["connectivity"]
    for key in result.distances:
        recomputed = result.distances[key]
        assert recomputed == pytest.approx(result.distances[key])
    assert diag.degenerate_exclusion_is_numerically_neutral is True


# --- VACUOUS_PARETO_SELECTION --------------------------------------------------

def test_vacuous_pareto_fires_with_a_single_pareto_member(params):
    candidates = [
        _candidate(5.0, 0.9, 0.9, 0.9, 0.9),   # dominates outright
        _candidate(5.5, 0.2, 0.2, 0.2, 0.2),
        _candidate(6.0, 0.1, 0.1, 0.1, 0.1),
    ]
    result = _select(candidates)
    verdicts = {c.key: _verdict(c.key) for c in candidates}

    diag = diagnose_selection(result, verdicts, params)
    assert result.pareto_members == [5.0]
    assert diag.vacuous_pareto is True
    assert "only non-dominated survivor" in diag.vacuous_reason


def test_vacuous_pareto_fires_when_every_objective_is_degenerate(params):
    candidates = [_candidate(r, 0.5, 0.5, 0.5, 0.5) for r in (5.0, 5.5, 6.0)]
    result = _select(candidates)
    verdicts = {c.key: _verdict(c.key) for c in candidates}

    diag = diagnose_selection(result, verdicts, params)
    assert set(diag.degenerate_objectives) == set(OBJECTIVES)
    assert diag.effective_objectives == []
    assert diag.vacuous_pareto is True
    assert "every objective" in diag.vacuous_reason
    assert diag.degenerate_exclusion_is_numerically_neutral is True


def test_vacuous_pareto_is_silent_on_a_genuine_trade_off(params):
    candidates = [
        _candidate(5.0, 0.9, 0.1, 0.5, 0.1),
        _candidate(5.5, 0.1, 0.9, 0.5, 0.5),
        _candidate(6.0, 0.5, 0.5, 0.9, 0.9),
    ]
    result = _select(candidates)
    verdicts = {c.key: _verdict(c.key) for c in candidates}

    diag = diagnose_selection(result, verdicts, params)
    assert len(result.pareto_members) > 1
    assert diag.vacuous_pareto is False


# --- ADMISSIBILITY_DOMINATES_SELECTION -----------------------------------------

def test_admissibility_dominates_fires_above_the_configured_fraction(params):
    assert params.admissibility_dominates_fraction == pytest.approx(0.50)
    candidates = [
        _candidate(5.0, 0.5, 0.5, 0.5, 0.5, admissible=False,
                  reason=f"{QC_F1}:coverage=0.9>0.5"),
        _candidate(5.5, 0.5, 0.5, 0.5, 0.5, admissible=False,
                  reason=f"{QC_F1}:coverage=0.9>0.5"),
        _candidate(6.0, 0.5, 0.5, 0.5, 0.5, admissible=False,
                  reason=f"{QC_F2}:bridging_index=0.9>0.5"),
        _candidate(6.5, 0.6, 0.6, 0.6, 0.6),
    ]
    result = _select(candidates)
    verdicts = {c.key: _verdict(c.key, admissible=c.admissible, reason=c.qc_failure_reason)
               for c in candidates}

    diag = diagnose_selection(result, verdicts, params)
    assert diag.inadmissible_fraction == pytest.approx(0.75)
    assert diag.admissibility_dominates is True
    assert diag.inadmissible_by_rule == {QC_F1: 2, QC_F2: 1}


def test_admissibility_dominates_is_silent_at_or_below_the_threshold(params):
    candidates = [
        _candidate(5.0, 0.5, 0.5, 0.5, 0.5, admissible=False,
                  reason=f"{QC_F3}:volume=0;mesh=extraction_failed"),
        _candidate(5.5, 0.6, 0.6, 0.6, 0.6),
    ]
    result = _select(candidates)
    verdicts = {c.key: _verdict(c.key, admissible=c.admissible, reason=c.qc_failure_reason)
               for c in candidates}

    diag = diagnose_selection(result, verdicts, params)
    assert diag.inadmissible_fraction == pytest.approx(0.50)
    assert diag.admissibility_dominates is False        # strictly MORE than half
    # QC_F3's own reason text has an internal ';' (volume=...;mesh=...); the rule
    # code must still be recovered exactly once, not split into a bogus fragment.
    assert diag.inadmissible_by_rule == {QC_F3: 1}


# --- wired into footprint_decision.json ----------------------------------------

def test_diagnostics_are_wired_into_footprint_decision_json(phase_c_run):
    ctx, _, result = phase_c_run
    payload = json.loads((ctx.full_results / "07_FOOTPRINT_RADIUS"
                         / "footprint_decision.json").read_text())
    for key in ("effective_objectives_in_distance",
               "degenerate_objectives_excluded_from_distance",
               "degenerate_exclusion_is_numerically_neutral",
               "VACUOUS_PARETO_SELECTION", "vacuous_pareto_reason",
               "ADMISSIBILITY_DOMINATES_SELECTION", "inadmissible_fraction",
               "inadmissible_by_rule", "significance_may_determine_admissibility"):
        assert key in payload, key
    assert payload["significance_may_determine_admissibility"] is False
    assert payload["degenerate_exclusion_is_numerically_neutral"] is True
