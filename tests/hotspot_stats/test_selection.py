"""II.6 — Pareto, normalization, distance-to-ideal, tie chain, boundary diagnostic."""
from __future__ import annotations

from types import SimpleNamespace

import pytest

from hotspot3d.hotspot.scan import OBJECTIVE_NAMES, RadiusResult
from hotspot3d.hotspot.selection import run_selection
from hotspot3d.utils.multiobjective import boundary_diagnostic

TIE_CHAIN = ["distance", "neighbor_stability", "perm_evidence_zg", "smaller_radius"]
SELECT = dict(tie_chain=TIE_CHAIN, w_k=1.0, metric="L2", near_tie_rel=0.05,
              near_tie_sep_steps=2, step=0.5, band_fraction=0.10, bw2_fraction=0.50,
              bw3_top_n=3, redundancy_rho=0.98)


def _result(radius, mcc, zg, fe, stab, admissible=True, reason="NA",
            diagnostics=()) -> RadiusResult:
    """One synthetic radius row.

    ``diagnostics`` carries the v2 §4 significance-conditioned diagnostics
    (QC_H1/QC_H2/QC_H4), which are recorded but must never affect admissibility;
    ``reason`` carries a genuine admissibility failure.
    """
    return RadiusResult(
        radius_A=radius, centers=[], center_rows=[],
        loo=SimpleNamespace(mcc=mcc), perm=SimpleNamespace(zg=zg),
        positional=SimpleNamespace(null_model="structure_aware_positional"),
        fold_enrichment=fe, fold_enrichment_defined=True,
        n_plp_in=0, n_labeled_in=0, coverage_fraction=0.1,
        n_singleton_centers=0, median_n_labeled=3.0,
        qc_failures=([] if admissible else [reason]),
        diagnostics_fired=list(diagnostics), diagnostic_values={},
        boundary_neutral_verified=True, bh_boundary_p=None, n_in_family=10,
        neighbor_stability=stab)


@pytest.mark.unit
def test_selection_is_not_argmax_mcc():
    """The max-MCC radius is recorded but holds no privileged status (F7/F8).

    r = 6.0 wins on LOO-MCC alone and loses on everything else; the balanced
    candidate r = 7.0 is nearer the utopia point and must be selected.
    """
    results = [
        _result(6.0, mcc=0.95, zg=0.10, fe=1.05, stab=0.05),
        _result(6.5, mcc=0.50, zg=1.50, fe=1.60, stab=0.50),
        _result(7.0, mcc=0.80, zg=3.00, fe=2.40, stab=0.90),
        _result(7.5, mcc=0.30, zg=2.00, fe=2.00, stab=0.60),
    ]
    sel = run_selection(results, [6.0, 6.5, 7.0, 7.5], **SELECT)

    assert sel.per_objective_argmax["loo_mcc"] == 6.0
    assert sel.selected == 7.0
    assert sel.selected != sel.per_objective_argmax["loo_mcc"]
    assert sel.independent_pareto_check


@pytest.mark.unit
def test_selection_is_not_the_pcf_peak_or_a_single_metric():
    """Each objective's argmax is a different radius; none of them decides alone."""
    results = [
        _result(5.0, mcc=0.9, zg=0.1, fe=1.0, stab=0.1),
        _result(5.5, mcc=0.1, zg=9.0, fe=1.0, stab=0.1),
        _result(6.0, mcc=0.1, zg=0.1, fe=9.0, stab=0.1),
        _result(6.5, mcc=0.1, zg=0.1, fe=1.0, stab=0.9),
        _result(7.0, mcc=0.7, zg=6.0, fe=6.0, stab=0.7),
    ]
    sel = run_selection(results, [5.0, 5.5, 6.0, 6.5, 7.0], **SELECT)
    argmaxes = set(sel.per_objective_argmax.values())
    assert len(argmaxes) == 4                       # four different single-metric winners
    assert sel.selected == 7.0                      # the balanced compromise wins


@pytest.mark.unit
def test_inadmissible_radii_are_excluded_from_selection_but_retained():
    results = [
        _result(5.0, 0.99, 9.0, 9.0, 0.99, admissible=False, reason="QC-H2"),
        _result(5.5, 0.40, 1.0, 1.2, 0.30),
        _result(6.0, 0.60, 2.0, 1.8, 0.50),
    ]
    sel = run_selection(results, [5.0, 5.5, 6.0], **SELECT)
    assert sel.selected == 6.0                      # the dominant radius is inadmissible
    assert sel.n_admissible == 2
    assert sel.result.rejected[5.0] == "QC-H2"
    assert 5.0 not in sel.result.normalized         # normalization is admissible-only
    assert 5.0 not in sel.result.distances


@pytest.mark.unit
def test_tie_chain_is_deterministic_and_total():
    """Identical objective vectors: distance ties, chain falls through to smaller r."""
    results = [_result(r, 0.5, 1.0, 1.5, 0.4) for r in (5.0, 5.5, 6.0)]
    first = run_selection(results, [5.0, 5.5, 6.0], **SELECT)
    second = run_selection(list(reversed(results)), [6.0, 5.5, 5.0], **SELECT)
    assert first.selected == second.selected == 5.0
    # Every objective is degenerate here, so all normalize to 1 and contribute 0.
    assert set(first.result.degenerate_objectives) == set(OBJECTIVE_NAMES)
    assert first.result.distances[5.0] == pytest.approx(0.0)


@pytest.mark.unit
def test_tie_chain_prefers_higher_stability_then_higher_zg():
    results = [
        _result(5.0, 0.5, 1.0, 1.0, 0.5),
        _result(6.0, 0.5, 2.0, 1.0, 0.5),           # same distance, higher Zg
    ]
    sel = run_selection(results, [5.0, 6.0], **SELECT)
    assert sel.selected == 6.0                      # Zg breaks it before smaller_radius


@pytest.mark.unit
def test_near_tie_is_flagged_and_non_blocking():
    # Two mutually non-dominated candidates five steps apart trade MCC against Zg and
    # land at the same distance-to-ideal.
    results = [
        _result(5.0, 1.0, 0.0, 1.0, 1.0),
        _result(5.5, 0.0, 0.0, 0.0, 0.0),
        _result(6.0, 0.0, 0.0, 0.0, 0.0),
        _result(7.5, 0.0, 1.0, 1.0, 1.0),
    ]
    sel = run_selection(results, [5.0, 5.5, 6.0, 7.5], **SELECT)
    assert sorted(sel.result.pareto_members) == [5.0, 7.5]
    assert sel.result.near_tie is True
    assert sel.selected is not None                  # non-blocking: a radius is still chosen
    assert sel.result.near_tie_detail["relative_gap"] < 0.05
    assert abs(sel.result.near_tie_detail["first"]
               - sel.result.near_tie_detail["second"]) > 2 * SELECT["step"]


@pytest.mark.unit
def test_objective_redundancy_is_detected():
    results = [_result(5.0 + 0.5 * i, 0.1 * i, 0.1 * i, 1.0 + i, 0.05 * i)
               for i in range(6)]
    sel = run_selection(results, [5.0 + 0.5 * i for i in range(6)], **SELECT)
    pairs = {(p["objective_a"], p["objective_b"]) for p in sel.redundant_pairs}
    assert ("loo_mcc", "perm_evidence_zg") in pairs


@pytest.mark.unit
def test_boundary_diagnostic_bw1_selected_in_band():
    grid = [5.0 + 0.5 * i for i in range(10)]
    distances = {r: (1.0 if r != grid[0] else 0.1) for r in grid}
    out = boundary_diagnostic(grid, grid[0], [grid[0], grid[5]], distances)
    assert out["BW1_selected_in_band"] is True
    assert out["DOMAIN_BOUNDARY_WARNING"] is True
    assert out["which_end"] == ["lower"]
    assert out["automatic_expansion"] is False


@pytest.mark.unit
def test_boundary_diagnostic_bw2_pareto_majority_in_band():
    grid = [5.0 + 0.5 * i for i in range(10)]
    distances = {r: 1.0 - 0.01 * i for i, r in enumerate(grid)}
    out = boundary_diagnostic(grid, grid[4], [grid[0], grid[9], grid[4]], distances)
    assert out["BW2_pareto_majority_in_band"] is True
    assert out["DOMAIN_BOUNDARY_WARNING"] is True


@pytest.mark.unit
def test_boundary_diagnostic_bw3_top_three_share_an_end():
    grid = [5.0 + 0.5 * i for i in range(30)]        # band width = ceil(0.1 * 30) = 3
    distances = {r: (0.1 if r in grid[-3:] else 5.0) for r in grid}
    out = boundary_diagnostic(grid, grid[10], [grid[10]], distances)
    assert out["BW3_top3_same_band"] is True
    assert out["which_end"] == ["upper"]


@pytest.mark.unit
def test_boundary_diagnostic_silent_when_the_optimum_is_interior():
    grid = [5.0 + 0.5 * i for i in range(20)]
    distances = {r: abs(r - 9.5) for r in grid}
    out = boundary_diagnostic(grid, 9.5, [9.5, 9.0, 10.0], distances)
    assert out["DOMAIN_BOUNDARY_WARNING"] is False
    assert out["which_end"] == ["none"]
