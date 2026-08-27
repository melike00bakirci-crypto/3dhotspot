"""II.10 QC/admissibility and multi-objective selection."""
from __future__ import annotations

import numpy as np
import pytest

from synthetic_footprint import lattice, synthetic_run

from hotspot3d.footprint import export
from hotspot3d.footprint.api import CenterSet, Universe, build_footprint
from hotspot3d.footprint.qc import coverage_constraint_failure, evaluate
from hotspot3d.footprint.stage import run_phase_c
from hotspot3d.utils.errors import EscalationRequired
from hotspot3d.utils.io import read_tsv

pytestmark = pytest.mark.unit


def _sublattice_run(n_side, spacing, every=2):
    """Centers on a coarser sublattice, which drives coverage upward on purpose."""
    ids, coords = lattice(n_side, spacing)
    selected = [i for i in range(len(coords))
                if all((c / spacing) % every < 0.5 for c in coords[i])]
    centers = CenterSet(
        ids=tuple(ids[i] for i in selected), coords=coords[selected],
        universe_rows=tuple(selected))
    return centers, Universe(ids=tuple(ids), coords=coords)


# --- QC-F1 -------------------------------------------------------------------

def test_qc_f1_marks_inadmissible_but_keeps_every_metric(params):
    """An excessive-coverage candidate is RETAINED with all metrics and a reason."""
    centers, universe = _sublattice_run(5, 4.0)
    solution = build_footprint(centers, universe, params)

    rows = export.build_scan_rows(solution)
    failing = [r for r in rows if r["excessive_footprint_coverage"]]
    passing = [r for r in rows if not r["excessive_footprint_coverage"]]
    assert failing and passing, "fixture must produce a mix of QC-F1 outcomes"

    for row in failing:
        assert row["qc_status"] == "FAIL"
        assert "QC_F1_EXCESSIVE_COVERAGE" in row["qc_failure_reason"]
        assert row["coverage_fraction"] > 0.50
        # nothing is dropped: the full metric set survives rejection
        for column in ("volume_A3", "surface_area_A2", "n_components",
                       "compactness", "convexity", "connectivity", "parsimony"):
            assert row[column] is not None
        # normalization is defined over the admissible set only (F15)
        assert row["compactness_normalized"] is None

    assert len(rows) == len(solution.sweep.rows), "no candidate is ever removed"


def test_qc_f1_threshold_is_read_from_config_and_is_strict(params):
    from hotspot3d.footprint.sweep import RadiusState

    def state(coverage):
        return RadiusState(
            rho=5.0, volume=100.0, surface_area=100.0, n_components=1,
            n_covered_residues=1, coverage=coverage, compactness=0.5,
            convexity=0.5, connectivity=1.0, n_bridges=0, bridging_index=0.0,
            n_components_eroded=1, covered=np.ones(1, dtype=bool),
            residue_labels=np.ones(1, dtype=int), center_labels=np.array([1]),
            mesh_ok=True, mesh_failure_reason="NA", mesh=None)

    class _Domain:
        rho_min, rho_max = 3.0, 8.0

    assert params.qc_f1_max_coverage == 0.50
    assert evaluate(state(0.50), params, _Domain()).admissible is True
    assert evaluate(state(0.500001), params, _Domain()).admissible is False


def test_coverage_constraint_failure_escalates_and_never_widens(tmp_path):
    """>= 95% QC-F1 failures -> report the curve and STOP. No rule is relaxed."""
    ctx, upstream = synthetic_run(
        tmp_path, n_side=4, spacing=6.0,
        center_ids=list(range(1, 65)), n_plp=20, n_blb=20)

    with pytest.raises(EscalationRequired) as excinfo:
        run_phase_c(ctx, upstream=upstream)

    assert "ALL_FOOTPRINT_CANDIDATES_EXCESSIVE_COVERAGE" in excinfo.value.ambiguity
    assert "STOP" in excinfo.value.recommendation

    # the coverage curve across the WHOLE domain is on disk before stopping
    scan = read_tsv(ctx.full_results / "07_FOOTPRINT_RADIUS" / "r_fp_scan.tsv")
    assert scan and all(r["excessive_footprint_coverage"] for r in scan)

    warnings = read_tsv(ctx.full_results / "07_FOOTPRINT_RADIUS"
                        / "warnings_07_FOOTPRINT_RADIUS.tsv")
    codes = {w["warning_code"]: w["severity"] for w in warnings}
    assert codes["FOOTPRINT_COVERAGE_CONSTRAINT_FAILURE"] == "BLOCKING"

    # 09_ROBUSTNESS exists with an explanatory status — never a mystery folder
    status = ctx.full_results / "09_ROBUSTNESS" / "stage_status.json"
    assert status.is_file()
    assert (ctx.full_results / "09_ROBUSTNESS" / "NOT_RUN.txt").is_file()


def test_coverage_constraint_helper_uses_the_configured_fraction(params):
    class _V:
        def __init__(self, bad):
            self.excessive_coverage = bad

    assert params.all_inadmissible_escalation_fraction == 0.95

    fired, fraction = coverage_constraint_failure([_V(True)] * 19 + [_V(False)],
                                                  params)
    assert (fired, fraction) == (True, 0.95)                # exactly at threshold

    fired, fraction = coverage_constraint_failure([_V(True)] * 18 + [_V(False)] * 2,
                                                  params)
    assert (fired, fraction) == (False, 0.90)

    fired, fraction = coverage_constraint_failure([_V(True)] * 10, params)
    assert (fired, fraction) == (True, 1.0)


# --- selection ---------------------------------------------------------------

def test_selection_is_not_argmax_compactness(phase_c_run):
    """A single-metric optimum is exactly what the Pareto machinery must not pick."""
    _, _, result = phase_c_run
    solution = result.solution
    admissible = [r for r in solution.sweep.rows
                  if solution.verdicts[r.rho].admissible]
    best_compactness = max(admissible, key=lambda r: r.compactness).rho
    assert solution.r_fp is not None
    assert solution.r_fp != best_compactness, (
        "compactness rises monotonically with rho for a union of balls; selecting "
        "its argmax would mean parsimony was not doing its job")


def test_parsimony_prevents_a_protein_sized_blob(phase_c_run):
    _, handoff, _ = phase_c_run
    assert handoff.payload["coverage"] <= 0.50


def test_selected_radius_is_on_grid_and_inside_the_domain(phase_c_run):
    _, _, result = phase_c_run
    domain = result.solution.domain
    assert result.r_fp in domain.grid
    assert domain.rho_min <= result.r_fp <= domain.rho_max


def test_dominated_and_inadmissible_candidates_are_retained(phase_c_run):
    import json

    ctx, _, _ = phase_c_run
    payload = json.loads((ctx.full_results / "07_FOOTPRINT_RADIUS"
                          / "footprint_pareto_dominated.json").read_text())
    for entry in payload["dominated"]:
        assert entry["dominated_by"], "a dominated candidate keeps its dominators"
        assert entry["objectives"]


def test_admissibility_table_is_a_verbatim_projection(phase_c_run):
    """P1: footprint_admissibility.tsv copies r_fp_scan.tsv; it recomputes nothing."""
    ctx, _, _ = phase_c_run
    stage = ctx.full_results / "07_FOOTPRINT_RADIUS"
    scan = {r["r_fp_A"]: r for r in read_tsv(stage / "r_fp_scan.tsv")}
    projection = read_tsv(stage / "footprint_admissibility.tsv")

    assert len(projection) == len(scan)
    for row in projection:
        source = scan[row["r_fp_A"]]
        for column in row:
            if column == "schema_version":
                continue
            assert row[column] == source[column], column


def test_objectives_are_the_four_frozen_ones(params):
    assert params.objective_names == ("compactness", "connectivity", "stability",
                                      "parsimony")
    assert params.w_k == 1.0
    assert params.distance_metric == "L2"
    assert params.tie_chain == ("distance", "stability", "compactness",
                                "smaller_radius")


def test_boundary_diagnostic_is_recorded_whether_or_not_it_fires(phase_c_run):
    import json

    ctx, _, _ = phase_c_run
    payload = json.loads((ctx.full_results / "07_FOOTPRINT_RADIUS"
                          / "footprint_domain_boundary_diagnostic.json").read_text())
    for key in ("BW1_selected_in_band", "BW2_pareto_majority_in_band",
                "BW3_top3_same_band", "DOMAIN_BOUNDARY_WARNING", "binding_edge"):
        assert key in payload
    assert payload["automatic_expansion"] is False
