"""II.10 QC/admissibility and multi-objective selection."""
from __future__ import annotations

import numpy as np
import pytest

from synthetic_footprint import lattice, synthetic_run

from hotspot3d.footprint import export
from hotspot3d.footprint.api import (CenterSet, Universe,
                                     build_footprint_contained)
from hotspot3d.footprint.inputs import EvidenceRecord, EvidenceSet
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


def _contained(centers, universe, params, binding_distance_A):
    """A contained solution whose ``r_fp`` is pinned by one evidence residue.

    [DECISION-RFP-CONTAINMENT-0001] ``r_contain`` is the largest nearest-center
    distance over the evidence set, so a single record fixes the radius exactly. No
    ``ReadOnlyHotspotRadius`` is involved: the evidence set is the whole input, and
    building one by hand is how a test states "the binding distance is D".
    """
    evidence = EvidenceSet(
        records=(EvidenceRecord(residue_id=int(universe.ids[0]),
                                nearest_center_id=int(centers.ids[0]),
                                distance_A=float(binding_distance_A),
                                in_evidence=True),),
        n_plp_mapped=1)
    return build_footprint_contained(centers, universe, params, evidence)


# --- QC-F1 -------------------------------------------------------------------

def test_qc_f1_marks_inadmissible_but_keeps_every_metric(params):
    """An excessive-coverage candidate is RETAINED with all metrics and a reason.

    Containment evaluates one radius, so the QC-F1 mix that used to come from a sweep
    now comes from two solutions with different binding distances. What is being
    protected is unchanged: failing QC must not delete the candidate or its metrics.
    """
    centers, universe = _sublattice_run(6, 5.5)
    failing = export.build_scan_rows(_contained(centers, universe, params, 14.0))
    passing = export.build_scan_rows(_contained(centers, universe, params, 3.0))

    assert len(failing) == len(passing) == 1, "containment evaluates one radius"
    assert failing[0]["excessive_footprint_coverage"]
    assert not passing[0]["excessive_footprint_coverage"]

    row = failing[0]
    assert row["qc_status"] == "FAIL"
    assert "QC_F1_EXCESSIVE_COVERAGE" in row["qc_failure_reason"]
    assert row["coverage_fraction"] > 0.50
    # nothing is dropped: the full metric set survives rejection
    for column in ("volume_A3", "surface_area_A2", "n_components",
                   "compactness", "convexity", "connectivity", "parsimony"):
        assert row[column] is not None
    # the retired search's columns are absent, not NA: it did not run
    for retired in ("compactness_normalized", "connectivity_normalized",
                    "stability_normalized", "parsimony_normalized",
                    "pareto_member", "distance_to_ideal", "boundary_warning"):
        assert retired not in row, retired


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
    scan = read_tsv(ctx.full_results / "07_FOOTPRINT_RADIUS"
                    / "footprint_descriptors.tsv")
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

def test_no_geometric_descriptor_participates_in_the_determination(params):
    """`footprint_rule.geometric_metrics_decide: false`, asserted on behaviour.

    The retired search picked the radius whose objective vector sat nearest the utopia
    point, so compactness and its three companions decided. Under containment the
    radius follows the binding distance and nothing else: the same centers and the same
    universe — hence the same geometry, the same compactness ordering — give three
    different radii purely because the evidence set says so.
    """
    centers, universe = _sublattice_run(6, 5.5)
    for distance_A, expected in ((6.2, 6.5), (8.1, 8.5), (9.0, 9.0)):
        solution = _contained(centers, universe, params, distance_A)
        assert solution.domain.grid == (expected,), distance_A
        assert solution.containment["r_contain_A"] == expected
        assert solution.containment["r_contain_exact_A"] == distance_A
        # the determination is upstream of QC, so it holds whether or not QC passes
        assert solution.domain.rho_min == solution.domain.rho_max == expected


def test_parsimony_prevents_a_protein_sized_blob(phase_c_run):
    _, handoff, _ = phase_c_run
    assert handoff.payload["coverage"] <= 0.50


def test_selected_radius_is_on_grid_and_inside_the_domain(phase_c_run):
    _, _, result = phase_c_run
    domain = result.solution.domain
    assert result.r_fp in domain.grid
    assert domain.rho_min <= result.r_fp <= domain.rho_max


def test_the_sensitivity_band_retains_every_radius_it_evaluated(phase_c_run):
    """The retention property, relocated to the artefact that now carries the trace.

    There is no Pareto front to keep dominated members of. What a reader needs instead
    is the band: the footprint at every radius from the resolution floor up to r_hot,
    with the determined one marked and the others kept rather than pruned to make the
    answer look inevitable.
    """
    ctx, _, result = phase_c_run
    band = read_tsv(ctx.full_results / "07_FOOTPRINT_RADIUS" / "sensitivity_band.tsv")

    assert band, "the band is written on every completed run"
    selected = [r for r in band if str(r["selected"]).upper() == "TRUE"]
    assert len(selected) == 1, "exactly one row is the determined radius"
    assert float(selected[0]["r_fp_A"]) == pytest.approx(result.r_fp)

    radii = [float(r["r_fp_A"]) for r in band]
    assert radii == sorted(radii), "ascending, so the trend is readable"
    assert len(set(radii)) == len(radii), "no radius is listed twice"
    for row in band:
        for column in ("n_residues", "n_components", "n_plp", "n_blb", "coverage"):
            assert row[column] is not None, column


def test_admissibility_table_is_a_verbatim_projection(phase_c_run):
    """P1: footprint_admissibility.tsv copies the descriptors; it recomputes nothing."""
    ctx, _, _ = phase_c_run
    stage = ctx.full_results / "07_FOOTPRINT_RADIUS"
    scan = {r["r_fp_A"]: r for r in read_tsv(stage / "footprint_descriptors.tsv")}
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


def test_the_containment_record_is_written_whether_or_not_its_flags_fire(phase_c_run):
    """The boundary diagnostic's replacement: there is no domain edge to sit near.

    A radius determined by containment cannot be "at the boundary" of anything, so the
    three BW clauses have no subject. What a reader must still be able to audit is the
    determination: which residue bound it, by how much the grid rounded it up, whether
    that residue is a distance outlier, and whether the radius fell below the scale at
    which the voxel grid can resolve anything. All of it is recorded unconditionally.
    """
    import json

    ctx, _, _ = phase_c_run
    payload = json.loads((ctx.full_results / "07_FOOTPRINT_RADIUS"
                          / "footprint_domain.json").read_text())
    record = payload["containment"]
    for key in ("r_contain_exact_A", "r_contain_A", "grid_slack_A", "binding_residue",
                "binding_residue_is_outlier", "outlier_fence_A",
                "n_plp_mapped", "n_evidence_set", "n_plp_outside_evidence",
                "fraction_plp_outside_evidence", "binding_nearest_center",
                "evidence_distance_median_A", "evidence_distance_iqr_A",
                "coverage_escalated", "detection_scale_below_resolution"):
        assert key in record, key
    assert record["r_contain_A"] >= record["r_contain_exact_A"], "rounds UP, never down"
    assert 0.0 <= record["grid_slack_A"] < payload["step_fp_A"]
    # the retired search's diagnostic is gone, not merely unread
    assert not (ctx.full_results / "07_FOOTPRINT_RADIUS"
                / "footprint_domain_boundary_diagnostic.json").exists()
