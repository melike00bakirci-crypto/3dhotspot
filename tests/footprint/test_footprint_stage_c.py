"""Stage C end to end: outputs, the freeze gate, F10 independence, determinism."""
from __future__ import annotations

import json

import pytest

from synthetic_footprint import synthetic_run, two_cluster_centers

from hotspot3d.footprint.freeze import load_freeze, verify_freeze
from hotspot3d.footprint.inputs import ReadOnlyHotspotRadius
from hotspot3d.footprint.stage import EXPECTED_07, EXPECTED_08, run_phase_c
from hotspot3d.utils.errors import BlockedError
from hotspot3d.utils.io import read_tsv

pytestmark = pytest.mark.unit


# --- canonical outputs -------------------------------------------------------

def test_every_expected_output_exists(phase_c_run):
    ctx, _, _ = phase_c_run
    for stage, expected in (("07_FOOTPRINT_RADIUS", EXPECTED_07),
                            ("08_FINAL_FOOTPRINT", EXPECTED_08)):
        for name, _description in expected:
            assert (ctx.full_results / stage / name).is_file(), f"{stage}/{name}"


def test_stage_manifest_declares_every_expected_file(phase_c_run):
    ctx, _, _ = phase_c_run
    rows = read_tsv(ctx.full_results / "07_FOOTPRINT_RADIUS" / "stage_manifest.tsv")
    declared = {r["file"]: r for r in rows}
    for name, _ in EXPECTED_07:
        assert declared[name]["status"] == "CREATED"
    assert all(r["status"] in ("CREATED", "NOT_CREATED") for r in rows)


def test_per_radius_component_membership_is_retained_for_every_radius(phase_c_run):
    ctx, _, result = phase_c_run
    directory = ctx.full_results / "07_FOOTPRINT_RADIUS" / "footprint_components"
    files = sorted(directory.glob("r_*.tsv"))
    assert len(files) == len(result.solution.domain.grid)


def test_per_radius_voxel_grids_are_not_stored_but_are_regenerable(phase_c_run):
    """P5: the grids are omitted and the exact regeneration recipe is recorded."""
    ctx, _, _ = phase_c_run
    stage = ctx.full_results / "07_FOOTPRINT_RADIUS"
    assert not list(stage.glob("**/*.npz"))
    domain = json.loads((stage / "footprint_domain.json").read_text())
    retention = domain["voxel_grid_retention"]
    assert retention["per_radius_voxel_grids_stored"] is False
    assert "occupancy_at(distance_field(spec" in retention["regeneration_recipe"]

    final = ctx.full_results / "08_FINAL_FOOTPRINT" / "footprint_occupancy.npz"
    assert final.is_file(), "the FINAL radius occupancy is retained"


def test_scan_records_the_full_decision_trace(phase_c_run):
    ctx, _, result = phase_c_run
    rows = read_tsv(ctx.full_results / "07_FOOTPRINT_RADIUS" / "footprint_descriptors.tsv")
    assert len(rows) == len(result.solution.domain.grid)
    assert sum(1 for r in rows if r["selected"]) == 1
    for row in rows:
        assert row["qc_status"] in ("PASS", "FAIL")
        assert row["volume_A3"] is not None


def test_merge_scales_are_reported_and_no_merge_event_can_exist(phase_c_run):
    """[DECISION-RFP-CONTAINMENT-0001] A merge event is a change between two radii.

    Containment evaluates ONE radius, so there is no previous state to merge against
    and merge_events.tsv is empty by construction, not by accident — asserted here so
    a future non-empty file is treated as a fault rather than a finding. What survives
    is the MST: its edge weights still say at which scale the center set would merge,
    and they are reported as geometry, never as a constraint on r_fp.
    """
    import json

    ctx, _, _ = phase_c_run
    rows = read_tsv(ctx.full_results / "07_FOOTPRINT_RADIUS" / "merge_events.tsv")
    assert rows == [], "one radius cannot produce a between-radius event"

    domain = json.loads((ctx.full_results / "07_FOOTPRINT_RADIUS"
                         / "footprint_domain.json").read_text())
    assert domain["n_grid_points"] == 1, "the domain is the determined radius alone"
    assert domain["merge_scales_A"], "the MST merge scales are still reported"
    for edge in domain["mst_edges"]:
        assert edge["weight_A"] > 0
        assert edge["merge_scale_A"] == pytest.approx(edge["weight_A"] / 2.0)
        assert "center_i" in edge and "center_j" in edge


def test_all_centers_are_contained_in_the_final_footprint(phase_c_run):
    """QC-F5 — a violation would be a computational fault, not a finding."""
    ctx, _, result = phase_c_run
    residues = {r["residue_index"] for r in read_tsv(
        ctx.full_results / "08_FINAL_FOOTPRINT" / "footprint_residues.tsv")}
    assert set(result.inputs.centers.ids) <= residues


def test_structures_and_view_scripts_are_written(phase_c_run):
    ctx, _, _ = phase_c_run
    structures = ctx.full_results / "08_FINAL_FOOTPRINT" / "structures"
    assert "ATOM" in (structures / "final_footprint.pdb").read_text()
    assert "_atom_site.Cartn_x" in (structures / "final_footprint.cif").read_text()
    assert "surface" in (structures / "view_footprint.cxc").read_text()
    assert "show surface, footprint" in (structures / "view_footprint.pml").read_text()


# --- F10: r_fp is independent of r_hot ---------------------------------------

def test_r_hot_cannot_be_used_as_a_number_in_a_footprint_context():
    r_hot = ReadOnlyHotspotRadius(12.0)
    assert r_hot.value_for_coverage_query == 12.0
    for operation in (float, lambda v: v + 1.0, lambda v: v > 3.0,
                      lambda v: v * 2):
        with pytest.raises(BlockedError):
            operation(r_hot)


def test_footprint_api_signature_never_mentions_r_hot():
    import inspect

    from hotspot3d.footprint import api

    for name in ("build_footprint", "reconstruct_at", "occupancy_on_grid"):
        parameters = inspect.signature(getattr(api, name)).parameters
        assert not any("hot" in p for p in parameters), name


def test_r_fp_is_bounded_above_by_r_hot_and_says_so(tmp_path):
    """L2, and the end of F10's claim of independence.

    The retired record asserted `r_fp_independent_of_r_hot: True` and
    `r_fp_bounded_by_r_hot: False`. Both are now false: r_hot defines the evidence
    set, so every evidence distance is <= r_hot and L2 makes the bound structural.
    The run must state that dependence rather than deny it — while still never
    letting the numeric r_hot reach the determination as anything but a boolean.
    """
    ctx, upstream = synthetic_run(tmp_path, center_ids=two_cluster_centers()[:4],
                                  r_hot=12.0)
    handoff, result = run_phase_c(ctx, upstream=upstream)

    assert handoff.payload["hotspot_radius"] == 12.0          # passthrough, unmodified
    assert handoff.payload["footprint_radius"] <= 12.0        # L2
    assert "r_fp_independent_of_r_hot" not in handoff.payload, (
        "the claim of independence is revoked, not restated")

    decision = json.loads((ctx.full_results / "07_FOOTPRINT_RADIUS"
                           / "footprint_decision.json").read_text())
    assert decision["rule"] == "DECISION-RFP-CONTAINMENT-0001"
    assert decision["selection_performed"] is False
    assert decision["r_hot_defines_the_evidence_set"] is True
    assert decision["r_fp_bounded_above_by_r_hot"] is True
    assert decision["blb_used_in_determination"] is False     # L5
    # the determination is reconstructible from the record alone
    assert decision["r_contain_A"] == pytest.approx(result.r_fp)
    assert decision["r_contain_A"] >= decision["r_contain_exact_A"]


# --- the freeze gate ---------------------------------------------------------

def test_freeze_records_every_scientific_artifact(phase_c_run):
    ctx, _, result = phase_c_run
    payload = load_freeze(ctx)

    assert payload["footprint_radius"] == result.r_fp
    assert payload["footprint_code_version"] == result.code_version
    frozen = payload["files_sha256"]
    assert any("footprint_descriptors.tsv" in k for k in frozen)
    assert any("footprint_residues.tsv" in k for k in frozen)
    assert any("footprint_occupancy.npz" in k for k in frozen)
    # files written after the gate are excluded, explicitly and auditably
    assert not any("footprint_freeze.json" in k for k in frozen)
    assert "handoff_03.json" in payload["excluded_written_after_gate"]

    verify_freeze(ctx, payload, when="test")                  # passes as written


def test_mutating_a_frozen_artifact_after_the_gate_is_blocking(tmp_path):
    ctx, upstream = synthetic_run(tmp_path, center_ids=two_cluster_centers()[:4])
    run_phase_c(ctx, upstream=upstream)
    payload = load_freeze(ctx)

    target = ctx.full_results / "07_FOOTPRINT_RADIUS" / "footprint_descriptors.tsv"
    target.write_text(target.read_text() + "# tampered\n")

    with pytest.raises(BlockedError) as excinfo:
        verify_freeze(ctx, payload, when="test")
    assert "BASELINE INTEGRITY FAILURE" in str(excinfo.value)


def test_stage_directories_must_be_empty_for_the_run(tmp_path):
    ctx, upstream = synthetic_run(tmp_path, center_ids=two_cluster_centers()[:4])
    run_phase_c(ctx, upstream=upstream)

    with pytest.raises(BlockedError) as excinfo:
        run_phase_c(ctx, upstream=upstream)
    assert "is not empty" in str(excinfo.value)


# --- determinism -------------------------------------------------------------

# [DECISION-RFP-CONTAINMENT-0001] footprint_pareto.json and
# footprint_objective_matrix.json are gone with the search. The artefacts that now
# carry the decision are the evidence set — which IS the determination — and the
# sensitivity band, so determinism is pinned on those instead.
DECISION_BEARING = [
    "07_FOOTPRINT_RADIUS/evidence_set.tsv",
    "07_FOOTPRINT_RADIUS/sensitivity_band.tsv",
    "07_FOOTPRINT_RADIUS/footprint_descriptors.tsv",
    "07_FOOTPRINT_RADIUS/footprint_admissibility.tsv",
    "07_FOOTPRINT_RADIUS/merge_events.tsv",
    "08_FINAL_FOOTPRINT/footprint_residues.tsv",
    "08_FINAL_FOOTPRINT/footprint_geometry.tsv",
    "08_FINAL_FOOTPRINT/footprint_components.tsv",
]


def test_phase_c_is_byte_identical_on_re_run(tmp_path):
    """Same inputs and derived seeds -> byte-identical decision-bearing outputs."""
    digests = []
    for name in ("a", "b"):
        ctx, upstream = synthetic_run(tmp_path / name,
                                      center_ids=two_cluster_centers()[:4])
        run_phase_c(ctx, upstream=upstream)
        digests.append({rel: (ctx.full_results / rel).read_bytes()
                        for rel in DECISION_BEARING})
    for rel in DECISION_BEARING:
        assert digests[0][rel] == digests[1][rel], rel


def test_decision_json_records_the_escalations_it_raised(phase_c_run):
    ctx, _, result = phase_c_run
    decision = json.loads((ctx.full_results / "07_FOOTPRINT_RADIUS"
                           / "footprint_decision.json").read_text())
    conditions = {e["condition"] for e in decision["escalations"]}
    assert conditions == {e["condition"] for e in result.escalations}
    for entry in decision["escalations"]:
        assert entry["recommendation"]
