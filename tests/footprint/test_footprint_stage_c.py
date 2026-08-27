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
    rows = read_tsv(ctx.full_results / "07_FOOTPRINT_RADIUS" / "r_fp_scan.tsv")
    assert len(rows) == len(result.solution.domain.grid)
    assert sum(1 for r in rows if r["selected"]) == 1
    for row in rows:
        assert row["qc_status"] in ("PASS", "FAIL")
        assert row["volume_A3"] is not None


def test_merge_events_carry_necks_and_a_connection_type(phase_c_run):
    ctx, _, _ = phase_c_run
    rows = read_tsv(ctx.full_results / "07_FOOTPRINT_RADIUS" / "merge_events.tsv")
    assert rows, "the fixture is built so that the two clusters merge in-domain"
    for row in rows:
        assert row["connection_type"] in ("natural", "artificial", "undetermined")
        assert row["neck_width_A"] is not None
        assert row["component_a_center_residues"]
        assert row["component_b_center_residues"]


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


def test_r_fp_below_r_hot_runs_end_to_end(tmp_path):
    """F10 permits r_fp < r_hot; the run must not notice or object."""
    ctx, upstream = synthetic_run(tmp_path, center_ids=two_cluster_centers()[:4],
                                  r_hot=12.0)
    handoff, result = run_phase_c(ctx, upstream=upstream)

    assert handoff.payload["hotspot_radius"] == 12.0          # passthrough, unmodified
    assert handoff.payload["footprint_radius"] < 12.0
    assert handoff.payload["r_fp_independent_of_r_hot"] is True

    decision = json.loads((ctx.full_results / "07_FOOTPRINT_RADIUS"
                           / "footprint_decision.json").read_text())
    assert decision["r_hot_used_in_selection"] is False
    assert decision["r_fp_bounded_by_r_hot"] is False


# --- the freeze gate ---------------------------------------------------------

def test_freeze_records_every_scientific_artifact(phase_c_run):
    ctx, _, result = phase_c_run
    payload = load_freeze(ctx)

    assert payload["footprint_radius"] == result.r_fp
    assert payload["footprint_code_version"] == result.code_version
    frozen = payload["files_sha256"]
    assert any("r_fp_scan.tsv" in k for k in frozen)
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

    target = ctx.full_results / "07_FOOTPRINT_RADIUS" / "r_fp_scan.tsv"
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

DECISION_BEARING = [
    "07_FOOTPRINT_RADIUS/r_fp_scan.tsv",
    "07_FOOTPRINT_RADIUS/footprint_admissibility.tsv",
    "07_FOOTPRINT_RADIUS/merge_events.tsv",
    "07_FOOTPRINT_RADIUS/footprint_pareto.json",
    "07_FOOTPRINT_RADIUS/footprint_objective_matrix.json",
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
