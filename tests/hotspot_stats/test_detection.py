"""II.7 — final detection: the three residue objects, regions and the center family."""
from __future__ import annotations

import numpy as np
import pytest

from hotspot3d.hotspot.detection import (
    ALL_CENTER_TEST_COLUMNS,
    final_detection,
    significant_subset_matches,
    verify_object_nesting,
)
from synthetic_hotspot import planted_cohort, synthetic_protein

Q = 0.05


@pytest.fixture
def detection(rng):
    coords = synthetic_protein(sphere_radius=16.0)
    positions, labels = planted_cohort(coords, n_cluster_plp=18, n_scatter_plp=4,
                                       n_blb=35)
    index = np.arange(1, len(coords) + 1, dtype=np.int64)
    return final_detection(7.0, coords, index, positions, coords[positions],
                           labels.astype(float), rng, B=1000, q=Q,
                           seed_context="final_detection"), coords, index, positions, labels


@pytest.mark.unit
def test_every_residue_of_u_struct_is_a_candidate_center(detection):
    det, coords, index, _, _ = detection
    assert len(det.center_rows) == len(coords)
    assert {r["center_residue_index"] for r in det.center_rows} == set(index.tolist())


@pytest.mark.unit
def test_uninformative_centers_leave_the_family_with_a_reason(detection):
    det, *_ = detection
    excluded = [r for r in det.center_rows if not r["in_test_family"]]
    for row in excluded:
        assert row["exclusion_reason"] in ("no_labeled_residue_in_sphere",
                                           "zero_permutation_variance")
        assert row["q_bh"] is None and row["significant"] is False
    assert det.counts["n_in_test_family"] + len(excluded) == len(det.center_rows)


@pytest.mark.unit
def test_significant_file_is_the_exact_significant_subset(detection):
    det, *_ = detection
    assert significant_subset_matches(det.center_rows, det.significant_rows)
    assert all(r["significant"] for r in det.significant_rows)
    assert set(det.significant_rows[0]) == set(ALL_CENTER_TEST_COLUMNS) - {"schema_version"}


@pytest.mark.unit
def test_three_objects_have_distinct_key_columns(detection):
    det, *_ = detection
    assert "center_residue_index" in det.significant_rows[0]
    assert "variant_residue_index" in det.classified_rows[0]
    assert "covered_residue_index" in det.covered_rows[0]
    # An accidental join on a shared key is impossible: the key names are disjoint.
    assert "center_residue_index" not in det.covered_rows[0]
    assert "covered_residue_index" not in det.classified_rows[0]


@pytest.mark.unit
def test_object_nesting_centers_subset_covered_subset_universe(detection):
    det, coords, index, positions, labels = detection
    verify_object_nesting(det)
    centers = {r["center_residue_index"] for r in det.significant_rows}
    covered = {r["covered_residue_index"] for r in det.covered_rows}
    classified = {r["variant_residue_index"] for r in det.classified_rows}
    labelled = set(index[positions].tolist())

    assert centers                                # the planted cluster is detected
    assert centers <= covered                     # a center is within r of itself
    assert classified <= covered
    assert classified <= labelled
    assert covered <= set(index.tolist())


@pytest.mark.unit
def test_centers_need_not_carry_a_clinvar_variant(detection):
    det, coords, index, positions, labels = detection
    labelled = set(index[positions].tolist())
    without = [r for r in det.significant_rows if not r["carries_clinvar_variant"]]
    assert det.n_centers_without_variant == len(without)
    assert all(r["center_residue_index"] not in labelled for r in without)
    assert det.counts["n_centers_without_variant"] == len(without)


@pytest.mark.unit
def test_regions_are_connected_components_under_twice_the_radius(detection):
    det, coords, index, *_ = detection
    radius = det.radius_A
    position_of = {int(v): i for i, v in enumerate(index)}
    for row in det.region_rows:
        members = [position_of[int(x)] for x in row["center_residue_indices"].split(";")]
        if len(members) > 1:
            d = np.linalg.norm(coords[members][:, None] - coords[members][None], axis=2)
            np.fill_diagonal(d, np.inf)
            assert d.min(axis=1).max() <= 2 * radius     # every member has a neighbour
    # Distinct regions are separated by more than 2r.
    ids = [r["hotspot_id"] for r in det.region_rows]
    assert len(ids) == len(set(ids))
    assert all(i.startswith("H") for i in ids)


@pytest.mark.unit
def test_hotspot_ids_are_consistent_across_the_three_objects(detection):
    det, *_ = detection
    region_ids = {r["hotspot_id"] for r in det.region_rows}
    for row in det.covered_rows:
        for hid in row["hotspot_ids"].split(";"):
            assert hid in region_ids
    for row in det.significant_rows:
        assert row["hotspot_id"] in region_ids


@pytest.mark.unit
def test_bh_table_matches_the_center_table(detection):
    det, *_ = detection
    rejected_in_table = {r["center_residue_index"] for r in det.bh_rows if r["reject_bh"]}
    significant = {r["center_residue_index"] for r in det.significant_rows}
    assert rejected_in_table == significant
    ranks = [r["rank"] for r in det.bh_rows]
    assert ranks == sorted(ranks)
    ps = [r["p_emp"] for r in det.bh_rows]
    assert ps == sorted(ps)
    if det.boundary_p is not None:
        assert any(r["is_boundary"] for r in det.bh_rows)


@pytest.mark.unit
def test_no_detection_on_a_cohort_without_structure(rng):
    """Labels spread uniformly at random produce no hotspot — a valid negative."""
    coords = synthetic_protein(sphere_radius=16.0)
    index = np.arange(1, len(coords) + 1, dtype=np.int64)
    picker = np.random.default_rng(2)
    positions = np.sort(picker.choice(len(coords), size=50, replace=False))
    labels = np.zeros(50)
    labels[picker.choice(50, size=20, replace=False)] = 1.0
    det = final_detection(7.0, coords, index, positions, coords[positions], labels, rng,
                          B=1000, q=Q, seed_context="final_detection")
    assert det.n_significant == 0
    assert det.significant_rows == [] and det.covered_rows == []
    assert det.region_rows == []
    verify_object_nesting(det)                    # nesting holds trivially, not by accident
