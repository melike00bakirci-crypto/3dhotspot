"""II.8 geometric core: the single EDT, the descriptors and their definitions."""
from __future__ import annotations

import numpy as np
import pytest
from scipy import ndimage

from hotspot3d.footprint import descriptors as desc
from hotspot3d.footprint.grid import (distance_field, make_grid, occupancy_at,
                                      seed_indices, seed_mask)

pytestmark = pytest.mark.unit


def _grid(points, pad=6.0, h=0.5):
    return make_grid(np.asarray(points, dtype=np.float64), pad=pad, h=h,
                     h_fallback=1.0, max_voxels=10 ** 8)


# --- one EDT, reused ---------------------------------------------------------

def test_single_edt_matches_naive_per_radius_recomputation():
    """The sublevel sets of ONE transform equal per-radius brute-force unions.

    This is the identity the whole sweep rests on: FP(rho) = {x : d(x, S) <= rho}.
    """
    points = np.array([[0.0, 0.0, 0.0], [4.0, 1.0, 0.0], [1.0, 5.0, 2.0]])
    spec = _grid(points, pad=5.0, h=0.5)
    dist = distance_field(spec, points)                     # computed exactly once

    idx = np.argwhere(np.ones(spec.shape, dtype=bool))
    centres = spec.centers_of(idx)
    seeds = spec.centers_of(seed_indices(spec, points))
    naive = np.linalg.norm(centres[:, None, :] - seeds[None, :, :],
                           axis=2).min(axis=1).reshape(spec.shape)

    assert np.allclose(dist, naive, atol=1e-9)
    for rho in (1.0, 2.5, 4.0, 6.5):
        assert np.array_equal(occupancy_at(dist, rho), naive <= rho)


def test_one_edt_call_serves_the_whole_sweep(monkeypatch):
    """A second transform per radius would be both wasteful and a drift source."""
    from hotspot3d.footprint import grid as grid_mod

    calls = {"n": 0}
    original = ndimage.distance_transform_edt

    def counting(*args, **kwargs):
        calls["n"] += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(grid_mod.ndimage, "distance_transform_edt", counting)
    points = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
    spec = _grid(points)
    distance_field(spec, points)
    assert calls["n"] == 1


def test_nested_sublevel_sets_give_monotone_volume_and_components():
    points = np.array([[0.0, 0, 0], [9.0, 0, 0], [0.0, 9.0, 0]])
    spec = _grid(points, pad=8.0)
    dist = distance_field(spec, points)

    volumes, counts = [], []
    for rho in np.arange(1.0, 8.0, 0.5):
        occ = occupancy_at(dist, rho)
        volumes.append(desc.volume(occ, spec.h))
        counts.append(desc.components(occ).n_components)

    assert all(b >= a for a, b in zip(volumes[:-1], volumes[1:]))
    assert all(b <= a for a, b in zip(counts[:-1], counts[1:]))


# --- descriptor definitions --------------------------------------------------

@pytest.mark.parametrize("radius", [4.0, 8.0])
def test_compactness_of_a_sphere_is_one_within_grid_tolerance(radius):
    spec = _grid(np.array([[0.0, 0.0, 0.0]]), pad=radius + 3.0, h=0.5)
    dist = distance_field(spec, np.array([[0.0, 0.0, 0.0]]))
    occ = occupancy_at(dist, radius)

    vol = desc.volume(occ, spec.h)
    area, mesh, reason = desc.surface_area(dist, radius, spec.h)
    assert mesh is not None, reason
    assert desc.compactness(vol, area) == pytest.approx(1.0, abs=0.03)


def test_compactness_penalises_an_elongated_shape():
    points = np.array([[0.0, 0, 0], [6.0, 0, 0], [12.0, 0, 0], [18.0, 0, 0]])
    spec = _grid(points, pad=6.0)
    dist = distance_field(spec, points)
    occ = occupancy_at(dist, 4.0)
    vol = desc.volume(occ, spec.h)
    area, _, _ = desc.surface_area(dist, 4.0, spec.h)
    assert desc.compactness(vol, area) < 0.85


def test_connectivity_is_the_largest_component_volume_fraction():
    points = np.array([[0.0, 0, 0], [40.0, 0, 0]])
    spec = _grid(points, pad=6.0)
    dist = distance_field(spec, points)

    two = desc.components(occupancy_at(dist, 4.0))
    assert two.n_components == 2
    assert desc.connectivity(two) == pytest.approx(0.5, abs=0.02)

    one = desc.components(occupancy_at(dist, 21.0))
    assert one.n_components == 1
    assert desc.connectivity(one) == 1.0


def test_parsimony_is_one_minus_coverage():
    from hotspot3d.footprint.sweep import RadiusState

    state = RadiusState(
        rho=5.0, volume=1.0, surface_area=1.0, n_components=1,
        n_covered_residues=30, coverage=0.3, compactness=0.5, convexity=0.5,
        connectivity=1.0, n_bridges=0, bridging_index=0.0, n_components_eroded=1,
        covered=np.zeros(100, dtype=bool), residue_labels=np.zeros(100, dtype=int),
        center_labels=np.array([1]), mesh_ok=True, mesh_failure_reason="NA",
        mesh=None)
    assert state.parsimony == pytest.approx(0.7)


# --- artificial bridging -----------------------------------------------------

def test_erosion_detects_a_deliberately_necked_shape():
    """Two lobes joined by a thin neck: erosion breaks the neck and counts a bridge."""
    occ = np.zeros((60, 40, 40), dtype=bool)
    occ[5:25, 12:28, 12:28] = True                 # lobe A
    occ[35:55, 12:28, 12:28] = True                # lobe B
    occ[25:35, 19:21, 19:21] = True                # 2-voxel-wide neck

    assert desc.components(occ).n_components == 1
    n_bridges, bi, n_eroded = desc.erosion_bridging(occ, h=1.0, depth=1.5)
    assert n_bridges == 1
    assert n_eroded == 2
    assert bi == pytest.approx(0.5)


def test_no_bridge_is_reported_for_a_solid_block():
    occ = np.zeros((40, 40, 40), dtype=bool)
    occ[8:32, 8:32, 8:32] = True
    n_bridges, bi, _ = desc.erosion_bridging(occ, h=1.0, depth=2.0)
    assert (n_bridges, bi) == (0, 0.0)


def test_bottleneck_radius_measures_the_neck_half_width():
    occ = np.zeros((60, 40, 40), dtype=bool)
    occ[5:25, 12:28, 12:28] = True
    occ[35:55, 12:28, 12:28] = True
    occ[25:35, 18:22, 18:22] = True                # 4-voxel-wide neck

    half = desc.bottleneck_radius(occ, h=1.0, a_idx=(10, 20, 20), b_idx=(50, 20, 20))
    assert 1.5 <= half <= 2.5                      # half-width of a 4-voxel neck


# --- ARI stability -----------------------------------------------------------

def test_ari_is_one_for_an_unchanged_partition_and_low_for_a_merge():
    labels_a = np.array([1, 1, 1, 2, 2, 2])
    labels_b = np.array([3, 3, 3, 4, 4, 4])        # relabelled, same partition
    merged = np.array([1, 1, 1, 1, 1, 1])
    mask = np.ones(6, dtype=bool)

    assert desc.adjusted_rand(labels_a, labels_b, mask) == pytest.approx(1.0)
    assert desc.adjusted_rand(labels_a, merged, mask) == pytest.approx(0.0)


def test_ari_with_no_shared_support_is_zero_by_convention():
    labels = np.array([1, 1, 2, 2])
    assert desc.adjusted_rand(labels, labels, np.zeros(4, dtype=bool)) == 0.0


# --- voxel resolution --------------------------------------------------------

def test_voxel_fallback_is_applied_and_recorded_above_the_threshold():
    points = np.array([[0.0, 0, 0], [50.0, 50.0, 50.0]])
    fine = make_grid(points, pad=10.0, h=0.5, h_fallback=1.0, max_voxels=10 ** 8)
    assert fine.h == 0.5 and fine.fallback_applied is False

    coarse = make_grid(points, pad=10.0, h=0.5, h_fallback=1.0, max_voxels=1000)
    assert coarse.h == 1.0 and coarse.fallback_applied is True
    assert coarse.as_dict()["voxel_resolution_fallback_applied"] is True


def test_every_center_voxel_is_a_seed_and_therefore_always_occupied():
    """QC-F5 holds by construction: a center is at distance 0 from itself."""
    points = np.array([[0.3, -1.7, 2.2], [11.0, 4.0, -3.0]])
    spec = _grid(points)
    mask = seed_mask(spec, points)
    assert mask.sum() == 2
    dist = distance_field(spec, points)
    idx = seed_indices(spec, points)
    assert np.all(dist[idx[:, 0], idx[:, 1], idx[:, 2]] == 0.0)
