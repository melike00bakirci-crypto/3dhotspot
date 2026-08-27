"""Bit-identical output pin for footprint/descriptors.py bottleneck_radius/
erosion_bridging, filed against bug brief QA-09-BOTTLENECK-PERF-0001.

Origin: profiling one real Stage D robustness iteration on KCNA2
(run 20260818T114125Z_f5c0984d_bb8ccb85) showed bottleneck_radius's manual,
per-voxel, dict-based union-find (its nested find() helper) accounting for
~36% of one iteration's wall time on that run's real occupancy grid — a
performance defect, not a correctness defect. The bug brief requires any
optimization of that union-find to reproduce CURRENT output bit-for-bit: same
erosion fraction semantics, same bridging detection, same numerical output.

This module is that reproducibility contract, pinned on a small fixed
synthetic fixture rather than the real KCNA2 grid (regression promotion
policy: smallest fixture that exercises the defect, not the originating
dataset). It does not touch, import, or exercise anything about N_CAP,
K_MAX, B, QC thresholds, or the 0.15*rho erosion-depth methodology constant
(``depth`` here is a literal test parameter passed directly to the function
under test, exactly as any caller would).

If footprint-robustness replaces the union-find with a faster equivalent
(e.g. a vectorized max-spanning-forest via scipy.sparse.csgraph), this test
must still pass unchanged. A value drifting here means the optimization
changed behavior, not just speed.
"""
from __future__ import annotations

import numpy as np
import pytest

from hotspot3d.footprint import descriptors as desc

pytestmark = pytest.mark.integration

H = 0.5          # voxel edge (A) -- matches the real run's voxel_h_A, a literal here
DEPTH = 0.75     # erosion depth (A) -- a literal test parameter, not derived from rho
A_IDX = (0, 1, 1)
B_IDX = (9, 1, 1)


def _bridging_fixture() -> np.ndarray:
    """Two 3x3x3 blobs joined by a single-voxel-wide, 4-voxel-long neck.

    58 occupied voxels total: deliberately small so this test is fast, but
    large enough to exercise erosion-based bridge detection AND the
    union-find join logic bottleneck_radius performs (unlike a trivial
    2-voxel fixture, which would trivially short-circuit the join loop).
    """
    occ = np.zeros((10, 3, 3), dtype=bool)
    occ[0:3, :, :] = True      # blob A: x=0..2, full 3x3 cross-section
    occ[7:10, :, :] = True     # blob B: x=7..9, full 3x3 cross-section
    occ[3:7, 1, 1] = True      # bridge: x=3..6, single voxel wide at (y=1, z=1)
    return occ


def test_fixture_is_one_connected_component_before_erosion():
    occ = _bridging_fixture()
    assert int(occ.sum()) == 58
    field = desc.components(occ)
    assert field.n_components == 1, (
        "fixture setup invariant: the neck must connect both blobs before "
        "erosion, or this fixture does not exercise bridging at all")


def test_erosion_bridging_reproduces_current_output_exactly():
    occ = _bridging_fixture()
    n_bridges, bridging_index, n_components_eroded = desc.erosion_bridging(
        occ, H, DEPTH)
    assert n_bridges == 1
    assert bridging_index == pytest.approx(0.5, abs=0.0)
    assert n_components_eroded == 2


def test_bottleneck_radius_reproduces_current_output_exactly():
    occ = _bridging_fixture()
    value = desc.bottleneck_radius(occ, H, A_IDX, B_IDX)
    assert value == pytest.approx(0.5, abs=0.0)


def test_bottleneck_radius_is_symmetric_in_its_two_endpoints():
    """Not part of the bug brief's evidence, but a cheap property any
    correct replacement must also keep: the widest-path bottleneck between
    a and b does not depend on which endpoint is passed first."""
    occ = _bridging_fixture()
    forward = desc.bottleneck_radius(occ, H, A_IDX, B_IDX)
    backward = desc.bottleneck_radius(occ, H, B_IDX, A_IDX)
    assert forward == pytest.approx(backward, abs=0.0)


def test_bottleneck_radius_nan_when_an_endpoint_is_unoccupied():
    occ = _bridging_fixture()
    empty_idx = (5, 0, 0)          # x=5 is bridge x-range but y=0,z=0 is background
    assert not occ[empty_idx]
    value = desc.bottleneck_radius(occ, H, empty_idx, B_IDX)
    assert np.isnan(value)
