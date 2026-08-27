"""P5 regeneration: an omitted per-radius grid must come back byte-identical."""
from __future__ import annotations

import numpy as np
import pytest

from synthetic_footprint import synthetic_run, two_cluster_centers

from hotspot3d.footprint import regenerate_grid
from hotspot3d.footprint.grid import distance_field, occupancy_at
from hotspot3d.footprint.stage import run_phase_c
from hotspot3d.utils.errors import BlockedError

pytestmark = pytest.mark.unit

CONFIG = regenerate_grid.DEFAULT_CONFIG


def _paths(ctx):
    return (ctx.full_results / "06_FINAL_HOTSPOTS" / "significant_hotspot_centers.tsv",
            ctx.full_results / "03_STRUCTURE_QC" / "residue_coordinates.tsv")


def test_regenerated_grid_matches_the_one_built_inside_the_sweep(phase_c_run):
    """The CLI reuses build_sweep_grid; a second builder could drift silently."""
    ctx, _, result = phase_c_run
    centers_tsv, coords_tsv = _paths(ctx)
    solution = result.solution

    for rho in (solution.domain.grid[0], result.r_fp, solution.domain.grid[-1]):
        spec, occupancy, center_ids = regenerate_grid.regenerate(
            centers_path=centers_tsv, coords_path=coords_tsv, radius=rho,
            config_path=CONFIG)

        assert spec == solution.spec, "same origin, shape, voxel edge and fallback"
        assert center_ids == result.inputs.centers.ids

        expected = occupancy_at(solution.sweep.dist, rho)
        assert np.array_equal(occupancy, expected)
        assert np.packbits(occupancy.ravel()).tobytes() == \
            np.packbits(expected.ravel()).tobytes()


def test_regenerated_final_grid_matches_the_stored_npz(phase_c_run):
    """It also reproduces the one grid that IS retained — the final radius."""
    ctx, _, result = phase_c_run
    centers_tsv, coords_tsv = _paths(ctx)

    stored = np.load(ctx.full_results / "08_FINAL_FOOTPRINT" / "footprint_occupancy.npz")
    _, occupancy, _ = regenerate_grid.regenerate(
        centers_path=centers_tsv, coords_path=coords_tsv, radius=result.r_fp,
        config_path=CONFIG)

    assert list(stored["shape"]) == list(occupancy.shape)
    assert np.allclose(stored["origin"], result.solution.spec.origin)
    assert stored["occupancy_packed"].tobytes() == \
        np.packbits(occupancy.ravel()).tobytes()


def test_the_cli_writes_a_readable_npz(tmp_path, phase_c_run):
    ctx, _, result = phase_c_run
    centers_tsv, coords_tsv = _paths(ctx)
    out = tmp_path / "grid_regen.npz"

    exit_code = regenerate_grid.main([
        "--centers", str(centers_tsv), "--coords", str(coords_tsv),
        "--radius", str(result.r_fp), "--voxel-h", str(result.solution.spec.h),
        "--config", str(CONFIG), "--out", str(out),
    ])
    assert exit_code == 0

    payload = np.load(out)
    occupancy = np.unpackbits(payload["occupancy_packed"])[
        :int(np.prod(payload["shape"]))].reshape(tuple(payload["shape"])).astype(bool)
    expected = occupancy_at(result.solution.sweep.dist, result.r_fp)
    assert np.array_equal(occupancy, expected)
    assert float(payload["voxel_h_A"][0]) == result.solution.spec.h
    assert float(payload["r_fp_A"][0]) == pytest.approx(result.r_fp)
    assert list(payload["center_residue_indices"]) == list(result.inputs.centers.ids)


def test_regeneration_is_repeatable(phase_c_run):
    ctx, _, result = phase_c_run
    centers_tsv, coords_tsv = _paths(ctx)
    runs = [regenerate_grid.regenerate(centers_path=centers_tsv,
                                       coords_path=coords_tsv, radius=result.r_fp,
                                       config_path=CONFIG)[1] for _ in range(2)]
    assert np.array_equal(runs[0], runs[1])


def test_a_radius_outside_the_domain_is_refused(phase_c_run):
    """The domain is never widened to accommodate a regeneration request."""
    ctx, _, result = phase_c_run
    centers_tsv, coords_tsv = _paths(ctx)
    outside = result.solution.domain.rho_max + 5.0

    with pytest.raises(BlockedError) as excinfo:
        regenerate_grid.regenerate(centers_path=centers_tsv, coords_path=coords_tsv,
                                   radius=outside, config_path=CONFIG)
    assert "outside the derived II.9 domain" in str(excinfo.value)


def test_an_explicit_universe_file_gives_the_same_grid(phase_c_run):
    """U_struct defaults to every usable-CA row, which is what U_struct means."""
    ctx, _, result = phase_c_run
    centers_tsv, coords_tsv = _paths(ctx)
    universe = ctx.full_results / "03_STRUCTURE_QC" / "positional_universe.tsv"

    implicit = regenerate_grid.regenerate(
        centers_path=centers_tsv, coords_path=coords_tsv, radius=result.r_fp,
        config_path=CONFIG)
    explicit = regenerate_grid.regenerate(
        centers_path=centers_tsv, coords_path=coords_tsv, radius=result.r_fp,
        universe_path=universe, config_path=CONFIG)

    assert implicit[0] == explicit[0]
    assert np.array_equal(implicit[1], explicit[1])


def test_regeneration_does_not_reimplement_the_grid_builder():
    """The audit requirement: exactly one grid-construction function exists."""
    source = (regenerate_grid.__file__)
    text = open(source, encoding="utf-8").read()
    assert "from .sweep import build_sweep_grid" in text
    assert "build_sweep_grid(" in text
    # no local bounding-box arithmetic that could drift from the sweep's
    for token in ("min(axis=0)", "max(axis=0)", "np.ceil(", "GridSpec("):
        assert token not in text, f"regenerate_grid.py builds its own grid: {token}"
