"""II.1 — global statistics under the structure-aware positional null."""
from __future__ import annotations

import numpy as np
import pytest

from hotspot3d.spatial.ripley import (
    build_grid,
    cohort_statistics,
    intensity,
    observed_curves,
    pcf_peaks,
    positional_null_draws,
)
from synthetic_hotspot import planted_cohort, synthetic_protein


@pytest.mark.unit
def test_positional_null_draws_without_replacement(rng):
    draws = positional_null_draws(rng, n_universe=50, n_points=12, B=200)
    assert draws.shape == (200, 12)
    for row in draws:
        assert len(set(row.tolist())) == 12          # no residue occupied twice
    assert draws.min() >= 0 and draws.max() < 50


@pytest.mark.unit
def test_positional_null_refuses_an_oversized_cohort(rng):
    with pytest.raises(ValueError):
        positional_null_draws(rng, n_universe=5, n_points=6, B=10)


@pytest.mark.unit
def test_null_draws_come_from_u_struct_only(rng):
    """The null relabels POSITIONS of U_struct; it never invents coordinates (F4)."""
    universe = np.arange(30)
    draws = positional_null_draws(rng, 30, 8, 50)
    assert set(np.unique(draws)).issubset(set(universe.tolist()))


@pytest.mark.unit
def test_grid_is_capped_at_half_the_diameter():
    coords = synthetic_protein(sphere_radius=12.0)
    grid, d_max, r_max = build_grid(coords, 1.0)
    assert r_max == int(np.floor(d_max / 2))
    assert grid[0] == 1.0 and grid[-1] <= r_max


@pytest.mark.unit
def test_clustered_cohort_shows_elevated_k_and_g():
    coords = synthetic_protein(sphere_radius=16.0)
    positions, labels = planted_cohort(coords)
    grid, d_max, _ = build_grid(coords, 1.0)
    _, volume, _ = intensity(coords)

    plp = coords[positions[labels == 1]]
    blb = coords[positions[labels == 0]]
    k_p, l_p, g_p = observed_curves(plp, grid, volume, 1.0)
    k_b, l_b, g_b = observed_curves(blb, grid, volume, 1.0)

    # Below the 5 A lattice spacing both curves are identically empty, so compare over
    # the range where pairs actually exist.
    short = (grid >= 5.0) & (grid <= 8.0)
    assert (l_p - grid)[short].max() > (l_b - grid)[short].max()
    assert k_p[short].max() > k_b[short].max()
    assert g_p[short].max() > 1.0


@pytest.mark.unit
def test_cohort_statistics_envelopes_and_global_test(rng):
    coords = synthetic_protein(sphere_radius=14.0)
    positions, labels = planted_cohort(coords, n_cluster_plp=12, n_scatter_plp=4,
                                       n_blb=25)
    grid, _, _ = build_grid(coords, 1.0)
    _, volume, _ = intensity(coords)
    stats = cohort_statistics("PLP", coords[positions[labels == 1]], coords, grid, volume,
                              1.0, 300, rng, "global_null|PLP")

    assert np.all(stats.k_env_lo <= stats.k_env_hi)
    assert np.all(stats.g_env_lo <= stats.g_env_hi)
    assert stats.t_obs >= 0
    assert 1 / 301 <= stats.p_global <= 1.0          # add-one estimator, never zero
    assert stats.p_global < 0.05                     # the planted cluster is detectable


@pytest.mark.unit
def test_pcf_peaks_are_recorded_but_never_used_as_r_hot():
    grid = np.arange(1.0, 11.0)
    g = np.array([0.5, 1.4, 2.2, 1.1, 0.8, 1.9, 1.2, 0.7, 0.6, 0.5])
    peaks = pcf_peaks(grid, g)
    assert peaks["USED_AS_R_HOT"] is False
    assert peaks["principal_peak"]["radius_A"] == 3.0
    assert {p["radius_A"] for p in peaks["peaks"]} == {3.0, 6.0}


@pytest.mark.unit
def test_intensity_falls_back_when_the_hull_is_degenerate():
    flat = np.array([[float(i), float(j), 0.0] for i in range(4) for j in range(4)])
    lam, volume, method = intensity(flat)
    assert method == "bounding_box"
    assert volume > 0 and lam > 0
