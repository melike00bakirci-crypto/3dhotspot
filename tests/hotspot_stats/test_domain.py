"""II.2 — the candidate radius domain, its clamps and all four fallback triggers."""
from __future__ import annotations

import numpy as np
import pytest

from hotspot3d.spatial.domain import (
    TRIGGER_DEFINITIONS,
    derive_domain,
    elevated_set,
    lattice,
    longest_contiguous_run,
)

GRID = np.arange(1.0, 21.0)
DEFAULTS = dict(r_floor=5.0, r_ceil=25.0, step=0.5, margin=1.0, n_min_pcf=10,
                min_run_length=2.0, min_grid_points=5)


def _pcf(elevated: list[float], grid: np.ndarray = GRID):
    """Observed g and its null 95th percentile, elevated exactly at ``elevated``."""
    g = np.full(len(grid), 0.8)
    p95 = np.full(len(grid), 1.0)
    for r in elevated:
        g[int(r) - 1] = 1.9
    return g, p95


@pytest.mark.unit
def test_elevated_set_requires_both_conditions():
    grid = np.array([1.0, 2.0, 3.0, 4.0])
    g = np.array([1.5, 0.9, 1.2, 2.0])
    p95 = np.array([1.6, 0.5, 1.0, 1.0])
    # r=1: above p95? no. r=2: > 1? no. r=3 and r=4: both conditions hold.
    assert list(elevated_set(grid, g, p95)) == [False, False, True, True]


@pytest.mark.unit
def test_longest_run_breaks_ties_by_the_smallest_radius():
    grid = np.arange(1.0, 9.0)
    mask = np.array([True, True, False, False, True, True, False, False])
    assert longest_contiguous_run(grid, mask) == (1.0, 2.0)
    assert longest_contiguous_run(grid, np.zeros(8, dtype=bool)) is None


@pytest.mark.unit
def test_lattice_stays_on_grid_and_never_exceeds_the_cap():
    assert list(lattice(5.0, 7.0, 0.5)) == [5.0, 5.5, 6.0, 6.5, 7.0]
    assert list(lattice(5.0, 6.3, 0.5)) == [5.0, 5.5, 6.0]
    assert len(lattice(9.0, 8.0, 0.5)) == 0


@pytest.mark.unit
def test_pcf_derived_branch():
    g, p95 = _pcf([6, 7, 8, 9, 10])
    d = derive_domain(GRID, g, p95, n_plp=30, d_max=80.0, **DEFAULTS)
    assert d.source == "pcf_derived" and d.fallback is False
    assert d.longest_run == (6.0, 10.0)
    assert (d.lo, d.hi) == (5.0, 11.0)                  # [6-1, 10+1], inside the clamps
    assert d.grid[0] == 5.0 and d.grid[-1] == 11.0
    assert d.trigger_id is None and d.admissible


@pytest.mark.unit
def test_clamping_to_r_floor_and_r_cap():
    g, p95 = _pcf([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    d = derive_domain(GRID, g, p95, n_plp=30, d_max=40.0, **DEFAULTS)
    assert d.r_cap == pytest.approx(10.0)               # min(25, 0.25 * 40)
    assert d.lo == 5.0                                   # clamped up to R_FLOOR
    assert d.hi == pytest.approx(10.0)                   # clamped down to r_cap
    assert d.detail["pre_clamp_interval"] == [1.0, 13.0]


@pytest.mark.unit
def test_ft1_empty_elevated_set_triggers_the_fallback():
    g, p95 = _pcf([])
    d = derive_domain(GRID, g, p95, n_plp=30, d_max=80.0, **DEFAULTS)
    assert d.trigger_id == "FT-1" and d.fallback is True
    assert d.source == "fallback_envelope"
    assert (d.lo, d.hi) == (5.0, 20.0)                   # [5, min(25, 0.25*80)]
    assert "FT-1" in TRIGGER_DEFINITIONS


@pytest.mark.unit
def test_ft2_short_run_triggers_the_fallback():
    g, p95 = _pcf([7, 8])                                # run length 1.0 A < 2.0 A
    d = derive_domain(GRID, g, p95, n_plp=30, d_max=80.0, **DEFAULTS)
    assert d.trigger_id == "FT-2" and d.fallback is True
    assert d.run_length_A == pytest.approx(1.0)


@pytest.mark.unit
def test_ft3_too_few_grid_points_triggers_the_fallback():
    """A long-enough run that the clamps squeeze below min_grid_points."""
    g, p95 = _pcf([1, 2, 3, 4])                          # run 1-4 -> [0,5] -> [5,5]
    d = derive_domain(GRID, g, p95, n_plp=30, d_max=80.0, **DEFAULTS)
    assert "FT-3" in d.triggers_fired and d.fallback is True
    assert d.trigger_id == "FT-3"


@pytest.mark.unit
def test_ft4_small_pathogenic_cohort_triggers_the_fallback():
    g, p95 = _pcf([6, 7, 8, 9, 10])
    d = derive_domain(GRID, g, p95, n_plp=9, d_max=80.0, **DEFAULTS)
    assert d.triggers_fired == ["FT-4"] and d.trigger_id == "FT-4"
    assert d.fallback is True and d.source == "fallback_envelope"


@pytest.mark.unit
def test_multiple_triggers_are_all_recorded():
    g, p95 = _pcf([])
    d = derive_domain(GRID, g, p95, n_plp=3, d_max=80.0, **DEFAULTS)
    assert d.triggers_fired == ["FT-1", "FT-4"]
    assert d.trigger_id == "FT-1"                        # first in canonical ID order
    assert d.as_json()["triggers_fired"] == ["FT-1", "FT-4"]


@pytest.mark.unit
def test_r_floor_above_r_cap_is_inadmissible():
    """A protein too small for the pre-registered floor: negative result, not a clamp."""
    g, p95 = _pcf([6, 7, 8])
    d = derive_domain(GRID, g, p95, n_plp=30, d_max=16.0, **DEFAULTS)   # r_cap = 4
    assert d.admissible is False and d.source == "inadmissible"
    assert "R_FLOOR" in d.detail["inadmissible_reason"]
    assert len(d.grid) == 0


@pytest.mark.unit
def test_the_pcf_peak_is_not_the_domain():
    """A single towering peak does not become the interval — F6/F14."""
    g = np.full(len(GRID), 0.5)
    p95 = np.full(len(GRID), 1.0)
    g[11] = 25.0                                          # a huge peak at r = 12
    d = derive_domain(GRID, g, p95, n_plp=30, d_max=80.0, **DEFAULTS)
    assert d.fallback is True                             # run length 0 -> FT-2
    assert d.as_json()["pcf_peak_used_as_r_hot"] is False
    assert not (d.lo <= 12.0 <= d.hi and len(d.grid) == 1)


@pytest.mark.unit
def test_json_records_the_full_audit_trail():
    g, p95 = _pcf([6, 7, 8, 9, 10])
    payload = derive_domain(GRID, g, p95, n_plp=30, d_max=80.0, **DEFAULTS).as_json()
    for key in ("branch_taken", "FALLBACK_RADIUS_DOMAIN", "trigger_id",
                "E_elevated_radii_A", "longest_contiguous_run", "clamped_interval",
                "r_cap", "grid_A", "automatic_domain_expansion"):
        assert key in payload
    assert payload["automatic_domain_expansion"] is False
