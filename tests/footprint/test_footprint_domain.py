"""II.9 domain derivation: MST merge scales, the rho_max rule, and its edge cases."""
from __future__ import annotations

import numpy as np
import pytest

from hotspot3d.footprint.domain import build_grid, derive_domain
from hotspot3d.utils.errors import NegativeResult
from hotspot3d.utils.geometry import mst_edges

pytestmark = pytest.mark.unit


def _universe(extent: float = 100.0, n: int = 40) -> np.ndarray:
    """A straight chain of CA positions with a controllable diameter."""
    return np.column_stack([np.linspace(0.0, extent, n),
                            np.zeros(n), np.zeros(n)])


def test_merge_scales_are_mst_edge_weights_halved(params):
    centers = np.array([[0.0, 0, 0], [6.0, 0, 0], [14.0, 0, 0]])
    domain = derive_domain(centers, _universe(), params)

    weights = sorted(w for _, _, w in mst_edges(centers))
    assert weights == pytest.approx([6.0, 8.0])
    assert sorted(domain.merge_scales) == pytest.approx([3.0, 4.0])
    assert domain.rho_all == pytest.approx(4.0)


def test_rho_max_is_the_minimum_of_the_two_frozen_terms(params):
    """DECISION-FOOTPRINT-DOMAIN-0001: rho_max = min(0.25 * D_max, 25.0).

    The MST-derived ``1.25 * rho_all`` term no longer participates in the bound.
    """
    # Tightly packed centers. 1.25 * rho_all = 3.75 would once have bound here and
    # emptied the domain against the floor — the exact defect that blocked 6 of 8
    # BLOCKED genes in the ACMG SF v3.2 batch. It must no longer bind at all.
    centers = np.array([[0.0, 0, 0], [6.0, 0, 0]])          # rho_all = 3.0
    domain = derive_domain(centers, _universe(extent=200.0), params)
    assert domain.rho_all == pytest.approx(3.0)             # still computed...
    assert "1.25*rho_all" not in dict(domain.rho_max_candidates)   # ...never a cap
    assert domain.rho_max == pytest.approx(25.0)            # 0.25 * 200 = 50 > 25
    assert domain.binding_constraint == "hard_ceiling"

    # widely separated centers in a small protein make 0.25 * D_max bind
    centers = np.array([[0.0, 0, 0], [40.0, 0, 0]])
    domain = derive_domain(centers, _universe(extent=30.0), params)
    assert domain.rho_max == pytest.approx(7.5)             # 0.25 * 30
    assert domain.binding_constraint == "0.25*D_max"

    # a large protein hits the 25 A hard ceiling
    centers = np.array([[0.0, 0, 0], [120.0, 0, 0]])        # rho_all = 60
    domain = derive_domain(centers, _universe(extent=400.0), params)
    assert domain.rho_max == pytest.approx(25.0)
    assert domain.binding_constraint == "hard_ceiling"


def test_single_center_drops_the_post_merge_term(params):
    """|S| = 1 -> the MST is empty -> rho_max = min(0.25 D_max, 25), recorded."""
    domain = derive_domain(np.array([[10.0, 0, 0]]), _universe(extent=60.0), params)

    assert domain.mst == ()
    assert domain.merge_scales == ()
    assert domain.rho_all is None
    assert domain.single_center_special_case is True
    assert domain.rho_max == pytest.approx(15.0)            # 0.25 * 60
    assert domain.as_dict()["single_center_special_case"] is True


def test_empty_domain_is_a_negative_result_not_a_crash(params):
    """rho_min > rho_max terminates the chain honestly (P4).

    Under DECISION-FOOTPRINT-DOMAIN-0001 centre packing can no longer empty the
    domain; only a protein too small to host a 5 A footprint can. D_max < 20 A
    means 0.25 * D_max < the 5.0 A floor.
    """
    centers = np.array([[0.0, 0, 0], [6.0, 0, 0]])
    with pytest.raises(NegativeResult) as excinfo:
        derive_domain(centers, _universe(extent=16.0), params)   # 0.25 * 16 = 4.0
    assert excinfo.value.condition == "NO_ADMISSIBLE_FOOTPRINT_DOMAIN"
    assert excinfo.value.context["binding_constraint"] == "0.25*D_max"


def test_grid_is_uniform_and_never_adaptive():
    grid = build_grid(3.0, 8.42, 0.5)
    assert grid == (3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0)
    steps = {round(b - a, 9) for a, b in zip(grid[:-1], grid[1:])}
    assert steps == {0.5}, "step_fp is uniform; adaptive grids would break Phase D"


def test_grid_endpoints_are_stable_keys():
    """Radii are dict keys across the sweep, selection and Phase D."""
    grid = build_grid(3.0, 10.0, 0.5)
    assert len(set(grid)) == len(grid)
    assert grid[-1] == 10.0


def test_r_hot_is_not_an_argument_to_the_domain(params):
    """F10 is structural: no signature in this module accepts r_hot."""
    import inspect

    from hotspot3d.footprint import domain as domain_mod

    for name, obj in vars(domain_mod).items():
        if inspect.isfunction(obj):
            assert "r_hot" not in inspect.signature(obj).parameters, name
