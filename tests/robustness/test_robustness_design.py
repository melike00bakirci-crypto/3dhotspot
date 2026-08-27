"""II.11 subset design: K_MAX, the budget arithmetic and combinadic unranking."""
from __future__ import annotations

from itertools import combinations
from math import comb

import pytest

from hotspot3d.robustness.design import (build_design, k_max, rank_subset,
                                         unrank_subset)
from hotspot3d.utils.errors import BlockedError
from hotspot3d.utils.seeds import SeedRegistry

pytestmark = pytest.mark.unit


def _seeds(run_id="TEST_RUN"):
    return SeedRegistry(master_seed=20250101, run_id=run_id)


# --- K_MAX -------------------------------------------------------------------

@pytest.mark.parametrize("n_s,expected", [
    (2, 1), (3, 1), (4, 2), (5, 2), (6, 3), (7, 3), (10, 5), (11, 5), (100, 50),
])
def test_k_max_rule(n_s, expected):
    """K_MAX = min(n_S - 1, floor(n_S / 2)) — frozen (A19)."""
    assert k_max(n_s) == expected


def test_k_max_leaves_at_least_one_center(n_s=9):
    assert n_s - k_max(n_s) >= 1


# --- combinadic unranking ----------------------------------------------------

@pytest.mark.parametrize("n,k", [(5, 2), (6, 3), (8, 4), (10, 1)])
def test_unranking_enumerates_every_subset_exactly_once(n, k):
    produced = [unrank_subset(rank, n, k) for rank in range(comb(n, k))]
    assert produced == sorted(combinations(range(n), k))
    assert len(set(produced)) == comb(n, k), "duplicate-free by construction"


@pytest.mark.parametrize("n,k", [(5, 2), (9, 4)])
def test_ranking_is_the_exact_inverse_of_unranking(n, k):
    for rank in range(comb(n, k)):
        assert rank_subset(unrank_subset(rank, n, k), n) == rank


def test_unranking_is_exact_for_astronomically_large_spaces():
    """Python integers, not float or int64: C(200, 100) has 59 digits."""
    n, k = 200, 100
    total = comb(n, k)
    assert total > 2 ** 63
    for rank in (0, 12345678901234567890, total - 1):
        subset = unrank_subset(rank, n, k)
        assert len(set(subset)) == k
        assert rank_subset(subset, n) == rank


def test_out_of_range_rank_is_rejected():
    with pytest.raises(ValueError):
        unrank_subset(comb(6, 3), 6, 3)


# --- budget allocation -------------------------------------------------------

def test_small_space_is_fully_exhaustive():
    design = build_design(6, n_cap=10000, seeds=_seeds(),
                          center_ids=tuple(range(1, 7)))
    assert design.K_MAX == 3
    assert design.total_subset_space == 6 + 15 + 20
    assert design.total_evaluated == design.total_subset_space
    assert all(level.mode == "exhaustive" for level in design.levels)
    assert all(level.sampling_fraction == 1.0 for level in design.levels)
    assert len(design.subsets) == design.total_evaluated


def test_budget_exhausts_ascending_then_splits_the_remainder_equally():
    n_s, cap = 40, 1000
    design = build_design(n_s, n_cap=cap, seeds=_seeds(),
                          center_ids=tuple(range(1, n_s + 1)))
    levels = {level.k: level for level in design.levels}

    assert design.K_MAX == 20
    assert levels[1].mode == "exhaustive" and levels[1].n_evaluated == 40
    assert levels[2].mode == "exhaustive" and levels[2].n_evaluated == comb(40, 2)
    # C(40,1) + C(40,2) = 820 <= 1000; C(40,3) would blow the cap
    assert levels[3].mode == "sampled"

    remaining = [levels[k] for k in range(3, 21)]
    budget = cap - 40 - comb(40, 2)
    assert sum(level.n_evaluated for level in remaining) == budget
    counts = {level.n_evaluated for level in remaining}
    assert max(counts) - min(counts) <= 1, "equal split up to integer remainder"
    assert design.total_evaluated <= cap


def test_a_saturated_sampled_level_becomes_exhaustive_and_the_surplus_moves():
    """A level whose allocation reaches T(k) is exhaustive; the surplus is reused."""
    n_s, cap = 12, 300
    design = build_design(n_s, n_cap=cap, seeds=_seeds(),
                          center_ids=tuple(range(1, n_s + 1)))
    levels = {level.k: level for level in design.levels}

    for level in design.levels:
        assert level.n_evaluated <= level.total_subsets
        if level.mode == "exhaustive":
            assert level.n_evaluated == level.total_subsets
    assert design.total_evaluated <= cap
    # k=1 and k=2 are tiny and fully enumerated; k=3 saturates from the split
    assert levels[1].mode == levels[2].mode == "exhaustive"
    assert levels[3].n_evaluated == levels[3].total_subsets


def test_sampled_levels_are_distinct_and_within_range():
    n_s, cap = 30, 500
    design = build_design(n_s, n_cap=cap, seeds=_seeds(),
                          center_ids=tuple(range(1, n_s + 1)))
    assert len(set(design.subsets)) == len(design.subsets)
    for subset in design.subsets:
        assert len(set(subset)) == len(subset)
        assert set(subset) <= set(range(1, n_s + 1))


def test_design_is_order_independent_and_reproducible():
    """Seeds come from context strings, never from execution order."""
    a = build_design(30, 500, _seeds(), tuple(range(1, 31)))
    b = build_design(30, 500, _seeds(), tuple(range(1, 31)))
    assert a.subsets == b.subsets

    # drawing a later level first must not change anything
    shuffled = _seeds()
    for level in (3, 1, 2):
        shuffled.rng(f"subset_design|k={level}")
    c = build_design(30, 500, shuffled, tuple(range(1, 31)))
    assert c.subsets == a.subsets


def test_a_different_run_id_gives_a_different_but_reproducible_sample():
    a = build_design(30, 500, _seeds("RUN_A"), tuple(range(1, 31)))
    b = build_design(30, 500, _seeds("RUN_B"), tuple(range(1, 31)))
    assert a.subsets != b.subsets
    assert a.total_evaluated == b.total_evaluated


def test_design_below_two_centers_is_refused_not_faked():
    with pytest.raises(BlockedError) as excinfo:
        build_design(1, 10000, _seeds(), (5,))
    assert "ROBUSTNESS_NOT_EVALUABLE" in str(excinfo.value)


def test_two_centers_give_the_minimal_honest_design():
    design = build_design(2, 10000, _seeds(), (7, 9))
    assert design.K_MAX == 1
    assert design.total_evaluated == 2
    assert design.subsets == ((7,), (9,))
