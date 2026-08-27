"""II.3 / II.5B — the label-permutation null and the p-value estimator."""
from __future__ import annotations

import numpy as np
import pytest

from hotspot3d.hotspot.permutation import permutation_pass, permute_labels
from hotspot3d.utils.geometry import membership_matrix


@pytest.mark.unit
def test_permutation_preserves_class_counts_exactly(rng):
    y = np.array([1, 1, 1, 0, 0, 0, 0, 0, 1, 0], dtype=np.int64)
    Y = permute_labels(rng, y, 500)
    assert Y.shape == (len(y), 500)
    assert np.all(Y.sum(axis=0) == y.sum())                    # N_P preserved, always
    assert np.all((Y == 0).sum(axis=0) == (y == 0).sum())      # N_B preserved, always
    assert set(np.unique(Y)) <= {0.0, 1.0}


@pytest.mark.unit
def test_permutation_actually_permutes(rng):
    y = np.array([1] * 10 + [0] * 10, dtype=np.int64)
    Y = permute_labels(rng, y, 200)
    assert not np.all(Y == y[:, None])                          # not the identity
    assert len({tuple(col) for col in Y.T}) > 1                 # not one repeated draw


@pytest.mark.unit
def test_p_emp_is_never_zero_and_never_below_the_floor(rng, small_cohort):
    """Phipson-Smyth add-one: p >= 1/(B+1) by construction (II.3)."""
    coords, y = small_cohort
    A = membership_matrix(coords, coords, 2.0)
    B = 500
    out = permutation_pass(A, y.astype(float), rng, B, 2.0, "test")
    assert out.p_emp.min() >= 1.0 / (B + 1) - 1e-15
    assert (out.p_emp > 0).all()
    assert out.p_emp.max() <= 1.0


@pytest.mark.unit
def test_centers_without_labelled_neighbours_leave_the_family(rng):
    coords = np.array([[0.0, 0, 0], [1.0, 0, 0], [200.0, 0, 0]])
    labelled = coords[:2]
    y = np.array([1.0, 0.0])
    A = membership_matrix(coords, labelled, 2.0)
    out = permutation_pass(A, y, rng, 200, 2.0, "test")
    assert out.n_labeled[2] == 0
    assert not out.in_family[2]
    assert out.exclusion_reason[2] == "no_labeled_residue_in_sphere"


@pytest.mark.unit
def test_zero_variance_centers_leave_the_family(rng):
    """A center containing the WHOLE cohort has a constant n_P — sigma_c = 0."""
    coords = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.0, 0, 0]])
    y = np.array([1.0, 0.0, 0.0])
    A = membership_matrix(coords, coords, 100.0)          # every center sees everything
    out = permutation_pass(A, y, rng, 200, 100.0, "test")
    assert np.all(out.null_sd == 0)
    assert not out.in_family.any()
    assert set(out.exclusion_reason) == {"zero_permutation_variance"}


@pytest.mark.unit
def test_observed_counts_match_a_direct_computation(rng, small_cohort):
    coords, y = small_cohort
    A = membership_matrix(coords, coords, 2.0)
    out = permutation_pass(A, y.astype(float), rng, 100, 2.0, "test")
    direct_l = (np.linalg.norm(coords[:, None] - coords[None], axis=2) <= 2.0).sum(axis=1)
    assert np.array_equal(out.n_labeled, direct_l.astype(float))
    assert out.n_plp_obs[0] == 4                        # itself plus three P/LP


@pytest.mark.unit
def test_ratio_statistic_and_zg(rng, small_cohort):
    """T(r) sums n_P(c)/n_L(c) over centers with at least one labelled neighbour."""
    # P/LP concentrated in the dense core, B/LB at the sparse rim: exactly the pattern
    # T(r) = sum_c n_P(c)/n_L(c) is built to reward.
    coords = np.array([[float(x), 0.0, 0.0] for x in range(8)])
    y = np.array([0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])
    A = membership_matrix(coords, coords, 2.5)
    out = permutation_pass(A, y, rng, 400, 2.5, "test")
    informative = out.n_labeled >= 1
    expected = float(np.sum(out.n_plp_obs[informative] / out.n_labeled[informative]))
    assert out.t_obs == pytest.approx(expected)
    assert out.t_obs > out.t_null_mean
    assert out.zg > 0                                    # the planted cluster is real
    assert 0 < out.p_ratio <= 1


@pytest.mark.unit
def test_degenerate_ratio_statistic_yields_zero_evidence(rng, small_cohort):
    """When the spheres partition the cohort, T(r) is constant: Zg = 0, p = 1.

    Float noise on a mathematically constant statistic must never be read as evidence.
    """
    coords, y = small_cohort
    A = membership_matrix(coords, coords, 2.0)          # disjoint blocks
    out = permutation_pass(A, y.astype(float), rng, 400, 2.0, "test")
    assert out.t_null_sd < 1e-12
    assert out.zg == 0.0
    assert out.p_ratio == 1.0


@pytest.mark.unit
def test_same_seed_gives_identical_permutation_results(small_cohort):
    coords, y = small_cohort
    A = membership_matrix(coords, coords, 2.0)
    a = permutation_pass(A, y.astype(float),
                         np.random.Generator(np.random.PCG64(99)), 300, 2.0, "c")
    b = permutation_pass(A, y.astype(float),
                         np.random.Generator(np.random.PCG64(99)), 300, 2.0, "c")
    assert np.array_equal(a.p_emp, b.p_emp)
    assert a.zg == b.zg


@pytest.mark.unit
def test_chunking_does_not_change_the_result(small_cohort):
    """The chunked matrix product is exact integer arithmetic — chunk size is invisible."""
    coords, y = small_cohort
    A = membership_matrix(coords, coords, 2.0)
    a = permutation_pass(A, y.astype(float),
                         np.random.Generator(np.random.PCG64(5)), 600, 2.0, "c", chunk=600)
    b = permutation_pass(A, y.astype(float),
                         np.random.Generator(np.random.PCG64(5)), 600, 2.0, "c", chunk=100)
    assert np.array_equal(a.ge_count, b.ge_count)
    assert np.array_equal(a.null_mean, b.null_mean)
