"""Benjamini-Hochberg at q = 0.05 — boundary ties, monotonicity, BY comparison."""
from __future__ import annotations

import numpy as np
import pytest

from hotspot3d.hotspot.fdr import benjamini_hochberg, bh_adjust, by_constant

Q = 0.05


@pytest.mark.unit
def test_bh_rejects_the_textbook_set():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216])
    out = benjamini_hochberg(p, Q)
    # Largest k with p_(k) <= k q / m is k = 4 (0.041 <= 0.02? no) -> k = 2 (0.008 <= 0.01)
    assert out.boundary_rank == 2
    assert out.boundary_p == pytest.approx(0.008)
    assert out.n_reject_bh == 2
    assert list(np.nonzero(out.reject_bh)[0]) == [0, 1]


@pytest.mark.unit
def test_all_boundary_ties_are_rejected():
    """Ties at the boundary are rejected TOGETHER — never an arbitrary subset."""
    p = np.array([0.004, 0.004, 0.004, 0.004, 0.9, 0.9, 0.9, 0.9])
    out = benjamini_hochberg(p, Q)
    assert out.n_reject_bh == 4
    assert np.array_equal(out.reject_bh, p <= out.boundary_p)
    # Sort order must not decide anything: reversing the input reverses only the mask.
    reversed_out = benjamini_hochberg(p[::-1], Q)
    assert reversed_out.n_reject_bh == out.n_reject_bh
    assert reversed_out.boundary_p == out.boundary_p


@pytest.mark.unit
def test_adjusted_p_values_are_monotone_and_capped():
    rng = np.random.default_rng(4)
    p = np.sort(rng.uniform(0, 1, 200))
    q_bh = bh_adjust(p)
    assert np.all(np.diff(q_bh) >= -1e-15)                # monotone in p
    assert q_bh.max() <= 1.0 and q_bh.min() >= 0.0
    assert np.all(q_bh >= p - 1e-15)                      # adjustment never shrinks p


@pytest.mark.unit
def test_threshold_rule_and_adjusted_rule_agree_exactly():
    """The two published columns can never disagree — asserted inside the estimator."""
    rng = np.random.default_rng(11)
    for _ in range(50):
        p = rng.uniform(0, 1, rng.integers(2, 120))
        out = benjamini_hochberg(p, Q)
        assert np.array_equal(out.reject_bh, out.q_bh <= Q)


@pytest.mark.unit
def test_by_is_more_conservative_and_reported_alongside():
    p = np.array([0.0001, 0.002, 0.01, 0.03, 0.2, 0.5, 0.7, 0.9])
    out = benjamini_hochberg(p, Q)
    assert out.by_constant == pytest.approx(by_constant(len(p)))
    assert out.by_constant > 1.0
    assert np.all(out.q_by >= out.q_bh - 1e-15)
    assert out.n_reject_by <= out.n_reject_bh


@pytest.mark.unit
def test_no_rejection_leaves_no_boundary():
    out = benjamini_hochberg(np.array([0.4, 0.5, 0.6, 0.9]), Q)
    assert out.n_reject_bh == 0
    assert out.boundary_p is None and out.boundary_rank is None


@pytest.mark.unit
def test_empty_family_is_handled_explicitly():
    out = benjamini_hochberg(np.array([]), Q)
    assert out.m == 0 and out.n_reject_bh == 0 and out.boundary_p is None


@pytest.mark.unit
def test_a_discrete_p_value_floor_can_block_every_rejection():
    """With m large and B small, even p at the floor cannot clear q/m — II.4's premise."""
    B, m = 200, 5000
    p = np.full(m, 1.0 / (B + 1))
    p[10:] = 0.6
    out = benjamini_hochberg(p, Q)
    assert out.n_reject_bh == 0                # 10 * q / m = 1e-4 < p_res = 4.98e-3
