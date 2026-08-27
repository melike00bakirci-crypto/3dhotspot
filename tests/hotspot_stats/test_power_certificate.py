"""Workflow v2 §5.4 — the power certificate.

The certificate decides whether a run is *allowed to make a claim*, so its
arithmetic is checked against closed forms rather than against itself, and both
floors are exercised separately: one is remediable by ``B`` and one is not, and
conflating them is how an unremediable design gets told to buy more permutations.
"""
from __future__ import annotations

from math import comb

import numpy as np
import pytest

from hotspot3d.hotspot.power import (
    BINDING_COMBINATORIAL,
    BINDING_PERMUTATION,
    LABEL_NULL,
    POSITIONAL_NULL,
    UnderpoweredResult,
    b_required_to_clear,
    best_attainable_floor,
    certify,
    domain_occupancies,
    label_null_floor,
    positional_null_floor,
    preflight,
)

pytestmark = pytest.mark.unit


# --- the two floors ----------------------------------------------------------

def test_label_null_floor_is_the_hypergeometric_tail():
    """C(N_P, n_L) / C(N_P + N_B, n_L) — the all-pathogenic sphere."""
    assert label_null_floor(38, 12, 20) == pytest.approx(comb(38, 20) / comb(50, 20))
    assert label_null_floor(10, 10, 3) == pytest.approx(comb(10, 3) / comb(20, 3))
    # n_L > N_P: the extreme is "every pathogenic residue is inside".
    assert label_null_floor(3, 10, 5) == pytest.approx(
        comb(3, 3) * comb(10, 2) / comb(13, 5))


def test_label_null_floor_falls_as_the_sphere_fills():
    floors = [label_null_floor(38, 12, n) for n in range(1, 20)]
    assert all(b < a for a, b in zip(floors, floors[1:]))


def test_positional_floor_is_orders_below_the_label_floor_at_the_same_cohort():
    """The §5.2 change is not cosmetic: it moves the floor by orders of magnitude."""
    label = label_null_floor(38, 12, 20)
    positional = positional_null_floor(M=499, n_U=20, N_P=38)
    assert positional < label / 1e15


def test_best_attainable_floor_minimises_over_realised_occupancies():
    occ_u = np.array([5, 20, 38, 120])
    best, detail = best_attainable_floor(
        POSITIONAL_NULL, N_P=38, N_B=12, M=499,
        occupancies_labeled=np.array([]), occupancies_universe=occ_u)
    assert detail["best_occupancy"] == 38
    assert detail["largest_occupancy"] == 120
    assert best < detail["floor_at_largest_occupancy"]
    assert best == pytest.approx(positional_null_floor(499, 38, 38))


def test_best_attainable_floor_for_the_label_null_is_the_largest_sphere():
    occ_l = np.array([3, 11, 20])
    best, detail = best_attainable_floor(
        LABEL_NULL, N_P=38, N_B=12, M=499,
        occupancies_labeled=occ_l, occupancies_universe=np.array([]))
    assert detail["best_occupancy"] == detail["largest_occupancy"] == 20
    assert best == pytest.approx(label_null_floor(38, 12, 20))


# --- the certificate ---------------------------------------------------------

def _cert(null_model, *, B, q, m, N_P, N_B, M, occ_l, occ_u, phase="posthoc"):
    return certify(phase, null_model=null_model, B=B, q=q, m=m, N_P=N_P, N_B=N_B,
                   n_universe=M, n_center_universe=M,
                   occupancies_labeled=np.asarray(occ_l),
                   occupancies_universe=np.asarray(occ_u))


def test_certificate_arithmetic():
    cert = _cert(POSITIONAL_NULL, B=10_000, q=0.05, m=452, N_P=38, N_B=12, M=499,
                 occ_l=[20], occ_u=[38])
    assert cert.p_res == pytest.approx(1 / 10_001)
    assert cert.c_1 == pytest.approx(0.05 / 452)
    assert cert.p_floor == max(cert.p_res, cert.p_comb_best)
    assert cert.ratio == pytest.approx(cert.p_floor / cert.c_1)
    assert cert.passes is (cert.p_floor <= cert.c_1)


def test_permutation_floor_binds_and_names_the_B_that_would_clear_c1():
    cert = _cert(POSITIONAL_NULL, B=200, q=0.05, m=1000, N_P=30, N_B=30, M=400,
                 occ_l=[10], occ_u=[30])
    assert cert.binding_floor == BINDING_PERMUTATION
    assert cert.increasing_B_is_futile is False
    assert cert.passes is False
    need = cert.B_to_clear_c_1
    assert 1 / (need + 1) <= cert.c_1 < 1 / need
    assert any(f"B = {need}" in r for r in cert.remedies)
    assert not any("FUTILE" in r for r in cert.remedies)


def test_combinatorial_floor_binds_and_forbids_recommending_more_permutations():
    cert = _cert(LABEL_NULL, B=10_000, q=0.05, m=452, N_P=38, N_B=12, M=499,
                 occ_l=[20], occ_u=[38])
    assert cert.binding_floor == BINDING_COMBINATORIAL
    assert cert.increasing_B_is_futile is True
    assert cert.passes is False
    assert any("FUTILE" in r for r in cert.remedies)
    assert any("benign cohort" in r for r in cert.remedies)
    assert any("different null" in r for r in cert.remedies)
    assert any("smaller test family" in r for r in cert.remedies)


def test_a_certificate_that_passes_recommends_nothing():
    cert = _cert(POSITIONAL_NULL, B=10_000, q=0.05, m=100, N_P=30, N_B=30, M=400,
                 occ_l=[10], occ_u=[30])
    assert cert.passes is True
    assert cert.remedies == []
    assert cert.as_json()["TEST_CANNOT_REJECT"] is False
    assert "uninformative" not in cert.as_json()["interpretation"]


def test_failed_certificate_json_carries_the_prohibited_inferences():
    cert = _cert(LABEL_NULL, B=10_000, q=0.05, m=452, N_P=38, N_B=12, M=499,
                 occ_l=[20], occ_u=[38])
    payload = cert.as_json()
    assert payload["TEST_CANNOT_REJECT"] is True
    assert "uninformative about the presence or absence" in payload["interpretation"]
    assert "No 3D hotspot is detectable in this gene." in \
        payload["prohibited_inferences_if_failed"]


def test_an_empty_test_family_never_reads_as_a_vacuous_pass():
    """m = 0 means no hypothesis exists; `p_floor <= inf` must not certify that."""
    cert = _cert(POSITIONAL_NULL, B=10_000, q=0.05, m=0, N_P=10, N_B=10, M=100,
                 occ_l=[], occ_u=[])
    assert cert.passes is False
    assert cert.binding_floor == "empty_test_family"
    assert cert.B_to_clear_c_1 is None
    assert cert.increasing_B_is_futile is True


def test_b_required_to_clear_is_the_smallest_sufficient_B():
    for c_1 in (0.05, 1e-3, 1.106e-4, 3.7e-5):
        b = b_required_to_clear(c_1)
        assert 1 / (b + 1) <= c_1
        assert 1 / b > c_1


# --- the pre-flight instance -------------------------------------------------

def _coords(n_side=6, spacing=4.0):
    grid = np.arange(n_side) * spacing
    return np.array([[x, y, z] for x in grid for y in grid for z in grid],
                    dtype=np.float64)


def test_domain_occupancies_count_spheres_at_every_radius():
    coords = _coords(4)
    labeled = np.array([0, 1, 2, 3], dtype=np.int64)
    occ = domain_occupancies(coords, labeled, [4.0, 8.0])
    assert [o.radius_A for o in occ] == [4.0, 8.0]
    assert (occ[1].n_universe >= occ[0].n_universe).all()   # monotone in r
    assert (occ[0].n_universe >= 1).all()                   # a center covers itself
    assert occ[0].m == int((occ[0].n_labeled >= 1).sum())


def test_preflight_representative_is_the_most_favourable_radius():
    """Terminating early must mean 'no radius could reject', not 'some could not'."""
    coords = _coords(6)
    labeled = np.arange(30, dtype=np.int64)
    occ = domain_occupancies(coords, labeled, [4.0, 6.0, 8.0, 12.0])
    best, rows = preflight(POSITIONAL_NULL, B=10_000, q=0.05, N_P=15, N_B=15,
                           n_universe=len(coords), n_center_universe=len(coords),
                           occupancies=occ)
    assert len(rows) == 4
    assert best.ratio == min(r["p_floor_over_c_1"] for r in rows)
    assert best.detail["n_radii_evaluated"] == 4
    assert best.detail["n_radii_passing"] == sum(1 for r in rows if r["passes"])


def test_underpowered_result_is_neither_a_negative_nor_a_blocked_error():
    """The terminal state must not be typeable as a finding or as a fault."""
    from hotspot3d.utils.errors import BlockedError, NegativeResult

    cert = _cert(LABEL_NULL, B=10_000, q=0.05, m=452, N_P=38, N_B=12, M=499,
                 occ_l=[20], occ_u=[38])
    exc = UnderpoweredResult(cert, "06_FINAL_HOTSPOTS")
    assert not isinstance(exc, NegativeResult)
    assert not isinstance(exc, BlockedError)
    assert exc.certificate is cert
    assert "TEST_CANNOT_REJECT" in str(exc)
