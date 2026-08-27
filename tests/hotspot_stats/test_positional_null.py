"""Workflow v2 §5.2 — the structure-aware positional null is the PRIMARY test.

Two things are checked here and nowhere else:

  1. the positional pass computes the null it claims to compute — the counts, the
     moments and the tail all agree with the exact hypergeometric law the null
     implies, so the p-values are not merely reproducible but *right*;
  2. the label-permutation null is genuinely SEPARATE — a different hypothesis, a
     different random stream, a different answer, and no path from it into the
     primary result.
"""
from __future__ import annotations

import numpy as np
import pytest

from hotspot3d.hotspot.detection import final_detection
from hotspot3d.hotspot.positional import (
    NULL_MODEL,
    hypergeom_sf,
    positional_combinatorial_floor,
    positional_pass,
)
from hotspot3d.utils.geometry import cross_distances

pytestmark = pytest.mark.unit


def _lattice(n_side: int = 6, spacing: float = 4.0) -> np.ndarray:
    grid = np.arange(n_side) * spacing
    return np.array([[x, y, z] for x in grid for y in grid for z in grid],
                    dtype=np.float64)


def _membership(coords: np.ndarray, radius: float) -> np.ndarray:
    return (cross_distances(coords, coords) <= radius).astype(np.uint8)


# --- the exact tail ----------------------------------------------------------

def test_hypergeometric_tail_matches_scipy():
    """The floor that decides whether a run may make a claim must be exact."""
    scipy_stats = pytest.importorskip("scipy.stats")
    rng = np.random.default_rng(0)
    for _ in range(200):
        M = int(rng.integers(5, 600))
        K = int(rng.integers(1, M))
        n = int(rng.integers(1, M))
        k = int(rng.integers(0, min(K, n) + 2))
        ours = hypergeom_sf(k, M, K, n)
        theirs = float(scipy_stats.hypergeom.sf(k - 1, M, K, n))
        assert ours == pytest.approx(theirs, rel=1e-9, abs=1e-300)


def test_hypergeometric_tail_edges():
    assert hypergeom_sf(0, 100, 10, 5) == 1.0          # every draw satisfies X >= 0
    assert hypergeom_sf(6, 100, 10, 5) == 0.0          # more successes than draws
    assert hypergeom_sf(5, 100, 10, 5) == pytest.approx(
        np.prod([(10 - i) / (100 - i) for i in range(5)]))


# --- the pass itself ---------------------------------------------------------

def test_observed_counts_are_in_sphere_pathogenic_counts():
    coords = _lattice()
    U = _membership(coords, 5.0)
    plp = np.array([0, 1, 6, 7], dtype=np.int64)
    out = positional_pass(U, plp, U[:, plp].sum(axis=1).astype(float),
                          np.random.default_rng(1), 200, 5.0, "t")

    expected = U[:, plp].sum(axis=1)
    assert np.array_equal(out.n_plp_obs.astype(int), expected)
    assert np.array_equal(out.n_universe, U.sum(axis=1))
    assert out.null_model == NULL_MODEL


def test_null_moments_match_the_hypergeometric_law():
    """The null is exactly Hypergeom(M, n_U(c), N_P) — check the mean it implies."""
    coords = _lattice()
    U = _membership(coords, 6.0)
    M = len(coords)
    plp = np.arange(20, dtype=np.int64)
    out = positional_pass(U, plp, U[:, plp].sum(axis=1).astype(float),
                          np.random.default_rng(7), 4000, 6.0, "t")

    analytic = len(plp) * out.n_universe / M
    assert np.allclose(out.null_mean, analytic, atol=0.15)
    assert np.allclose(out.expected_plp, analytic)


def test_empirical_p_tracks_the_exact_tail():
    """p_emp is a Monte-Carlo estimate of the exact tail; the two must agree."""
    coords = _lattice()
    U = _membership(coords, 6.0)
    plp = np.arange(0, 40, 2, dtype=np.int64)
    out = positional_pass(U, plp, U[:, plp].sum(axis=1).astype(float),
                          np.random.default_rng(3), 8000, 6.0, "t")

    # Only where the tail is resolvable at this B; below the floor p_emp saturates.
    resolvable = out.p_exact > 5e-3
    assert resolvable.sum() > 20
    assert np.allclose(out.p_emp[resolvable], out.p_exact[resolvable], atol=0.03)


def test_p_emp_never_goes_below_the_permutation_floor():
    coords = _lattice()
    U = _membership(coords, 5.0)
    plp = np.array([0, 1, 2, 6, 7, 8], dtype=np.int64)
    B = 500
    out = positional_pass(U, plp, U[:, plp].sum(axis=1).astype(float),
                          np.random.default_rng(5), B, 5.0, "t")
    assert out.p_emp.min() >= 1.0 / (B + 1) - 1e-15
    assert (out.p_emp > 0).all()


def test_the_family_is_the_centers_holding_a_classified_residue():
    """v2 §5.1 — m is the number of centers actually tested and sets the bar."""
    coords = _lattice()
    U = _membership(coords, 4.5)
    plp = np.array([0], dtype=np.int64)
    n_labeled = np.zeros(len(coords))
    n_labeled[:10] = 2.0                       # only the first ten see any label
    out = positional_pass(U, plp, n_labeled, np.random.default_rng(2), 200, 4.5, "t")

    assert out.in_family[:10].all() or out.exclusion_reason[0] == "zero_permutation_variance"
    assert not out.in_family[10:].any()
    assert set(np.asarray(out.exclusion_reason)[10:]) == {"no_labeled_residue_in_sphere"}


def test_combinatorial_floor_is_not_monotone_in_sphere_size():
    """The most favourable occupancy is n_U = N_P, not the largest sphere.

    This is why the certificate minimises over the occupancies that occur rather
    than evaluating "the biggest sphere" and calling it the best case.
    """
    M, N_P = 499, 38
    floors = {n_u: float(positional_combinatorial_floor(np.array([n_u]), M, N_P)[0])
              for n_u in (5, 20, 38, 60, 120)}
    assert floors[38] == min(floors.values())
    assert floors[120] > floors[38]
    assert floors[5] > floors[38]


def test_same_seed_reproduces_the_pass_bit_for_bit():
    coords = _lattice()
    U = _membership(coords, 5.0)
    plp = np.array([0, 3, 9, 27], dtype=np.int64)
    n_l = U[:, plp].sum(axis=1).astype(float)
    a = positional_pass(U, plp, n_l, np.random.default_rng(11), 300, 5.0, "t")
    b = positional_pass(U, plp, n_l, np.random.default_rng(11), 300, 5.0, "t")
    assert np.array_equal(a.p_emp, b.p_emp)
    assert np.array_equal(a.ge_count, b.ge_count)


# --- separation from the secondary null --------------------------------------

def _detection(label: bool, seed: int = 4):
    """A clustered cohort inside a much larger universe."""
    coords = _lattice(7, 4.0)
    universe_index = np.arange(1, len(coords) + 1)
    # 12 P/LP packed into one corner, 12 B/LB scattered far away.
    order = np.argsort(np.linalg.norm(coords - coords[0], axis=1))
    plp_positions = order[:12]
    blb_positions = order[-12:]
    positions = np.concatenate([plp_positions, blb_positions]).astype(np.int64)
    y = np.array([1.0] * 12 + [0.0] * 12)
    rng = np.random.default_rng(seed)
    return final_detection(
        6.0, coords, universe_index, positions, coords[positions], y, rng,
        B=2000, q=0.05, seed_context="primary",
        label_rng=(np.random.default_rng(seed + 1) if label else None),
        label_seed_context="secondary")


def test_the_two_nulls_are_computed_separately_and_reported_side_by_side():
    det = _detection(label=True)
    assert det.counts["primary_null"] == NULL_MODEL
    assert det.counts["secondary_null"] == "label_permutation"
    assert det.secondary["computed"] is True
    assert det.secondary["IS_SECONDARY_NEVER_PRIMARY"] is True
    assert det.secondary["hypothesis"] != det.secondary["primary_hypothesis_for_contrast"]

    # significance / q_bh / the residue objects come from the PRIMARY null only.
    primary_sig = {r["center_residue_index"] for r in det.significant_rows}
    assert primary_sig == {r["center_residue_index"] for r in det.center_rows
                           if r["significant"]}
    assert all(r["null_model"] == NULL_MODEL for r in det.center_rows)

    # A center the secondary null likes but the primary does not must NOT appear in
    # any published object — that is what "never mixed" has to mean operationally.
    secondary_only = set(det.secondary["significant_center_residue_indices"]) - primary_sig
    covered = {r["covered_residue_index"] for r in det.covered_rows}
    region_members = {int(i) for r in det.region_rows
                      for i in r["center_residue_indices"].split(";")}
    assert not (secondary_only & region_members)
    assert primary_sig <= covered
    for row in det.center_rows:
        if row["center_residue_index"] in secondary_only:
            assert row["significant"] is False
            assert row["significant_label_null"] is True


def test_the_secondary_null_never_changes_the_primary_result():
    """Computing the label null must leave every published primary number identical."""
    with_label = _detection(label=True)
    without_label = _detection(label=False)

    assert [r["center_residue_index"] for r in with_label.significant_rows] == \
        [r["center_residue_index"] for r in without_label.significant_rows]
    for a, b in zip(with_label.center_rows, without_label.center_rows):
        for key in ("p_emp", "q_bh", "significant", "z_score", "p_exact_hypergeom"):
            assert a[key] == b[key], key
    assert without_label.secondary["computed"] is False


def test_the_positional_null_detects_what_the_label_null_cannot():
    """The nulls answer different questions, so their answers may differ.

    A tight P/LP cluster inside a large universe is extreme under uniform placement
    over U_struct; whether it is extreme under label permutation depends only on how
    the labels sit among the variant-bearing positions.
    """
    det = _detection(label=True)
    assert det.n_significant > 0, "the positional null must reject a planted cluster"
    assert det.n_significant >= det.secondary["n_significant_bh"]
