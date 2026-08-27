"""Workflow v2 Appendix B — the KCNA2 regression case.

Run ``20260817T155124Z_df8ec59f_f30b0d52`` reported zero significant hotspot centers
as a scientific negative. It was not one: under the configuration the agent derived
from v1, the minimum p-value the per-center test could return (7.125e-04) exceeded
the rank-1 Benjamini-Hochberg critical value (1.106e-04), so no center could have
been declared significant regardless of the data.

| quantity                        | value      |
|---------------------------------|------------|
| N_P / N_B                       | 38 / 12    |
| \\|U_struct\\|                   | 499        |
| m (family size, as executed)    | 452        |
| null used for the per-center test | label permutation |
| best observed configuration     | center 335, r = 21 A, 20 of 20 classified P/LP |
| combinatorial floor C(38,20)/C(50,20) | 7.125e-04 |
| rank-1 critical value q/m       | 1.106e-04  |
| verdict under §5.4              | p_floor > c_1 by a factor of 6.4 -> UNDERPOWERED |

**These are regression assertions on the executed configuration, not targets.**
Nothing anywhere is tuned to make KCNA2 positive; the second half of this module
asserts only that the §5.2-conforming null has rejection capability *in principle*
at this cohort size, which is what makes the first half a statement about the design
rather than about the gene.
"""
from __future__ import annotations

from math import comb

import numpy as np
import pytest

from hotspot3d.hotspot.positional import hypergeom_sf, positional_pass
from hotspot3d.hotspot.power import (
    BINDING_COMBINATORIAL,
    BINDING_PERMUTATION,
    LABEL_NULL,
    POSITIONAL_NULL,
    certify,
)
from hotspot3d.utils.geometry import cross_distances

pytestmark = pytest.mark.unit

# Appendix B, verbatim.
N_P, N_B = 38, 12
U_STRUCT = 499
M_FAMILY = 452
B = 10_000
Q = 0.05
BEST_SPHERE_N_L = 20            # center 335 at r = 21 A: 20 of 20 classified P/LP


def _cert(null_model, occ_l, occ_u, m=M_FAMILY, b=B):
    return certify("posthoc", null_model=null_model, B=b, q=Q, m=m, N_P=N_P, N_B=N_B,
                   n_universe=U_STRUCT, n_center_universe=U_STRUCT,
                   occupancies_labeled=np.asarray(occ_l),
                   occupancies_universe=np.asarray(occ_u), radius_A=21.0)


# --- regression 1: the executed configuration is UNDERPOWERED ----------------

def test_appendix_b_numbers_are_reproduced_exactly():
    cert = _cert(LABEL_NULL, occ_l=[BEST_SPHERE_N_L], occ_u=[])
    assert cert.p_comb_best == pytest.approx(comb(38, 20) / comb(50, 20))
    assert cert.p_comb_best == pytest.approx(7.125e-04, rel=1e-3)
    assert cert.c_1 == pytest.approx(1.106e-04, rel=1e-3)
    assert cert.p_floor == pytest.approx(cert.p_comb_best)
    assert cert.ratio == pytest.approx(6.4, rel=0.02)


def test_the_executed_kcna2_configuration_returns_underpowered():
    """The v2 verdict on the run that motivated the specification revision."""
    cert = _cert(LABEL_NULL, occ_l=[BEST_SPHERE_N_L], occ_u=[])
    assert cert.passes is False
    assert cert.as_json()["TEST_CANNOT_REJECT"] is True
    assert cert.binding_floor == BINDING_COMBINATORIAL


def test_kcna2_is_never_reportable_as_a_valid_negative():
    """v2 Appendix A — the certificate failed, so no absence claim is licensed."""
    payload = _cert(LABEL_NULL, occ_l=[BEST_SPHERE_N_L], occ_u=[]).as_json()
    assert "uninformative about the presence or absence" in payload["interpretation"]
    assert payload["prohibited_inferences_if_failed"] == [
        "No 3D hotspot is detectable in this gene.",
        "This is a valid, complete scientific negative.",
        "This result licenses no modification of the method.",
    ]


def test_raising_B_could_never_have_rescued_the_kcna2_run():
    """The binding floor is combinatorial, so more permutations buy nothing."""
    for b in (1e4, 1e5, 1e6, 1e7, 1e9):
        cert = _cert(LABEL_NULL, occ_l=[BEST_SPHERE_N_L], occ_u=[], b=int(b))
        assert cert.passes is False
        assert cert.increasing_B_is_futile is True
        assert cert.binding_floor == BINDING_COMBINATORIAL
        assert not any("Re-run at B" in r for r in cert.remedies)


def test_no_occupancy_reached_by_the_kcna2_run_could_have_cleared_the_bar():
    """Not just the best observed sphere: no occupancy the run reached clears c_1.

    The executed configuration's largest sphere held 20 classified residues, so the
    certificate is evaluated over ``n_L = 1..20``. Every one of them fails.
    """
    cert = _cert(LABEL_NULL, occ_l=list(range(1, BEST_SPHERE_N_L + 1)), occ_u=[])
    assert cert.passes is False
    assert cert.p_floor > cert.c_1
    assert cert.detail["combinatorial_floor_detail"]["best_occupancy"] == BEST_SPHERE_N_L


def test_the_floor_is_evaluated_over_occupancies_that_actually_occur():
    """A sphere absorbing almost the whole cohort would have a tiny floor.

    ``n_L = 40`` of 50 classified residues gives a formally minuscule tail, so the
    certificate would "pass" if it were allowed to invent occupancies. It is not:
    the floor is minimised over the occupancies the family actually realises, and a
    hypothetical sphere that never occurred can never be used to claim power.
    """
    hypothetical = _cert(LABEL_NULL, occ_l=[40], occ_u=[])
    assert hypothetical.passes is True                # arithmetically, in isolation
    reached = _cert(LABEL_NULL, occ_l=list(range(1, BEST_SPHERE_N_L + 1)), occ_u=[])
    assert reached.passes is False                    # what the run could actually do


# --- regression 2: the §5.2-conforming null CAN reject -----------------------

def test_the_positional_null_clears_c1_by_orders_of_magnitude():
    """Same cohort, same family, the §5.2 null: the certificate passes.

    Appendix B records p = 1.2e-07 at r = 10 A around center 404 — three orders below
    c_1. Under the positional null the attainable tail is far below that, so the only
    remaining floor is the permutation resolution 1/(B+1), which sits BELOW c_1 at
    the frozen B = 10000.
    """
    cert = _cert(POSITIONAL_NULL, occ_l=[], occ_u=[20, 30, 38, 60])
    assert cert.passes is True
    assert cert.binding_floor == BINDING_PERMUTATION
    assert cert.p_comb_best < 1e-7 / 1000
    assert cert.p_floor == pytest.approx(1 / (B + 1))
    assert cert.p_floor < cert.c_1


def test_an_exact_positional_tail_at_kcna2_scale_reaches_appendix_b_magnitude():
    """A sphere of ~40 residues holding 14 of the 38 P/LP gives p ~ 1.2e-07.

    This is the arithmetic behind Appendix B's r = 10 A / center 404 row: it shows
    the null has rejection capability at THIS cohort size, which is the claim the
    regression needs. It is not a re-analysis of KCNA2 and no real data is read.
    """
    c_1 = Q / M_FAMILY
    p = hypergeom_sf(14, U_STRUCT, 40, N_P)
    assert p == pytest.approx(1.2e-07, rel=0.2)      # the value Appendix B records
    assert c_1 / p > 900                             # "three orders of magnitude"


@pytest.mark.integration
def test_the_positional_null_rejects_a_planted_cluster_at_kcna2_scale():
    """End to end at Appendix B's dimensions: N_P = 38, N_B = 12, |U_struct| = 499.

    The point is capability, not KCNA2: with a genuine cluster present, the primary
    null returns p-values that clear the same bar the executed configuration could
    never have cleared.
    """
    rng = np.random.default_rng(20250101)
    coords = rng.normal(0.0, 22.0, size=(U_STRUCT, 3))
    centre = coords[0]
    order = np.argsort(np.linalg.norm(coords - centre, axis=1))
    plp = order[:24].tolist() + rng.choice(order[200:], size=14, replace=False).tolist()
    blb = rng.choice(order[300:], size=N_B, replace=False).tolist()
    plp = np.array(sorted(set(plp))[:N_P], dtype=np.int64)
    labeled = np.array(sorted(set(plp.tolist()) | set(int(b) for b in blb)),
                       dtype=np.int64)

    radius = 12.0
    U = (cross_distances(coords, coords) <= radius).astype(np.uint8)
    n_labeled = U[:, labeled].sum(axis=1).astype(float)
    out = positional_pass(U, plp, n_labeled, np.random.default_rng(7), B, radius, "reg")

    family = out.in_family
    assert family.sum() > 0
    c_1 = Q / int(family.sum())
    assert out.p_emp[family].min() <= c_1, \
        "the primary null must be able to clear the rank-1 bar it is judged against"
    assert out.p_exact[family].min() < c_1 / 1000
