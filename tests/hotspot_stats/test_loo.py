"""II.5A — leave-one-out MCC: hold-out completeness, kappa neutrality, conventions."""
from __future__ import annotations

import numpy as np
import pytest

from hotspot3d.hotspot.loo import kappa_is_boundary_neutral, leave_one_out
from hotspot3d.utils.geometry import mcc_from_confusion

KAPPA = 2.0


@pytest.mark.unit
def test_holdout_excludes_residue_from_its_own_neighbourhood(small_cohort):
    """Residue i never counts itself among its neighbours."""
    coords, y = small_cohort
    out = leave_one_out(coords, y, radius=2.0, kappa=KAPPA)

    # The four cluster residues sit within 2 A of each other but not of themselves.
    assert out.n_labeled[0] == 3
    assert out.n_plp[0] == 3                       # the other three P/LP, not itself
    # A brute-force recomputation with i explicitly removed must agree exactly.
    for i in range(len(y)):
        mask = np.ones(len(y), dtype=bool)
        mask[i] = False
        d = np.linalg.norm(coords[mask] - coords[i], axis=1)
        assert out.n_labeled[i] == int((d <= 2.0).sum())
        assert out.n_plp[i] == int(((d <= 2.0) & (y[mask] == 1)).sum())


@pytest.mark.unit
def test_holdout_excludes_residue_from_the_background_rate(small_cohort):
    """The prevalence is held out too: pi_-i = (N_P - 1[y_i=P/LP]) / (N - 1)."""
    coords, y = small_cohort
    n, n_p = len(y), int(y.sum())
    out = leave_one_out(coords, y, radius=2.0, kappa=KAPPA)

    for i in range(n):
        expected = (n_p - int(y[i])) / (n - 1)
        assert out.pi_minus[i] == pytest.approx(expected)
    # A P/LP residue is judged against a LOWER prevalence than a B/LB residue: this is
    # exactly the leakage the two-sided hold-out removes.
    assert out.pi_minus[y == 1].max() < out.pi_minus[y == 0].min()
    assert not np.allclose(out.pi_minus, n_p / n)


@pytest.mark.unit
@pytest.mark.parametrize("radius", [1.5, 2.0, 5.0, 45.0])
def test_kappa_is_boundary_neutral(small_cohort, radius):
    """f_i > pi_-i  <=>  n_P(i) > pi_-i * n_L(i): kappa cancels exactly."""
    coords, y = small_cohort
    out = leave_one_out(coords, y, radius=radius, kappa=KAPPA)
    assert kappa_is_boundary_neutral(out.n_plp, out.n_labeled, out.pi_minus, KAPPA)

    direct = out.n_plp > out.pi_minus * out.n_labeled
    predicted = out.predictions.astype(bool)
    # Zero-neighbour residues are forced to B/LB regardless, so compare where n_L > 0.
    informative = out.n_labeled > 0
    assert np.array_equal(predicted[informative], direct[informative])


@pytest.mark.unit
@pytest.mark.parametrize("kappa", [0.5, 1.0, 2.0, 8.0])
def test_predictions_do_not_depend_on_kappa(small_cohort, kappa):
    """Smoothing can never be used to shift the decision boundary (F13)."""
    coords, y = small_cohort
    reference = leave_one_out(coords, y, radius=5.0, kappa=2.0)
    other = leave_one_out(coords, y, radius=5.0, kappa=kappa)
    assert np.array_equal(reference.predictions, other.predictions)


@pytest.mark.unit
def test_zero_neighbour_resolves_to_blb(small_cohort):
    coords, y = small_cohort
    out = leave_one_out(coords, y, radius=0.5, kappa=KAPPA)     # nobody has a neighbour
    assert out.n_zero_neighbour == len(y)
    assert out.predictions.sum() == 0                            # all predicted B/LB
    assert bool(out.tie.all())                                   # f_i == pi_-i exactly


@pytest.mark.unit
def test_exact_tie_resolves_to_blb():
    """|f_i - pi_-i| <= 1e-12 predicts B/LB, never P/LP."""
    # Two P/LP and two B/LB in one tight group: every residue sees a neighbourhood whose
    # pathogenic fraction can be made to sit exactly on its held-out prevalence.
    coords = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.0, 0, 0], [3.0, 0, 0]])
    y = np.array([1, 0, 1, 0])
    out = leave_one_out(coords, y, radius=10.0, kappa=KAPPA)
    # residue 0: n_L = 3, n_P = 1, pi_-0 = (2-1)/3 = 1/3  ->  n_P == pi * n_L exactly
    assert out.n_labeled[0] == 3 and out.n_plp[0] == 1
    assert out.pi_minus[0] == pytest.approx(1 / 3)
    assert bool(out.tie[0])
    assert out.predictions[0] == 0


@pytest.mark.unit
def test_no_residue_is_dropped_for_sparsity(small_cohort):
    """n_L < 3 is a REPORTING cut-off; every residue stays in the confusion matrix."""
    coords, y = small_cohort
    out = leave_one_out(coords, y, radius=2.0, kappa=KAPPA, sparse_cutoff=3)
    m = out.metrics
    assert m["tp"] + m["fp"] + m["tn"] + m["fn"] == len(y)
    assert out.n_sparse > 0                       # sparse residues exist ...
    assert len(out.predictions) == len(y)         # ... and none was removed


@pytest.mark.unit
def test_sparse_diagnostics_are_split_by_true_and_predicted_class(small_cohort):
    coords, y = small_cohort
    out = leave_one_out(coords, y, radius=2.0, kappa=KAPPA)
    scopes = {row["scope"] for row in out.diagnostics}
    assert scopes == {"overall", "true_PLP", "true_BLB", "pred_PLP", "pred_BLB"}
    overall = next(r for r in out.diagnostics if r["scope"] == "overall")
    assert overall["n_residues"] == len(y)
    assert overall["n_l_median"] is not None


@pytest.mark.unit
def test_mcc_zero_denominator_convention():
    """A degenerate confusion matrix yields the recorded convention value, not NaN."""
    assert mcc_from_confusion(0, 0, 8, 0, undefined_value=0.0) == 0.0
    assert mcc_from_confusion(4, 0, 0, 0, undefined_value=0.0) == 0.0
    assert mcc_from_confusion(2, 0, 2, 0) == pytest.approx(1.0)


@pytest.mark.unit
def test_mcc_is_perfect_on_a_separable_cohort(small_cohort):
    coords, y = small_cohort
    out = leave_one_out(coords, y, radius=2.0, kappa=KAPPA)
    assert out.mcc > 0.0
    assert out.metrics["tp"] == 4                 # the planted cluster is recovered
