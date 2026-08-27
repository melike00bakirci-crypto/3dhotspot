"""METHOD_SPEC II.5A — leave-one-out MCC.

**This is a RADIUS-SELECTION OBJECTIVE AND NOTHING ELSE (F7/F13).** It is never
hotspot validation, it is never reported as validation, and there is no post-hoc
hotspot LOO stage anywhere in this architecture (F9). ``config.loo_mcc.role`` is
asserted to be ``radius_selection_only`` before this module is used.

The hold-out is **complete and two-sided**: residue ``i`` is removed from its own
neighbourhood counts *and* from the background rate against which it is judged.

    pi_-i = (N_P - 1[y_i = P/LP]) / (N - 1)
    f_i   = (n_P(i) + kappa * pi_-i) / (n_L(i) + kappa)          kappa = 2, FROZEN
    y_i   = P/LP  iff  f_i > pi_-i   (strictly)

**No classification threshold is fitted, per radius or at all.** The decision
boundary is the held-out prevalence itself, and the smoothing is *boundary-neutral*
by construction:

    f_i > pi_-i
      <=>  n_P(i) + kappa*pi_-i > pi_-i * (n_L(i) + kappa)
      <=>  n_P(i) > pi_-i * n_L(i)

so ``kappa`` cancels exactly and can never be used to shift the boundary. That
consequence is reported, not hidden — see :func:`kappa_is_boundary_neutral`.

Zero-neighbour residues and exact ties (``|f_i - pi_-i| <= 1e-12``) both resolve to
B/LB. **No residue is ever dropped for sparsity**: ``n_L < 3`` is a *reporting*
cut-off that appears in the diagnostics and never in a filter.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..utils.geometry import classification_metrics, pairwise_distances

SPARSE_COLUMNS = [
    "radius_A", "scope", "n_residues", "n_zero_neighbour", "prop_zero_neighbour",
    "n_sparse", "prop_sparse", "n_tie", "prop_tie", "n_l_min", "n_l_q1",
    "n_l_median", "n_l_q3", "n_l_max", "n_l_mean",
]


@dataclass
class LooOutcome:
    """LOO results at one radius, with the mandatory sparse-evidence reporting."""

    radius_A: float
    kappa: float
    mcc: float
    metrics: dict
    predictions: np.ndarray          # 1 = predicted P/LP
    n_labeled: np.ndarray            # n_L(i), self excluded
    n_plp: np.ndarray                # n_P(i), self excluded
    pi_minus: np.ndarray
    f_hat: np.ndarray
    tie: np.ndarray
    zero_neighbour: np.ndarray
    sparse: np.ndarray
    diagnostics: list[dict] = field(default_factory=list)

    @property
    def n_zero_neighbour(self) -> int:
        return int(self.zero_neighbour.sum())

    @property
    def n_sparse(self) -> int:
        return int(self.sparse.sum())

    @property
    def n_tie(self) -> int:
        return int(self.tie.sum())

    @property
    def prop_zero_neighbour(self) -> float:
        return float(self.zero_neighbour.mean()) if len(self.zero_neighbour) else 0.0

    @property
    def prop_sparse(self) -> float:
        return float(self.sparse.mean()) if len(self.sparse) else 0.0


def kappa_is_boundary_neutral(n_plp: np.ndarray, n_labeled: np.ndarray,
                              pi_minus: np.ndarray, kappa: float,
                              tol: float = 1e-12) -> bool:
    """Verify ``f_i > pi_-i  <=>  n_P(i) > pi_-i * n_L(i)`` at the given kappa.

    Called by the scan so the claim is *checked at run time*, not merely asserted
    in prose. A violation would mean the smoothing had become a tunable decision
    boundary, which F13 forbids.
    """
    f = (n_plp + kappa * pi_minus) / (n_labeled + kappa)
    smoothed = (f - pi_minus) > tol
    direct = (n_plp - pi_minus * n_labeled) > tol * (n_labeled + kappa)
    return bool(np.array_equal(smoothed, direct))


def _distribution(values: np.ndarray) -> dict:
    if len(values) == 0:
        return {k: None for k in ("n_l_min", "n_l_q1", "n_l_median", "n_l_q3",
                                  "n_l_max", "n_l_mean")}
    return {
        "n_l_min": float(np.min(values)), "n_l_q1": float(np.percentile(values, 25)),
        "n_l_median": float(np.median(values)), "n_l_q3": float(np.percentile(values, 75)),
        "n_l_max": float(np.max(values)), "n_l_mean": float(np.mean(values)),
    }


def leave_one_out(coords: np.ndarray, y: np.ndarray, radius: float, kappa: float,
                  tie_tolerance: float = 1e-12, sparse_cutoff: int = 3,
                  undefined_mcc: float = 0.0) -> LooOutcome:
    """Complete two-sided leave-one-out over the labelled cohort ``L`` at one radius."""
    coords = np.asarray(coords, dtype=np.float64)
    y = np.asarray(y, dtype=np.int64)
    n = len(y)
    n_p = int(y.sum())

    adj = pairwise_distances(coords) <= radius
    np.fill_diagonal(adj, False)                  # the hold-out: i never counts itself
    n_labeled = adj.sum(axis=1).astype(np.float64)
    n_plp = (adj.astype(np.float64) @ y.astype(np.float64))

    # ... and the background rate is held out too: residue i is removed from the
    # numerator when it is itself P/LP, and from the denominator always.
    pi_minus = (n_p - y) / (n - 1)

    f_hat = (n_plp + kappa * pi_minus) / (n_labeled + kappa)
    delta = f_hat - pi_minus
    tie = np.abs(delta) <= tie_tolerance
    zero_neighbour = n_labeled == 0
    sparse = n_labeled < sparse_cutoff

    predictions = (delta > 0) & ~tie
    predictions[zero_neighbour] = False           # FROZEN: zero-neighbour -> B/LB
    predictions = predictions.astype(np.int64)

    tp = int(((predictions == 1) & (y == 1)).sum())
    fp = int(((predictions == 1) & (y == 0)).sum())
    tn = int(((predictions == 0) & (y == 0)).sum())
    fn = int(((predictions == 0) & (y == 1)).sum())
    metrics = classification_metrics(tp, fp, tn, fn, undefined_mcc)

    outcome = LooOutcome(
        radius_A=float(radius), kappa=float(kappa), mcc=float(metrics["mcc"]),
        metrics=metrics, predictions=predictions, n_labeled=n_labeled, n_plp=n_plp,
        pi_minus=pi_minus, f_hat=f_hat, tie=tie, zero_neighbour=zero_neighbour,
        sparse=sparse,
    )
    outcome.diagnostics = sparse_evidence_rows(outcome, y)
    return outcome


def sparse_evidence_rows(o: LooOutcome, y: np.ndarray) -> list[dict]:
    """Mandatory per-radius sparse-evidence reporting: overall and BY CLASS.

    Split by *true* class and by *predicted* class, because sparse evidence that is
    concentrated in one class is a different problem from sparse evidence spread
    evenly. No residue is excluded on the basis of any of these counts.
    """
    scopes = {
        "overall": np.ones(len(y), dtype=bool),
        "true_PLP": y == 1,
        "true_BLB": y == 0,
        "pred_PLP": o.predictions == 1,
        "pred_BLB": o.predictions == 0,
    }
    rows = []
    for scope, mask in scopes.items():
        k = int(mask.sum())
        row = {
            "radius_A": o.radius_A, "scope": scope, "n_residues": k,
            "n_zero_neighbour": int(o.zero_neighbour[mask].sum()),
            "prop_zero_neighbour": float(o.zero_neighbour[mask].mean()) if k else 0.0,
            "n_sparse": int(o.sparse[mask].sum()),
            "prop_sparse": float(o.sparse[mask].mean()) if k else 0.0,
            "n_tie": int(o.tie[mask].sum()),
            "prop_tie": float(o.tie[mask].mean()) if k else 0.0,
        }
        row.update(_distribution(o.n_labeled[mask]))
        rows.append(row)
    return rows
