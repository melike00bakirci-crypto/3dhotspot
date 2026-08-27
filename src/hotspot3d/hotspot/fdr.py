"""Multiple-testing control — Benjamini-Hochberg at ``q = 0.05`` (FROZEN, A10).

BH is the **primary and only decision procedure**. Benjamini-Yekutieli is computed
and reported *alongside* it for comparison and is never a substitute; a material
disagreement between the two is escalated, never silently resolved.

Boundary ties are ALL rejected: rejection is expressed as ``p <= p_(k)``, so any
center sharing the boundary p-value is rejected with it. This is required because
permutation p-values are discrete and ties at the boundary are common; rejecting an
arbitrary subset of tied hypotheses would make the result depend on sort order.

``q``, ``B`` and the QC thresholds are fixed before any p-value is seen. ``B`` is
never increased to chase significance and the method is never switched to a laxer
correction — see the permutation-resolution diagnostic (II.4) for the honest
response to a resolution-limited family.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class FDRResult:
    p: np.ndarray
    q_bh: np.ndarray
    q_by: np.ndarray
    reject_bh: np.ndarray
    reject_by: np.ndarray
    boundary_p: float | None
    boundary_rank: int | None
    m: int
    q: float
    by_constant: float

    @property
    def n_reject_bh(self) -> int:
        return int(self.reject_bh.sum())

    @property
    def n_reject_by(self) -> int:
        return int(self.reject_by.sum())


def bh_adjust(p: np.ndarray) -> np.ndarray:
    """Monotone BH-adjusted p-values ``q_(i) = min_{j>=i} (m/j) p_(j)``, capped at 1."""
    p = np.asarray(p, dtype=np.float64)
    m = len(p)
    if m == 0:
        return p.copy()
    order = np.argsort(p, kind="stable")
    ranks = np.arange(1, m + 1, dtype=np.float64)
    scaled = m * p[order] / ranks
    adjusted = np.minimum.accumulate(scaled[::-1])[::-1]
    out = np.empty(m, dtype=np.float64)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def by_constant(m: int) -> float:
    """``c(m) = sum_{i=1..m} 1/i`` — the Benjamini-Yekutieli dependency penalty."""
    return float(np.sum(1.0 / np.arange(1, m + 1))) if m > 0 else 1.0


def benjamini_hochberg(p: np.ndarray, q: float) -> FDRResult:
    """BH step-up at level ``q`` with all boundary ties rejected, plus BY alongside."""
    p = np.asarray(p, dtype=np.float64)
    m = len(p)
    if m == 0:
        empty = np.zeros(0, dtype=np.float64)
        return FDRResult(p=empty, q_bh=empty.copy(), q_by=empty.copy(),
                         reject_bh=np.zeros(0, dtype=bool),
                         reject_by=np.zeros(0, dtype=bool), boundary_p=None,
                         boundary_rank=None, m=0, q=float(q), by_constant=1.0)

    order = np.argsort(p, kind="stable")
    ranks = np.arange(1, m + 1, dtype=np.float64)
    critical = ranks * q / m
    passing = np.nonzero(p[order] <= critical)[0]

    if len(passing):
        k = int(passing[-1]) + 1
        boundary = float(p[order][k - 1])
        # "<= boundary" rejects EVERY center sharing the boundary p-value (FROZEN).
        reject_bh = p <= boundary
    else:
        k, boundary = None, None
        reject_bh = np.zeros(m, dtype=bool)

    q_bh = bh_adjust(p)
    c_m = by_constant(m)
    q_by = np.minimum(q_bh * c_m, 1.0)
    reject_by = q_by <= q

    # The threshold rule and the adjusted-p rule are equivalent by construction;
    # asserting it here means the two published columns can never disagree.
    if not np.array_equal(reject_bh, q_bh <= q):
        raise AssertionError(
            "BH inconsistency: the p <= p_(k) rejection set differs from q_bh <= q. "
            "This must never happen and indicates a defect in the FDR implementation."
        )

    return FDRResult(p=p, q_bh=q_bh, q_by=q_by, reject_bh=reject_bh, reject_by=reject_by,
                     boundary_p=boundary, boundary_rank=k, m=m, q=float(q),
                     by_constant=c_m)
