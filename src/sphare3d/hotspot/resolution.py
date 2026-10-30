"""METHOD_SPEC II.4 — the permutation-resolution diagnostic.

A permutation p-value cannot go below ``p_res = 1/(B+1)``. When the test family is
large, the smallest BH critical value ``q/m`` can sit *below* that floor, and then
**no center can be significant no matter how strong the signal is**. An empty result
under those conditions is a statement about ``B``, not about the protein.

This diagnostic runs twice — **pre-flight** (condition C1 on ``m = |U_struct|``,
before any scan compute is spent) and **post-hoc** on the realized family — and it:

  * **never changes ``B``**;
  * **never changes the FDR method**;
  * never converts an empty result into a positive one.

It reports, and the Lead decides. ``B`` is never increased to chase significance;
a re-run at a larger ``B`` is a Lead decision under a NEW ``RUN_ID``.

    p_res   = 1 / (B + 1)
    R       = p_res * m / q     hypotheses that must sit at the floor together
                                before BH can reject any of them
    B_req   = ceil(s * m / (k_target * q)) - 1        s = 10 (safety factor, FROZEN)
    B_rec   = smallest ladder value >= B_req over {1e4, 1e5, 1e6, 1e7}

``k_target`` is the number of true centers one wants to remain detectable. It is read
from the frozen config (``permutation_resolution_diagnostic.k_target = 1``, the most
demanding choice) and recorded in the output. It affects only the *recommendation*;
it can never affect a p-value, a rejection or any published result.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# Fallback for direct library callers only. Stage B always passes the config value.
K_TARGET = 1

CONDITION_DEFINITIONS = {
    "C0": ("p_res > q — the p-value floor exceeds the FDR level, so no center can "
           "reach significance even before multiplicity correction."),
    "C1": ("R = p_res * m / q > 1 — the rank-1 BH critical value q/m lies below the "
           "p-value floor, so a single isolated center can never be rejected. "
           "Computable PRE-FLIGHT from m = |U_struct|."),
    "C2": ("0 < n_floor < R — centers sit exactly at the p-value floor but there are "
           "too few of them for BH to reject any."),
    "C3": ("n_floor > 0 and at least one floor center was NOT rejected — the "
           "rejection decision is bounded by the permutation resolution itself."),
}


@dataclass
class ResolutionDiagnostic:
    phase: str                       # preflight | posthoc
    B: int
    q: float
    m: int
    p_res: float
    R: float
    n_floor: int | None
    conditions: dict[str, bool | None]
    limited: bool
    k_target: int
    safety_factor_s: int
    B_req: int
    B_rec: int | None
    escalate_above_ladder: bool
    detail: dict = field(default_factory=dict)

    def as_json(self) -> dict:
        return {
            "phase": self.phase, "B": self.B, "q": self.q, "m_family_size": self.m,
            "p_res": self.p_res, "R": self.R, "n_floor": self.n_floor,
            "conditions": self.conditions,
            "condition_definitions": CONDITION_DEFINITIONS,
            "PERMUTATION_RESOLUTION_LIMITED": self.limited,
            "k_target": self.k_target,
            "k_target_source": "config:permutation_resolution_diagnostic.k_target",
            "k_target_note": ("Number of true centers required to remain detectable. "
                              "Affects only the B recommendation, never any p-value "
                              "or rejection."),
            "safety_factor_s": self.safety_factor_s,
            "B_req": self.B_req, "B_rec": self.B_rec,
            "B_req_exceeds_ladder": self.escalate_above_ladder,
            "B_was_changed": False,
            "fdr_method_was_changed": False,
            "note": ("This diagnostic never changes B and never changes the FDR "
                     "method. If the flag is set and the center set is empty, the "
                     "run is reported as RESOLUTION-LIMITED and escalated — it is "
                     "explicitly NOT evidence of absence."),
            **self.detail,
        }


def b_recommendation(m: int, q: float, s: int, ladder: list[int], max_b: int,
                     k_target: int = K_TARGET) -> tuple[int, int | None, bool]:
    """``(B_req, B_rec, exceeds_ladder)`` — a recommendation only, never applied."""
    b_req = int(np.ceil(s * m / (k_target * q))) - 1
    if b_req > max_b:
        return b_req, None, True
    candidates = [b for b in sorted(ladder) if b >= b_req]
    return b_req, (candidates[0] if candidates else None), False


def evaluate(phase: str, B: int, q: float, m: int, p_emp: np.ndarray | None,
             rejected: np.ndarray | None, s: int, ladder: list[int],
             max_b: int, k_target: int = K_TARGET) -> ResolutionDiagnostic:
    """Evaluate C0-C3. ``p_emp``/``rejected`` are ``None`` in the pre-flight phase."""
    p_res = 1.0 / (B + 1)
    R = p_res * m / q if q > 0 else float("inf")

    c0 = bool(p_res > q)
    c1 = bool(R > 1.0)
    if p_emp is None:
        n_floor, c2, c3 = None, None, None
    else:
        at_floor = np.isclose(np.asarray(p_emp, dtype=np.float64), p_res, rtol=0, atol=1e-15)
        n_floor = int(at_floor.sum())
        c2 = bool(0 < n_floor < R)
        c3 = bool(n_floor > 0 and rejected is not None and
                  not np.all(np.asarray(rejected, dtype=bool)[at_floor]))

    limited = bool(c0 or c1 or (c2 is True) or (c3 is True))
    b_req, b_rec, exceeds = b_recommendation(m, q, s, ladder, max_b, k_target)

    return ResolutionDiagnostic(
        phase=phase, B=int(B), q=float(q), m=int(m), p_res=float(p_res), R=float(R),
        n_floor=n_floor, conditions={"C0": c0, "C1": c1, "C2": c2, "C3": c3},
        limited=limited, k_target=int(k_target), safety_factor_s=int(s), B_req=b_req,
        B_rec=b_rec, escalate_above_ladder=exceeds,
        detail={"ladder": sorted(ladder), "max_recommendable_B": max_b},
    )
