"""II.10 — footprint QC and admissibility.

Constraints, not objectives. An inadmissible candidate keeps **all** of its
computed metrics plus a rejection reason and stays in ``r_fp_scan.tsv``: the
decision trace has to show what was rejected and why, not just what won.

QC-F1's 0.50 threshold is never relaxed and its denominator is always
``|U_struct|``. QC-F4 and QC-F5 are assertions — a violation is a computational
fault, not a candidate property.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..utils.errors import BlockedError
from .domain import FootprintDomain
from .params import FootprintParams
from .sweep import RadiusState

QC_F1 = "QC_F1_EXCESSIVE_COVERAGE"
QC_F2 = "QC_F2_EXCESSIVE_BRIDGING"
QC_F3 = "QC_F3_GEOMETRY_UNAVAILABLE"


@dataclass(frozen=True)
class QCVerdict:
    rho: float
    admissible: bool
    reasons: tuple[str, ...]
    excessive_coverage: bool

    @property
    def status(self) -> str:
        return "PASS" if self.admissible else "FAIL"

    @property
    def reason_text(self) -> str:
        return ";".join(self.reasons) if self.reasons else "NA"


def assert_qc_f4(rho: float, domain: FootprintDomain) -> None:
    """QC-F4: the evaluated radius lies inside the derived domain."""
    if not (domain.rho_min - 1e-9 <= rho <= domain.rho_max + 1e-9):
        raise BlockedError(
            f"QC-F4 assertion failed: r_fp = {rho} is outside the derived domain "
            f"[{domain.rho_min}, {domain.rho_max}]. This is a computational fault."
        )


def assert_qc_f5(state: RadiusState, center_row_index: np.ndarray) -> None:
    """QC-F5: every significant center is contained in the footprint.

    Containment holds by construction (a center is at distance 0 from itself), so
    a failure here means the grid, the mask or the center mapping is wrong. It is
    escalated, never recorded as a fact about the protein.
    """
    if len(center_row_index) == 0:
        return
    missing = [int(i) for i in center_row_index if not bool(state.covered[i])]
    if missing:
        raise BlockedError(
            f"QC-F5 assertion failed at r_fp = {state.rho}: significant centers at "
            f"U_struct rows {missing[:10]} are not contained in the footprint. A "
            f"final footprint that excludes a significant center is a computational "
            f"fault and is escalated to the Lead."
        )
    if (state.center_labels == 0).any():
        raise BlockedError(
            f"QC-F5 assertion failed at r_fp = {state.rho}: a center's voxel is not "
            f"part of any occupancy component."
        )


def evaluate(state: RadiusState, params: FootprintParams,
             domain: FootprintDomain) -> QCVerdict:
    """Apply QC-F1..QC-F3 to one candidate radius (F4/F5 are assertions)."""
    assert_qc_f4(state.rho, domain)

    reasons: list[str] = []
    excessive = bool(state.coverage > params.qc_f1_max_coverage)
    if excessive:
        reasons.append(
            f"{QC_F1}:coverage={state.coverage:.6g}>"
            f"{params.qc_f1_max_coverage:g}")
    if not np.isfinite(state.bridging_index) or \
            state.bridging_index > params.qc_f2_max_bridging_index:
        reasons.append(
            f"{QC_F2}:bridging_index={state.bridging_index:.6g}>"
            f"{params.qc_f2_max_bridging_index:g}")
    if state.volume <= 0 or not state.mesh_ok:
        reasons.append(f"{QC_F3}:volume={state.volume:.6g};"
                       f"mesh={state.mesh_failure_reason}")

    return QCVerdict(rho=state.rho, admissible=not reasons,
                     reasons=tuple(reasons), excessive_coverage=excessive)


def coverage_constraint_failure(verdicts: list[QCVerdict],
                                params: FootprintParams) -> tuple[bool, float]:
    """``ALL_FOOTPRINT_CANDIDATES_EXCESSIVE_COVERAGE`` (II.10).

    Fires when all — or at least ``all_inadmissible_escalation_fraction`` — of the
    candidates fail QC-F1. The response is to report the coverage curve and stop.
    Widening the domain or weakening the 0.50 rule is prohibited.
    """
    if not verdicts:
        return False, 0.0
    fraction = sum(1 for v in verdicts if v.excessive_coverage) / len(verdicts)
    return bool(fraction >= params.all_inadmissible_escalation_fraction), float(fraction)
