"""Canonical exception and outcome types.

The Output Contract (P4) requires that a *scientific negative* is never conflated
with a *technical failure*. These types make that distinction structural.
"""
from __future__ import annotations


class Hotspot3DError(Exception):
    """Base class for all pipeline errors."""


class BlockedError(Hotspot3DError):
    """A blocking precondition failed. Stage status BLOCKED / TECHNICAL_FAILURE.

    Raised for: hash mismatch, qc_status=FAIL upstream, missing config parameter,
    forbidden column loaded, empty class, unversioned source, code_version drift.
    Never raised for a scientifically empty result.
    """


class EscalationRequired(Hotspot3DError):
    """A decision belongs to the Lead. The agent stops and reports.

    Carries the ambiguity, the options, the consequence of each and a
    recommendation. The agent never resolves it itself.
    """

    def __init__(self, ambiguity: str, options: list[str], consequences: list[str],
                 recommendation: str):
        self.ambiguity = ambiguity
        self.options = options
        self.consequences = consequences
        self.recommendation = recommendation
        super().__init__(f"ESCALATION: {ambiguity} | recommend: {recommendation}")


class NegativeResult(Hotspot3DError):
    """A valid, complete scientific negative. Status COMPLETED_NEGATIVE.

    Examples: no admissible radius, S = empty at final r_hot, no admissible
    footprint, cohort too small. Terminates the chain honestly and licenses no
    method modification.
    """

    def __init__(self, condition: str, detail: str = "", **context):
        self.condition = condition
        self.detail = detail
        self.context = context
        super().__init__(f"NEGATIVE_RESULT[{condition}]: {detail}")


class LeakageError(Hotspot3DError):
    """An information barrier was violated. Always fatal, never recoverable.

    Raised when a stage loads a forbidden column, reads a downstream directory,
    or attempts to write outside its owned paths.
    """


class ContaminationEvent(LeakageError):
    """Information from a barred stage reached this stage. Stop and report."""
