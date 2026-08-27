"""Regression pin for a real defect found running Stage E against real TP53RK
data for the first time (2026-08-19, run 20260819T082215Z_4ec9d6a8_ad8f32d6).

Stage D (footprint-robustness) legitimately terminates 09_ROBUSTNESS as
COMPLETED_NEGATIVE / NOT_APPLICABLE with reason ROBUSTNESS_NOT_EVALUABLE
whenever n_S < config's robustness.min_n_S_for_evaluation (fewer than 2
significant hotspot centers -> no non-empty proper subset exists to
perturb). handoff_04.json correctly carries this as
``robustness_profile.evaluable = False`` plus a human-readable reason, and
``center_sensitivity.tsv`` correctly has zero data rows (no per-center
measurement was ever computed, because none could be).

Stage E's ``_per_center_influence`` (annotation/stage.py) only ever consults
``data.center_sensitivity``. When that list is empty -- which happens both for
a benign "no data survived for this one residue" case AND for the structural
"robustness was never evaluable" case -- it collapses both to the literal
string "NA" for ``per_center_influence_summary`` and ``center_recurrence_rate``.
``assert_mandatory_context`` (annotation/schema.py) then rejects the row: its
own docstring establishes that "NA" is legal for a measurement genuinely
unavailable, but is explicitly disallowed for these two fields, precisely so a
genuinely-not-evaluated value can never look indistinguishable from a silently
dropped one. Net effect: any gene whose Stage B run finds fewer than 2
significant hotspot centers can never complete Stage E, even though every
upstream stage terminated in an internally consistent, designed state
(TP53RK: BlockedError: SCHEMA violation in hotspot_annotation.tsv: hotspot
'H1' has empty/NA mandatory context ['per_center_influence_summary',
'center_recurrence_rate']).

This is pure cross-stage schema wiring -- it never touched N_CAP, K_MAX, B,
r_hot/r_fp selection, any FDR threshold or any cohort definition.
"""
from __future__ import annotations

from types import SimpleNamespace

import pytest

from hotspot3d.annotation.schema import MANDATORY_ROW_CONTEXT, assert_mandatory_context
from hotspot3d.annotation.stage import _per_center_influence

pytestmark = pytest.mark.unit


def _not_evaluable_data(centers: dict) -> SimpleNamespace:
    """The exact shape handoff_04 produced for TP53RK: 1 significant center,
    zero robustness iterations possible, zero rows in center_sensitivity.tsv."""
    return SimpleNamespace(
        centers=centers,
        center_sensitivity=[],
        robustness_profile={
            "result_type": "geometric_footprint_robustness_profile",
            "evaluable": False,
            "reason": (
                "ROBUSTNESS_NOT_EVALUABLE: n_S = 1 significant hotspot center(s). "
                "Perturbation removes geometric centers from S; with fewer than 2 "
                "centers no non-empty proper subset exists, so no iteration is "
                "possible. No surrogate analysis is substituted."
            ),
            "n_S": 1,
        },
    )


def test_per_center_influence_does_not_silently_collapse_to_na_when_not_evaluable():
    """The originating defect: a single significant center with zero robustness
    iterations produced the literal string "NA" for both fields, indistinguishable
    from a silently dropped value."""
    data = _not_evaluable_data({"H1": [{"residue_index": 42}]})
    influence = _per_center_influence(data)
    assert influence["H1"]["influence"] != "NA", (
        "per_center_influence_summary must distinguish 'robustness was not "
        "evaluable upstream' (a designed Stage D terminal sub-state) from a "
        "silently dropped value"
    )
    assert influence["H1"]["recurrence"] != "NA", (
        "center_recurrence_rate must distinguish 'robustness was not evaluable "
        "upstream' (a designed Stage D terminal sub-state) from a silently "
        "dropped value"
    )


def test_a_not_evaluable_row_can_still_be_written():
    """The end-to-end failure mode: TP53RK run 20260819T082215Z_4ec9d6a8_ad8f32d6
    BLOCKED with `BlockedError: SCHEMA violation in hotspot_annotation.tsv:
    hotspot 'H1' has empty/NA mandatory context ['per_center_influence_summary',
    'center_recurrence_rate']`. Every upstream stage (01-09) was internally
    consistent and reported its own designed terminal state correctly; Stage E
    must still be able to write a schema-compliant row for it."""
    data = _not_evaluable_data({"H1": [{"residue_index": 42}]})
    influence = _per_center_influence(data)
    row = {key: "placeholder" for key in MANDATORY_ROW_CONTEXT}
    row["hotspot_id"] = "H1"
    row["center_recurrence_rate"] = influence["H1"]["recurrence"]
    row["per_center_influence_summary"] = influence["H1"]["influence"]
    assert_mandatory_context([row])  # must not raise BlockedError
