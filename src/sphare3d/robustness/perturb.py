"""One perturbation iteration — fixed-radius footprint reconstruction only (II.11).

What is perturbed: **geometric center positions**. ``S_iter = S \\ T``.

What is NOT touched, in this module or anywhere in Phase D: ClinVar records, class
labels, the cohort ``L``, ``U_struct``, ``r_hot``, hotspot significance. Nothing in
here imports ``hotspot_stats``; Ripley's K, pair correlation, the radius scan,
``r_hot`` selection, permutation testing and FDR are not recomputed, not at stage
level and not inside an iteration. The withdrawn full-pipeline re-execution design
(A7/F11) is not implemented under this or any other name.

**[DECISION-STAGE-D-FIXED-RFP-0001, AUTHORIZED_BY=user]** Every iteration
reconstructs ``S_iter`` at the FROZEN Stage C ``r_fp`` — no per-iteration domain
re-derivation for selection purposes, no multi-scale sweep, no QC-driven
admissibility, no Pareto or distance-to-ideal re-selection. Radius selection is
Stage C's job, performed once on the full data with the multi-scale sweep, ARI
neighboring-radius stability and Pareto criteria; Phase D answers only "does the
SELECTED footprint survive loss of its supporting centers," never "would Stage C
have chosen a different radius here." This function calls
:func:`sphare3d.footprint.api.reconstruct_at` and never
:func:`sphare3d.footprint.api.build_footprint` on a perturbed center set — see
:mod:`sphare3d.robustness.stage`'s ``assert_no_radius_optimization_in_perturb``
for the mechanical pin, and ``docs/decisions/DECISION-STAGE-D-FIXED-RFP-0001.md``
for the full rationale. Prior to this decision, each iteration additionally ran
the full frozen methodology (its own domain re-derivation, sweep, QC, Pareto and
distance-to-ideal selection) to obtain a **primary** reconstruction at a
re-selected ``r_fp^iter``, alongside this **fixed-radius** reconstruction — see
the decision record for what that comparison showed and why it was removed.

**[v2 §8] "UNDERPOWERED iteration" investigated and found INAPPLICABLE under the
frozen design.** v2 §8's final [NEW] sentence reads: "Re-running the analysis
inside a validation iteration must re-run §5.4. If removing a subset of centers
causes the power certificate to fail, that iteration is recorded as
``UNDERPOWERED``... separately from genuine non-reproductions." §5.4's power
certificate (``p_floor`` vs ``c_1 = q/m``) is a property of the per-center
significance TEST — it needs ``N_P``, ``N_B``, ``m`` and the correction level,
all recomputed only when the primary hotspot detection (Ripley's K / positional or
label-permutation null / BH) is re-run. This module, and Phase D as a whole, never
re-runs that test: ``robustness.rerun_hotspot_pipeline: false`` (FROZEN, F11/A7)
and ``robustness.perturbation_universe: SIGNIFICANT_HOTSPOT_CENTERS`` (FROZEN,
A21) in ``config/pipeline.yaml`` fix Phase D to pure footprint-geometry
reconstruction on an already-frozen center subset — no p-value is ever computed
here, so there is no power certificate to pass or fail, and therefore nothing an
``UNDERPOWERED`` verdict could report. The antecedent action the [NEW] sentence
presupposes ("re-run §5.4 inside a validation iteration") is exactly the withdrawn
full-pipeline re-execution design that F11/A7 prohibits under any name; the
consequent instruction has no operational referent while that prohibition holds.
This is asserted mechanically, not just by convention: ``sphare3d.hotspot``
(which contains ``hotspot/power.py``, the §5.4 implementation) is in
:data:`FORBIDDEN_IMPORTS <sphare3d.robustness.stage.FORBIDDEN_IMPORTS>` and
checked by :func:`assert_no_hotspot_statistics_import
<sphare3d.robustness.stage.assert_no_hotspot_statistics_import>`.

This is NOT the same claim as "an existing FAILED_* status secretly means
UNDERPOWERED". ``STATUS_GEOMETRY_FAULT`` below is a GEOMETRIC failure (a
computational fault in the fixed-radius reconstruction) — a categorically
different thing from "the design could not have rejected any center". Reusing
it, or emitting the registered ``ROBUSTNESS_ITERATION_UNDERPOWERED`` warning
code, for a signal this design cannot compute would misrepresent what actually
failed and would require exactly the recomputation A7/F11 forbids. Reported to the
Lead for a ruling rather than resolved unilaterally; RULING (Lead): "genuinely
inapplicable, documented and pinned" is correct. No code path in this package
emits ``ROBUSTNESS_ITERATION_UNDERPOWERED`` — pinned by
``tests/robustness/test_robustness_v2_underpowered_investigation.py``. The finding
is additionally surfaced machine-readably, unconditionally, in every
``robustness_profile.json`` (and, through it, ``handoff_04.json``) as
``v2_section8_underpowered_iteration_applicable`` / ``_note`` — see
:data:`V2_SECTION8_UNDERPOWERED_NOTE <sphare3d.robustness.stage.V2_SECTION8_UNDERPOWERED_NOTE>`.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from ..footprint.api import (CenterSet, Universe, reconstruct_at)
from ..footprint.params import FootprintParams
from ..utils.errors import BlockedError
from . import classify
from .compare import GeometryComparison, compare

STATUS_OK = "OK"
STATUS_GEOMETRY_FAULT = "FAILED_GEOMETRY_FAULT"

FAILURE_COLUMNS = ["iteration_index", "k", "removed_center_residues", "status",
                   "outcome_type", "failure_reason"]


@dataclass
class IterationResult:
    index: int
    k: int
    removed: tuple[int, ...]
    enumeration_mode: str
    seed_context: str
    status: str
    outcome_type: str
    failure_reason: str
    robustness_r_fp: float | None
    wall_seconds: float
    row: dict = field(default_factory=dict)
    residues: frozenset = frozenset()          # the sole reconstruction's coverage

    @property
    def failed(self) -> bool:
        return self.status != STATUS_OK


CONTAINMENT_DIAGNOSTIC_FIELDS = ["robustness_r_contain",
                                 "robustness_r_fp_would_be",
                                 "robustness_r_contain_binding_residue"]


def containment_diagnostic(evidence_residues: tuple[int, ...], universe: Universe,
                           iter_centers: CenterSet, step_fp_A: float) -> dict:
    """What ``r_contain`` would have been on this iteration's center set.

    [DECISION-RFP-CONTAINMENT-0001] Purely diagnostic, and it must stay that way:
    DECISION-STAGE-D-FIXED-RFP-0001 holds the radius at Stage C's ``r_fp`` for every
    iteration, so this number is reported and never reconstructed at.

    The evidence set is the FROZEN one from Stage C — the same residues, re-measured
    against the surviving centers. Re-deriving membership here would need ``r_hot``
    and would conflate two different effects: residues leaving the evidence set, and
    the surviving ones getting further from a center. Only the second is what losing a
    center does to the radius, so only the second is measured.

    Removing centers can only increase a nearest-center distance, so
    ``r_contain_iter >= r_contain_original`` always. The gap is how much the radius
    would have grown — a direct read on how much of ``r_fp`` rests on one center.
    """
    import math

    if not evidence_residues or len(iter_centers) == 0:
        return {"robustness_r_contain": None, "robustness_r_fp_would_be": None,
                "robustness_r_contain_binding_residue": None}

    row_of = {int(rid): k for k, rid in enumerate(universe.ids)}
    rows = [row_of[int(r)] for r in evidence_residues if int(r) in row_of]
    if not rows:
        return {"robustness_r_contain": None, "robustness_r_fp_would_be": None,
                "robustness_r_contain_binding_residue": None}

    uc = np.asarray(universe.coords, dtype=np.float64)[rows]
    cc = np.atleast_2d(np.asarray(iter_centers.coords, dtype=np.float64))
    dmin = np.linalg.norm(uc[:, None, :] - cc[None, :, :], axis=-1).min(axis=1)
    worst = int(dmin.argmax())
    exact = float(dmin[worst])
    step = float(step_fp_A)
    return {
        "robustness_r_contain": round(exact, 6),
        "robustness_r_fp_would_be": round(math.ceil(exact / step - 1e-9) * step, 6),
        "robustness_r_contain_binding_residue": int(
            [r for r in evidence_residues if int(r) in row_of][worst]),
    }


def run_iteration(index: int, removed: tuple[int, ...], enumeration_mode: str,
                  seed_context: str, *, centers: CenterSet, universe: Universe,
                  params: FootprintParams, original_state, original_residues:
                  frozenset, r_fp_original: float, plp: frozenset, blb: frozenset,
                  r_hot_value: float, undefined_mcc: float,
                  evidence_residues: tuple[int, ...] = ()) -> IterationResult:
    t0 = time.perf_counter()
    iter_centers = centers.without(removed)

    row: dict = {
        "iteration_index": index,
        "k": len(removed),
        "removed_center_residues": list(removed),
        "n_removed": len(removed),
        "n_remaining": len(iter_centers),
        "enumeration_mode": enumeration_mode,
        "seed_context": seed_context,
    }

    # --- the center-loss family that does not depend on reconstruction -----
    row.update(classify.hotspot_coverage_family(
        iter_centers, universe, r_hot_value, plp, blb, undefined_mcc))

    # --- diagnostic only: the radius this iteration would have determined ---
    row.update(containment_diagnostic(evidence_residues, universe, iter_centers,
                                      params.step_fp_A))

    # --- the sole reconstruction: fixed at the frozen Stage C r_fp ----------
    status, reason = STATUS_OK, "NA"
    try:
        recon = reconstruct_at(iter_centers, universe, params, r_fp_original)
    except BlockedError as exc:
        status, reason = STATUS_GEOMETRY_FAULT, str(exc)
        row.update(GeometryComparison.failed())
        row.update(classify.failed_row())
        row["robustness_r_fp"] = None
        row["status"] = status
        row["outcome_type"] = "TECHNICAL_FAILURE"
        row["failure_reason"] = reason
        return IterationResult(
            index=index, k=len(removed), removed=removed,
            enumeration_mode=enumeration_mode, seed_context=seed_context,
            status=status, outcome_type="TECHNICAL_FAILURE", failure_reason=reason,
            robustness_r_fp=None, wall_seconds=time.perf_counter() - t0, row=row)

    # DECISION-STAGE-D-FIXED-RFP-0001's invariant, asserted at its origin: the
    # reconstruction was never asked to select anything, so this can only fail if
    # a future edit accidentally wires in a different radius at the call site.
    if float(recon.rho) != float(r_fp_original):
        raise BlockedError(
            f"BLOCKED — Phase D reconstructed iteration {index} at rho="
            f"{recon.rho}, not the frozen r_fp={r_fp_original}. Per "
            f"DECISION-STAGE-D-FIXED-RFP-0001 every iteration must reconstruct "
            f"at exactly the frozen Stage C radius."
        )

    iter_residues = frozenset(int(universe.ids[i])
                              for i, flag in enumerate(recon.state.covered) if flag)
    geometry = compare(original_state, original_residues, centers, r_fp_original,
                       recon.state, iter_residues, iter_centers,
                       r_fp_original, universe, params)
    row.update(geometry.as_dict())
    row.update(classify.metrics_for(iter_residues, plp, blb, undefined_mcc))
    row["robustness_r_fp"] = float(recon.rho)
    row["status"] = STATUS_OK
    row["outcome_type"] = "COMPLETED"
    row["failure_reason"] = "NA"

    return IterationResult(
        index=index, k=len(removed), removed=removed,
        enumeration_mode=enumeration_mode, seed_context=seed_context,
        status=STATUS_OK, outcome_type="COMPLETED", failure_reason="NA",
        robustness_r_fp=float(recon.rho), wall_seconds=time.perf_counter() - t0,
        row=row, residues=iter_residues)
