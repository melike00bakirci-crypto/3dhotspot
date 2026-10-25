"""Stage C+D entry point: ``run_stage_cd``.

    Phase C  ->  FREEZE GATE  ->  Phase D

The order is structural. Phase D imports the *frozen* footprint API and can never
write back into it; the gate between them is a hash set that both phases verify.
Phase C is implemented here; Phase D lives in :mod:`hotspot3d.robustness.stage`
and is imported lazily at the gate so that the dependency runs one way only.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path

from ..orchestration.contracts import (QC_PASS, QC_PASS_WITH_WARNINGS, Handoff,
                                       assert_handoff_shape)
from ..utils.errors import BlockedError, EscalationRequired, NegativeResult
from ..utils.hashing import manifest_for
from ..utils.io import write_json, write_text
from ..utils.runctx import RunContext
from ..utils.status import (OutcomeType, Severity, StageStatus, Status,
                            WarningCollector, utc_now)
from . import export, figures, freeze as freeze_mod, inputs as inputs_mod
from . import report as report_mod
from .domain import derive_containment_domain
from .api import (FOOTPRINT_API_VERSION, FootprintSolution,
                  build_footprint_contained,
                  footprint_code_version)
from .inputs import UpstreamInputs
from .params import FootprintParams
from .staging import stage_files, write_stage_manifest

AGENT = "footprint-robustness"
STAGE_07 = "07_FOOTPRINT_RADIUS"
STAGE_08 = "08_FINAL_FOOTPRINT"
STAGE_09 = "09_ROBUSTNESS"

# [DECISION-RFP-CONTAINMENT-0001] The objective matrix, the Pareto files and the
# boundary diagnostic are gone with the search. evidence_set.tsv is the new decision
# trace — r_fp is the largest in-evidence distance in it — and sensitivity_band.tsv
# is the mitigation for r_fp inheriting r_hot's margin.
EXPECTED_07 = [
    ("evidence_set.tsv", "THE decision trace: every mapped P/LP, its nearest center, "
                         "the containment verdict, and which residue binds r_fp"),
    ("sensitivity_band.tsv", "the footprint at every radius from the resolution floor "
                             "to r_hot, with the selected row marked"),
    ("footprint_descriptors.tsv", "geometric descriptors of the determined radius"),
    ("footprint_domain.json", "containment record: r_contain, evidence counts, "
                              "binding residue, outlier fence, MST merge scales"),
    ("footprint_admissibility.tsv", "verbatim QC projection (P1)"),
    ("merge_events.tsv", "merge events with participating components and neck widths"),
    ("footprint_decision.json", "the determination record"),
    ("stage_status.json", "stage status (always present)"),
    (f"warnings_{STAGE_07}.tsv", "stage warnings"),
    ("stage_manifest.tsv", "this manifest"),
    ("figures/F7_footprint_multiscale_evolution.png", "F7 (300 dpi)"),
    ("figures/F7_footprint_multiscale_evolution.svg", "F7 (vector)"),
    ("figures/F8_footprint_pareto_selection.png", "F8 (300 dpi)"),
    ("figures/F8_footprint_pareto_selection.svg", "F8 (vector)"),
]

EXPECTED_08 = [
    ("footprint_residues.tsv", "FP_original residues with component membership"),
    ("footprint_geometry.tsv", "FP_original geometry at the selected radius"),
    ("footprint_components.tsv", "per-component summary of FP_original"),
    ("footprint_occupancy.npz", "voxel occupancy at the final radius only (P5)"),
    ("footprint_surface.ply", "marching-cubes surface mesh"),
    ("structures/final_footprint.pdb", "CA-only footprint model (PDB)"),
    ("structures/final_footprint.cif", "CA-only footprint model (mmCIF)"),
    ("structures/view_footprint.cxc", "ChimeraX view script"),
    ("structures/view_footprint.pml", "PyMOL view script"),
    ("footprint_freeze.json", "FREEZE GATE hash record"),
    ("handoff_03.json", "Phase C -> Phase D handoff"),
    ("stage_c_report.md", "narrative Phase C report"),
    ("stage_status.json", "stage status (always present)"),
    (f"warnings_{STAGE_08}.tsv", "stage warnings"),
    ("stage_manifest.tsv", "this manifest"),
]


@dataclass
class PhaseCResult:
    """Everything Phase D is allowed to know about Phase C — and nothing more."""

    inputs: UpstreamInputs
    params: FootprintParams
    solution: FootprintSolution
    freeze_payload: dict
    code_version: str
    r_fp: float
    warnings_07: WarningCollector
    escalations: list[dict] = field(default_factory=list)


# --------------------------------------------------------------------------
# entry point
# --------------------------------------------------------------------------

def run_stage_cd(ctx: RunContext, *, upstream: Handoff) -> tuple[Handoff, Handoff]:
    """Run Phase C, seal the freeze gate, then run Phase D. Returns (h3, h4)."""
    handoff_03, phase_c = run_phase_c(ctx, upstream=upstream)

    if phase_c is None:                       # scientific negative — chain ends here
        handoff_04 = _negative_handoff_04(ctx, handoff_03)
        return handoff_03, handoff_04

    from ..robustness.stage import run_phase_d      # deferred: one-way dependency
    handoff_04 = run_phase_d(ctx, phase_c=phase_c, handoff_03=handoff_03)
    return handoff_03, handoff_04


# --------------------------------------------------------------------------
# Phase C
# --------------------------------------------------------------------------

def run_phase_c(ctx: RunContext, *, upstream: Handoff
                ) -> tuple[Handoff, PhaseCResult | None]:
    started = utc_now()
    t0 = time.perf_counter()
    dir_07 = ctx.stage_dir(STAGE_07)
    dir_08 = ctx.stage_dir(STAGE_08)
    ctx.stage_dir(STAGE_09)                   # created now so it is never mysterious
    warn07 = WarningCollector(STAGE_07, AGENT)
    warn08 = WarningCollector(STAGE_08, AGENT)
    escalations: list[dict] = []

    # -- 1. validate ------------------------------------------------------
    inputs_mod.assert_stage_dirs_empty(ctx, (STAGE_07, STAGE_08, STAGE_09))
    params = FootprintParams.from_config(ctx.config)
    up = inputs_mod.load(ctx, upstream)
    code_ver = footprint_code_version()

    if ctx.synthetic:
        warn07.add("SYNTHETIC_INPUT_MODE",
                   "Stage C ran on synthetic inputs; results are not biological.",
                   status="ACCEPTED_BY_DESIGN")

    # -- 2-6. build the solution ------------------------------------------
    try:
        solution = build_footprint_contained(up.centers, up.universe, params,
                                             up.evidence)
    except NegativeResult as neg:
        return _finish_negative(ctx, dir_07, dir_08, warn07, warn08, up, params,
                                neg, started, t0), None

    # -- 7. write 07_* (always, including for a negative outcome) ----------
    scan_rows = export.build_scan_rows(solution)
    export.write_scan(dir_07, scan_rows)
    export.write_admissibility(dir_07, scan_rows)
    export.write_merge_events(dir_07, solution, up.centers.ids)
    export.write_component_memberships(dir_07, solution, up.universe, up.centers.ids)
    h_sensitivity = (export.voxel_h_sensitivity(solution, up.centers.coords, params)
                     if solution.r_fp is not None else None)
    export.write_domain(dir_07, solution, up.universe, up.centers.ids, params,
                        h_sensitivity=h_sensitivity)
    # [DECISION-RFP-CONTAINMENT-0001] No objective matrix, Pareto front or boundary
    # diagnostic: with one candidate there is no domination, no utopia distance and
    # no domain edge. The containment record and the evidence set replace them.
    export.write_evidence_set(dir_07, up.evidence, solution)
    export.write_sensitivity_band(dir_07, up, params, solution)

    # [DECISION-RFP-CONTAINMENT-0001] The objective/constraint separation diagnostic
    # is gone with the objectives: a one-candidate set cannot have a degenerate
    # objective change a distance, because no distance decides anything.
    _emit_containment_warnings(warn07, solution, up, params, escalations,
                               h_sensitivity)

    figures.figure_f7(dir_07 / "figures", solution, params)
    figures.figure_f8(dir_07 / "figures", solution, params)
    export.write_decision(dir_07, solution, params, escalations)

    # -- escalation: coverage constraint failure --------------------------
    if solution.coverage_constraint_failure:
        return _finish_coverage_escalation(ctx, dir_07, dir_08, warn07, warn08,
                                           solution, params, up, started, t0)

    # -- negative: no admissible candidate --------------------------------
    if solution.r_fp is None:
        neg = NegativeResult(
            "NO_ADMISSIBLE_FOOTPRINT",
            "The containment radius failed QC-F1/F2/F3. There is no alternative "
            "candidate to fall back to: DECISION-RFP-CONTAINMENT-0001 determines a "
            "single radius and never shrinks it to pass a gate. The radius and its "
            "metrics are retained in footprint_descriptors.tsv.",
            n_candidates=len(solution.sweep.rows),
            excessive_coverage_fraction=solution.excessive_coverage_fraction,
        )
        return _finish_negative(ctx, dir_07, dir_08, warn07, warn08, up, params,
                                neg, started, t0), None

    # -- 8. write 08_* ----------------------------------------------------
    export.write_final_residues(dir_08, solution, up.universe, up.centers.ids,
                                up.cohort_plp, up.cohort_blb)
    export.write_final_geometry(dir_08, solution)
    export.write_final_components(dir_08, solution, up.universe, up.centers.ids)
    export.write_occupancy(dir_08, solution)
    _, mesh_reason = export.write_surface_ply(dir_08, solution)
    if mesh_reason != "NA":
        warn08.add("RENDERING_UNAVAILABLE",
                   f"surface mesh could not be extracted at r_fp: {mesh_reason}",
                   affected_output="footprint_surface.ply")
    export.write_structures(dir_08 / "structures", solution, up.universe,
                            up.centers.ids, up.residue_names, ctx.gene)
    report_mod.write_stage_c_report(
        dir_08 / "stage_c_report.md", solution, up.universe, params,
        up.centers.ids, warn07.items + warn08.items, escalations, code_ver,
        ctx.run_id, ctx.gene)

    # -- 9. FREEZE GATE ---------------------------------------------------
    freeze = freeze_mod.write_freeze(
        ctx, solution.r_fp, code_ver, FOOTPRINT_API_VERSION,
        extra={"center_set_sha256": up.centers.digest(),
               "upstream_sha256": up.hashes})

    # -- 10. handoff, statuses, manifests, provenance ---------------------
    handoff = _build_handoff_03(ctx, solution, up, params, code_ver, warn07, warn08)
    _write_status(dir_07, STAGE_07, Status.COMPLETED, OutcomeType.COMPLETED,
                  "footprint radius selected", started, t0, warn07, upstream,
                  len(EXPECTED_07))
    _write_status(dir_08, STAGE_08, Status.COMPLETED, OutcomeType.COMPLETED,
                  "final footprint exported and frozen", started, t0, warn08,
                  upstream, len(EXPECTED_08))
    warn07.write(dir_07)
    warn08.write(dir_08)
    write_stage_manifest(dir_07, EXPECTED_07,
                         extra_created=stage_files(dir_07))
    write_stage_manifest(dir_08, EXPECTED_08, extra_created=stage_files(dir_08))
    _record_provenance(ctx, STAGE_07, up, params, solution, code_ver, warn07)
    _record_provenance(ctx, STAGE_08, up, params, solution, code_ver, warn08)

    return handoff, PhaseCResult(
        inputs=up, params=params, solution=solution,
        freeze_payload=freeze["payload"], code_version=code_ver,
        r_fp=float(solution.r_fp), warnings_07=warn07, escalations=escalations,
    )


# --------------------------------------------------------------------------
# warnings and escalations
# --------------------------------------------------------------------------

# [DECISION-RFP-CONTAINMENT-0001] _emit_selection_warnings and
# _emit_objective_constraint_warnings are gone. They warned about NEAR_TIE_RADIUS_SELECTION,
# OBJECTIVE_REDUNDANCY and BOUNDARY_OPTIMUM_WARNING — all properties of a candidate set.
# Two of their warnings were NOT about selection and were kept rather than dropped with
# them: VOXEL_RESOLUTION_FALLBACK moved into _emit_containment_warnings, and the per-row
# EXCESSIVE_FOOTPRINT_COVERAGE is subsumed by FOOTPRINT_COVERAGE_CONSTRAINT_FAILURE,
# which fires at BLOCKING severity because with one candidate "this radius exceeds
# coverage" and "every candidate exceeds coverage" are the same event.
def _finish_coverage_escalation(ctx, dir_07: Path, dir_08: Path,
                                warn07: WarningCollector, warn08: WarningCollector,
                                solution: FootprintSolution, params: FootprintParams,
                                up: UpstreamInputs, started: str, t0: float):
    fraction = solution.excessive_coverage_fraction
    warn07.add("FOOTPRINT_COVERAGE_CONSTRAINT_FAILURE",
               f"{fraction:.1%} of candidate radii exceed the frozen QC-F1 coverage "
               f"limit of {params.qc_f1_max_coverage:g} "
               f"(escalation threshold {params.all_inadmissible_escalation_fraction:.0%})",
               potential_consequence="no footprint can be selected without weakening a "
                                     "frozen rule",
               recommended_action="Lead decision. The domain is NOT widened and the "
                                  "0.50 rule is NOT relaxed.",
               affected_output="footprint_descriptors.tsv")
    reason = (
        f"ALL_FOOTPRINT_CANDIDATES_EXCESSIVE_COVERAGE = TRUE: {fraction:.1%} of the "
        f"{len(solution.sweep.rows)} candidate radii in "
        f"[{solution.domain.rho_min:g}, {solution.domain.rho_max:g}] A cover more than "
        f"{params.qc_f1_max_coverage:g} of U_struct. The full coverage curve is in "
        f"footprint_descriptors.tsv."
    )
    _write_status(dir_07, STAGE_07, Status.BLOCKED, OutcomeType.SCIENTIFIC_NEGATIVE,
                  reason, started, t0, warn07, None, len(EXPECTED_07),
                  recommended_action="Lead decision required; see stage report.")
    _write_status(dir_08, STAGE_08, Status.NOT_RUN, OutcomeType.NOT_APPLICABLE,
                  "Phase C stopped at the QC-F1 escalation; no final footprint exists.",
                  started, t0, warn08, None, len(EXPECTED_08))
    _write_status(ctx.stage_dir(STAGE_09), STAGE_09, Status.NOT_RUN,
                  OutcomeType.NOT_APPLICABLE,
                  "Phase D not evaluable: Phase C produced no FP_original "
                  "(FOOTPRINT_COVERAGE_CONSTRAINT_FAILURE).",
                  started, t0, WarningCollector(STAGE_09, AGENT), None, 0)
    warn07.write(dir_07)
    warn08.write(dir_08)
    write_stage_manifest(dir_07, EXPECTED_07, extra_created=stage_files(dir_07))
    write_stage_manifest(dir_08, EXPECTED_08, extra_created=stage_files(dir_08))
    raise EscalationRequired(
        ambiguity=reason,
        options=[
            "STOP and report the negative result as it stands",
            "Lead-authorized re-run under a NEW RUN_ID with a pre-registered, "
            "explicitly different center set or domain",
        ],
        consequences=[
            "no footprint claim is made; the coverage curve is fully documented",
            "any change to the domain or the 0.50 rule made here would be domain "
            "shopping and would invalidate the pre-registration",
        ],
        recommendation="STOP. Report the coverage curve across the whole domain. "
                       "Do not widen the domain and do not relax QC-F1.",
    )


def _finish_negative(ctx: RunContext, dir_07: Path, dir_08: Path,
                     warn07: WarningCollector, warn08: WarningCollector,
                     up: UpstreamInputs, params: FootprintParams,
                     neg: NegativeResult, started: str, t0: float) -> Handoff:
    """A valid, complete, terminating scientific negative (P4)."""
    # [DECISION-REPORTING-CONTEXT-0001] The orchestrator's declared-negative hook keys
    # on this flag. Stage C does not raise NegativeResult to the orchestrator — it
    # records the negative, writes NOT_RUN for 08/09 and returns — so without the flag
    # the hook skips the run, Stage E is invoked, blocks on the absent
    # footprint_residues.tsv, and the run is typed BLOCKED/TECHNICAL_FAILURE with
    # negative_result: null. Found live in the v3 batch: CACNA1S and CASQ2 both ended
    # 07/08 = COMPLETED_NEGATIVE and then came back BLOCKED. Third instance of one
    # seam — Stage B always set it, Stage A was fixed earlier in this decision, and
    # Stage C is this.
    detail = {"condition": neg.condition, "detail": neg.detail,
              "settled_negative": True, **neg.context}
    if not (dir_07 / "footprint_domain.json").is_file():
        write_json(dir_07 / "footprint_domain.json",
                   {"negative_result": detail,
                    "note": "the II.9 domain could not be populated; see NOT_RUN.txt"})

    _write_status(dir_07, STAGE_07, Status.COMPLETED_NEGATIVE,
                  OutcomeType.SCIENTIFIC_NEGATIVE, neg.detail, started, t0, warn07,
                  None, len(EXPECTED_07), negative=detail)
    _write_status(dir_08, STAGE_08, Status.COMPLETED_NEGATIVE,
                  OutcomeType.SCIENTIFIC_NEGATIVE,
                  f"No final footprint exists: {neg.condition}.", started, t0, warn08,
                  None, len(EXPECTED_08), negative=detail)
    _write_status(ctx.stage_dir(STAGE_09), STAGE_09, Status.NOT_RUN,
                  OutcomeType.NOT_APPLICABLE,
                  f"Phase D not applicable: Phase C ended in the scientific negative "
                  f"{neg.condition}. No FP_original exists to perturb.",
                  started, t0, WarningCollector(STAGE_09, AGENT), None, 0)
    write_text(dir_08 / "stage_c_report.md", _negative_report(neg, up, params))
    warn07.write(dir_07)
    warn08.write(dir_08)
    write_stage_manifest(dir_07, EXPECTED_07, extra_created=stage_files(dir_07))
    write_stage_manifest(dir_08, EXPECTED_08, extra_created=stage_files(dir_08))

    payload = {
        "hotspot_radius": up.r_hot.value_for_coverage_query,
        "footprint_radius": None, "n_components": None, "volume_A3": None,
        "surface_area_A2": None, "n_footprint_residues": None, "coverage": None,
        "voxel_h_A": params.voxel_h_A, "voxel_h": params.voxel_h_A,
        "r_contain_A": None, "r_contain_exact_A": None, "grid_slack_A": None,
        "n_plp_mapped": int(up.evidence.n_plp_mapped),
        "n_evidence_set": len(up.evidence.evidence),
        "n_plp_outside_evidence": len(up.evidence.outside),
        "fraction_plp_outside_evidence": (
            len(up.evidence.outside) / up.evidence.n_plp_mapped
            if up.evidence.n_plp_mapped else None),
        "binding_residue": None, "binding_residue_is_outlier": None,
        "detection_scale_below_resolution": None, "coverage_escalated": None,
        "r_hot_runner_up_A": (up.r_hot_runner_up.value_for_coverage_query
                              if up.r_hot_runner_up is not None else None),
        "r_hot_runner_up_margin_rel": up.r_hot_runner_up_margin_rel,
        "r_contain_at_runner_up_A": None,
        "code_version": footprint_code_version(),
        "footprint_api_version": FOOTPRINT_API_VERSION,
    }
    handoff = Handoff(
        name="handoff_03", run_id=ctx.run_id, config_sha256=ctx.config.sha256,
        qc_status=QC_PASS_WITH_WARNINGS,
        manifest=manifest_for(stage_files(dir_07) + stage_files(dir_08), ctx.run_root),
        payload=payload, negative_result=detail,
    )
    assert_handoff_shape(handoff)
    write_json(ctx.handoff_path("handoff_03"), handoff.as_dict())
    return handoff


def _negative_report(neg: NegativeResult, up: UpstreamInputs,
                     params: FootprintParams) -> str:
    return "\n".join([
        "# Stage C — NEGATIVE RESULT",
        "",
        f"**Condition:** `{neg.condition}`",
        "",
        neg.detail,
        "",
        "## How to read this",
        "",
        "This is a TRUE SCIENTIFIC NEGATIVE, not a technical failure. The footprint "
        "construction ran correctly and the result is that no admissible footprint "
        "exists under the frozen methodology.",
        "",
        "No threshold was relaxed, no domain was widened, no center was dropped and "
        "the voxel grid was not coarsened to produce a positive result.",
        "",
        f"- |S| = {len(up.centers)} significant hotspot centers (immutable input)",
        f"- |U_struct| = {len(up.universe)}",
        f"- frozen QC-F1 coverage limit: {params.qc_f1_max_coverage:g}",
        f"- context: {neg.context}",
        "",
        "`09_ROBUSTNESS/` exists with its own `stage_status.json` recording "
        "NOT_RUN / NOT_APPLICABLE and naming this stage as the cause — no "
        "mysteriously empty folder is ever left behind.",
        "",
    ])


def _negative_handoff_04(ctx: RunContext, handoff_03: Handoff) -> Handoff:
    """handoff_04 for the case where Phase C terminated the chain."""
    dir_09 = ctx.stage_dir(STAGE_09)
    empty = {"n": 0, "mean": None, "sd": None, "median": None, "q1": None,
             "q3": None, "iqr_lo": None, "iqr_hi": None, "iqr": None,
             "min": None, "max": None}
    profile = {
        "result_type": "geometric_footprint_robustness_profile",
        "evaluable": False,
        "reason": "Phase C produced no FP_original; there is no footprint to perturb.",
        "upstream_condition": (handoff_03.negative_result or {}).get("condition", "NA"),
        "categorical_verdict_emitted": False,
        # Contract keys present but null, so a consumer reads "not evaluable"
        # rather than tripping over a missing key or reading 0.0 as a low score.
        "jaccard_residues": dict(empty),
        "dice_residues": dict(empty),
        "preservation_freq_050": None,
        "preservation_freq_070": None,
        "iterations": {"planned": 0, "completed": 0, "failed": 0},
        "clinical_dataset_unmodified": True,
        "geometric_footprint_robustness_only": True,
    }
    payload = {
        "footprint_code_version": footprint_code_version(),
        "robustness_profile": profile, "n_S": None, "K_MAX": None,
        "total_subset_space": 0, "total_evaluated": 0, "per_level": [],
        "robustness_r_fp_invariant": {},
        "clinical_dataset_unmodified": True,
        "geometric_footprint_robustness_only": True,
    }
    handoff = Handoff(
        name="handoff_04", run_id=ctx.run_id, config_sha256=ctx.config.sha256,
        qc_status=QC_PASS_WITH_WARNINGS,
        manifest=manifest_for(stage_files(dir_09), ctx.run_root),
        payload=payload, negative_result=handoff_03.negative_result,
    )
    assert_handoff_shape(handoff)
    write_json(ctx.handoff_path("handoff_04"), handoff.as_dict())
    return handoff


# --------------------------------------------------------------------------
# handoff / status / provenance helpers
# --------------------------------------------------------------------------

def _build_handoff_03(ctx: RunContext, solution: FootprintSolution,
                      up: UpstreamInputs, params: FootprintParams, code_ver: str,
                      warn07: WarningCollector, warn08: WarningCollector) -> Handoff:
    state = solution.selected
    counts = warn07.counts()
    qc_status = QC_PASS
    if counts[Severity.MAJOR.value] or warn08.counts()[Severity.MAJOR.value] \
            or counts[Severity.ADVISORY.value]:
        qc_status = QC_PASS_WITH_WARNINGS

    payload = {
        # passthrough, unmodified — recorded, never used in the r_fp decision
        "hotspot_radius": up.r_hot.value_for_coverage_query,
        "footprint_radius": float(solution.r_fp),
        "n_components": int(state.n_components),
        "volume_A3": float(state.volume),
        "surface_area_A2": float(state.surface_area),
        "n_footprint_residues": int(state.n_covered_residues),
        "coverage": float(state.coverage),
        "voxel_h": float(solution.spec.h),
        "voxel_h_A": float(solution.spec.h),
        "voxel_resolution_fallback_applied": bool(solution.spec.fallback_applied),
        "code_version": code_ver,
        "footprint_api_version": FOOTPRINT_API_VERSION,
        "n_significant_centers": len(up.centers),
        "center_set_sha256": up.centers.digest(),
        # [DECISION-RFP-CONTAINMENT-0001] the containment record
        **{k: solution.containment[k] for k in (
            "r_contain_A", "r_contain_exact_A", "grid_slack_A",
            "n_plp_mapped", "n_evidence_set", "n_plp_outside_evidence",
            "fraction_plp_outside_evidence", "binding_residue",
            "binding_residue_is_outlier", "detection_scale_below_resolution",
            "coverage_escalated")},
        "r_hot_runner_up_A": (up.r_hot_runner_up.value_for_coverage_query
                              if up.r_hot_runner_up is not None else None),
        "r_hot_runner_up_margin_rel": up.r_hot_runner_up_margin_rel,
        "r_contain_at_runner_up_A": _r_contain_at_runner_up(up, params),
    }
    handoff = Handoff(
        name="handoff_03", run_id=ctx.run_id, config_sha256=ctx.config.sha256,
        qc_status=qc_status,
        manifest=manifest_for(
            stage_files(ctx.full_results / STAGE_07)
            + stage_files(ctx.full_results / STAGE_08), ctx.run_root),
        payload=payload,
    )
    assert_handoff_shape(handoff)
    write_json(ctx.handoff_path("handoff_03"), handoff.as_dict())
    return handoff


def _write_status(stage_dir: Path, stage: str, status: Status,
                  outcome: OutcomeType, reason: str, started: str, t0: float,
                  warn: WarningCollector, upstream: Handoff | None,
                  n_expected: int, negative: dict | None = None,
                  recommended_action: str = "") -> None:
    created = sum(1 for _ in stage_files(stage_dir))
    StageStatus(
        stage=stage, agent_owner=AGENT, status=status, outcome_type=outcome,
        reason=reason, started_utc=started, ended_utc=utc_now(),
        wall_seconds=time.perf_counter() - t0,
        upstream_handoff="handoff_02" if upstream is not None else "NA",
        upstream_qc_status=upstream.qc_status if upstream is not None else "NA",
        n_outputs_expected=n_expected, n_outputs_created=created,
        warnings=warn.counts(), recommended_action=recommended_action,
        negative_result=negative,
    ).write(stage_dir)


def _record_provenance(ctx: RunContext, stage: str, up: UpstreamInputs,
                       params: FootprintParams, solution: FootprintSolution,
                       code_ver: str, warn: WarningCollector) -> None:
    from ..utils.runctx import software_versions

    ctx.record_provenance(
        stage, AGENT,
        inputs={"upstream_files_sha256": up.hashes,
                "upstream_paths": {k: str(v) for k, v in up.paths.items()},
                "center_set_sha256": up.centers.digest(),
                "n_significant_centers": len(up.centers),
                "n_universe_residues": len(up.universe)},
        outputs={"selected_r_fp_A": solution.r_fp,
                 "n_candidate_radii": len(solution.sweep.rows),
                 "n_admissible": solution.n_admissible},
        parameters={**params.as_dict(),
                    "footprint_code_version": code_ver,
                    "footprint_api_version": FOOTPRINT_API_VERSION,
                    "r_hot_recorded_readonly": up.r_hot.value_for_coverage_query,
                    "r_hot_used_in_r_fp_selection": False,
                    "grid_spec": solution.spec.as_dict(),
                    "edt_calls_for_sweep": 1,
                    # V and A are computed on DIFFERENT representations by design.
                    # Marching cubes on a binary mask carries a ~8% staircase bias
                    # that never converges, which would make compactness of a
                    # sphere 0.92 at every resolution; taking the isosurface of the
                    # same single distance field makes Psi = 1 for a sphere true
                    # and convergent (A18 freezes marching cubes, not the array it
                    # is handed).
                    "volume_basis": "voxel count x h^3 on the occupancy grid",
                    "surface_area_basis": (
                        "marching-cubes isosurface at level rho of the SAME single "
                        "Euclidean distance transform that defines the occupancy"),
                    "compactness_definition": "Psi = (36 pi V^2)^(1/3) / A",
                    "compactness_of_a_sphere": (
                        "1.00 +/- 0.02 at h = 0.5 A, verified by test"),
                    "software_versions": software_versions()},
        commands=[f"hotspot3d.footprint.stage.run_phase_c(run_id={ctx.run_id})"],
        warnings=[w.as_row() for w in warn.items],
    )


def _r_contain_at_runner_up(up: UpstreamInputs, params: FootprintParams):
    """``r_contain`` recomputed against the runner-up ``r_hot``. Reporting only.

    [DECISION-RFP-CONTAINMENT-0001, weakness 1] ``r_fp`` is a function of ``r_hot``,
    and ``r_hot`` is still chosen by the multi-objective machinery this decision
    rejects for ``r_fp``. Measured on KCNA2, whose ``r_hot`` has flipped between 9.5
    and 10.0 across runs at a 0.18-1.58 % margin, that propagates to 117 vs 143
    footprint residues and 2 vs 1 components. Reporting both radii makes the
    fragility visible; it changes nothing.
    """
    if up.evidence_runner_up is None:
        return None
    try:
        _, record = derive_containment_domain(
            up.centers.coords, up.universe.coords, up.evidence_runner_up, params)
    except NegativeResult:
        return None                      # runner-up admits no evidence set
    return float(record["r_contain_A"])


def _emit_containment_warnings(warn: WarningCollector, solution: FootprintSolution,
                               up: UpstreamInputs, params: FootprintParams,
                               escalations: list[dict],
                               h_sensitivity: dict | None = None) -> None:
    """The reporting the rule owes, per DECISION-RFP-CONTAINMENT-0001.

    None of these change the radius. The decision is explicit that the coverage rule
    escalates rather than shrinking ``r_fp``, and that the resolution flag never
    raises it.
    """
    rc = solution.containment

    if rc["detection_scale_below_resolution"]:
        warn.add(
            "DETECTION_SCALE_BELOW_3D_RESOLUTION",
            f"r_contain = {rc['r_contain_exact_A']:.3f} A is below the "
            f"{params_detection_floor():g} A flag: the detection scale is at or "
            f"under what Ca-level geometry resolves. The radius is NOT raised — "
            f"raising it would report a volume the evidence does not support.",
            severity=Severity.ADVISORY,
            potential_consequence=("Footprint geometry at this scale is dominated by "
                                   "the Ca spacing of adjacent residues."),
            recommended_action="None. Report the flag beside the footprint.",
            affected_output="footprint_descriptors.tsv")

    if rc["binding_residue_is_outlier"]:
        warn.add(
            "CONTAINMENT_BINDING_RESIDUE_IS_OUTLIER",
            f"residue {rc['binding_residue']} sets r_contain at "
            f"{rc['r_contain_exact_A']:.3f} A while the evidence-set distances have "
            f"median {rc['evidence_distance_median_A']:.3f} A and an outlier fence "
            f"of {rc['outlier_fence_A']:.3f} A. The radius is extrapolated from the "
            f"bulk to one point.",
            severity=Severity.ADVISORY,
            potential_consequence=("Reclassifying that single residue would move the "
                                   "footprint boundary."),
            recommended_action="None. Report the binding residue and the fence.",
            affected_output="evidence_set.tsv")

    frac = rc["fraction_plp_outside_evidence"]
    if frac is not None and frac > 0.25:
        warn.add(
            "PLP_OUTSIDE_EVIDENCE_SET",
            f"{rc['n_plp_outside_evidence']} of {rc['n_plp_mapped']} mapped P/LP "
            f"residues ({100 * frac:.1f} %) lie outside the detection neighbourhood "
            f"of every significant center. No center-anchored footprint can contain "
            f"them; the radius is not grown to reach them.",
            severity=Severity.ADVISORY,
            potential_consequence=("The footprint covers a minority of this gene's "
                                   "pathogenic variation and must not be described "
                                   "as its pathogenic region."),
            recommended_action=("State the covered fraction wherever the footprint "
                                "is reported."),
            affected_output="evidence_set.tsv")

    if up.r_hot_runner_up_margin_rel is not None and up.r_hot_runner_up_margin_rel < 0.05:
        warn.add(
            "R_FP_INHERITS_R_HOT_NEAR_TIE",
            f"r_hot beat its runner-up by "
            f"{100 * up.r_hot_runner_up_margin_rel:.3f} %. r_fp is derived from the "
            f"evidence set that r_hot defines, so the footprint boundary rests on "
            f"that margin.",
            severity=Severity.ADVISORY,
            potential_consequence=("Component count and residue count can change "
                                   "with a re-seeded run."),
            recommended_action=("Report r_contain at both radii; the sensitivity "
                                "band is in sensitivity_band.tsv."),
            affected_output="sensitivity_band.tsv")


    # [DECISION-RFP-CONTAINMENT-0001] Carried over from the retired
    # _emit_selection_warnings: this one was never about selection. If the box
    # exceeded the voxel budget the edge coarsens, and every geometric descriptor
    # below is then measured on a coarser lattice. Dropping it with the rest of that
    # block would have made the condition unreportable. The pad fix in this same
    # change makes the fallback rarer, not impossible.
    if solution.spec.fallback_applied:
        warn.add("VOXEL_RESOLUTION_FALLBACK",
                 f"voxel edge fell back from {params.voxel_h_A} A to "
                 f"{solution.spec.h} A above "
                 f"{params.voxel_count_fallback_threshold:,} voxels",
                 potential_consequence="absolute volumes are coarser; comparisons "
                                       "within this run remain internally consistent",
                 affected_output="footprint_domain.json")
        relative = (h_sensitivity or {}).get("relative_difference")
        if relative is not None and relative > 0.10:
            escalations.append({
                "condition": "VOXEL_RESOLUTION_FALLBACK_MATERIAL",
                "detail": (f"the resolution fallback changes the reported footprint "
                           f"volume by {relative:.1%} "
                           f"({h_sensitivity['volume_used_A3']:.4g} A^3 at h="
                           f"{h_sensitivity['voxel_h_used_A']} vs "
                           f"{h_sensitivity['volume_alternative_A3']:.4g} A^3 at h="
                           f"{h_sensitivity['voxel_h_alternative_A']})"),
                "options": ["accept the coarser grid and its recorded volumes",
                            "Lead-authorized re-run under a new RUN_ID with more "
                            "memory so the fine grid fits"],
                "recommendation": ("report both volumes; the radius itself is "
                                   "unaffected because containment is measured on "
                                   "C-alpha distances, not on the voxel grid. The "
                                   "voxel edge is never tuned toward a volume."),
            })


def params_detection_floor() -> float:
    from .api import DETECTION_SCALE_FLOOR_A
    return DETECTION_SCALE_FLOOR_A

