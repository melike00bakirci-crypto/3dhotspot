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
from . import selection as sel
from .api import (FOOTPRINT_API_VERSION, FootprintSolution, build_footprint,
                  footprint_code_version)
from .inputs import UpstreamInputs
from .params import FootprintParams
from .staging import stage_files, write_stage_manifest

AGENT = "footprint-robustness"
STAGE_07 = "07_FOOTPRINT_RADIUS"
STAGE_08 = "08_FINAL_FOOTPRINT"
STAGE_09 = "09_ROBUSTNESS"

EXPECTED_07 = [
    ("r_fp_scan.tsv", "the single decision trace: every candidate radius, every metric"),
    ("footprint_domain.json", "II.9 domain: MST, merge scales, rho_all, grid, voxel h"),
    ("footprint_admissibility.tsv", "verbatim QC projection of r_fp_scan.tsv (P1)"),
    ("merge_events.tsv", "merge events with participating components and neck widths"),
    ("footprint_objective_matrix.json", "objective matrix and correlation diagnostic"),
    ("footprint_pareto.json", "Pareto members with normalized values and distances"),
    ("footprint_pareto_dominated.json", "dominated and inadmissible candidates retained"),
    ("footprint_decision.json", "the selection record"),
    ("footprint_domain_boundary_diagnostic.json", "BW1/BW2/BW3 boundary diagnostic"),
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
        solution = build_footprint(up.centers, up.universe, params)
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
    redundancy = sel.redundancy_flags(solution.selection, params)
    export.write_objective_matrix(dir_07, solution, params, redundancy)
    export.write_pareto(dir_07, solution)
    export.write_boundary_diagnostic(dir_07, solution)

    # -- v2 §4/§7: objective/constraint separation diagnostics, mirrored from
    # Stage B (hotspot/selection.py). Post-hoc over the already-selected result and
    # QC verdicts — never a change to build_footprint(), which Phase D re-executes
    # per iteration (see footprint/selection.py module docstring).
    diagnostics = sel.diagnose_selection(solution.selection, solution.verdicts, params)
    if not diagnostics.degenerate_exclusion_is_numerically_neutral:
        raise BlockedError(
            "BLOCKED — excluding the degenerate footprint objectives from "
            "distance-to-utopia changes a distance. Under the frozen min-max rule a "
            "degenerate objective normalizes to 1.0 everywhere and must contribute "
            "exactly 0; a discrepancy means the selection machinery has drifted."
        )

    _emit_selection_warnings(warn07, solution, params, redundancy, escalations,
                             scan_rows, h_sensitivity, diagnostics)

    figures.figure_f7(dir_07 / "figures", solution, params)
    figures.figure_f8(dir_07 / "figures", solution, params)
    export.write_decision(dir_07, solution, params, escalations, diagnostics)

    # -- escalation: coverage constraint failure --------------------------
    if solution.coverage_constraint_failure:
        return _finish_coverage_escalation(ctx, dir_07, dir_08, warn07, warn08,
                                           solution, params, up, started, t0)

    # -- negative: no admissible candidate --------------------------------
    if solution.r_fp is None:
        neg = NegativeResult(
            "NO_ADMISSIBLE_FOOTPRINT",
            "No candidate radius in the derived domain satisfied QC-F1/F2/F3. "
            "Every candidate is retained in r_fp_scan.tsv with its metrics and "
            "rejection reason.",
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
        ctx.run_id, ctx.gene, diagnostics)

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

def _emit_selection_warnings(warn: WarningCollector, solution: FootprintSolution,
                             params: FootprintParams, redundancy: list[dict],
                             escalations: list[dict], scan_rows: list[dict],
                             h_sensitivity: dict | None = None,
                             diagnostics: "sel.SelectionDiagnostics | None" = None
                             ) -> None:
    if diagnostics is not None:
        _emit_objective_constraint_warnings(warn, diagnostics, params, escalations)

    for row in scan_rows:
        if row["excessive_footprint_coverage"]:
            warn.add("EXCESSIVE_FOOTPRINT_COVERAGE",
                     f"candidate r_fp = {row['r_fp_A']:g} A covers "
                     f"{row['coverage_fraction']:.4g} of U_struct, above the frozen "
                     f"QC-F1 limit {params.qc_f1_max_coverage:g}",
                     potential_consequence="candidate is inadmissible",
                     recommended_action="none — the 0.50 threshold is never relaxed",
                     affected_output="r_fp_scan.tsv", status="ACCEPTED_BY_DESIGN")

    if solution.selection.near_tie:
        detail = solution.selection.near_tie_detail
        warn.add("NEAR_TIE_RADIUS_SELECTION",
                 f"r_fp = {detail.get('first')} A and {detail.get('second')} A differ "
                 f"by {detail.get('relative_gap'):.3g} in distance-to-ideal while "
                 f"being more than {params.near_tie_radius_separation_steps} steps "
                 f"apart",
                 potential_consequence="a materially different radius was nearly "
                                       "selected",
                 affected_output="footprint_decision.json")
        escalations.append({
            "condition": "NEAR_TIE_RADIUS_SELECTION",
            "detail": f"near-tie between {detail.get('first')} A and "
                      f"{detail.get('second')} A (relative gap "
                      f"{detail.get('relative_gap'):.3g})",
            "options": ["accept the tie-chain outcome",
                        "Lead-authorized re-run under a new RUN_ID with a "
                        "pre-registered finer grid"],
            "recommendation": "accept the tie-chain outcome; it is deterministic and "
                              "pre-registered. Reported for visibility.",
        })

    if redundancy:
        pairs = ", ".join(f"{r['objective_a']}~{r['objective_b']} "
                          f"(r={r['pearson_r']:.3f})" for r in redundancy)
        warn.add("OBJECTIVE_REDUNDANCY",
                 f"objective pairs exceeding |rho| = "
                 f"{params.objective_redundancy_abs_rho}: {pairs}",
                 potential_consequence="the effective dimensionality of the Pareto "
                                       "problem is lower than four",
                 affected_output="footprint_objective_matrix.json")
        escalations.append({
            "condition": "OBJECTIVE_REDUNDANCY",
            "detail": f"objective correlation exceeds the configured threshold: {pairs}",
            "options": ["proceed — the objective set is FROZEN and equally weighted",
                        "Lead-authorized methodological revision under a new RUN_ID"],
            "recommendation": "proceed and report. The objective set and w_k = 1 are "
                              "frozen; changing them after seeing the correlation "
                              "would be outcome-dependent weighting.",
        })

    if solution.boundary.get("DOMAIN_BOUNDARY_WARNING"):
        warn.add("BOUNDARY_OPTIMUM_WARNING",
                 f"the r_fp optimum sits in the "
                 f"{solution.boundary.get('which_end')} boundary band "
                 f"({solution.boundary.get('band_width_points')} grid points); "
                 f"binding constraint: {solution.boundary.get('binding_edge')}",
                 potential_consequence="the true optimum may lie outside the derived "
                                       "domain",
                 recommended_action="Lead-authorized re-run under a NEW RUN_ID with an "
                                    "explicitly widened, pre-registered domain; "
                                    "automatic widening is prohibited",
                 affected_output="footprint_domain_boundary_diagnostic.json")
        escalations.append({
            "condition": "DOMAIN_BOUNDARY_WARNING",
            "detail": f"band = {solution.boundary.get('band_members')}, binding edge "
                      f"= {solution.boundary.get('binding_edge')}",
            "options": ["accept the selection as-is",
                        "Lead-authorized re-run with a widened, pre-registered domain "
                        "under a new RUN_ID"],
            "recommendation": "report the band and the binding constraint; recommend a "
                              "Lead-authorized widened re-run. The domain is never "
                              "widened here.",
        })

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
                "recommendation": ("report both volumes; the selection itself is "
                                   "unaffected because every candidate was measured "
                                   "on the same grid. The voxel edge is never tuned "
                                   "toward a preferred volume."),
            })


def _emit_objective_constraint_warnings(warn: WarningCollector,
                                        diagnostics: "sel.SelectionDiagnostics",
                                        params: FootprintParams,
                                        escalations: list[dict]) -> None:
    """v2 §4/§7 objective/constraint separation, mirrored from Stage B's
    ``hotspot/stage.py::select_radius`` block and adapted to footprint's own
    objective names and QC_F1/F2/F3 rule codes."""
    if diagnostics.degenerate_objectives:
        warn.add(
            "DEGENERATE_OBJECTIVE",
            f"Objective(s) {diagnostics.degenerate_objectives} take the identical "
            f"value at every admissible r_fp candidate and discriminate nothing. "
            f"Excluded from distance-to-utopia; the effective objective set is "
            f"{diagnostics.effective_objectives}.",
            potential_consequence=(
                "Fewer objectives are actually deciding r_fp than the four the "
                "design declares, so 'multi-objective' overstates the decision."),
            recommended_action=("Report the effective objective set alongside "
                                "r_fp. The objective set is frozen for this run."),
            affected_output="footprint_decision.json")
    if diagnostics.vacuous_pareto:
        warn.add(
            "VACUOUS_PARETO_SELECTION",
            f"VACUOUS_PARETO_SELECTION: {diagnostics.vacuous_reason}.",
            potential_consequence=(
                "r_fp was not chosen by trading objectives off against each "
                "other; the geometry was the only non-dominated survivor, not a "
                "trade-off outcome."),
            recommended_action=("Report the Pareto set size and the effective "
                                "objective set prominently beside r_fp."),
            affected_output="footprint_pareto.json")
        escalations.append({
            "condition": "VACUOUS_PARETO_SELECTION",
            "detail": f"VACUOUS_PARETO_SELECTION at r_fp: {diagnostics.vacuous_reason}.",
            "options": ["accept r_fp and report the vacuity prominently "
                        "(current behaviour)",
                        "Lead review of whether the objective set discriminates at "
                        "all for this footprint, under a NEW RUN_ID"],
            "recommendation": ("accept and report; the objective set is frozen "
                               "for this run. The reported radius carries no "
                               "trade-off information; changing the objective set "
                               "after seeing results would be outcome-dependent and "
                               "is prohibited."),
        })
    if diagnostics.admissibility_dominates:
        warn.add(
            "ADMISSIBILITY_DOMINATES_SELECTION",
            f"{diagnostics.inadmissible_fraction:.1%} of scanned r_fp candidates "
            f"were ruled inadmissible (threshold "
            f"{params.admissibility_dominates_fraction:.0%}); by rule: "
            f"{diagnostics.inadmissible_by_rule}.",
            potential_consequence=("Admissibility, not the objectives, decided the "
                                   "outcome."),
            recommended_action=("Report which QC_F rule did the eliminating. "
                                "Footprint admissibility is purely geometric "
                                "(coverage, bridging, undefined geometry) and is "
                                "never conditioned on significant-center count."),
            affected_output="r_fp_scan.tsv")


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
               affected_output="r_fp_scan.tsv")
    reason = (
        f"ALL_FOOTPRINT_CANDIDATES_EXCESSIVE_COVERAGE = TRUE: {fraction:.1%} of the "
        f"{len(solution.sweep.rows)} candidate radii in "
        f"[{solution.domain.rho_min:g}, {solution.domain.rho_max:g}] A cover more than "
        f"{params.qc_f1_max_coverage:g} of U_struct. The full coverage curve is in "
        f"r_fp_scan.tsv."
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
    detail = {"condition": neg.condition, "detail": neg.detail, **neg.context}
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
        "near_tie_flag": False, "domain_boundary_warning": False,
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
        "near_tie_flag": bool(solution.selection.near_tie),
        "domain_boundary_warning": bool(
            solution.boundary.get("DOMAIN_BOUNDARY_WARNING")),
        "code_version": code_ver,
        "footprint_api_version": FOOTPRINT_API_VERSION,
        "n_significant_centers": len(up.centers),
        "center_set_sha256": up.centers.digest(),
        "r_fp_independent_of_r_hot": True,
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
