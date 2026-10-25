"""``stage_c_report.md`` — the narrative record of the footprint decision.

Reports what happened, including the parts that are unflattering: rejected
candidates with their reasons, artificial bridging present in the final solution,
and any boundary or near-tie condition. An awkward footprint is a finding.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..utils.io import write_text
from .api import FootprintSolution, Universe
from .params import FootprintParams


def _fmt(value, spec: str = ".4g") -> str:
    if value is None:
        return "NA"
    if isinstance(value, float) and not np.isfinite(value):
        return "NA"
    return format(value, spec) if isinstance(value, (int, float)) else str(value)


def write_stage_c_report(path: Path, solution: FootprintSolution,
                         universe: Universe, params: FootprintParams,
                         center_ids: tuple[int, ...], warnings: list,
                         escalations: list[dict], code_version: str,
                         run_id: str, gene: str) -> Path:
    domain = solution.domain
    result = solution.selection
    state = solution.selected
    rows = solution.sweep.rows

    lines: list[str] = []
    add = lines.append

    add(f"# Stage C — footprint radius and final footprint ({gene})")
    add("")
    add(f"- run_id: `{run_id}`")
    add(f"- footprint code_version: `{code_version}`")
    add(f"- |S| (significant hotspot centers, immutable input): **{len(center_ids)}**")
    add(f"- |U_struct|: **{len(universe)}**")
    add("")

    add("## SCIENTIFIC DECISIONS")
    add("")
    add("### Derived search domain (II.9)")
    add("")
    add(f"- `rho_min` = {_fmt(domain.rho_min)} A (frozen)")
    add(f"- `rho_max` = {_fmt(domain.rho_max)} A, binding constraint: "
        f"**{domain.binding_constraint}**")
    add(f"- `rho_all` (max MST merge scale) = {_fmt(domain.rho_all)} A")
    add(f"- `D_max` (U_struct diameter) = {_fmt(domain.d_max)} A")
    add(f"- uniform step = {_fmt(domain.step)} A, grid points G = {len(domain.grid)}")
    if domain.single_center_special_case:
        add("- **special case**: |S| = 1, the MST is empty, so the post-merge term "
            "drops out of `rho_max` (recorded, not worked around).")
    if domain.merge_scales:
        driving = sorted(domain.merge_scales, reverse=True)[:5]
        add(f"- driving merge scales (largest first, A): "
            f"{', '.join(_fmt(v) for v in driving)}")
    add("")
    add("Merge radii informed `rho_max` and are reported below as events. They never "
        "determined `r_fp` (F10). `r_hot` is not an argument to the footprint API and "
        "was used nowhere in this decision.")
    add("")

    add("### Selected footprint radius")
    add("")
    # [DECISION-RFP-CONTAINMENT-0001] This section used to report the normalized
    # objective vector, the L2 distance to utopia (1,1,1,1), the Pareto membership
    # count, the tie chain and the v2 degeneracy diagnostics. None of those happened:
    # r_fp is DETERMINED by containment, and a report describing a selection that did
    # not run is how a reader ends up citing the wrong methodology.
    record = solution.containment or {}
    if solution.r_fp is None:
        add("**No footprint radius survived QC.** This is a complete scientific "
            "result, not a failure. The radius was determined — the determination is "
            "in `evidence_set.tsv` and `footprint_decision.json` — and then rejected "
            "on geometry; it is never shrunk to pass a gate, because that would drop "
            "an evidence residue. Every metric and the rejection reason are retained "
            "in `footprint_descriptors.tsv`, and `sensitivity_band.tsv` shows what the "
            "footprint would have been at every radius in the plausible range.")
    else:
        add(f"- **`r_fp` = {_fmt(solution.r_fp)} A** — determined, not selected")
        add(f"- `r_contain` = {_fmt(record.get('r_contain_exact_A'))} A exactly, "
            f"rounded UP to the {_fmt(params.step_fp_A)} A grid "
            f"(slack {_fmt(record.get('grid_slack_A'))} A). Rounding up is a "
            f"correctness condition: rounding down would drop the binding residue.")
        add(f"- bound by residue **{record.get('binding_residue')}**, whose nearest "
            f"significant center is {record.get('binding_nearest_center')}")
        add(f"- evidence set: {record.get('n_evidence_set')} of "
            f"{record.get('n_plp_mapped')} mapped P/LP residues lie within `r_hot` of "
            f"some significant center; {record.get('n_plp_outside_evidence')} "
            f"({_fmt(record.get('fraction_plp_outside_evidence'))}) lie outside and no "
            f"center-anchored footprint can reach them")
        add(f"- evidence distances: median {_fmt(record.get('evidence_distance_median_A'))} A, "
            f"IQR {_fmt(record.get('evidence_distance_iqr_A'))} A "
            f"(Q1 {_fmt(record.get('evidence_distance_q1_A'))}, "
            f"Q3 {_fmt(record.get('evidence_distance_q3_A'))})")
        add("")
        if record.get("binding_residue_is_outlier"):
            add(f"> **CONTAINMENT_BINDING_RESIDUE_IS_OUTLIER** — residue "
                f"{record.get('binding_residue')} sits beyond the "
                f"{_fmt(record.get('outlier_fence_A'))} A upper fence of the evidence "
                f"distance distribution. `r_fp` therefore rests on one unusually "
                f"distant residue rather than on the body of the evidence. The radius "
                f"is NOT adjusted for this — dropping the residue would break "
                f"containment — but the dependence is reported so it can be weighed.")
            add("")
        if record.get("detection_scale_below_resolution"):
            add(f"> **DETECTION_SCALE_BELOW_3D_RESOLUTION** — the determined radius is "
                f"below the scale at which this voxel grid resolves structure. The "
                f"geometric descriptors below are reported but should not be read as "
                f"shape measurements.")
            add("")
        add("No objective vector, utopia point, Pareto front or tie chain appears "
            "here: with one candidate there is nothing to trade off and nothing to "
            "break a tie between. The geometric descriptors in "
            "`footprint_descriptors.tsv` are reported and decide nothing.")
        add("")
        add("#### Sensitivity to `r_hot`")
        add("")
        add("`r_fp` is a function of the evidence set, which `r_hot` defines, and "
            "`r_hot` is still chosen by a multi-objective procedure that can turn on a "
            "sub-percent margin. `sensitivity_band.tsv` is the mitigation: it lists "
            "the footprint at every radius from the resolution floor up to `r_hot`, "
            "with the determined one marked, so how much of the result rests on the "
            "exact radius is readable rather than argued.")
        add("")

    add("### QC verdict on the determined radius")
    add("")
    inadmissible = [(r.rho, solution.verdicts[r.rho].reason_text)
                    for r in rows if not solution.verdicts[r.rho].admissible]
    add(f"- candidates evaluated: {len(rows)} (containment determines one radius)")
    add(f"- inadmissible (QC): {len(inadmissible)}")
    for rho, reason in inadmissible:
        add(f"  - {_fmt(rho)} A — {reason}")
    add("")
    add("There is no Pareto-dominated category: domination is a relation between "
        "candidates and there is only one. A QC failure escalates and is reported; "
        "it never moves the radius. Every metric and the rejection reason are "
        "retained in `footprint_descriptors.tsv`.")
    add("")

    add("### Artificial bridging in the final solution")
    add("")
    if state is None:
        add("- not applicable (no radius selected)")
    else:
        final_events = [e for e in solution.sweep.merge_events
                        if abs(e.rho - state.rho) < 1e-9]
        artificial = [e for e in final_events if e.connection_type == "artificial"]
        add(f"- bridges detected by relative erosion at "
            f"{_fmt(params.erosion_fraction * state.rho)} A "
            f"(= {params.erosion_fraction:g} x r_fp): **{state.n_bridges}**")
        add(f"- bridging index BI = **{_fmt(state.bridging_index)}** "
            f"(QC-F2 limit {params.qc_f2_max_bridging_index:g})")
        add(f"- merge events at this radius: {len(final_events)}, of which "
            f"{len(artificial)} classified **artificial**")
        for e in artificial:
            add(f"  - neck width {_fmt(e.neck_width_A)} A joining centers "
                f"{[int(center_ids[i]) for i in e.component_a_centers]} and "
                f"{[int(center_ids[i]) for i in e.component_b_centers]}")
        if artificial:
            add("  - bridges are measured and reported, never hand-removed; QC-F2 "
                "decides admissibility.")
    add("")

    add("## QC RESULTS")
    add("")
    if state is not None:
        add(f"- coverage at `r_fp` = **{_fmt(state.coverage)}** "
            f"(QC-F1 limit {params.qc_f1_max_coverage:g}, denominator |U_struct| = "
            f"{len(universe)})")
        add(f"- bridging index at `r_fp` = **{_fmt(state.bridging_index)}** "
            f"(QC-F2 limit {params.qc_f2_max_bridging_index:g})")
        add(f"- connected components at `r_fp` = **{state.n_components}** "
            f"(a measurement, not a target)")
        add(f"- QC-F5: all {len(center_ids)} significant centers contained — asserted")
    add(f"- fraction of candidates failing QC-F1: "
        f"{_fmt(solution.excessive_coverage_fraction)}")
    add(f"- voxel edge h = {_fmt(solution.spec.h)} A"
        + (" (**fallback applied**)" if solution.spec.fallback_applied else ""))
    add(f"- grid = {list(solution.spec.shape)} voxels "
        f"({solution.spec.n_voxels:,} total); one Euclidean distance transform was "
        f"computed and reused across all {len(domain.grid)} radii")
    violations = solution.sweep.monotonicity_violations
    add(f"- monotonicity checks (V non-decreasing, n_comp non-increasing): "
        f"{'PASS' if not violations else f'{len(violations)} VIOLATION(S) — recorded as numerical faults, never smoothed'}")
    add("")

    add("## WARNINGS")
    add("")
    if not warnings:
        add("- none")
    for w in warnings:
        add(f"- **{w.severity.value}** `{w.warning_code}` — {w.message}")
    add("")
    if escalations:
        add("## ESCALATIONS TO THE LEAD")
        add("")
        for esc in escalations:
            add(f"- **{esc['condition']}** — {esc['detail']}")
            add(f"  - recommendation: {esc['recommendation']}")
        add("")

    add("## Boundary diagnostic")
    add("")
    add(f"- `DOMAIN_BOUNDARY_WARNING` = "
        f"{solution.boundary.get('DOMAIN_BOUNDARY_WARNING')} "
        f"(BW1={solution.boundary.get('BW1_selected_in_band')}, "
        f"BW2={solution.boundary.get('BW2_pareto_majority_in_band')}, "
        f"BW3={solution.boundary.get('BW3_top3_same_band')})")
    add(f"- band width b = {solution.boundary.get('band_width_points')} grid points; "
        f"binding edge: {solution.boundary.get('binding_edge')}")
    add("- the diagnostic is evaluated and recorded whether or not it fires; the "
        "domain is never widened automatically.")
    add("")
    return write_text(path, "\n".join(lines))
