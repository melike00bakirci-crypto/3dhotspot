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
                         run_id: str, gene: str, diagnostics=None) -> Path:
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
    if solution.r_fp is None:
        add("**No admissible footprint radius was identified.** This is a complete "
            "scientific result, not a failure; the domain, every candidate and every "
            "rejection reason are retained in `r_fp_scan.tsv`.")
    else:
        norm = result.normalized.get(solution.r_fp, {})
        add(f"- **`r_fp` = {_fmt(solution.r_fp)} A**")
        add(f"- normalized objective vector "
            f"(compactness, connectivity, stability, parsimony) = "
            f"({_fmt(norm.get('compactness'))}, {_fmt(norm.get('connectivity'))}, "
            f"{_fmt(norm.get('stability'))}, {_fmt(norm.get('parsimony'))})")
        add(f"- L2 distance to utopia (1,1,1,1) = "
            f"{_fmt(result.distances.get(solution.r_fp))}")
        add(f"- Pareto members: {len(result.pareto_members)} of "
            f"{result.n_admissible} admissible candidates")
        add(f"- tie chain applied: {' -> '.join(result.tie_chain_applied)}")
        add("")
        if diagnostics is not None:
            add(f"- effective objectives entering distance-to-utopia: "
                f"`{diagnostics.effective_objectives}`. Degenerate (identical at "
                f"every admissible candidate, therefore non-informative and "
                f"excluded): `{diagnostics.degenerate_objectives}`.")
            add("")
            if diagnostics.vacuous_pareto:
                add(f"> **VACUOUS_PARETO_SELECTION** — {diagnostics.vacuous_reason}. "
                    f"`r_fp` was not selected by trading objectives off; it was the "
                    f"only non-dominated survivor.")
                add("")
            if diagnostics.admissibility_dominates:
                add(f"> **ADMISSIBILITY_DOMINATES_SELECTION** — "
                    f"{_fmt(diagnostics.inadmissible_fraction)} of scanned "
                    f"candidates were inadmissible, by rule "
                    f"`{diagnostics.inadmissible_by_rule}`. Admissibility, not the "
                    f"objectives, decided the outcome. Significance is never one "
                    f"of the eliminating rules (QC_F1/F2/F3 are purely geometric).")
                add("")
        add("#### Stability behaviour around the selection")
        add("")
        add("| rho (A) | stability | ARI(prev) | ARI(next) | n_comp | coverage |")
        add("|---|---|---|---|---|---|")
        centre = [i for i, r in enumerate(rows) if abs(r.rho - solution.r_fp) < 1e-9]
        span = range(max(0, centre[0] - 2), min(len(rows), centre[0] + 3)) if centre \
            else range(len(rows))
        for i in span:
            r = rows[i]
            mark = " **<-**" if abs(r.rho - solution.r_fp) < 1e-9 else ""
            add(f"| {_fmt(r.rho)}{mark} | {_fmt(r.stability)} | {_fmt(r.ari_prev)} | "
                f"{_fmt(r.ari_next)} | {r.n_components} | {_fmt(r.coverage)} |")
        add("")

    add("### Rejected candidates")
    add("")
    inadmissible = [(r.rho, solution.verdicts[r.rho].reason_text)
                    for r in rows if not solution.verdicts[r.rho].admissible]
    dominated = sorted(result.dominators)
    add(f"- inadmissible (QC): {len(inadmissible)} of {len(rows)}")
    for rho, reason in inadmissible[:15]:
        add(f"  - {_fmt(rho)} A — {reason}")
    if len(inadmissible) > 15:
        add(f"  - ... {len(inadmissible) - 15} further inadmissible radii, all "
            f"retained with full metrics in `r_fp_scan.tsv`")
    add(f"- Pareto-dominated (admissible but dominated): {len(dominated)}")
    add("")
    add("Every rejected candidate keeps all of its computed metrics and its rejection "
        "reason in `r_fp_scan.tsv`; `footprint_pareto_dominated.json` records the "
        "dominating vector for each dominated candidate.")
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
