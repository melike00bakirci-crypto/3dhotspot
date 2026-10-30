"""``stage_b_report.md`` — the narrative Stage B record.

Follows the agent's reporting protocol exactly: STATUS / INPUTS USED / METHODS
EXECUTED / OUTPUTS GENERATED / QC RESULTS / SCIENTIFIC DECISIONS / WARNINGS / FAILED
OR REJECTED ANALYSES / UNRESOLVED ISSUES / HANDOFF.

The report is written for a reader who did not run the pipeline, so it states the
things that are easy to get wrong: that a hotspot *center* need not carry a ClinVar
variant, that LOO-MCC is a radius-selection objective and not validation, that the
PCF peak is not ``r_hot``, and that the sensitivity analyses cannot redefine anything.

**[v2 §6] Every report opens with a one-page decision summary** giving, in order:
terminal state; ``N_P``, ``N_B``, ``|U_struct|``, ``|U_center|``, ``m``; the null used
for the primary test; ``p_floor``, ``c_1`` and their ratio; the selected radius and
whether the Pareto set was vacuous; and the number of significant centers. A reader
must be able to decide from that page alone whether a negative result is
interpretable at all.
"""
from __future__ import annotations


def _fmt(value, digits: int = 6) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return str(value)


def _table(headers: list[str], rows: list[list]) -> str:
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join("---" for _ in headers) + "|"]
    for row in rows:
        out.append("| " + " | ".join(_fmt(c) for c in row) + " |")
    return "\n".join(out)


PROHIBITED_INFERENCES = (
    "> **Prohibited inferences (v2 Appendix A).** None of the following may be "
    "written about this run: *\"No 3D hotspot is detectable in this gene\"*, "
    "*\"This is a valid, complete scientific negative\"*, *\"This result licenses no "
    "modification of the method\"*. The correct statement is: **Under the "
    "configuration executed, the test could not have rejected any center regardless "
    "of the data; the analysis is uninformative about the presence or absence of "
    "hotspots in this gene.**"
)


def _decision_summary(s: dict, status: str) -> list[str]:
    """v2 §6 — the one page that decides whether the rest is interpretable."""
    inputs = s.get("inputs", {})
    detection = s.get("detection", {})
    selection = s.get("selection")
    cert_block = s.get("power_certificate") or {}
    cert = (cert_block.get("binding") or cert_block.get("posthoc")
            or cert_block.get("preflight") or {})

    L = ["## DECISION SUMMARY", "",
         f"**Terminal state: {status}**", ""]
    L.append(_table(["quantity", "value"], [
        ["terminal state", status],
        ["N_P (P/LP residues)", inputs.get("N_P")],
        ["N_B (B/LB residues)", inputs.get("N_B")],
        ["|U_struct| (positional reference)", inputs.get("M")],
        ["|U_center| (test family before restriction)",
         cert.get("n_center_universe_U_center",
                  inputs.get("n_center_universe", inputs.get("M")))],
        ["m (family actually tested)",
         cert.get("m_test_family_size", detection.get("n_in_test_family"))],
        ["null used for the PRIMARY per-center test",
         s.get("primary_null") or cert.get("null_model")],
        ["p_floor", cert.get("p_floor")],
        ["c_1 = q/m", cert.get("c_1_rank1_critical_value")],
        ["p_floor / c_1", cert.get("p_floor_over_c_1")],
        ["binding floor", cert.get("binding_floor")],
        ["power certificate", "PASSED" if cert.get("PASSES") else "FAILED"],
        ["selected r_hot (A)", (selection or {}).get("r_hot")],
        ["Pareto set vacuous", (selection or {}).get("vacuous_pareto")],
        ["significant hotspot centers", detection.get("n_significant_bh", 0)],
    ]))
    L.append("")
    if cert and not cert.get("PASSES"):
        L.append(PROHIBITED_INFERENCES)
        L.append("")
        if cert.get("increasing_B_is_futile"):
            L.append("> The binding floor is **combinatorial**, so increasing `B` is "
                     "futile and is not recommended. The remedies are a larger benign "
                     "cohort, a different null (v2 §5.2), or a smaller test family "
                     "(v2 §2/§5.1).")
        else:
            L.append(f"> The binding floor is the **permutation resolution**. "
                     f"`B = {cert.get('B_to_clear_c_1')}` would clear `c_1`; the run "
                     f"must be re-executed at that `B` under the RUN_ID policy. The "
                     f"agent never raises `B` itself.")
        L.append("")
    return L


def render_report(s: dict) -> str:
    negative = s.get("negative_result")
    selection = s.get("selection")
    detection = s.get("detection", {})
    domain = s.get("domain", {})
    inputs = s.get("inputs", {})
    underpowered = bool(s.get("underpowered"))

    if underpowered:
        status = "UNDERPOWERED — TEST CANNOT REJECT (uninformative, NOT a negative)"
    elif negative and not selection:
        status = f"COMPLETED_NEGATIVE — {negative['condition']}"
    elif negative:
        status = ("COMPLETED_NEGATIVE — no center survives FDR"
                  if negative.get("settled_negative")
                  else "BLOCKED — RESOLUTION-LIMITED empty result (NOT evidence of absence)")
    else:
        majors = sum(1 for w in s.get("warnings", []) if w["severity"] == "MAJOR")
        status = "COMPLETED" + (f" (with {majors} MAJOR warning(s))" if majors else "")

    L: list[str] = []
    L.append(f"# Stage B — Primary spatial hotspot discovery ({s.get('gene', 'NA')})")
    L.append("")
    L.append(f"*Agent:* `{s.get('agent')}` *Run:* `{s.get('run_id')}` "
             f"*Code version:* `{str(s.get('code_version'))[:12]}` "
             f"*Config:* `{str(s.get('config_sha256'))[:12]}`")
    L.append("")

    L.extend(_decision_summary(s, status))

    L.append("## STATUS")
    L.append("")
    L.append(f"**{status}**")
    if s.get("synthetic"):
        L.append("")
        L.append("> Synthetic input mode — these numbers describe generated data, "
                 "not a real gene.")
    L.append("")

    L.append("## INPUTS USED")
    L.append("")
    L.append(_table(["quantity", "value"], [
        ["upstream handoff", inputs.get("upstream_handoff", "handoff_01")],
        ["upstream qc_status", inputs.get("upstream_qc")],
        ["|U_struct| (M)", inputs.get("M")],
        ["|U_center|", inputs.get("n_center_universe")],
        ["|L| (N)", inputs.get("N")],
        ["N_P (P/LP residues)", inputs.get("N_P")],
        ["N_B (B/LB residues)", inputs.get("N_B")],
        ["D_max", inputs.get("D_max_A")],
        ["r_max = floor(D_max/2)", inputs.get("r_max_A")],
        ["B (permutations)", s.get("B")],
        ["q (BH level)", s.get("q")],
        ["kappa (LOO smoothing)", s.get("kappa")],
    ]))
    L.append("")
    L.append(f"Column allowlist actually loaded: `{inputs.get('allowlist', [])}`. "
             f"Forbidden downstream columns are refused at load time; review-star "
             f"metadata reaches only the separate, non-redefining sensitivity channel.")
    L.append("")

    L.append("## METHODS EXECUTED")
    L.append("")
    L.append("1. **II.1** Ripley-type K and pair correlation, separately for P/LP and "
             "B/LB, under the structure-aware positional null (uniform "
             "without-replacement subsets of `U_struct`), with 2.5/97.5 envelopes and "
             "the global test `T = max_r |Z(r)|`.")
    L.append("2. **II.2** Deterministic candidate radius domain from the P/LP PCF; the "
             "four fallback triggers evaluated; the PCF peak used for nothing.")
    L.append("3. **v2 §5.4** Power certificate, pre-flight — the minimum p-value the "
             "test could return versus the rank-1 BH critical value `c_1 = q/m`, per "
             "candidate radius, before any permutation compute was spent.")
    L.append("4. **II.4** Permutation-resolution diagnostic, pre-flight and post-hoc.")
    L.append("5. **II.5** Radius scan — LOO-MCC, permutation evidence `Zg`, fold "
             "enrichment and neighbouring-radius stability, each radius evaluated "
             "independently under its own derived seed.")
    L.append("6. **v2 §4** Admissibility as constraints on *undefined quantities only*. "
             "The number of significant centers never removes a radius from the Pareto "
             "set; `QC_H1`/`QC_H2`/`QC_H4` are recorded diagnostics. Inadmissible radii "
             "are retained with their rule and reason.")
    L.append("7. **II.6** Pareto -> min-max normalization over the admissible set -> "
             "utopia (1,1,1,1) -> equally weighted L2 distance over the *effective* "
             "(non-degenerate) objectives -> total tie chain.")
    L.append("8. **II.7 / v2 §5.2** Final residue-centered detection over every residue "
             "of `U_center` (v2 §2 pLDDT restriction of `U_struct`) under the "
             "**structure-aware positional null**, placed over the full `U_struct` (the "
             "same null as the global analysis), BH at q with all boundary ties rejected, "
             "BY reported alongside, three residue objects and hotspot regions. The "
             "label-permutation null is computed and reported **beside** it as a "
             "secondary analysis of a different hypothesis.")
    L.append("9. **v2 §5.4** Power certificate, post-hoc, on the realized family — "
             "binding.")
    L.append("10. **II.7** Three post-primary, explicitly **non-redefining** sensitivity "
             "analyses at the frozen `r_hot`.")
    L.append("")
    if underpowered and s.get("stopped_at_stage"):
        L.append(f"> The list above is the full procedure. **This run stopped at "
                 f"`{s['stopped_at_stage']}`** because the v2 §5.4 power certificate "
                 f"failed; the steps that actually executed are timed below, and every "
                 f"later step was deliberately not run rather than run and discarded.")
        L.append("")
    if s.get("steps"):
        L.append(_table(["step", "wall seconds"],
                        [[st["step"], st["wall_seconds"]] for st in s["steps"]]))
        L.append("")

    L.append("## OUTPUTS GENERATED")
    L.append("")
    L.append("- `04_GLOBAL_CLUSTERING/` — K and PCF curves for both cohorts, envelopes, "
             "the global test, PCF peaks (descriptive only) and "
             "`candidate_radius_domain.json`.")
    L.append("- `05_HOTSPOT_RADIUS/` — `r_hot_scan.tsv` (the single decision trace), the "
             "per-radius center sets, LOO diagnostics, the objective matrix, the Pareto "
             "front and dominated set, `radius_decision.json` and the boundary "
             "diagnostic.")
    L.append("- `06_FINAL_HOTSPOTS/` — `all_residue_center_tests.tsv`, the three residue "
             "objects, hotspot regions, the BH table, the permutation-resolution "
             "diagnostic, structures and `handoff_02.json`.")
    L.append("- `11_SENSITIVITY/` — the three non-redefining analyses and "
             "`SENSITIVITY_IS_NON_REDEFINING.txt`.")
    L.append("")

    # ---- QC RESULTS ---------------------------------------------------------
    L.append("## QC RESULTS")
    L.append("")
    if selection:
        res = s.get("resolution", {}).get("posthoc", {})
        L.append(_table(["quantity", "value"], [
            ["significant hotspot centers (S)", detection.get("n_significant_bh")],
            ["significant centers carrying NO ClinVar variant",
             detection.get("n_centers_without_variant")],
            ["BH boundary p-value", detection.get("bh_boundary_p")],
            ["test family size m", detection.get("n_in_test_family")],
            ["centers excluded (n_L = 0)", detection.get("n_excluded_no_labeled")],
            ["centers excluded (sigma_c = 0)", detection.get("n_excluded_zero_variance")],
            ["coverage at r_hot", detection.get("coverage_at_r_hot")],
            ["LOO sparse proportion at r_hot", detection.get("loo_prop_sparse")],
            ["LOO zero-neighbour proportion at r_hot",
             detection.get("loo_prop_zero_neighbour")],
            ["LOO exact ties at r_hot", detection.get("loo_n_tie")],
            ["hotspot regions", detection.get("n_regions")],
            ["covered residues", detection.get("n_covered_residues")],
            ["classified variants inside hotspots",
             detection.get("n_classified_in_hotspots")],
            ["BY comparison (rejections)", detection.get("n_significant_by")],
            ["BY constant c(m)", detection.get("by_constant")],
            ["provisional centers at r_hot during the scan",
             detection.get("n_provisional_centers_at_r_hot")],
            ["scan vs final disagreement",
             detection.get("scan_final_disagreement")],
        ]))
        L.append("")
        if detection.get("scan_final_disagreement"):
            L.append("> " + (detection.get("scan_final_disagreement_note") or ""))
            L.append("")
        L.append("**A significant hotspot center is a geometric test position and need "
                 "not carry a ClinVar variant.** Every residue of `U_center` "
                 "(pLDDT >= threshold, a v2 §2 restriction of `U_struct`) was a "
                 "candidate center; the positional-null placement population "
                 "remains the full `U_struct`.")
        L.append("")
        boundary = selection.get("boundary", {})
        L.append(f"**DOMAIN_BOUNDARY_WARNING:** {_fmt(boundary.get('DOMAIN_BOUNDARY_WARNING'))} "
                 f"(BW1={_fmt(boundary.get('BW1_selected_in_band'))}, "
                 f"BW2={_fmt(boundary.get('BW2_pareto_majority_in_band'))}, "
                 f"BW3={_fmt(boundary.get('BW3_top3_same_band'))}); band "
                 f"{boundary.get('band_members')}, end {boundary.get('which_end')}. "
                 f"Binding constraint: "
                 f"{'R_FLOOR' if 'lower' in (boundary.get('which_end') or []) else 'r_cap = min(R_CEIL, 0.25*D_max)'}. "
                 f"Automatic domain expansion is prohibited.")
        L.append("")
        margin = detection.get("bh_margin") or {}
        if margin:
            L.append("**How close the evidence came to the BH bar.** The rank-1 critical "
                     "value `q/m` is the hardest bar in the procedure — a single isolated "
                     "center must clear it alone.")
            L.append("")
            L.append(_table(["quantity", "value"], [
                ["test family size m", margin.get("m_family")],
                ["rank-1 BH critical value q/m", margin.get("rank1_critical_value")],
                ["smallest observed p", margin.get("smallest_observed_p")],
                ["permutation resolution floor 1/(B+1)",
                 margin.get("p_resolution_floor")],
                ["margin ratio (smallest p / rank-1 critical value)",
                 margin.get("margin_ratio")],
            ]))
            L.append("")
            L.append(_table(["rank", "observed p", "BH critical value k*q/m", "rejected"],
                            [[r["rank"], r["p_emp"], r["bh_critical_value"],
                              r["reject_bh"]] for r in margin.get("top_ranks", [])]))
            L.append("")
            L.append(f"*{margin.get('margin_ratio_note', '')}*")
            L.append("")

        cert_block = s.get("power_certificate") or {}
        post = cert_block.get("posthoc") or {}
        pre = cert_block.get("preflight") or {}
        if post or pre:
            L.append("**Power certificate (v2 §5.4)** — the mandatory precondition for "
                     "any negative result. A negative may be reported only when this "
                     "has been computed and has passed.")
            L.append("")
            L.append(_table(["quantity", "pre-flight", "post-hoc"], [
                ["null model", pre.get("null_model"), post.get("null_model")],
                ["m (family size)", pre.get("m_test_family_size"),
                 post.get("m_test_family_size")],
                ["p_res = 1/(B+1)", pre.get("p_res"), post.get("p_res")],
                ["p_comb_best", pre.get("p_comb_best"), post.get("p_comb_best")],
                ["p_floor", pre.get("p_floor"), post.get("p_floor")],
                ["c_1 = q/m", pre.get("c_1_rank1_critical_value"),
                 post.get("c_1_rank1_critical_value")],
                ["p_floor / c_1", pre.get("p_floor_over_c_1"),
                 post.get("p_floor_over_c_1")],
                ["binding floor", pre.get("binding_floor"), post.get("binding_floor")],
                ["PASSES", pre.get("PASSES"), post.get("PASSES")],
            ]))
            L.append("")

        secondary = s.get("secondary_null") or {}
        if secondary.get("computed"):
            L.append("**Secondary null, reported side by side (v2 §5.2).** The two "
                     "nulls answer different questions and are never mixed:")
            L.append("")
            L.append(_table(["null", "hypothesis", "significant centers"], [
                [s.get("primary_null"), secondary.get("primary_hypothesis_for_contrast"),
                 detection.get("n_significant_bh")],
                [secondary.get("null_model"), secondary.get("hypothesis"),
                 secondary.get("n_significant_bh")],
            ]))
            L.append("")
            L.append("Only the **primary** row decides anything. Nothing in the "
                     "secondary row enters `significant`, `q_bh`, the three residue "
                     "objects, the hotspot regions or any downstream stage.")
            L.append("")

        L.append("**Permutation-resolution diagnostic** (never changes `B`, never "
                 "changes the FDR method):")
        L.append("")
        L.append(_table(["quantity", "pre-flight", "post-hoc"], [
            ["m", s["resolution"]["preflight"].get("m_family_size"),
             res.get("m_family_size")],
            ["p_res", s["resolution"]["preflight"].get("p_res"), res.get("p_res")],
            ["R", s["resolution"]["preflight"].get("R"), res.get("R")],
            ["n_floor", s["resolution"]["preflight"].get("n_floor"), res.get("n_floor")],
            ["C0", s["resolution"]["preflight"]["conditions"].get("C0"),
             res.get("conditions", {}).get("C0")],
            ["C1", s["resolution"]["preflight"]["conditions"].get("C1"),
             res.get("conditions", {}).get("C1")],
            ["C2", s["resolution"]["preflight"]["conditions"].get("C2"),
             res.get("conditions", {}).get("C2")],
            ["C3", s["resolution"]["preflight"]["conditions"].get("C3"),
             res.get("conditions", {}).get("C3")],
            ["B_req", s["resolution"]["preflight"].get("B_req"), res.get("B_req")],
            ["B_rec", s["resolution"]["preflight"].get("B_rec"), res.get("B_rec")],
        ]))
        L.append("")
        L.append(f"`PERMUTATION_RESOLUTION_LIMITED = "
                 f"{_fmt(s['resolution'].get('limited'))}`")
    else:
        L.append(_table(["quantity", "value"], [
            ["|U_struct|", inputs.get("M")], ["|L|", inputs.get("N")],
            ["N_P", inputs.get("N_P")], ["N_B", inputs.get("N_B")],
            ["outcome", negative.get("condition") if negative else "NA"],
        ]))
    L.append("")

    # ---- SCIENTIFIC DECISIONS ----------------------------------------------
    L.append("## SCIENTIFIC DECISIONS")
    L.append("")
    L.append(f"**Domain branch:** `{domain.get('branch_taken')}`; "
             f"`FALLBACK_RADIUS_DOMAIN = {_fmt(domain.get('FALLBACK_RADIUS_DOMAIN'))}`"
             + (f", trigger `{domain.get('trigger_id')}` "
                f"(all fired: {domain.get('triggers_fired')})"
                if domain.get("FALLBACK_RADIUS_DOMAIN") else "")
             + f". Interval `[{_fmt(domain.get('clamped_interval', [None, None])[0])}, "
               f"{_fmt(domain.get('clamped_interval', [None, None])[1])}]` A on a "
               f"{_fmt(domain.get('step_hot_A'))} A lattice with "
               f"{domain.get('n_grid_points')} grid points.")
    L.append("")
    principal = s.get("global", {}).get("principal_pcf_peak")
    if principal:
        L.append(f"The principal PCF peak is at {_fmt(principal['radius_A'])} A. "
                 f"**It was not used as `r_hot` and did not define the interval alone** "
                 f"(F6/F14).")
        L.append("")
    if selection:
        L.append(f"**Selected `r_hot` = {_fmt(selection['r_hot'])} A**, "
                 f"distance-to-ideal {_fmt(selection['distance_to_ideal'])}, from a "
                 f"Pareto set of {len(selection['pareto_members'])} member(s) over "
                 f"{selection['n_admissible']} admissible of "
                 f"{selection['n_candidates']} scanned radii.")
        L.append("")
        L.append(f"Effective objectives entering distance-to-utopia: "
                 f"`{selection.get('effective_objectives')}`. Degenerate (identical at "
                 f"every admissible radius, therefore non-informative and excluded): "
                 f"`{selection.get('degenerate_objectives')}`.")
        L.append("")
        if selection.get("vacuous_pareto"):
            L.append(f"> **VACUOUS_PARETO_SELECTION** — {selection.get('vacuous_reason')}. "
                     f"The radius was not selected by trading objectives off; combined "
                     f"with a boundary optimum this is a symptom of a power failure "
                     f"(v2 §4/§5.4), not a property of the radius domain.")
            L.append("")
        if selection.get("admissibility_dominates"):
            L.append(f"> **ADMISSIBILITY_DOMINATES_SELECTION** — "
                     f"{_fmt(selection.get('inadmissible_fraction'))} of scanned radii "
                     f"were inadmissible, by rule "
                     f"`{selection.get('inadmissible_by_rule')}`. Admissibility, not the "
                     f"objectives, decided the outcome. Significance is never one of "
                     f"these rules (v2 §4).")
            L.append("")
        L.append("Normalized objective vector at the selected radius (utopia = 1 on "
                 "every axis, equal weights `w_k = 1`):")
        L.append("")
        L.append(_table(["objective", "normalized", "raw"], [
            [k, v, selection["raw_objectives"][f"{selection['r_hot']:g}"][k]]
            for k, v in selection["normalized"].items()]))
        L.append("")
        zg = detection.get("perm_evidence_zg_at_r_hot")
        if zg is not None:
            L.append(f"**Reading `perm_evidence_zg` = {_fmt(zg)}.** "
                     + (s.get("perm_evidence_zg_direction") or ""))
            L.append("")
            if zg < 0:
                L.append("> `perm_evidence_zg` is negative at the selected radius. Under "
                         "this frozen statistic that is the **expected** reading when "
                         "P/LP clustering is present — it is **not** evidence against a "
                         "hotspot. `Zg` is derived from the SECONDARY label-permutation "
                         "null and is a selection objective only; significance comes "
                         "solely from the per-center **structure-aware positional** "
                         "test with BH at `q`, reported above (v2 §5.2).")
                L.append("")
        L.append(f"Per-objective argmax (recorded so the reader can verify that the "
                 f"max-MCC radius holds **no** privileged status, F7/F8): "
                 f"`{selection['per_objective_argmax']}`.")
        L.append("")
        L.append(f"Near-tie flag: `{_fmt(selection['near_tie'])}`"
                 + (f" — {selection['near_tie_detail']}" if selection["near_tie"] else "")
                 + f". Tie chain applied: `{selection['tie_chain']}`.")
        L.append("")
        L.append("Top rejected alternatives:")
        L.append("")
        L.append(_table(["radius (A)", "reason"],
                        [[r, reason] for r, reason in selection.get("top_rejected", [])]))
        L.append("")
        L.append("LOO-MCC is a **radius-selection objective only**. It is never hotspot "
                 "validation and no post-hoc hotspot LOO stage exists (F9/F13). No "
                 "classification threshold was fitted, per radius or at all; the "
                 "`kappa = 2` smoothing is boundary-neutral by construction and this was "
                 "verified numerically at every radius.")
    L.append("")

    # ---- sensitivity --------------------------------------------------------
    L.append("### Post-primary sensitivity analyses (NON-REDEFINING)")
    L.append("")
    sens = s.get("sensitivity", {})
    if sens:
        L.append(_table(["analysis", "status", "n centers", "overlap with primary S",
                         "Jaccard"],
                        [[k, v.get("status"), v.get("n_centers"),
                          v.get("overlap_fraction_of_primary"), v.get("jaccard")]
                         for k, v in sens.items()]))
        L.append("")
    L.append(f"Global-clustering caveat flag: `{s.get('global', {}).get('flag', 'NA')}`. "
             f"A non-significant global result does **not** terminate local discovery "
             f"(F5); it is carried as a caveat by every downstream claim.")
    L.append("")
    L.append("None of these analyses may add or remove a center, change `r_hot`, `q`, "
             "`B` or the FDR method. See `11_SENSITIVITY/SENSITIVITY_IS_NON_REDEFINING.txt`.")
    L.append("")

    # ---- warnings -----------------------------------------------------------
    L.append("## WARNINGS")
    L.append("")
    warnings = s.get("warnings", [])
    if warnings:
        L.append(_table(["severity", "code", "stage", "message"],
                        [[w["severity"], w["warning_code"], w["stage"], w["message"]]
                         for w in warnings]))
    else:
        L.append("None.")
    L.append("")

    # ---- rejected -----------------------------------------------------------
    L.append("## FAILED OR REJECTED ANALYSES")
    L.append("")
    if selection:
        inadmissible = selection.get("inadmissible", {})
        if inadmissible:
            L.append("Radii retained but ruled inadmissible by II.5E (constraints, not "
                     "objectives — none was dropped from the trace):")
            L.append("")
            L.append(_table(["radius (A)", "failing rule(s)"],
                            [[k, v] for k, v in sorted(inadmissible.items(),
                                                       key=lambda kv: float(kv[0]))]))
        else:
            L.append("Every scanned radius was admissible.")
    elif underpowered:
        cert = ((s.get("power_certificate") or {}).get("binding")) or {}
        L.append("The run terminated **UNDERPOWERED**: the per-center test could not "
                 "have rejected any center regardless of the data.")
        L.append("")
        L.append(_table(["quantity", "value"], [
            ["null model", cert.get("null_model")],
            ["N_P", cert.get("N_P")], ["N_B", cert.get("N_B")],
            ["|U_struct|", cert.get("n_universe_U_struct")],
            ["|U_center|", cert.get("n_center_universe_U_center")],
            ["m", cert.get("m_test_family_size")],
            ["B", cert.get("B")],
            ["p_res", cert.get("p_res")],
            ["p_comb_best", cert.get("p_comb_best")],
            ["p_floor", cert.get("p_floor")],
            ["c_1 = q/m", cert.get("c_1_rank1_critical_value")],
            ["p_floor / c_1", cert.get("p_floor_over_c_1")],
            ["binding floor", cert.get("binding_floor")],
            ["increasing B is futile", cert.get("increasing_B_is_futile")],
            ["B that would clear c_1", cert.get("B_to_clear_c_1")],
        ]))
        L.append("")
        for remedy in cert.get("remedies", []):
            L.append(f"- {remedy}")
        L.append("")
        L.append(PROHIBITED_INFERENCES)
    elif negative:
        L.append(f"`{negative['condition']}` — {negative['detail']}")
    L.append("")
    if negative:
        L.append("### Negative result")
        L.append("")
        L.append(f"- **Condition:** `{negative['condition']}`")
        L.append(f"- **Settled negative:** {_fmt(negative.get('settled_negative'))}")
        L.append(f"- **Interpretation:** {negative.get('interpretation')}")
        L.append("- This licenses **no** method modification: `B` is not increased, the "
                 "correction is not loosened, the domain is not widened, no residue is "
                 "dropped and no re-seeding is performed.")
        L.append("")

    # ---- unresolved ---------------------------------------------------------
    L.append("## UNRESOLVED ISSUES")
    L.append("")
    escalations = s.get("escalations", [])
    if escalations:
        for e in escalations:
            L.append(f"- **{e['ambiguity']}**")
            for opt, cons in zip(e["options"], e["consequences"]):
                L.append(f"  - option: {opt} — consequence: {cons}")
            L.append(f"  - *recommendation:* {e['recommendation']}")
            L.append("  - *resolved by the agent:* NO — this is a Lead decision.")
    else:
        L.append("None.")
    L.append("")

    # ---- handoff ------------------------------------------------------------
    L.append("## HANDOFF")
    L.append("")
    if underpowered:
        L.append("`handoff_02.json` carries `terminal_state = UNDERPOWERED`, "
                 "`qc_status = FAIL` and **no** `negative_result` block — because this "
                 "is not a negative result. Every consumer contract blocks on "
                 "`qc_status = FAIL`, so the chain stops here without any downstream "
                 "stage having to interpret the science.")
    elif selection:
        L.append(f"`handoff_02.json` -> `footprint-robustness`: `r_hot = "
                 f"{_fmt(selection['r_hot'])} A`, "
                 f"`n_significant_centers = {detection.get('n_significant_bh')}`, "
                 f"`q = {s.get('q')}`, `fdr_method = BH`, `B = {s.get('B')}`, "
                 f"primary null `{s.get('primary_null')}`, power certificate "
                 f"`{'PASSED' if (s.get('power_certificate') or {}).get('passes') else 'FAILED'}`.")
    else:
        L.append("`handoff_02.json` carries a negative result: "
                 "`n_significant_centers = 0`, which terminates the chain before "
                 "footprint construction. The v2 §5.4 power certificate is attached as "
                 "the evidence that the design could have found something.")
    L.append("")
    L.append("The consumer must verify every hash, that each center has a usable CA, "
             "that all centers lie at or below the BH boundary, that "
             "`n_significant_centers > 0` and that `qc_status != FAIL`, and must "
             "propagate `PERMUTATION_RESOLUTION_LIMITED` into every downstream report. "
             "The center set, `r_hot` and everything under `04_*`, `05_*`, `06_*` and "
             "`11_*` are immutable downstream.")
    L.append("")
    L.append("**Stage B stops here.** Footprint radius, footprint geometry, robustness "
             "and biological annotation are other agents' work.")
    L.append("")
    L.append("### Derived seeds")
    L.append("")
    L.append(_table(["context", "seed"],
                    [[k, v] for k, v in (s.get("seeds") or {}).items()]))
    L.append("")
    return "\n".join(L)
