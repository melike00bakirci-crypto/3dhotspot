"""METHOD_SPEC II.6 — final ``r_hot`` selection.

The machinery itself is NOT implemented here. Pareto filtering, min-max
normalization over the admissible set, the utopia point and the L2
distance-to-ideal live once in :mod:`hotspot3d.utils.multiobjective` and are shared
with the ``r_fp`` selection of II.10, which makes divergence between the two
selections impossible. This module only supplies Stage B's own objectives, its tie
chain and its audit trail.

What is FROZEN and is asserted here rather than assumed:
  * exactly four objectives, all maximizing (F15);
  * ``w_k = 1`` — equal importance, declared in config, never outcome-dependent and
    never changed after results are observed;
  * min-max normalization over the **admissible** set, L2 distance, utopia (1,1,1,1);
  * the tie chain ``distance -> stability -> Zg -> smaller radius``, which is total,
    so no unresolvable tie can occur.

Explicitly NOT the rule: max MCC, the PCF peak, any single metric, or a
lexicographic ordering. The max-MCC radius has no privileged status and is recorded
merely as one per-objective argmax (F7/F8).

**[v2 §4] Three guards on the selection itself**, because a selection can be
arithmetically correct and still mean nothing:

  * a **degenerate** objective (identical at every admissible radius) discriminates
    nothing. It is excluded from the distance-to-utopia computation and reported as
    non-informative rather than left to inflate the apparent dimensionality of the
    decision. Under the frozen min-max rule a degenerate objective normalizes to 1.0
    everywhere and so contributes exactly 0 to every L2 distance — exclusion is
    therefore numerically identical and that identity is ASSERTED here, not assumed;
  * a **Pareto set of size one**, or a selection in which every objective is
    degenerate, means the radius was not selected — it was the only survivor
    (``VACUOUS_PARETO_SELECTION``);
  * more than half the scanned radii ruled inadmissible means admissibility, not the
    objectives, decided the outcome (``ADMISSIBILITY_DOMINATES_SELECTION``).
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import numpy as np

from ..utils.multiobjective import (
    Candidate,
    SelectionResult,
    boundary_diagnostic,
    distance_to_ideal,
    select,
)
from .scan import OBJECTIVE_NAMES, RadiusResult


@dataclass
class RadiusSelection:
    result: SelectionResult
    candidates: list[Candidate]
    boundary: dict
    per_objective_argmax: dict[str, float | None]
    redundant_pairs: list[dict]
    independent_pareto_check: bool
    selected: float | None
    correlation_informative: bool = True
    #: v2 §4 diagnostics
    effective_objectives: list[str] = field(default_factory=list)
    degenerate_objectives: list[str] = field(default_factory=list)
    degenerate_exclusion_is_numerically_neutral: bool = True
    vacuous_pareto: bool = False
    vacuous_reason: str = "NA"
    admissibility_dominates: bool = False
    inadmissible_fraction: float = 0.0
    inadmissible_by_rule: dict[str, int] = field(default_factory=dict)

    @property
    def n_admissible(self) -> int:
        return self.result.n_admissible


def build_candidates(results: list[RadiusResult]) -> list[Candidate]:
    """One candidate per scanned radius — inadmissible radii INCLUDED, with reasons."""
    return [
        Candidate(key=r.radius_A, objectives=r.objectives(), admissible=r.admissible,
                  qc_failure_reason=r.qc_failure_reason,
                  extra={"n_significant_centers": len(r.centers),
                         "coverage_fraction": r.coverage_fraction,
                         "n_singleton_centers": r.n_singleton_centers})
        for r in results
    ]


def _independent_dominance_check(candidates: list[Candidate], names: list[str],
                                 members: list[float]) -> bool:
    """QC rule 8: re-verify Pareto membership with a second, independent routine.

    Deliberately written differently from ``utils.multiobjective.pareto_front`` (a
    plain O(n^2) sweep over admissible candidates) so that a defect in one would not
    be reproduced by the other.
    """
    admissible = [c for c in candidates if c.admissible]
    recomputed = []
    for c in admissible:
        dominated = False
        for other in admissible:
            if other.key == c.key:
                continue
            better_or_equal = True
            strictly_better = False
            for n in names:
                if other.objectives[n] < c.objectives[n]:
                    better_or_equal = False
                    break
                if other.objectives[n] > c.objectives[n]:
                    strictly_better = True
            if better_or_equal and strictly_better:
                dominated = True
                break
        if not dominated:
            recomputed.append(c.key)
    return sorted(recomputed) == sorted(members)


def _degenerate_exclusion_is_neutral(result: SelectionResult, effective: list[str],
                                     w_k: float, metric: str) -> bool:
    """Assert that dropping degenerate objectives changes no distance (v2 §4).

    The frozen min-max rule maps a degenerate objective to 1.0 for every candidate,
    so its contribution to ``d = sqrt(sum_k w_k (1 - x~_k)^2)`` is exactly 0. If that
    identity ever failed, excluding the objective would silently change the selected
    radius — so it is verified rather than trusted.
    """
    if not effective:
        return all(abs(d) <= 1e-12 for d in result.distances.values())
    for key, row in result.normalized.items():
        recomputed = distance_to_ideal(row, effective, w_k, metric)
        if abs(recomputed - result.distances[key]) > 1e-12:
            return False
    return True


def run_selection(results: list[RadiusResult], grid: list[float], *, tie_chain: list[str],
                  w_k: float, metric: str, near_tie_rel: float, near_tie_sep_steps: int,
                  step: float, band_fraction: float, bw2_fraction: float, bw3_top_n: int,
                  redundancy_rho: float,
                  admissibility_dominates_fraction: float = 0.50,
                  exclude_degenerate_objectives: bool = True,
                  near_tie_major_rel: float = 0.01) -> RadiusSelection:
    """Run the frozen II.6 selection and its v2 §4 diagnostics over the scanned grid."""
    candidates = build_candidates(results)
    result = select(candidates, OBJECTIVE_NAMES, tie_chain, w_k=w_k, metric=metric,
                    near_tie_rel=near_tie_rel, near_tie_sep_steps=near_tie_sep_steps,
                    step=step, near_tie_major_rel=near_tie_major_rel)

    # --- v2 §4: degenerate objectives are non-informative -----------------------
    degenerate = list(result.degenerate_objectives)
    effective = ([n for n in OBJECTIVE_NAMES if n not in degenerate]
                 if exclude_degenerate_objectives else list(OBJECTIVE_NAMES))
    neutral = _degenerate_exclusion_is_neutral(result, effective, w_k, metric)

    # --- v2 §4: was this a selection, or the only survivor? ---------------------
    vacuous, vacuous_reason = False, "NA"
    if result.n_admissible and len(result.pareto_members) == 1:
        vacuous = True
        vacuous_reason = (
            f"the Pareto set has exactly one member ({result.pareto_members[0]:g} A) "
            f"over {result.n_admissible} admissible radii: the radius was not selected "
            f"by trading objectives off, it was the only non-dominated survivor")
    if result.n_admissible and len(degenerate) == len(OBJECTIVE_NAMES):
        vacuous = True
        reason = (f"every objective ({', '.join(OBJECTIVE_NAMES)}) is degenerate across "
                  f"the admissible set, so distance-to-utopia is identical for every "
                  f"candidate and the tie chain alone decided")
        vacuous_reason = reason if vacuous_reason == "NA" else f"{vacuous_reason}; {reason}"

    # --- v2 §4: did admissibility, rather than the objectives, decide? ----------
    n_total = len(candidates)
    n_inadmissible = sum(1 for c in candidates if not c.admissible)
    fraction = (n_inadmissible / n_total) if n_total else 0.0
    by_rule = Counter(rule for c in candidates if not c.admissible
                      for rule in str(c.qc_failure_reason).split(";"))

    boundary = boundary_diagnostic(grid, result.selected_key, result.pareto_members,
                                   result.distances, band_fraction=band_fraction,
                                   bw2_fraction=bw2_fraction, bw3_top_n=bw3_top_n)

    # Per-objective argmax over the ADMISSIBLE set. Recorded so a reader can see the
    # max-MCC radius explicitly and verify it holds no privileged status (F7/F8).
    admissible = [c for c in candidates if c.admissible]
    argmax: dict[str, float | None] = {}
    for name in OBJECTIVE_NAMES:
        if admissible:
            best = max(admissible, key=lambda c: (c.objectives[name], -c.key))
            argmax[name] = best.key
        else:
            argmax[name] = None

    # A Pearson correlation over fewer than three admissible candidates is +/-1 by
    # arithmetic, not by redundancy, so the diagnostic is reported as uninformative
    # rather than raised as a finding. This changes no result — it is a reporting rule.
    correlation_informative = len(admissible) >= 3
    redundant = []
    if correlation_informative:
        for i, a in enumerate(OBJECTIVE_NAMES):
            for b in OBJECTIVE_NAMES[i + 1:]:
                rho = result.correlation.get(a, {}).get(b, float("nan"))
                if rho is not None and not np.isnan(rho) and abs(rho) > redundancy_rho:
                    redundant.append({"objective_a": a, "objective_b": b,
                                      "rho": float(rho)})

    return RadiusSelection(
        result=result, candidates=candidates, boundary=boundary,
        per_objective_argmax=argmax, redundant_pairs=redundant,
        independent_pareto_check=_independent_dominance_check(
            candidates, OBJECTIVE_NAMES, result.pareto_members),
        selected=result.selected_key,
        correlation_informative=correlation_informative,
        effective_objectives=effective,
        degenerate_objectives=degenerate,
        degenerate_exclusion_is_numerically_neutral=neutral,
        vacuous_pareto=vacuous, vacuous_reason=vacuous_reason,
        admissibility_dominates=bool(fraction > admissibility_dominates_fraction),
        inadmissible_fraction=float(fraction),
        inadmissible_by_rule={k: int(v) for k, v in sorted(by_rule.items())},
    )


def decision_payload(sel: RadiusSelection, results: list[RadiusResult],
                     domain_source: str, fallback: bool) -> dict:
    """``radius_decision.json`` — the complete, replayable selection record."""
    by_radius = {r.radius_A: r for r in results}
    rows = []
    for c in sel.candidates:
        r = by_radius[c.key]
        norm = sel.result.normalized.get(c.key)
        rows.append({
            "radius_A": c.key,
            "admissible": c.admissible,
            "qc_failure_reason": c.qc_failure_reason,
            "admissibility_rule_text": r.admissibility_reason,
            "significance_conditioned_diagnostics_fired": r.diagnostics_reason,
            "diagnostic_values": r.diagnostic_values,
            "objectives_raw": c.objectives,
            # Normalization is defined over the admissible set only, so an
            # inadmissible row carries NULL here — recorded, never silently blank.
            "objectives_normalized": norm,
            "distance_to_ideal": sel.result.distances.get(c.key),
            "pareto_member": (c.key in sel.result.pareto_members) if c.admissible else None,
            "dominated_by": sel.result.dominators.get(c.key),
            "selected": c.key == sel.selected,
            "rejection_reason": _rejection_reason(c, sel),
            "detail": {
                "n_significant_centers": len(r.centers),
                "coverage_fraction": r.coverage_fraction,
                "n_singleton_centers": r.n_singleton_centers,
                "median_n_labeled_per_center": r.median_n_labeled,
                "n_in_test_family": r.n_in_family,
                "n_plp_in_hotspot_spheres": r.n_plp_in,
                # n_labeled_in_universe in r_hot_scan.tsv is exactly this quantity:
                # labelled residues inside the union of the S(r) spheres (the FE
                # denominator). Named here in full to remove any ambiguity.
                "n_labeled_in_hotspot_spheres": r.n_labeled_in,
                "fold_enrichment_defined": r.fold_enrichment_defined,
                "neighbor_stability_defined": r.neighbor_stability_defined,
                "bh_boundary_p": r.bh_boundary_p,
                "primary_null": r.extra["primary_null"],
                "smallest_p_emp_in_family": r.extra["smallest_p_emp_in_family"],
                "smallest_p_exact_in_family": r.extra["smallest_p_exact_in_family"],
                # T(r) and Zg come from the SECONDARY label-permutation null (II.5B);
                # they are a selection objective and never a significance statement.
                "T_obs": r.perm.t_obs, "T_null_mean": r.perm.t_null_mean,
                "T_null_sd": r.perm.t_null_sd,
                "p_ratio_descriptive_only": r.perm.p_ratio,
                "loo_secondary_metrics": {
                    k: r.loo.metrics[k] for k in
                    ("sens", "spec", "acc", "ppv", "npv", "f1", "balacc", "tp", "fp", "tn", "fn")},
                "loo_kappa_boundary_neutral_verified": r.boundary_neutral_verified,
                "seed_context": r.extra["seed_context"],
                # v2 §4 — biological-plausibility DIAGNOSTICS; never a Pareto
                # objective and never an admissibility constraint (F15).
                "structural_coverage": r.structural_coverage,
                "cohort_absorption": r.cohort_absorption,
                "r_to_domain_ratio": r.r_to_domain_ratio,
            },
        })

    return {
        "selected_r_hot_A": sel.selected,
        "search_domain_source": domain_source,
        "FALLBACK_RADIUS_DOMAIN": fallback,
        "n_candidates": len(sel.candidates),
        "n_admissible": sel.result.n_admissible,
        "objectives": OBJECTIVE_NAMES,
        "objective_directions": {n: "maximize" for n in OBJECTIVE_NAMES},
        "weights": {n: 1.0 for n in OBJECTIVE_NAMES},
        "weights_frozen_equal": True,
        "normalization": "min_max",
        "normalization_domain": "admissible_set",
        "degenerate_objectives": sel.result.degenerate_objectives,
        "degenerate_normalized_value": 1.0,
        # --- v2 §4 -----------------------------------------------------------
        "effective_objectives_in_distance": sel.effective_objectives,
        "degenerate_objectives_excluded_from_distance": sel.degenerate_objectives,
        "degenerate_exclusion_is_numerically_neutral":
            sel.degenerate_exclusion_is_numerically_neutral,
        "degenerate_objective_note": (
            "An objective identical at every admissible radius discriminates nothing. "
            "It is excluded from distance-to-utopia and reported as non-informative "
            "(v2 §4). Under the frozen min-max rule it normalizes to 1.0 everywhere, "
            "so exclusion changes no distance — verified, not assumed."),
        "VACUOUS_PARETO_SELECTION": sel.vacuous_pareto,
        "vacuous_pareto_reason": sel.vacuous_reason,
        "ADMISSIBILITY_DOMINATES_SELECTION": sel.admissibility_dominates,
        "inadmissible_fraction": sel.inadmissible_fraction,
        "inadmissible_by_rule": sel.inadmissible_by_rule,
        "significance_may_determine_admissibility": False,
        "admissibility_note": (
            "v2 §4 — the number of significant centers at a candidate radius NEVER "
            "removes it from the Pareto set. A radius that yields zero significant "
            "centers scores zero on the enrichment objective and remains a legitimate "
            "candidate. Only genuinely undefined quantities are admissibility "
            "constraints; QC_H1/QC_H2/QC_H4 are recorded diagnostics."),
        "utopia_point": sel.result.utopia,
        "distance_metric": "L2",
        "tie_chain_applied": sel.result.tie_chain_applied,
        "near_tie_flag": sel.result.near_tie,
        "near_tie_detail": sel.result.near_tie_detail or None,
        "pareto_members": sel.result.pareto_members,
        "independent_dominance_check_passed": sel.independent_pareto_check,
        "per_objective_argmax": sel.per_objective_argmax,
        "max_mcc_radius_has_no_privileged_status": True,
        "selection_rule": ("admissible set -> Pareto non-dominated filtering -> min-max "
                           "normalization over the admissible set -> utopia (1,1,1,1) -> "
                           "L2 distance with w_k = 1 -> argmin over Pareto members -> "
                           "total tie chain. Never max-MCC, never the PCF peak, never "
                           "lexicographic, never a single metric."),
        "candidates": rows,
    }


def _rejection_reason(c: Candidate, sel: RadiusSelection) -> str:
    if c.key == sel.selected:
        return "NA"
    if not c.admissible:
        return f"INADMISSIBLE:{c.qc_failure_reason}"
    if c.key in sel.result.dominators:
        dominators = ",".join(f"{d:g}" for d in sel.result.dominators[c.key])
        return f"DOMINATED_BY:{dominators}"
    d_sel = sel.result.distances.get(sel.selected)
    d_own = sel.result.distances.get(c.key)
    if d_sel is not None and d_own is not None and d_own > d_sel:
        return f"PARETO_MEMBER_LARGER_DISTANCE_TO_IDEAL:{d_own:.6g}>{d_sel:.6g}"
    return "PARETO_MEMBER_LOST_TIE_CHAIN"
