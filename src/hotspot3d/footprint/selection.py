"""II.10 selection — a thin adapter onto the ONE shared decision machinery.

Pareto filtering, min-max normalization over the admissible set, the all-ones
utopia point, ``w_k = 1`` L2 distance-to-ideal, the tie chain and the boundary
diagnostic all live in :mod:`hotspot3d.utils.multiobjective`, which II.6 (``r_hot``)
and II.10 (``r_fp``) share. This module only supplies the *geometric* objective
vector and reports which edge of the II.9 domain binds. Nothing is recomputed here.

**[v2 §4/§7] Objective/constraint separation diagnostics**, mirrored from
``hotspot/selection.py``'s ``run_selection``/``RadiusSelection`` pattern and
adapted to footprint's own four objectives and its own QC_F1/F2/F3 rule names.
Footprint admissibility (:mod:`hotspot3d.footprint.qc`) is already purely
geometric — coverage, bridging, undefined geometry — and was never conditioned on
significant-center count, since ``r_fp`` operates on the already-frozen
``SIGNIFICANT_HOTSPOT_CENTERS`` set and never re-tests significance. What v2 §7
adds on top is the DIAGNOSTIC layer: which objectives are degenerate and excluded
from distance-to-utopia (numerically verified, not assumed), whether the Pareto
set is vacuous, and whether admissibility rather than the objectives decided the
outcome.

``diagnose_selection`` is deliberately a POST-HOC function over the already
returned :class:`~hotspot3d.utils.multiobjective.SelectionResult` and the QC
verdicts, not a change to ``run_selection`` or :func:`build_footprint
<hotspot3d.footprint.api.build_footprint>` themselves. ``build_footprint`` is the
one frozen API Phase D calls on every perturbation iteration (tens of thousands of
times); keeping this diagnostic layer outside it means Phase D's per-iteration
reconstruction pays nothing for it and emits none of these warnings — only Phase
C's one-shot ``r_fp`` selection does.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import numpy as np

from ..utils.multiobjective import (Candidate, SelectionResult, boundary_diagnostic,
                                    distance_to_ideal, select)
from .domain import FootprintDomain, binding_constraint_detail
from .params import FootprintParams
from .qc import QC_F1, QC_F2, QC_F3, QCVerdict
from .sweep import RadiusState

#: The three FROZEN, enumerable QC_F rule codes (II.10). None is a substring of
#: another, so presence in ``reason_text`` is checked directly rather than by
#: splitting on ';' — QC_F3's own reason text embeds an internal ';'
#: (``volume=...;mesh=...``), which would otherwise misparse.
_QC_F_RULE_CODES = (QC_F1, QC_F2, QC_F3)


def objective_vector(state: RadiusState) -> dict[str, float]:
    """The four FROZEN maximizing objectives (II.10).

    ``parsimony = 1 - coverage`` is the anti-inflation term: compactness,
    connectivity and stability would all happily reward a protein-sized blob.
    """
    return {
        "compactness": _finite(state.compactness),
        "connectivity": _finite(state.connectivity),
        "stability": _finite(state.stability),
        "parsimony": _finite(state.parsimony),
    }


def _finite(value: float) -> float:
    """A non-finite objective is the worst possible value, never a free pass."""
    return float(value) if np.isfinite(value) else 0.0


def build_candidates(rows: list[RadiusState],
                     verdicts: dict[float, QCVerdict]) -> list[Candidate]:
    return [
        Candidate(
            key=row.rho,
            objectives=objective_vector(row),
            admissible=verdicts[row.rho].admissible,
            qc_failure_reason=verdicts[row.rho].reason_text,
            extra={"n_components": row.n_components, "volume_A3": row.volume},
        )
        for row in rows
    ]


def run_selection(candidates: list[Candidate],
                  params: FootprintParams, step: float) -> SelectionResult:
    return select(
        candidates,
        objective_names=params.objective_names,
        tie_chain=params.tie_chain,
        w_k=params.w_k,
        metric=params.distance_metric,
        near_tie_rel=params.near_tie_relative_threshold,
        near_tie_sep_steps=params.near_tie_radius_separation_steps,
        step=step,
    )


def run_boundary_diagnostic(domain: FootprintDomain, result: SelectionResult,
                            params: FootprintParams) -> dict:
    """BW-1/BW-2/BW-3 over the full ``[rho_min, rho_max]`` grid, plus which edge binds."""
    diagnostic = boundary_diagnostic(
        grid=list(domain.grid),
        selected=result.selected_key,
        pareto_members=result.pareto_members,
        distances=result.distances,
        band_fraction=params.band_fraction,
        bw2_fraction=params.bw2_pareto_fraction,
        bw3_top_n=params.bw3_top_n,
    )
    diagnostic["domain_constraints"] = binding_constraint_detail(domain, params)
    diagnostic["binding_edge"] = _binding_edge(diagnostic, domain)
    return diagnostic


def _binding_edge(diagnostic: dict, domain: FootprintDomain) -> str:
    """Name the constraint responsible for the band the optimum sits in."""
    if not diagnostic.get("DOMAIN_BOUNDARY_WARNING"):
        return "NA"
    ends = diagnostic.get("which_end", ["none"])
    if "lower" in ends:
        return "rho_min"
    if "upper" in ends:
        return domain.binding_constraint
    return "NA"


def redundancy_flags(result: SelectionResult, params: FootprintParams) -> list[dict]:
    """Objective pairs whose |rho| exceeds the configured redundancy threshold."""
    out: list[dict] = []
    names = list(params.objective_names)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            value = result.correlation.get(a, {}).get(b, float("nan"))
            if np.isfinite(value) and abs(value) > params.objective_redundancy_abs_rho:
                out.append({"objective_a": a, "objective_b": b,
                            "pearson_r": float(value),
                            "threshold": params.objective_redundancy_abs_rho})
    return out


@dataclass
class SelectionDiagnostics:
    """v2 §4/§7 objective/constraint-separation diagnostics for ``r_fp`` selection."""

    effective_objectives: list[str] = field(default_factory=list)
    degenerate_objectives: list[str] = field(default_factory=list)
    degenerate_exclusion_is_numerically_neutral: bool = True
    vacuous_pareto: bool = False
    vacuous_reason: str = "NA"
    admissibility_dominates: bool = False
    inadmissible_fraction: float = 0.0
    inadmissible_by_rule: dict[str, int] = field(default_factory=dict)


def _degenerate_exclusion_is_neutral(result: SelectionResult, effective: list[str],
                                     w_k: float, metric: str) -> bool:
    """Assert that dropping degenerate objectives changes no distance (v2 §4).

    The frozen min-max rule maps a degenerate objective to 1.0 for every candidate,
    so its contribution to ``d = sqrt(sum_k w_k (1 - x~_k)^2)`` is exactly 0. If that
    identity ever failed, excluding the objective would silently change the selected
    radius — so it is verified rather than trusted (same pattern as Stage B).
    """
    if not effective:
        return all(abs(d) <= 1e-12 for d in result.distances.values())
    for key, row in result.normalized.items():
        recomputed = distance_to_ideal(row, effective, w_k, metric)
        if abs(recomputed - result.distances[key]) > 1e-12:
            return False
    return True


def diagnose_selection(result: SelectionResult, verdicts: dict[float, QCVerdict],
                       params: FootprintParams) -> SelectionDiagnostics:
    """Compute the v2 §7 diagnostics for one completed ``r_fp`` selection.

    Deliberately a function of the already-computed ``SelectionResult`` and QC
    verdicts, never of ``run_selection`` internals or of ``build_footprint`` itself
    — see the module docstring for why this must stay outside the frozen API that
    Phase D re-executes per iteration.
    """
    names = list(params.objective_names)

    # --- degenerate objectives are non-informative --------------------------
    degenerate = list(result.degenerate_objectives)
    effective = ([n for n in names if n not in degenerate]
                if params.exclude_degenerate_objectives else list(names))
    neutral = _degenerate_exclusion_is_neutral(result, effective, params.w_k,
                                               params.distance_metric)

    # --- was this a selection, or the only survivor? -------------------------
    vacuous, vacuous_reason = False, "NA"
    if result.n_admissible and len(result.pareto_members) == 1:
        vacuous = True
        vacuous_reason = (
            f"the Pareto set has exactly one member ({result.pareto_members[0]:g} A) "
            f"over {result.n_admissible} admissible radii: r_fp was not selected by "
            f"trading objectives off, it was the only non-dominated survivor")
    if result.n_admissible and len(degenerate) == len(names):
        reason = (f"every objective ({', '.join(names)}) is degenerate across the "
                  f"admissible set, so distance-to-utopia is identical for every "
                  f"candidate and the tie chain alone decided")
        vacuous = True
        vacuous_reason = reason if vacuous_reason == "NA" else f"{vacuous_reason}; {reason}"

    # --- did admissibility, rather than the objectives, decide? --------------
    n_total = len(verdicts)
    n_inadmissible = sum(1 for v in verdicts.values() if not v.admissible)
    fraction = (n_inadmissible / n_total) if n_total else 0.0
    # reason_text embeds numeric detail and, for QC_F3, an internal ';'
    # ("QC_F3_GEOMETRY_UNAVAILABLE:volume=0;mesh=..."), so rule codes are matched
    # by direct substring check against the three known QC_F codes rather than by
    # splitting the string — "which constraint did the eliminating" needs the rule
    # code alone, and a naive ';'/':' split would misparse QC_F3 and make every
    # key near-unique from the embedded numeric detail.
    by_rule = Counter(code for v in verdicts.values() if not v.admissible
                      for code in _QC_F_RULE_CODES if code in v.reason_text)

    return SelectionDiagnostics(
        effective_objectives=effective,
        degenerate_objectives=degenerate,
        degenerate_exclusion_is_numerically_neutral=neutral,
        vacuous_pareto=vacuous, vacuous_reason=vacuous_reason,
        admissibility_dominates=bool(fraction > params.admissibility_dominates_fraction),
        inadmissible_fraction=float(fraction),
        inadmissible_by_rule={k: int(v) for k, v in sorted(by_rule.items())},
    )
