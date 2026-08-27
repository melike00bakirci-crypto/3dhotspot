"""Pareto -> normalization -> distance-to-ideal — ONE implementation, two callers.

METHOD_SPEC II.6 (r_hot) and II.10 (r_fp) specify *identical machinery* over
*different objective sets*. Implementing it once here makes divergence between the
two selections impossible; each caller supplies only its own objectives and tie
chain.

FROZEN (F15):
  * equal objective importance, ``w_k = 1``, declared in config, never
    outcome-dependent and never changed after results are observed;
  * deterministic min-max normalization computed over the **admissible** set;
  * degenerate objective (max == min) -> normalized value 1.0 for all, contributing
    0 to every distance, and recorded as non-discriminating;
  * utopia point = (1, 1, ..., 1) by construction;
  * L2 distance as the primary metric;
  * selection restricted to Pareto members of the admissible set.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np


@dataclass
class Candidate:
    """One scanned radius with its raw objective values and QC verdict."""
    key: float                                  # the radius (r_hot or r_fp)
    objectives: dict[str, float]
    admissible: bool
    qc_failure_reason: str = "NA"
    extra: dict = field(default_factory=dict)


@dataclass
class SelectionResult:
    selected_key: float | None
    normalized: dict[float, dict[str, float]]
    distances: dict[float, float]
    pareto_members: list[float]
    dominators: dict[float, list[float]]
    degenerate_objectives: list[str]
    utopia: dict[str, float]
    tie_chain_applied: list[str]
    near_tie: bool
    near_tie_detail: dict
    correlation: dict[str, dict[str, float]]
    n_admissible: int
    rejected: dict[float, str]


def _dominates(a: dict[str, float], b: dict[str, float], names: Sequence[str]) -> bool:
    """``a`` dominates ``b`` iff a >= b on every objective and > on at least one.

    All objectives are maximizing by construction (II.6 / II.10).
    """
    ge_all = all(a[n] >= b[n] for n in names)
    gt_any = any(a[n] > b[n] for n in names)
    return ge_all and gt_any


def pareto_front(candidates: Sequence[Candidate],
                 names: Sequence[str]) -> tuple[list[float], dict[float, list[float]]]:
    """Non-dominated filtering. Dominated candidates are retained with dominators."""
    members: list[float] = []
    dominators: dict[float, list[float]] = {}
    for c in candidates:
        doms = [o.key for o in candidates
                if o.key != c.key and _dominates(o.objectives, c.objectives, names)]
        if doms:
            dominators[c.key] = sorted(doms)
        else:
            members.append(c.key)
    return sorted(members), dominators


def min_max_normalize(candidates: Sequence[Candidate],
                      names: Sequence[str]) -> tuple[dict[float, dict[str, float]], list[str]]:
    """FROZEN min-max normalization over the ADMISSIBLE set (F15).

    Returns ``(normalized, degenerate_objective_names)``.
    """
    normalized: dict[float, dict[str, float]] = {c.key: {} for c in candidates}
    degenerate: list[str] = []
    for name in names:
        values = np.array([c.objectives[name] for c in candidates], dtype=np.float64)
        lo, hi = float(values.min()), float(values.max())
        if hi == lo:
            degenerate.append(name)
            for c in candidates:
                normalized[c.key][name] = 1.0
        else:
            for c in candidates:
                normalized[c.key][name] = float((c.objectives[name] - lo) / (hi - lo))
    return normalized, degenerate


def distance_to_ideal(normalized_row: dict[str, float], names: Sequence[str],
                      w_k: float = 1.0, metric: str = "L2") -> float:
    """``d = sqrt(sum_k w_k (1 - x~_k)^2)`` with the utopia point at all-ones."""
    gaps = np.array([1.0 - normalized_row[n] for n in names], dtype=np.float64)
    if metric == "L2":
        return float(np.sqrt(np.sum(w_k * gaps ** 2)))
    if metric == "L1":
        return float(np.sum(w_k * np.abs(gaps)))
    if metric == "Linf":
        return float(np.max(w_k * np.abs(gaps)))
    raise ValueError(f"unsupported distance metric: {metric!r}")


def objective_correlation(candidates: Sequence[Candidate],
                          names: Sequence[str]) -> dict[str, dict[str, float]]:
    """Pairwise Pearson correlation — the redundancy diagnostic (II.6)."""
    out: dict[str, dict[str, float]] = {a: {} for a in names}
    for a in names:
        va = np.array([c.objectives[a] for c in candidates], dtype=np.float64)
        for b in names:
            vb = np.array([c.objectives[b] for c in candidates], dtype=np.float64)
            if len(va) < 2 or np.std(va) == 0 or np.std(vb) == 0:
                out[a][b] = float("nan")
            else:
                out[a][b] = float(np.corrcoef(va, vb)[0, 1])
    return out


def select(candidates: Sequence[Candidate],
           objective_names: Sequence[str],
           tie_chain: Sequence[str],
           w_k: float = 1.0,
           metric: str = "L2",
           near_tie_rel: float = 0.05,
           near_tie_sep_steps: int = 2,
           step: float = 0.5) -> SelectionResult:
    """Run the frozen selection machinery end to end.

    ``tie_chain`` names objectives (after the implicit leading ``distance``) to break
    exact ties, always ending in ``smaller_radius`` — the radius is unique so the
    chain always terminates and no unresolvable tie can occur by construction.
    """
    admissible = [c for c in candidates if c.admissible]
    rejected = {c.key: c.qc_failure_reason for c in candidates if not c.admissible}

    if not admissible:
        return SelectionResult(
            selected_key=None, normalized={}, distances={}, pareto_members=[],
            dominators={}, degenerate_objectives=[], utopia={},
            tie_chain_applied=[], near_tie=False, near_tie_detail={},
            correlation={}, n_admissible=0, rejected=rejected,
        )

    normalized, degenerate = min_max_normalize(admissible, objective_names)
    members, dominators = pareto_front(admissible, objective_names)
    distances = {
        c.key: distance_to_ideal(normalized[c.key], objective_names, w_k, metric)
        for c in admissible
    }

    ranked = sorted(
        members,
        key=lambda k: _tie_key(k, distances, normalized, tie_chain, objective_names),
    )
    selected = ranked[0]

    near_tie, detail = False, {}
    if len(ranked) >= 2:
        d1, d2 = distances[ranked[0]], distances[ranked[1]]
        rel = abs(d1 - d2) / d1 if d1 > 0 else 0.0
        if rel < near_tie_rel and abs(ranked[0] - ranked[1]) > near_tie_sep_steps * step:
            near_tie = True
            detail = {"first": ranked[0], "second": ranked[1], "d_first": d1,
                      "d_second": d2, "relative_gap": rel}

    return SelectionResult(
        selected_key=selected,
        normalized=normalized,
        distances=distances,
        pareto_members=members,
        dominators=dominators,
        degenerate_objectives=degenerate,
        utopia={n: 1.0 for n in objective_names},
        # The chain always leads with `distance`; config chains already name it, so
        # de-duplicate rather than reporting "distance -> distance -> ...".
        tie_chain_applied=["distance", *[s for s in tie_chain if s != "distance"]],
        near_tie=near_tie,
        near_tie_detail=detail,
        correlation=objective_correlation(admissible, objective_names),
        n_admissible=len(admissible),
        rejected=rejected,
    )


def _tie_key(key: float, distances, normalized, tie_chain, objective_names):
    """Deterministic TOTAL ordering: distance, then the chain, then smaller radius."""
    parts: list[float] = [distances[key]]
    for step_name in tie_chain:
        if step_name in ("distance",):
            continue
        if step_name == "smaller_radius":
            parts.append(key)
        elif step_name in objective_names:
            parts.append(-normalized[key][step_name])   # higher is better
        else:
            raise ValueError(f"tie-chain entry {step_name!r} is not a known objective")
    parts.append(key)      # guarantees totality
    return tuple(parts)


def boundary_diagnostic(grid: Sequence[float], selected: float | None,
                        pareto_members: Sequence[float],
                        distances: dict[float, float],
                        band_fraction: float = 0.10,
                        bw2_fraction: float = 0.50,
                        bw3_top_n: int = 3) -> dict:
    """DOMAIN_BOUNDARY_WARNING — computed identically on both grids (II.6 / II.10).

    Non-blocking and never self-correcting: it reports the band, its members and the
    implied direction, then recommends a Lead-authorized re-run under an explicitly
    widened, pre-registered domain. Automatic domain expansion is prohibited.
    """
    grid = sorted(grid)
    G = len(grid)
    if G == 0:
        return {"DOMAIN_BOUNDARY_WARNING": False, "reason": "empty_grid"}
    b = max(1, int(np.ceil(band_fraction * G)))
    lower, upper = set(grid[:b]), set(grid[-b:])
    band = lower | upper

    bw1 = selected is not None and selected in band
    bw2 = bool(pareto_members) and (
        sum(1 for m in pareto_members if m in band) / len(pareto_members) >= bw2_fraction)
    top = sorted(distances, key=lambda k: (distances[k], k))[:bw3_top_n]
    bw3 = len(top) == bw3_top_n and (
        all(t in lower for t in top) or all(t in upper for t in top))

    which = []
    if selected is not None and selected in lower:
        which.append("lower")
    if selected is not None and selected in upper:
        which.append("upper")
    if bw3 and not which:
        which.append("lower" if all(t in lower for t in top) else "upper")

    return {
        "band_width_points": b, "band_members": sorted(band),
        "lower_band": sorted(lower), "upper_band": sorted(upper),
        "BW1_selected_in_band": bool(bw1),
        "BW2_pareto_majority_in_band": bool(bw2),
        "BW3_top3_same_band": bool(bw3),
        "DOMAIN_BOUNDARY_WARNING": bool(bw1 or bw2 or bw3),
        "which_end": which or ["none"],
        "top_candidates": top,
        "automatic_expansion": False,
        "recommended_action": (
            "Lead decision: consider a re-run with an explicitly widened, "
            "pre-registered domain under a NEW RUN_ID. Automatic domain expansion "
            "is prohibited (it would be domain shopping)."),
    }
