"""II.9 — the candidate ``r_fp`` search domain, derived from center geometry alone.

    merge scales : {w_k / 2} over the edges of the Euclidean MST of S
    rho_all      : max_k w_k / 2
[REVISED — DECISION-RFP-CONTAINMENT-0001] ``r_fp`` is no longer searched for; it is
DETERMINED by containment:

    r_fp = max over the evidence set of (distance to nearest significant center),
           rounded UP to the step grid

The evidence set arrives as an :class:`~hotspot3d.footprint.inputs.EvidenceSet` —
per-residue distances and an ``in_evidence`` verdict. It carries no radius.

``r_hot`` STILL plays no part in this module: it is not an argument to any function
here, which is the mechanical form of the guarantee this docstring has always made.
The ``r_hot`` comparison lives in exactly one place, the single sanctioned
``ReadOnlyHotspotRadius.derive_evidence_set``, and its numeric value never leaves
that method.

``derive_domain`` and its search envelope are retained to supply ``D_max`` and the
MST merge scales, which the decision keeps as reported events. They bound nothing.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..utils.errors import NegativeResult
from ..utils.geometry import mst_edges, protein_diameter
from .params import FootprintParams

# Grid values are snapped to this many decimals so that a radius is a stable dict
# key across the sweep, the selection and Phase D.
_GRID_DECIMALS = 6


@dataclass(frozen=True)
class FootprintDomain:
    """The derived domain plus every quantity that produced it."""

    rho_min: float
    rho_max: float
    rho_all: float | None
    step: float
    grid: tuple[float, ...]
    d_max: float
    mst: tuple[tuple[int, int, float], ...]
    merge_scales: tuple[float, ...]
    single_center_special_case: bool
    binding_constraint: str
    rho_max_candidates: tuple[tuple[str, float], ...]

    def as_dict(self, center_ids: list[int] | None = None) -> dict:
        edges = []
        for i, j, w in self.mst:
            edge = {"i": int(i), "j": int(j), "weight_A": float(w),
                    "merge_scale_A": float(w) / 2.0}
            if center_ids is not None:
                edge["center_i"] = int(center_ids[i])
                edge["center_j"] = int(center_ids[j])
            edges.append(edge)
        return {
            "rho_min_A": self.rho_min,
            "rho_max_A": self.rho_max,
            "rho_all_A": self.rho_all,
            "step_fp_A": self.step,
            "n_grid_points": len(self.grid),
            "grid_A": list(self.grid),
            "D_max_A": self.d_max,
            "mst_edges": edges,
            "merge_scales_A": list(self.merge_scales),
            "single_center_special_case": self.single_center_special_case,
            "rho_max_binding_constraint": self.binding_constraint,
            "rho_max_candidates_A": {k: v for k, v in self.rho_max_candidates},
            "note": ("[DECISION-RFP-CONTAINMENT-0001] r_fp is DETERMINED by "
                     "containment, so rho_min == rho_max == r_fp and the grid has one "
                     "member. Merge scales are reported as events and determine "
                     "nothing. r_hot is not an argument to any function in this "
                     "module; it reaches the stage only as the per-residue boolean in "
                     "the evidence set. See the `containment` block for the "
                     "determination itself."),
        }


def derive_domain(center_coords: np.ndarray, universe_coords: np.ndarray,
                  params: FootprintParams) -> FootprintDomain:
    """Derive the II.9 domain. Raises :class:`NegativeResult` when it is empty."""
    centers = np.atleast_2d(np.asarray(center_coords, dtype=np.float64))
    n_s = len(centers)

    d_max = float(protein_diameter(np.asarray(universe_coords, dtype=np.float64)))
    cap_dmax = params.dmax_fraction_cap * d_max
    ceiling = params.hard_ceiling_A

    edges = tuple((int(i), int(j), float(w)) for i, j, w in mst_edges(centers))
    merge_scales = tuple(round(w / 2.0, _GRID_DECIMALS) for _, _, w in edges)

    candidates: dict[str, float] = {
        "0.25*D_max": cap_dmax,
        "hard_ceiling": ceiling,
    }
    if n_s == 1:
        # |S| = 1 -> MST empty -> rho_all undefined; the post-merge term drops out.
        rho_all = None
        single_case = True
    else:
        rho_all = float(max(merge_scales))
        # [REVISED — DECISION-FOOTPRINT-DOMAIN-0001] rho_all is still computed and
        # still reported as merge events; it no longer CAPS rho_max. It measures how
        # tightly the centers are packed, which is not what r_fp is for: with
        # sequence-adjacent centers 1.25*rho_all falls below the floor and empties
        # the domain.
        if params.post_merge_factor_caps_rho_max:
            candidates["1.25*rho_all"] = params.post_merge_factor * rho_all
        single_case = False

    binding = min(candidates, key=lambda k: (candidates[k], k))
    rho_max = float(min(candidates.values()))
    rho_min = float(params.rho_min_A)

    if rho_min > rho_max:
        raise NegativeResult(
            "NO_ADMISSIBLE_FOOTPRINT_DOMAIN",
            f"rho_min ({rho_min:.3f} A) exceeds rho_max ({rho_max:.3f} A); the "
            f"II.9 domain is empty, so no footprint radius can be evaluated. The "
            f"binding upper constraint was {binding}.",
            rho_min=rho_min, rho_max=rho_max, binding_constraint=binding,
            d_max=d_max, rho_all=rho_all, n_centers=int(n_s),
        )

    grid = build_grid(rho_min, rho_max, params.step_fp_A)
    if not grid:
        raise NegativeResult(
            "NO_ADMISSIBLE_FOOTPRINT_DOMAIN",
            f"the uniform {params.step_fp_A} A grid over [{rho_min}, {rho_max}] is "
            f"empty.",
            rho_min=rho_min, rho_max=rho_max, binding_constraint=binding,
        )

    return FootprintDomain(
        rho_min=rho_min, rho_max=rho_max, rho_all=rho_all,
        step=float(params.step_fp_A), grid=grid, d_max=d_max,
        mst=edges, merge_scales=merge_scales,
        single_center_special_case=single_case, binding_constraint=binding,
        rho_max_candidates=tuple(sorted((k, float(v)) for k, v in candidates.items())),
    )


def build_grid(rho_min: float, rho_max: float, step: float) -> tuple[float, ...]:
    """Uniform ascending grid on ``[rho_min, rho_max]`` — never adaptive (II.9).

    Built by integer multiples of ``step`` from ``rho_min`` so that Phase D
    reproduces Phase C's grid exactly rather than approximately.
    """
    if step <= 0:
        raise ValueError("step_fp must be positive")
    n = int(np.floor(round((rho_max - rho_min) / step, 9))) + 1
    values = [round(rho_min + i * step, _GRID_DECIMALS) for i in range(max(n, 0))]
    return tuple(v for v in values if v <= rho_max + 1e-9)


def binding_constraint_detail(domain: FootprintDomain,
                              params: FootprintParams) -> dict:
    """Which edge of the domain binds — required by the boundary diagnostic."""
    detail = {
        "rho_min_A": domain.rho_min,
        "cap_0.25_D_max_A": params.dmax_fraction_cap * domain.d_max,
        "hard_ceiling_A": params.hard_ceiling_A,
        "post_merge_1.25_rho_all_A": (
            params.post_merge_factor * domain.rho_all
            if domain.rho_all is not None else None),
        "binding_upper_constraint": domain.binding_constraint,
        "single_center_special_case": domain.single_center_special_case,
    }
    return detail


def derive_containment_domain(center_coords: np.ndarray, universe_coords: np.ndarray,
                              evidence, params: FootprintParams
                              ) -> tuple[FootprintDomain, dict]:
    """DECISION-RFP-CONTAINMENT-0001 — ``r_fp`` as the containment minimum.

    Takes an ``EvidenceSet``, never a radius. Returns a ONE-POINT domain at ``r_fp``
    plus the record that the reporting layer consumes. The grid has a single member
    by construction: the radius is the unique solution of a constraint (L4), so there
    is nothing to select among. Downstream QC, export and reporting therefore see a
    sweep of exactly one candidate and keep working unchanged — the geometric
    descriptors are computed and reported, and cannot decide anything because no
    alternative exists to prefer.

    Raises :class:`NegativeResult` when the evidence set is empty. By L1 that cannot
    happen for a non-empty center set: a center is significant because its sphere is
    P/LP-enriched, and a sphere holding no P/LP cannot reject. The guard is here so a
    violation of the lemma fails loudly instead of yielding a degenerate radius.
    """
    centers = np.atleast_2d(np.asarray(center_coords, dtype=np.float64))
    universe = np.asarray(universe_coords, dtype=np.float64)
    d_max = float(protein_diameter(universe))

    # Merge scales are still computed: the decision keeps them as reported events.
    edges = tuple((int(i), int(j), float(w)) for i, j, w in mst_edges(centers))
    merge_scales = tuple(round(w / 2.0, _GRID_DECIMALS) for _, _, w in edges)
    rho_all = float(max(merge_scales)) if merge_scales else None

    inside = list(evidence.evidence)
    if not inside:
        raise NegativeResult(
            "NO_EVIDENCE_SET",
            "no P/LP residue lies inside the detection neighbourhood of any "
            "significant center, so the containment radius is undefined. L1 of "
            "DECISION-RFP-CONTAINMENT-0001 says this cannot occur for a non-empty "
            "center set; reaching it means the significance decision and the "
            "classified cohort disagree.",
            n_centers=int(len(centers)), n_plp_mapped=int(evidence.n_plp_mapped),
        )

    binding = max(inside, key=lambda r: r.distance_A)
    r_contain_exact = float(binding.distance_A)

    # Rounded UP to the step grid. The direction is a correctness condition, not a
    # style choice: rounding DOWN would drop the binding residue out of the footprint
    # and break L3, which is the entire content of the rule. Measured cost over six
    # archived runs: 0-2 extra unclassified residues, with P/LP, B/LB and component
    # counts unchanged in all six. The exact value is recorded beside it so the slack
    # stays auditable.
    step = float(params.step_fp_A)
    r_fp = math.ceil(r_contain_exact / step) * step

    dists = sorted(r.distance_A for r in inside)
    q1, med, q3 = (float(np.percentile(dists, p)) for p in (25, 50, 75))
    iqr = q3 - q1
    fence = q3 + 1.5 * iqr

    n_out = evidence.n_plp_mapped - len(inside)
    record = {
        "rule": "DECISION-RFP-CONTAINMENT-0001",
        "selector": "containment",
        "r_contain_exact_A": r_contain_exact,
        "r_contain_A": r_fp,
        "grid_slack_A": r_fp - r_contain_exact,
        "n_plp_mapped": int(evidence.n_plp_mapped),
        "n_evidence_set": len(inside),
        "n_plp_outside_evidence": int(n_out),
        "fraction_plp_outside_evidence": (
            n_out / evidence.n_plp_mapped if evidence.n_plp_mapped else None),
        "binding_residue": int(binding.residue_id),
        "binding_nearest_center": int(binding.nearest_center_id),
        "evidence_distance_median_A": med,
        "evidence_distance_q1_A": q1,
        "evidence_distance_q3_A": q3,
        "evidence_distance_iqr_A": iqr,
        "binding_residue_is_outlier": bool(r_contain_exact > fence),
        "outlier_fence_A": fence,
    }

    domain = FootprintDomain(
        rho_min=r_fp, rho_max=r_fp, rho_all=rho_all, step=step,
        grid=(r_fp,), d_max=d_max, mst=edges, merge_scales=merge_scales,
        single_center_special_case=(len(centers) == 1),
        binding_constraint="containment",
        rho_max_candidates=(("containment", r_fp),),
    )
    return domain, record
