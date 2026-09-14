"""II.9 — the candidate ``r_fp`` search domain, derived from center geometry alone.

    merge scales : {w_k / 2} over the edges of the Euclidean MST of S
    rho_all      : max_k w_k / 2
    rho_min      : 5.0 A                                   (DECISION-FOOTPRINT-DOMAIN-0001)
    rho_max      : min(0.25 * D_max, 25.0)                 (DECISION-FOOTPRINT-DOMAIN-0001)
    step_fp      : 0.5 A, uniform and never adaptive

The merge scales **inform the upper bound and are reported as events; they never
determine ``r_fp``** (F10). ``r_hot`` plays no part here whatsoever: it is not an
argument to any function in this module, which is the mechanical form of the
guarantee that ``r_fp`` is neither set equal to nor bounded below by ``r_hot``.
"""
from __future__ import annotations

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
            "note": ("Merge scales inform rho_max and are reported as events; they "
                     "never determine r_fp (F10). r_hot is not an input to this "
                     "derivation."),
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
