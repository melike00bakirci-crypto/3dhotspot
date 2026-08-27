"""METHOD_SPEC II.7 — the three mandated post-primary sensitivity analyses.

**NON-REDEFINING, without exception.** Each analysis re-runs the final detection at
the *frozen* ``r_hot`` on a restricted input and reports the overlap with the primary
center set ``S``. None of them may add a center, remove a center, move ``r_hot``,
change ``q``, change ``B``, or alter any published number. They are read *after* the
primary result, never before it.

  1. **pLDDT >= 70** — structural-confidence sensitivity (F2: pLDDT is recorded and
     is NEVER a primary filter).
  2. **Review status >= 1 star and >= 2 stars** (F12: stars are NEVER an inclusion
     criterion). Review metadata is read through the separate, explicitly
     non-redefining channel declared in ``handoff_01`` — never through the Stage B
     primary-path allowlist. If no such channel is declared the analysis is recorded
     as NOT_EVALUABLE with the reason; it is never silently skipped and the forbidden
     columns are never read from the primary tables.
  3. **Global-clustering caveat** — the F5 flag carried by every downstream claim.

No pre-registered numeric threshold exists for "the result is sensitive", so none is
invented here. The MAJOR warning fires only on a *qualitative* flip (the primary
result was non-empty and the restricted analysis is empty); the continuous overlap
statistics are always reported so a reader can judge for themselves.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..utils.geometry import jaccard
from .detection import final_detection

SENSITIVITY_OVERLAP_COLUMNS = [
    "schema_version", "analysis", "status", "radius_A", "n_labeled_used", "n_plp_used",
    "n_blb_used", "n_centers_universe", "n_centers_primary", "n_centers_sensitivity",
    "n_overlap", "overlap_fraction_of_primary", "jaccard", "conclusion_preserved",
    "is_non_redefining", "reason",
]

NON_REDEFINING_TEXT = """SENSITIVITY ANALYSES ARE NON-REDEFINING
=======================================

Everything in FULL_RESULTS/11_SENSITIVITY/ was computed AFTER the primary result was
final, at the FROZEN hotspot radius r_hot, and none of it may change the primary
result in any way.

Specifically, no file in this directory may be used to:

  * add a residue to, or remove a residue from, SIGNIFICANT_HOTSPOT_CENTERS;
  * change r_hot, the search domain, q, B, or the FDR method;
  * change any p-value, q-value, fold enrichment or hotspot region;
  * re-run, re-seed or re-select any part of the primary analysis.

The primary analysis deliberately applies NO pLDDT filter (F2) and NO ClinVar review
star filter (F12). These analyses exist to show how the primary result would look if
such filters HAD been applied, which is a statement about robustness — not a competing
result and not a better result.

A disagreement between a sensitivity analysis and the primary result is reported as a
caveat on the primary result. It is never resolved by adopting the sensitivity result.
"""


@dataclass
class SensitivityOutcome:
    name: str
    status: str                       # COMPLETED | NOT_EVALUABLE | NOT_APPLICABLE
    reason: str
    centers: list[int] = field(default_factory=list)
    detail: dict = field(default_factory=dict)


def _overlap(primary: set[int], sens: set[int]) -> dict:
    inter = primary & sens
    return {
        "n_centers_primary": len(primary),
        "n_centers_sensitivity": len(sens),
        "n_overlap": len(inter),
        "overlap_fraction_of_primary": (len(inter) / len(primary)) if primary else None,
        "jaccard": jaccard(primary, sens),
        "centers_lost": sorted(primary - sens),
        "centers_gained_not_adopted": sorted(sens - primary),
        "conclusion_preserved": bool(sens) if primary else None,
    }


def run_restricted_detection(name: str, radius: float, universe_coords: np.ndarray,
                             universe_index: np.ndarray, universe_keep: np.ndarray,
                             labeled_positions: np.ndarray, labeled_coords: np.ndarray,
                             y: np.ndarray, labeled_keep: np.ndarray,
                             rng: np.random.Generator, *, B: int, q: float,
                             seed_context: str, primary_centers: set[int],
                             chunk: int = 2000) -> SensitivityOutcome:
    """Re-run II.7 at the FROZEN radius on a restricted universe/cohort."""
    keep_u = np.asarray(universe_keep, dtype=bool)
    keep_l = np.asarray(labeled_keep, dtype=bool)
    y_sub = np.asarray(y)[keep_l]

    if keep_u.sum() == 0 or keep_l.sum() == 0:
        return SensitivityOutcome(name, "NOT_EVALUABLE",
                                  "restriction leaves an empty universe or cohort")
    if y_sub.sum() == 0 or (len(y_sub) - y_sub.sum()) == 0:
        return SensitivityOutcome(
            name, "NOT_EVALUABLE",
            f"restriction leaves one class empty (N_P={int(y_sub.sum())}, "
            f"N_B={int(len(y_sub) - y_sub.sum())}); the LOO prevalence and the secondary "
            f"label-permutation null are both undefined, and a one-class cohort cannot "
            f"reproduce the primary contrast")

    # Positions of the retained labelled residues within the retained universe.
    index_map = {int(idx): i for i, idx in enumerate(universe_index[keep_u])}
    kept_labeled_index = universe_index[np.asarray(labeled_positions)[keep_l]]
    if not all(int(i) in index_map for i in kept_labeled_index):
        return SensitivityOutcome(
            name, "NOT_EVALUABLE",
            "a retained labelled residue is absent from the retained universe")
    positions = np.array([index_map[int(i)] for i in kept_labeled_index], dtype=np.int64)

    detection = final_detection(
        radius, universe_coords[keep_u], universe_index[keep_u], positions,
        np.asarray(labeled_coords)[keep_l], y_sub, rng, B=B, q=q,
        seed_context=seed_context, chunk=chunk)

    centers = sorted(r["center_residue_index"] for r in detection.significant_rows)
    detail = _overlap(primary_centers, set(centers))
    detail.update({
        "radius_A": float(radius),
        "n_labeled_used": int(keep_l.sum()),
        "n_plp_used": int(y_sub.sum()),
        "n_blb_used": int(len(y_sub) - y_sub.sum()),
        "n_centers_universe": int(keep_u.sum()),
        "bh_boundary_p": detection.boundary_p,
        "n_centers_without_variant": detection.n_centers_without_variant,
        "seed_context": seed_context, "B": int(B), "q": float(q),
        "IS_NON_REDEFINING": True,
    })
    return SensitivityOutcome(name, "COMPLETED", "NA", centers=centers, detail=detail)


def overlap_row(o: SensitivityOutcome) -> dict:
    d = o.detail
    return {
        "analysis": o.name, "status": o.status,
        "radius_A": d.get("radius_A"), "n_labeled_used": d.get("n_labeled_used"),
        "n_plp_used": d.get("n_plp_used"), "n_blb_used": d.get("n_blb_used"),
        "n_centers_universe": d.get("n_centers_universe"),
        "n_centers_primary": d.get("n_centers_primary"),
        "n_centers_sensitivity": d.get("n_centers_sensitivity"),
        "n_overlap": d.get("n_overlap"),
        "overlap_fraction_of_primary": d.get("overlap_fraction_of_primary"),
        "jaccard": d.get("jaccard"),
        "conclusion_preserved": d.get("conclusion_preserved"),
        "is_non_redefining": True,
        "reason": o.reason,
    }


def conclusion_flipped(o: SensitivityOutcome) -> bool:
    """True when a non-empty primary center set becomes empty under the restriction."""
    if o.status != "COMPLETED":
        return False
    primary = o.detail.get("n_centers_primary") or 0
    return bool(primary > 0 and (o.detail.get("n_centers_sensitivity") or 0) == 0)
