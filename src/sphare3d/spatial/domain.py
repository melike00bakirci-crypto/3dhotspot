"""METHOD_SPEC II.2 — the candidate ``r_hot`` search domain.

The domain is derived **deterministically, before any radius is scored**, from the
P/LP pair correlation function and its positional null:

    E     = { r : g_PLP(r) > 1  AND  g_PLP_obs(r) > 95th percentile of its null }
    [a,b] = the longest contiguous run of E on the PCF grid
    D_hot = [a - margin, b + margin]  clamped to [R_FLOOR, r_cap]
    r_cap = min(R_CEIL, 0.25 * D_max)
    grid  = D_hot on a 0.5 A lattice

Four pre-registered fallback triggers force the frozen fallback envelope
``[5.0, min(25.0, 0.25 * D_max)]`` instead:

    FT-1  E is empty
    FT-2  the longest run is shorter than ``min_run_length_A``
    FT-3  the derived domain carries fewer than ``min_grid_points`` grid points
    FT-4  ``N_P < N_MIN_PCF`` — too few pathogenic residues for a usable PCF

Every trigger sets ``FALLBACK_RADIUS_DOMAIN = TRUE`` with the triggering condition
recorded and surfaced; the switch is never silent. ``R_FLOOR > r_cap`` makes the
domain inadmissible — a negative result with escalation, never a quiet widening.

**No domain shopping.** The branch and the clamps are fixed here and are never
re-derived because a result looked unattractive; automatic domain expansion is
prohibited by config (``radius_domain.automatic_domain_expansion: false``).
**The PCF peak is never r_hot and never defines the interval alone** (F6/F14).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

TRIGGER_DEFINITIONS = {
    "FT-1": "E is empty: no radius has g_PLP(r) > 1 AND g_PLP(r) above its null 95th percentile.",
    "FT-2": "The longest contiguous run of E spans less than radius_domain.min_run_length_A.",
    "FT-3": "The derived domain carries fewer than radius_domain.min_grid_points grid points.",
    "FT-4": "N_P < radius_domain.N_MIN_PCF — too few P/LP residues for a usable PCF.",
}


@dataclass
class RadiusDomain:
    """The candidate radius domain plus the complete audit trail of its derivation."""

    lo: float
    hi: float
    grid: np.ndarray
    source: str                       # pcf_derived | fallback_envelope | inadmissible
    fallback: bool
    trigger_id: str | None
    triggers_fired: list[str]
    elevated_radii: list[float]       # the set E
    longest_run: tuple[float, float] | None
    run_length_A: float
    r_floor: float
    r_cap: float
    d_max: float
    step: float
    admissible: bool
    detail: dict = field(default_factory=dict)

    def as_json(self) -> dict:
        return {
            "branch_taken": self.source,
            "FALLBACK_RADIUS_DOMAIN": self.fallback,
            "trigger_id": self.trigger_id or "NA",
            "triggers_fired": self.triggers_fired,
            "trigger_definitions": TRIGGER_DEFINITIONS,
            "E_elevated_radii_A": self.elevated_radii,
            "E_criterion": ("g_PLP(r) > 1 AND g_PLP_obs(r) > 95th percentile of the "
                            "structure-aware positional null at r"),
            "longest_contiguous_run": (list(self.longest_run) if self.longest_run else None),
            "run_length_A": self.run_length_A,
            "margin_A": self.detail.get("margin_A"),
            "pre_clamp_interval": self.detail.get("pre_clamp_interval"),
            "clamped_interval": [self.lo, self.hi],
            "R_FLOOR": self.r_floor,
            "R_CEIL": self.detail.get("R_CEIL"),
            "r_cap": self.r_cap,
            "r_cap_rule": "min(R_CEIL, 0.25 * D_max)",
            "D_max_A": self.d_max,
            "step_hot_A": self.step,
            "grid_A": [float(x) for x in self.grid],
            "n_grid_points": int(len(self.grid)),
            "domain_admissible": self.admissible,
            "inadmissible_reason": self.detail.get("inadmissible_reason", "NA"),
            "automatic_domain_expansion": False,
            "pcf_peak_used_as_r_hot": False,
            "note": ("Derived before any radius was scored. The PCF peak is never r_hot "
                     "and never defines this interval alone (F6/F14). Automatic domain "
                     "expansion is prohibited; widening requires a Lead-authorized, "
                     "pre-registered domain under a NEW RUN_ID."),
        }


def elevated_set(grid: np.ndarray, g_obs: np.ndarray, g_null_p95: np.ndarray) -> np.ndarray:
    """Boolean mask of ``E`` — both conditions, evaluated jointly (II.2)."""
    return (g_obs > 1.0) & (g_obs > g_null_p95)


def longest_contiguous_run(grid: np.ndarray, mask: np.ndarray) -> tuple[float, float] | None:
    """Longest run of consecutive ``True`` grid points.

    Ties in length are broken by the SMALLEST starting radius — a fixed rule, so the
    interval can never be chosen after the fact.
    """
    best: tuple[int, int] | None = None
    start: int | None = None
    for i, flag in enumerate(list(mask) + [False]):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            run = (start, i - 1)
            if best is None or (run[1] - run[0]) > (best[1] - best[0]):
                best = run
            start = None
    if best is None:
        return None
    return float(grid[best[0]]), float(grid[best[1]])


def lattice(lo: float, hi: float, step: float) -> np.ndarray:
    """Grid points ``lo, lo+step, ...`` not exceeding ``hi`` (deterministic, on-lattice)."""
    if hi < lo:
        return np.zeros(0, dtype=np.float64)
    n = int(np.floor(round((hi - lo) / step, 9))) + 1
    return np.round(lo + step * np.arange(n, dtype=np.float64), 6)


def derive_domain(grid: np.ndarray, g_obs: np.ndarray, g_null_p95: np.ndarray,
                  n_plp: int, d_max: float, *, r_floor: float, r_ceil: float,
                  step: float, margin: float, n_min_pcf: int,
                  min_run_length: float, min_grid_points: int) -> RadiusDomain:
    """Execute II.2 exactly, recording every branch, clamp and trigger."""
    r_cap = float(min(r_ceil, 0.25 * d_max))
    detail: dict = {"margin_A": margin, "R_CEIL": r_ceil,
                    "N_MIN_PCF": n_min_pcf, "min_run_length_A": min_run_length,
                    "min_grid_points": min_grid_points}

    mask = elevated_set(grid, g_obs, g_null_p95)
    e_radii = [float(r) for r in grid[mask]]
    run = longest_contiguous_run(grid, mask)
    run_length = float(run[1] - run[0]) if run else 0.0

    triggers: list[str] = []
    if run is None:
        triggers.append("FT-1")
    elif run_length < min_run_length:
        triggers.append("FT-2")

    derived_lo = derived_hi = None
    derived_grid = np.zeros(0, dtype=np.float64)
    if run is not None:
        derived_lo = max(r_floor, run[0] - margin)
        derived_hi = min(r_cap, run[1] + margin)
        detail["pre_clamp_interval"] = [run[0] - margin, run[1] + margin]
        derived_grid = lattice(derived_lo, derived_hi, step)
        if len(derived_grid) < min_grid_points:
            triggers.append("FT-3")
    if n_plp < n_min_pcf:
        triggers.append("FT-4")

    # R_FLOOR > r_cap makes every branch inadmissible — the protein is too small for
    # the pre-registered floor. Negative result plus escalation; never a silent clamp.
    if r_floor > r_cap:
        detail["inadmissible_reason"] = (
            f"R_FLOOR ({r_floor} A) exceeds r_cap ({r_cap:.6g} A = min(R_CEIL, 0.25*D_max) "
            f"with D_max = {d_max:.6g} A): the pre-registered radius floor does not fit "
            f"inside this structure.")
        return RadiusDomain(
            lo=float("nan"), hi=float("nan"), grid=np.zeros(0), source="inadmissible",
            fallback=bool(triggers), trigger_id=(sorted(triggers)[0] if triggers else None),
            triggers_fired=sorted(triggers), elevated_radii=e_radii, longest_run=run,
            run_length_A=run_length, r_floor=r_floor, r_cap=r_cap, d_max=d_max,
            step=step, admissible=False, detail=detail)

    if triggers:
        lo, hi = float(r_floor), float(min(r_ceil, 0.25 * d_max))
        grid_out = lattice(lo, hi, step)
        detail["fallback_envelope_rule"] = "[R_FLOOR, min(R_CEIL, 0.25 * D_max)]"
        if len(grid_out) < min_grid_points:
            detail["fallback_grid_below_min_points"] = True
        return RadiusDomain(
            lo=lo, hi=hi, grid=grid_out, source="fallback_envelope", fallback=True,
            trigger_id=sorted(triggers)[0], triggers_fired=sorted(triggers),
            elevated_radii=e_radii, longest_run=run, run_length_A=run_length,
            r_floor=r_floor, r_cap=r_cap, d_max=d_max, step=step, admissible=True,
            detail=detail)

    return RadiusDomain(
        lo=float(derived_lo), hi=float(derived_hi), grid=derived_grid,
        source="pcf_derived", fallback=False, trigger_id=None, triggers_fired=[],
        elevated_radii=e_radii, longest_run=run, run_length_A=run_length,
        r_floor=r_floor, r_cap=r_cap, d_max=d_max, step=step, admissible=True,
        detail=detail)
