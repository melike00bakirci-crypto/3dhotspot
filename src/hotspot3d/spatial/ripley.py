"""METHOD_SPEC II.1 — global spatial statistics under the structure-aware positional null.

Ripley-type ``K`` and the pair correlation function ``g(r)`` are computed
**separately for P/LP and B/LB**, on a 1.0 A grid up to ``r_max = floor(D_max/2)``,
with pointwise 2.5/97.5 envelopes and the global test ``T = max_r |Z(r)|`` carrying
its own permutation distribution.

FROZEN facts implemented here:
  * the null is **positional** (F4): uniform *without-replacement* subsets of
    ``U_struct`` of the same size as the observed cohort. The protein geometry is
    therefore held fixed and only the *positions* of the cohort are randomized;
  * this is a **global assessment only**. Ripley's K never selects ``r_hot`` and a
    non-significant global result does not terminate local discovery (F5) — it sets
    the ``global_clustering`` caveat flag carried by every downstream claim;
  * the PCF is used solely to derive the *candidate radius domain* of II.2. A PCF
    peak is never ``r_hot`` (F6/F14).

Estimator conventions (recorded, never re-derived per run):
  ``Kbar(r) = (1/n) * #{ordered pairs (i,j), i != j, d_ij <= r}``
  ``K(r)    = Kbar(r) / lambda_C``   with ``lambda_C = n_C / V_struct``
  ``L(r)    = (3 K(r) / (4 pi))^(1/3)``          so ``L(r) - r == 0`` under CSR
  ``g(r)    = shellbar(r) / (lambda_C * 4 pi r^2 * w)``  shell ``(r-w/2, r+w/2]``

The intensity is the **cohort's own** intensity ``n_C / V_struct``, not the intensity
of ``U_struct``. That is what makes ``L(r) - r = 0`` and ``g(r) = 1`` the correct CSR
reference lines for a cohort of ``n_C`` points spread over the protein volume — and
it is what makes the II.2 criterion ``g_PLP(r) > 1`` meaningful rather than a
statement about how small the cohort is. The observed and null cohorts share ``n_C``,
so the choice cancels in every comparison and affects only the absolute scale.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..utils.geometry import pairwise_distances, protein_diameter, r_max_for_estimators

# Z is defined as 0 where the null standard deviation is exactly 0 (a degenerate
# radius at which every null replicate produced the identical count). Recorded
# convention, not a tuning choice.
_ZERO_SD_Z = 0.0


@dataclass
class CohortStatistics:
    """K / L / g for one cohort with its positional-null envelopes."""

    cohort: str
    n_points: int
    grid: np.ndarray
    k_obs: np.ndarray
    k_null_mean: np.ndarray
    k_null_sd: np.ndarray
    k_env_lo: np.ndarray
    k_env_hi: np.ndarray
    l_obs: np.ndarray
    l_null_mean: np.ndarray
    l_env_lo: np.ndarray
    l_env_hi: np.ndarray
    g_obs: np.ndarray
    g_null_mean: np.ndarray
    g_null_sd: np.ndarray
    g_env_lo: np.ndarray
    g_env_hi: np.ndarray
    g_null_p95: np.ndarray
    k_z: np.ndarray
    g_z: np.ndarray
    t_obs: float
    t_null: np.ndarray
    p_global: float
    B: int
    seed_context: str
    extra: dict = field(default_factory=dict)

    @property
    def significant(self) -> bool:
        return bool(self.p_global <= self.extra.get("alpha", 0.05))


# --- intensity ---------------------------------------------------------------

def structural_volume(coords: np.ndarray) -> tuple[float, str]:
    """``V_struct`` for the intensity ``lambda = M / V_struct``.

    Convex hull of ``U_struct``; axis-aligned bounding box only if the hull is
    degenerate (collinear/coplanar synthetic inputs). The method actually used is
    recorded — the value is a shared constant and cancels in every comparison.
    """
    coords = np.asarray(coords, dtype=np.float64)
    if len(coords) >= 4:
        try:
            from scipy.spatial import ConvexHull

            hull = ConvexHull(coords)
            if hull.volume > 0:
                return float(hull.volume), "convex_hull"
        except Exception:                                     # pragma: no cover
            pass
    span = coords.max(axis=0) - coords.min(axis=0)
    span = np.where(span > 0, span, 1.0)
    return float(np.prod(span)), "bounding_box"


def intensity(coords_universe: np.ndarray) -> tuple[float, float, str]:
    """``(lambda, V_struct, method)`` from the positional universe."""
    volume, method = structural_volume(coords_universe)
    return float(len(coords_universe) / volume), volume, method


# --- estimators --------------------------------------------------------------

def _pair_distances_sorted(dist: np.ndarray, idx: np.ndarray | None = None) -> np.ndarray:
    """Sorted upper-triangular distances of the (sub)set ``idx``."""
    sub = dist if idx is None else dist[np.ix_(idx, idx)]
    iu = np.triu_indices(sub.shape[0], k=1)
    return np.sort(sub[iu])


def _counts(dsorted: np.ndarray, grid: np.ndarray, shell_width: float) -> tuple[np.ndarray, np.ndarray]:
    """Exact cumulative (``d <= r``) and shell (``r-w/2 < d <= r+w/2``) pair counts."""
    cum = np.searchsorted(dsorted, grid, side="right").astype(np.float64)
    hi = np.searchsorted(dsorted, grid + shell_width / 2.0, side="right")
    lo = np.searchsorted(dsorted, grid - shell_width / 2.0, side="right")
    return cum, (hi - lo).astype(np.float64)


def k_l_g_from_counts(cum: np.ndarray, shell: np.ndarray, n: int, grid: np.ndarray,
                      volume: float, shell_width: float
                      ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert unordered pair counts to ``(K, L, g)`` under the recorded conventions.

    ``volume`` is ``V_struct``; the intensity used is the cohort's own ``n / V_struct``.
    """
    if n < 2:
        zeros = np.zeros_like(grid, dtype=np.float64)
        return zeros, zeros.copy(), zeros.copy()
    lam = n / volume
    k = (2.0 * cum / n) / lam
    l = np.cbrt(3.0 * k / (4.0 * np.pi))
    shell_volume = 4.0 * np.pi * grid ** 2 * shell_width
    g = (2.0 * shell / n) / (lam * shell_volume)
    return k, l, g


def observed_curves(coords: np.ndarray, grid: np.ndarray, volume: float,
                    shell_width: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dist = pairwise_distances(coords)
    cum, shell = _counts(_pair_distances_sorted(dist), grid, shell_width)
    return k_l_g_from_counts(cum, shell, len(coords), grid, volume, shell_width)


# --- the structure-aware positional null -------------------------------------

def positional_null_draws(rng: np.random.Generator, n_universe: int, n_points: int,
                          B: int) -> np.ndarray:
    """``B`` uniform **without-replacement** subsets of ``U_struct`` of size ``n_points``.

    Without replacement is the whole point of F4: the null cohort is a *set of
    residues*, exactly like the observed cohort, so no residue can be occupied twice.
    """
    if n_points > n_universe:
        raise ValueError(f"cohort size {n_points} exceeds |U_struct| = {n_universe}")
    return np.array([rng.choice(n_universe, size=n_points, replace=False)
                     for _ in range(B)], dtype=np.int64)


def _z(obs: np.ndarray, mean: np.ndarray, sd: np.ndarray) -> np.ndarray:
    out = np.full_like(obs, _ZERO_SD_Z, dtype=np.float64)
    ok = sd > 0
    out[ok] = (obs[ok] - mean[ok]) / sd[ok]
    return out


def cohort_statistics(cohort: str, cohort_coords: np.ndarray, universe_coords: np.ndarray,
                      grid: np.ndarray, volume: float, shell_width: float, B: int,
                      rng: np.random.Generator, seed_context: str,
                      envelope_percentiles: tuple[float, float] = (2.5, 97.5),
                      alpha: float = 0.05) -> CohortStatistics:
    """Full II.1 analysis for one cohort: curves, envelopes and the global test.

    The global test statistic is ``T = max_r |Z_K(r)|``. Its null distribution is
    obtained by recomputing ``T`` on each null replicate against the pooled null
    mean and sd (a fixed, documented convention), giving
    ``p = (1 + #{b : T_b >= T_obs}) / (1 + B)`` — the same add-one estimator used
    everywhere else in the pipeline, so a global p-value is never exactly zero.
    """
    universe_coords = np.asarray(universe_coords, dtype=np.float64)
    n = len(cohort_coords)
    dist_u = pairwise_distances(universe_coords)

    k_obs, l_obs, g_obs = observed_curves(cohort_coords, grid, volume, shell_width)

    draws = positional_null_draws(rng, len(universe_coords), n, B)
    k_null = np.empty((B, len(grid)), dtype=np.float64)
    g_null = np.empty((B, len(grid)), dtype=np.float64)
    l_null = np.empty((B, len(grid)), dtype=np.float64)
    for b, idx in enumerate(draws):
        cum, shell = _counts(_pair_distances_sorted(dist_u, idx), grid, shell_width)
        kb, lb, gb = k_l_g_from_counts(cum, shell, n, grid, volume, shell_width)
        k_null[b], l_null[b], g_null[b] = kb, lb, gb

    lo_p, hi_p = envelope_percentiles
    k_mean, k_sd = k_null.mean(axis=0), k_null.std(axis=0, ddof=1)
    g_mean, g_sd = g_null.mean(axis=0), g_null.std(axis=0, ddof=1)

    k_z = _z(k_obs, k_mean, k_sd)
    g_z = _z(g_obs, g_mean, g_sd)

    # Global test on K, pooled standardization (recorded convention).
    z_null = np.zeros_like(k_null)
    ok = k_sd > 0
    z_null[:, ok] = (k_null[:, ok] - k_mean[ok]) / k_sd[ok]
    t_null = np.abs(z_null).max(axis=1)
    t_obs = float(np.abs(k_z).max()) if len(grid) else 0.0
    p_global = float((1 + int((t_null >= t_obs).sum())) / (1 + B))

    return CohortStatistics(
        cohort=cohort, n_points=n, grid=grid,
        k_obs=k_obs, k_null_mean=k_mean, k_null_sd=k_sd,
        k_env_lo=np.percentile(k_null, lo_p, axis=0),
        k_env_hi=np.percentile(k_null, hi_p, axis=0),
        l_obs=l_obs, l_null_mean=l_null.mean(axis=0),
        l_env_lo=np.percentile(l_null, lo_p, axis=0),
        l_env_hi=np.percentile(l_null, hi_p, axis=0),
        g_obs=g_obs, g_null_mean=g_mean, g_null_sd=g_sd,
        g_env_lo=np.percentile(g_null, lo_p, axis=0),
        g_env_hi=np.percentile(g_null, hi_p, axis=0),
        g_null_p95=np.percentile(g_null, 95.0, axis=0),
        k_z=k_z, g_z=g_z, t_obs=t_obs, t_null=t_null, p_global=p_global,
        B=B, seed_context=seed_context,
        extra={"alpha": alpha, "envelope_percentiles": list(envelope_percentiles),
               "V_struct_A3": volume, "lambda_cohort": (n / volume if volume else None),
               "shell_width_A": shell_width,
               "null_model": "structure_aware_positional",
               "global_test_statistic": "max_abs_z_on_K"},
    )


def build_grid(universe_coords: np.ndarray, grid_step: float) -> tuple[np.ndarray, float, int]:
    """``(grid, D_max, r_max)`` with ``r_max = floor(D_max / 2)`` (II.1)."""
    d_max = protein_diameter(universe_coords)
    r_max = r_max_for_estimators(d_max)
    grid = np.arange(1, int(np.floor(r_max / grid_step)) + 1, dtype=np.float64) * grid_step
    return grid, float(d_max), int(r_max)


def pcf_peaks(grid: np.ndarray, g_obs: np.ndarray) -> dict:
    """Local maxima of the observed PCF — DESCRIPTIVE ONLY.

    Recorded so a reader can see what the classical "pick the PCF peak" shortcut
    would have chosen, and see that the pipeline did not do that (F6/F14).
    """
    peaks: list[dict] = []
    for i in range(len(grid)):
        left = g_obs[i - 1] if i > 0 else -np.inf
        right = g_obs[i + 1] if i < len(grid) - 1 else -np.inf
        if g_obs[i] > left and g_obs[i] >= right:
            peaks.append({"radius_A": float(grid[i]), "g": float(g_obs[i])})
    principal = max(peaks, key=lambda p: (p["g"], -p["radius_A"]), default=None)
    return {
        "peaks": peaks,
        "principal_peak": principal,
        "USED_AS_R_HOT": False,
        "note": ("The pair-correlation peak is NEVER r_hot and never defines the "
                 "search interval on its own (F6/F14). Peaks are recorded here only "
                 "so the reader can verify that the selected radius was not taken "
                 "from this list."),
    }
