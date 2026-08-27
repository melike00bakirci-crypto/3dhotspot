"""Stage B global spatial statistics (METHOD_SPEC II.1 / II.2).

Frozen, deterministic and importable by other stages; never edited by them.
"""
from .domain import (
    TRIGGER_DEFINITIONS,
    RadiusDomain,
    derive_domain,
    elevated_set,
    lattice,
    longest_contiguous_run,
)
from .ripley import (
    CohortStatistics,
    build_grid,
    cohort_statistics,
    intensity,
    k_l_g_from_counts,
    observed_curves,
    pcf_peaks,
    positional_null_draws,
    structural_volume,
)

__all__ = [
    "CohortStatistics", "build_grid", "cohort_statistics", "intensity",
    "k_l_g_from_counts", "observed_curves", "pcf_peaks", "positional_null_draws",
    "structural_volume", "RadiusDomain", "derive_domain", "elevated_set",
    "lattice", "longest_contiguous_run", "TRIGGER_DEFINITIONS",
]
