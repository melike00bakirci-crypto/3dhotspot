"""Stage D — geometric footprint robustness (METHOD_SPEC II.11).

Perturbs the SIGNIFICANT_HOTSPOT_CENTERS set and rebuilds the footprint with the
frozen Phase C methodology only. Never re-runs hotspot discovery, never modifies
ClinVar records or class labels, and never emits a categorical robustness verdict.
"""
from .design import SubsetDesign, build_design, k_max, rank_subset, unrank_subset
from .params import RobustnessParams
from .stage import run_phase_d

__all__ = [
    "RobustnessParams",
    "SubsetDesign",
    "build_design",
    "k_max",
    "rank_subset",
    "run_phase_d",
    "unrank_subset",
]
