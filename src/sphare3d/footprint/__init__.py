"""Stage C — footprint construction and ``r_fp`` selection (METHOD_SPEC II.8-II.10).

The public surface is deliberately small. Phase D imports :func:`build_footprint`,
:func:`reconstruct_at` and :data:`FOOTPRINT_API_VERSION` and nothing else — it
consumes this package, it never edits it.
"""
from .api import (FOOTPRINT_API_VERSION, CenterSet, FixedReconstruction,
                  FootprintSolution, Universe, build_footprint,
                  footprint_code_version, occupancy_on_grid, reconstruct_at)
from .params import FootprintParams
from .stage import run_stage_cd

__all__ = [
    "FOOTPRINT_API_VERSION",
    "CenterSet",
    "FixedReconstruction",
    "FootprintParams",
    "FootprintSolution",
    "Universe",
    "build_footprint",
    "footprint_code_version",
    "occupancy_on_grid",
    "reconstruct_at",
    "run_stage_cd",
]
