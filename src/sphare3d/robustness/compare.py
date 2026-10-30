"""¶47 comparison metrics: ``FP_original`` versus ``FP_iteration``.

Covered-residue Jaccard and Dice are the **primary** measures because they are
resolution-independent — they do not move when the voxel edge changes. Voxel
overlap is reported alongside them, but only ever after both footprints have been
re-evaluated on one **common reference grid**: each iteration derives its own
bounding box, so voxel sets from two different grids are not comparable and any
"overlap" between them would be an artifact of the frames.

Every metric of ¶47 is reported, including the unflattering ones. Reporting Dice
while suppressing Jaccard is prohibited (§8).
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from ..footprint.api import CenterSet, Universe, occupancy_on_grid
from ..footprint.grid import make_grid
from ..footprint.params import FootprintParams
from ..utils.geometry import centroid, dice, hausdorff95, jaccard

GEOMETRY_FIELDS = [
    "jaccard_residues", "dice_residues", "jaccard_voxels", "dice_voxels",
    "n_footprint_residues", "coverage_fraction", "volume_A3", "volume_ratio",
    "volume_delta_A3", "surface_area_A2", "surface_ratio", "n_components",
    "component_preserved", "continuity", "centroid_shift_A", "hausdorff95_A",
]


@dataclass
class GeometryComparison:
    jaccard_residues: float
    dice_residues: float
    jaccard_voxels: float
    dice_voxels: float
    n_footprint_residues: int
    coverage_fraction: float
    volume_A3: float
    volume_ratio: float
    volume_delta_A3: float
    surface_area_A2: float
    surface_ratio: float
    n_components: int
    component_preserved: bool
    continuity: bool
    centroid_shift_A: float
    hausdorff95_A: float

    def as_dict(self, suffix: str = "") -> dict:
        return {f"{k}{suffix}": v for k, v in asdict(self).items()}

    @staticmethod
    def failed(suffix: str = "") -> dict:
        """A failed iteration contributes J = 0 and stays in every denominator."""
        row = {k: None for k in GEOMETRY_FIELDS}
        row["jaccard_residues"] = 0.0
        row["dice_residues"] = 0.0
        row["jaccard_voxels"] = 0.0
        row["dice_voxels"] = 0.0
        row["component_preserved"] = False
        row["continuity"] = False
        return {f"{k}{suffix}": v for k, v in row.items()}


def comparison_grid(original_centers: CenterSet, params: FootprintParams,
                    radii: list[float]):
    """One frame both footprints are evaluated in.

    ``S_iter`` is a subset of ``S``, so the original bounding box already contains
    every iteration's centers; padding by the larger of the two radii guarantees
    it also contains both footprints.
    """
    return make_grid(original_centers.coords,
                     pad=max(radii) + params.bbox_padding_A,
                     h=params.voxel_h_A, h_fallback=params.voxel_h_fallback_A,
                     max_voxels=params.voxel_count_fallback_threshold)


def compare(original_state, original_residues: frozenset,
            original_centers: CenterSet, r_fp_original: float,
            iter_state, iter_residues: frozenset, iter_centers: CenterSet,
            r_fp_iter: float, universe: Universe,
            params: FootprintParams) -> GeometryComparison:
    """All ¶47 geometry metrics for one reconstruction."""
    spec = comparison_grid(original_centers, params,
                           [float(r_fp_original), float(r_fp_iter)])
    occ_orig = occupancy_on_grid(spec, original_centers, float(r_fp_original))
    occ_iter = occupancy_on_grid(spec, iter_centers, float(r_fp_iter))
    inter = int(np.count_nonzero(occ_orig & occ_iter))
    union = int(np.count_nonzero(occ_orig | occ_iter))
    n_orig = int(np.count_nonzero(occ_orig))
    n_iter = int(np.count_nonzero(occ_iter))

    coords_by_id = {rid: universe.coords[i] for i, rid in enumerate(universe.ids)}
    orig_xyz = np.array([coords_by_id[r] for r in sorted(original_residues)]) \
        if original_residues else np.empty((0, 3))
    iter_xyz = np.array([coords_by_id[r] for r in sorted(iter_residues)]) \
        if iter_residues else np.empty((0, 3))

    shift = float("nan")
    if len(orig_xyz) and len(iter_xyz):
        shift = float(np.linalg.norm(centroid(iter_xyz) - centroid(orig_xyz)))

    return GeometryComparison(
        jaccard_residues=jaccard(set(original_residues), set(iter_residues)),
        dice_residues=dice(set(original_residues), set(iter_residues)),
        jaccard_voxels=float(inter / union) if union else 0.0,
        dice_voxels=float(2 * inter / (n_orig + n_iter)) if (n_orig + n_iter) else 0.0,
        n_footprint_residues=len(iter_residues),
        coverage_fraction=float(len(iter_residues) / len(universe)),
        volume_A3=float(iter_state.volume),
        volume_ratio=(float(iter_state.volume / original_state.volume)
                      if original_state.volume > 0 else float("nan")),
        volume_delta_A3=float(iter_state.volume - original_state.volume),
        surface_area_A2=float(iter_state.surface_area),
        surface_ratio=(float(iter_state.surface_area / original_state.surface_area)
                       if original_state.surface_area > 0 else float("nan")),
        n_components=int(iter_state.n_components),
        component_preserved=bool(iter_state.n_components
                                 == original_state.n_components),
        continuity=bool(len(set(original_residues) & set(iter_residues)) > 0),
        centroid_shift_A=shift,
        hausdorff95_A=hausdorff95(orig_xyz, iter_xyz),
    )
