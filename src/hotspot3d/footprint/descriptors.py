"""II.8 geometric descriptors evaluated on one occupancy grid.

Every descriptor here is a *measurement*. Only four of them (compactness,
connectivity, stability, parsimony) are later promoted to Pareto objectives;
convexity, expansion rate, component count, bridge count and merge events are
descriptive by construction (II.10) and are never optimized. ``n_comp`` in
particular is a measurement, not a target — "one clean blob" is not a goal.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage

# 6-connectivity == face adjacency (A18). Fixed here so no caller can widen it.
_STRUCT_6 = ndimage.generate_binary_structure(3, 1)


def crop_slices(mask: np.ndarray, margin: int = 2
                ) -> tuple[tuple[slice, slice, slice], np.ndarray]:
    """Tight bounding box of ``mask`` plus a background ``margin``.

    Every per-radius descriptor is local to the occupied region, so cropping is an
    exact optimization, not an approximation: the nearest background voxel to any
    occupied voxel lies inside the margin, and the level-``rho`` isosurface lies
    within one voxel of the occupancy boundary. At small radii this is the
    difference between scanning the padded protein box and scanning a few spheres.
    """
    idx = np.argwhere(mask)
    if idx.size == 0:
        zero = (slice(0, 0), slice(0, 0), slice(0, 0))
        return zero, np.zeros(3, dtype=np.int64)
    lo = np.maximum(idx.min(axis=0) - margin, 0)
    hi = np.minimum(idx.max(axis=0) + margin + 1, np.asarray(mask.shape))
    return tuple(slice(int(lo[d]), int(hi[d])) for d in range(3)), lo.astype(np.int64)


@dataclass(frozen=True)
class ComponentField:
    labels: np.ndarray          # int label per voxel, 0 = background
    n_components: int
    sizes: np.ndarray           # voxel count per label, index 0 unused


def components(occ: np.ndarray) -> ComponentField:
    """6-connectivity connected components of the occupancy grid."""
    labels, n = ndimage.label(occ, structure=_STRUCT_6)
    if n == 0:
        return ComponentField(labels=labels, n_components=0,
                              sizes=np.zeros(1, dtype=np.int64))
    sizes = np.bincount(labels.ravel(), minlength=n + 1).astype(np.int64)
    return ComponentField(labels=labels, n_components=int(n), sizes=sizes)


def volume(occ: np.ndarray, h: float) -> float:
    return float(int(occ.sum()) * (h ** 3))


def surface_area(dist: np.ndarray, rho: float, h: float,
                 occupied: bool = True) -> tuple[float, dict | None, str]:
    """Marching-cubes surface area of ``FP(rho)`` (A18).

    The isosurface is extracted at level ``rho`` from **the same single distance
    field** that defines the occupancy, so the surface and the volume describe one
    object and the whole sweep still rests on one EDT. Extracting instead from the
    binary mask would carry the well-known ~8% staircase bias of marching cubes on
    indicator data and would put ``Psi`` of a perfect sphere at 0.92 rather than 1.

    Returns ``(area, mesh, failure_reason)``; ``area`` is NaN and ``mesh`` is None
    when extraction fails, which QC-F3 turns into an inadmissible candidate rather
    than a crash.
    """
    from skimage import measure

    if not occupied:
        return float("nan"), None, "EMPTY_VOLUME"
    outside = float(np.max(dist)) + 10.0 * h
    padded = np.pad(dist, 1, mode="constant", constant_values=outside)
    try:
        verts, faces, normals, values = measure.marching_cubes(
            padded, level=float(rho), spacing=(h, h, h))
    except (ValueError, RuntimeError) as exc:       # no isosurface / degenerate
        return float("nan"), None, f"MESH_EXTRACTION_FAILED:{type(exc).__name__}"
    area = float(measure.mesh_surface_area(verts, faces))
    if not np.isfinite(area) or area <= 0:
        return float("nan"), None, "MESH_SURFACE_AREA_NON_POSITIVE"
    # Sample [i,j,k] of ``dist`` sits at the CENTER of voxel (i,j,k), i.e. at
    # origin + (idx + 0.5) * h. Undo the 1-voxel pad (-h) and shift to the voxel
    # center (+0.5h) so that world coordinates are simply ``origin + verts``.
    mesh = {"verts": verts - 0.5 * h, "faces": faces, "normals": normals}
    return area, mesh, ""


def compactness(vol: float, area: float) -> float:
    """``Psi = (36 pi V^2)^(1/3) / A`` — 1.0 for a perfect sphere, < 1 otherwise."""
    if not np.isfinite(area) or area <= 0 or vol <= 0:
        return float("nan")
    return float((36.0 * np.pi * vol ** 2) ** (1.0 / 3.0) / area)


def surface_voxels(occ: np.ndarray) -> np.ndarray:
    """Occupied voxels with at least one background face neighbour."""
    interior = ndimage.binary_erosion(occ, structure=_STRUCT_6, border_value=0)
    return occ & ~interior


def convexity(occ: np.ndarray, h: float, origin: tuple[float, float, float],
              vol: float) -> float:
    """``X = V / V_hull`` over the occupied voxel centers (descriptive only).

    Only boundary voxels are passed to Qhull: the convex hull of a solid equals
    the convex hull of its boundary, so this is exact and orders of magnitude
    cheaper than hulling the full interior.
    """
    from scipy.spatial import ConvexHull, QhullError

    idx = np.argwhere(surface_voxels(occ))
    if len(idx) < 4:
        return float("nan")
    pts = np.asarray(origin) + (idx + 0.5) * h
    try:
        hull = ConvexHull(pts)
    except (QhullError, ValueError):
        return float("nan")
    if hull.volume <= 0:
        return float("nan")
    return float(vol / hull.volume)


def connectivity(field: ComponentField) -> float:
    """``C = V_largest / V_total``."""
    if field.n_components == 0:
        return float("nan")
    total = float(field.sizes[1:].sum())
    if total <= 0:
        return float("nan")
    return float(field.sizes[1:].max() / total)


def erosion_bridging(occ: np.ndarray, h: float, depth: float) -> tuple[int, float, int]:
    """Artificial bridging by relative erosion at ``depth = 0.15 * rho`` (A16).

    A connection that vanishes when the footprint is eroded by ``depth`` was held
    together by a neck thinner than ``2 * depth``: eroding therefore *counts* the
    thin connections instead of removing them. Bridges are measured and reported;
    they are never hand-removed from the geometry.

    Returns ``(n_bridges, bridging_index, n_components_eroded)``.
    """
    field = components(occ)
    if field.n_components == 0:
        return 0, 0.0, 0
    if depth <= 0:
        return 0, 0.0, field.n_components
    # distance from each occupied voxel to the nearest background voxel
    d_bg = ndimage.distance_transform_edt(occ, sampling=(h, h, h))
    eroded = occ & (d_bg >= depth)
    eroded_field = components(eroded)
    n_bridges = max(0, eroded_field.n_components - field.n_components)
    denom = max(1, eroded_field.n_components)
    return int(n_bridges), float(n_bridges / denom), int(eroded_field.n_components)


def bottleneck_radius(occ: np.ndarray, h: float, a_idx: tuple[int, int, int],
                      b_idx: tuple[int, int, int]) -> float:
    """Widest-path bottleneck between two voxels — the half-width of the neck.

    Returns ``max over paths of min over path voxels of d_bg``: the largest
    threshold ``tau`` such that ``a`` and ``b`` remain connected within
    ``occ & (d_bg >= tau)``. That predicate is monotone in ``tau``, so the value
    is found by binary search over the distinct ``d_bg`` levels at occupied
    voxels, each level probed with one compiled ``ndimage.label`` connectivity
    check — O(log n) compiled calls in place of an O(n) per-voxel Python
    union-find (QA-09-BOTTLENECK-PERF-0001). The two are exactly equivalent, not
    just empirically close: within any single ``d_bg`` level every voxel shares
    the same value, so which voxel first joins ``a`` and ``b`` never affects the
    returned value — only the (here unused) traversal order does.
    """
    if not occ[a_idx] or not occ[b_idx]:
        return float("nan")
    sl, offset = crop_slices(occ, margin=2)
    occ = occ[sl]
    a_idx = tuple(int(a_idx[d] - offset[d]) for d in range(3))
    b_idx = tuple(int(b_idx[d] - offset[d]) for d in range(3))
    d_bg = ndimage.distance_transform_edt(occ, sampling=(h, h, h))

    def connected_at(tau: float) -> bool:
        mask = occ & (d_bg >= tau)
        if not mask[a_idx] or not mask[b_idx]:
            return False
        labels, _ = ndimage.label(mask, structure=_STRUCT_6)
        return bool(labels[a_idx] == labels[b_idx])

    levels = np.unique(d_bg[occ])            # ascending; distinct d_bg values only
    if levels.size == 0 or not connected_at(levels[0]):
        return float("nan")
    lo, hi = 0, levels.size - 1              # connected_at(levels[lo]) holds
    while lo < hi:                            # find the largest index that still holds
        mid = (lo + hi + 1) // 2
        if connected_at(levels[mid]):
            lo = mid
        else:
            hi = mid - 1
    return float(levels[lo])


def center_component_labels(field: ComponentField, seed_idx: np.ndarray) -> np.ndarray:
    """Component label of every center's own voxel.

    A center's voxel is a seed, so its distance to the center set is 0 and it is
    occupied at every radius — the label is always defined, which is what makes
    QC-F5 true by construction and the center<->component map total.
    """
    return field.labels[seed_idx[:, 0], seed_idx[:, 1], seed_idx[:, 2]]


def residue_component_labels(centers: np.ndarray, residue_coords: np.ndarray,
                             covered: np.ndarray, rho: float,
                             center_labels: np.ndarray) -> np.ndarray:
    """Assign each covered residue the component of its NEAREST covering center.

    Residue-level membership is resolution-independent, which is exactly why the
    ARI stability term and the per-component residue lists are defined on it
    rather than on voxels.  Uncovered residues get label 0.
    """
    from ..utils.geometry import cross_distances

    labels = np.zeros(len(residue_coords), dtype=np.int64)
    if not covered.any() or len(centers) == 0:
        return labels
    d = cross_distances(residue_coords[covered], centers)
    d = np.where(d <= rho, d, np.inf)
    nearest = np.argmin(d, axis=1)
    labels[covered] = center_labels[nearest]
    return labels


def adjusted_rand(labels_a: np.ndarray, labels_b: np.ndarray,
                  mask: np.ndarray) -> float:
    """ARI between two residue partitions on their shared support.

    Fewer than two residues in common -> 0.0 by recorded convention: with no
    shared support there is no evidence of stability, and inventing 1.0 would
    reward a radius for having nothing to compare.
    """
    from sklearn.metrics import adjusted_rand_score

    if mask.sum() < 2:
        return 0.0
    return float(adjusted_rand_score(labels_a[mask], labels_b[mask]))
