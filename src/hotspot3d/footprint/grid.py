"""Occupancy grid and the SINGLE Euclidean distance transform (METHOD_SPEC II.8).

The whole multi-scale sweep rests on one observation: the footprint at radius rho
is exactly the rho-sublevel set of the distance field to the hotspot centers,

    FP(rho) = { x : d(x, S) <= rho }

so **one** EDT computed against the seed voxels serves every radius in the domain.
Recomputing a transform per radius would be both wasteful and a source of drift
between radii; the sweep is an evolution of one object, not a stack of unrelated
reconstructions (¶23).

The grid is deterministic given (points, pad, h): origin, shape and voxel centers
are fixed by the bounding box, so two runs — and Phase C versus Phase D — place
the same geometry in the same cells.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage

from ..utils.errors import BlockedError


@dataclass(frozen=True)
class GridSpec:
    """A deterministic axis-aligned voxel grid.

    Voxel ``(i, j, k)`` spans ``[origin + idx*h, origin + (idx+1)*h)`` and its
    center sits at ``origin + (idx + 0.5) * h``.
    """

    origin: tuple[float, float, float]
    shape: tuple[int, int, int]
    h: float
    h_requested: float
    fallback_applied: bool

    @property
    def n_voxels(self) -> int:
        return int(self.shape[0]) * int(self.shape[1]) * int(self.shape[2])

    @property
    def voxel_volume(self) -> float:
        return float(self.h ** 3)

    def index_of(self, points: np.ndarray) -> np.ndarray:
        """Voxel indices containing ``points``, clipped into the grid."""
        pts = np.atleast_2d(np.asarray(points, dtype=np.float64))
        idx = np.floor((pts - np.asarray(self.origin)) / self.h).astype(np.int64)
        upper = np.asarray(self.shape, dtype=np.int64) - 1
        return np.clip(idx, 0, upper)

    def centers_of(self, idx: np.ndarray) -> np.ndarray:
        """World coordinates of the centers of the given voxel indices."""
        idx = np.atleast_2d(np.asarray(idx, dtype=np.float64))
        return np.asarray(self.origin) + (idx + 0.5) * self.h

    def as_dict(self) -> dict:
        return {
            "origin_A": [float(v) for v in self.origin],
            "shape": [int(v) for v in self.shape],
            "voxel_h_A": self.h,
            "voxel_h_requested_A": self.h_requested,
            "voxel_resolution_fallback_applied": self.fallback_applied,
            "n_voxels": self.n_voxels,
            "voxel_volume_A3": self.voxel_volume,
        }


def make_grid(points: np.ndarray, pad: float, h: float, h_fallback: float,
              max_voxels: int) -> GridSpec:
    """Bounding box of ``points`` expanded by ``pad``, voxelized at ``h``.

    Above ``max_voxels`` the edge falls back to ``h_fallback`` automatically and
    the fact is recorded (A18) — a silent resolution change would make volumes
    incomparable between runs.
    """
    pts = np.atleast_2d(np.asarray(points, dtype=np.float64))
    if pts.size == 0:
        raise BlockedError("cannot build an occupancy grid from an empty point set")
    if not np.all(np.isfinite(pts)):
        raise BlockedError("non-finite coordinate reached the occupancy grid")

    lo = pts.min(axis=0) - pad
    hi = pts.max(axis=0) + pad

    def _shape(edge: float) -> tuple[int, int, int]:
        extent = hi - lo
        return tuple(int(np.ceil(extent[d] / edge)) + 1 for d in range(3))

    shape = _shape(h)
    used, fallback = h, False
    if int(np.prod(shape, dtype=np.int64)) > int(max_voxels):
        used, fallback = h_fallback, True
        shape = _shape(used)

    return GridSpec(origin=(float(lo[0]), float(lo[1]), float(lo[2])),
                    shape=shape, h=float(used), h_requested=float(h),
                    fallback_applied=fallback)


def seed_mask(spec: GridSpec, points: np.ndarray) -> np.ndarray:
    """Boolean grid marking the voxel that contains each point."""
    mask = np.zeros(spec.shape, dtype=bool)
    idx = spec.index_of(points)
    mask[idx[:, 0], idx[:, 1], idx[:, 2]] = True
    return mask


def distance_field(spec: GridSpec, points: np.ndarray) -> np.ndarray:
    """THE distance transform: physical distance from every voxel to nearest seed.

    Computed exactly once per center set and reused for every radius in the sweep.
    """
    seeds = seed_mask(spec, points)
    if not seeds.any():
        raise BlockedError("distance field requested with no seed voxels")
    return ndimage.distance_transform_edt(~seeds, sampling=(spec.h, spec.h, spec.h))


def occupancy_at(dist: np.ndarray, rho: float) -> np.ndarray:
    """``FP(rho)`` as the rho-sublevel set of the (single) distance field."""
    return dist <= rho


def seed_indices(spec: GridSpec, points: np.ndarray) -> np.ndarray:
    """Voxel indices of the seeds, in input order (used for center<->component)."""
    return spec.index_of(points)
