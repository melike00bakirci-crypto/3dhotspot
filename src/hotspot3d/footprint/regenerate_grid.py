"""Regenerate an omitted per-radius voxel grid (Output Contract P5).

Per-radius occupancy grids are not stored: they are a deterministic function of
the center set, the radius and the voxel edge, so storing forty of them would add
bulk without adding evidence. What the decision rests on — per-radius component
membership and every derived metric — is retained in full.

This module is the recorded way to get any of them back:

    python -m hotspot3d.footprint.regenerate_grid \\
        --centers FULL_RESULTS/06_FINAL_HOTSPOTS/significant_hotspot_centers.tsv \\
        --coords  FULL_RESULTS/03_STRUCTURE_QC/residue_coordinates.tsv \\
        --radius <RHO> --voxel-h <H> --out grid_<RHO>.npz

It builds nothing of its own. The bounding box, the padding, the voxel edge and
the fallback rule all come from :func:`hotspot3d.footprint.sweep.build_sweep_grid`
— the very function the sweep called — and the occupancy is the rho-sublevel set
of the same single distance transform. A second grid builder here would be a
duplicate scientific implementation and the two could drift apart silently.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from ..utils.config import load_config
from ..utils.errors import BlockedError
from ..utils.io import read_tsv, read_tsv_columns
from .domain import derive_domain
from .grid import GridSpec, distance_field, occupancy_at
from .params import FootprintParams
from .sweep import build_sweep_grid

DEFAULT_CONFIG = "config/pipeline.yaml"

_CENTER_ID_ALIASES = ("center_residue_index", "residue_index", "center_residue")


def load_coordinates(coords_path: Path, universe_path: Path | None
                     ) -> tuple[dict[int, np.ndarray], np.ndarray]:
    """CA coordinates, and ``U_struct`` as the set that defines ``D_max``.

    Without an explicit ``--universe`` the universe is every row with a usable CA,
    which is what ``U_struct`` means (METHOD_SPEC II.0). Passing the file makes the
    choice explicit when a run needs it to be.
    """
    rows = read_tsv_columns(coords_path,
                            ["residue_index", "x_ca", "y_ca", "z_ca", "ca_usable"])
    coord_by_id: dict[int, np.ndarray] = {}
    for row in rows:
        rid, xyz = row.get("residue_index"), (row.get("x_ca"), row.get("y_ca"),
                                              row.get("z_ca"))
        if rid is None or any(v is None for v in xyz):
            continue
        flag = row.get("ca_usable")
        if flag is False:
            continue
        coord_by_id[int(rid)] = np.array([float(v) for v in xyz], dtype=np.float64)
    if not coord_by_id:
        raise BlockedError(f"no usable CA coordinates in {coords_path}")

    if universe_path is None:
        ids = sorted(coord_by_id)
    else:
        ids = sorted({int(r["residue_index"])
                      for r in read_tsv_columns(universe_path, ["residue_index"])
                      if r.get("residue_index") is not None})
        missing = [i for i in ids if i not in coord_by_id]
        if missing:
            raise BlockedError(
                f"{len(missing)} universe residues have no CA coordinate "
                f"(first: {missing[:5]})")
    return coord_by_id, np.vstack([coord_by_id[i] for i in ids])


def load_centers(centers_path: Path,
                 coord_by_id: dict[int, np.ndarray]) -> tuple[tuple[int, ...],
                                                              np.ndarray]:
    rows = read_tsv(centers_path)
    if not rows:
        raise BlockedError(f"{centers_path} contains no centers")
    column = next((name for name in _CENTER_ID_ALIASES if name in rows[0]), None)
    if column is None:
        raise BlockedError(
            f"{centers_path} has no center identifier column "
            f"(accepted: {list(_CENTER_ID_ALIASES)})")
    ids = sorted({int(row[column]) for row in rows if row.get(column) is not None})
    missing = [i for i in ids if i not in coord_by_id]
    if missing:
        raise BlockedError(f"centers {missing[:5]} have no CA coordinate")
    return tuple(ids), np.vstack([coord_by_id[i] for i in ids])


def regenerate(centers_path: Path, coords_path: Path, radius: float,
               voxel_h: float | None = None, universe_path: Path | None = None,
               config_path: Path | str = DEFAULT_CONFIG,
               ) -> tuple[GridSpec, np.ndarray, tuple[int, ...]]:
    """Rebuild ``FP(radius)`` exactly as the sweep built it.

    Returns ``(spec, occupancy, center_ids)``.
    """
    cfg = load_config(config_path)
    params = FootprintParams.from_config(cfg)
    if voxel_h is not None and float(voxel_h) != params.voxel_h_A:
        # An explicit edge is honoured, but the fallback rule and the padding stay
        # exactly as configured, so the geometry is still the sweep's geometry.
        params = replace_voxel_h(params, float(voxel_h))

    coord_by_id, universe_coords = load_coordinates(coords_path, universe_path)
    center_ids, center_coords = load_centers(centers_path, coord_by_id)

    domain = derive_domain(center_coords, universe_coords, params)
    if not domain.rho_min - 1e-9 <= float(radius) <= domain.rho_max + 1e-9:
        raise BlockedError(
            f"radius {radius} lies outside the derived II.9 domain "
            f"[{domain.rho_min}, {domain.rho_max}]. The domain is never widened to "
            f"accommodate a regeneration request."
        )

    spec = build_sweep_grid(center_coords, domain, params)      # THE grid builder
    occupancy = occupancy_at(distance_field(spec, center_coords), float(radius))
    return spec, occupancy, center_ids


def replace_voxel_h(params: FootprintParams, voxel_h: float) -> FootprintParams:
    import dataclasses

    return dataclasses.replace(params, voxel_h_A=voxel_h)


def write_npz(path: Path, spec: GridSpec, occupancy: np.ndarray, radius: float,
              center_ids: tuple[int, ...]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        occupancy_packed=np.packbits(occupancy.ravel()),
        shape=np.asarray(occupancy.shape, dtype=np.int64),
        origin=np.asarray(spec.origin, dtype=np.float64),
        voxel_h_A=np.asarray([spec.h], dtype=np.float64),
        r_fp_A=np.asarray([float(radius)], dtype=np.float64),
        center_residue_indices=np.asarray(center_ids, dtype=np.int64),
    )
    return path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m hotspot3d.footprint.regenerate_grid",
        description=("Regenerate a per-radius footprint occupancy grid that was "
                     "omitted under Output Contract P5."))
    parser.add_argument("--centers", required=True, type=Path,
                        help="06_FINAL_HOTSPOTS/significant_hotspot_centers.tsv")
    parser.add_argument("--coords", required=True, type=Path,
                        help="03_STRUCTURE_QC/residue_coordinates.tsv")
    parser.add_argument("--radius", required=True, type=float,
                        help="the footprint radius RHO in Angstrom")
    parser.add_argument("--voxel-h", type=float, default=None,
                        help="voxel edge in Angstrom (default: the configured value)")
    parser.add_argument("--universe", type=Path, default=None,
                        help="03_STRUCTURE_QC/positional_universe.tsv; defaults to "
                             "every coordinate row with a usable CA")
    parser.add_argument("--config", type=Path, default=Path(DEFAULT_CONFIG),
                        help="frozen pipeline config (padding, fallback rule)")
    parser.add_argument("--out", required=True, type=Path,
                        help="output .npz path")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    spec, occupancy, center_ids = regenerate(
        centers_path=args.centers, coords_path=args.coords, radius=args.radius,
        voxel_h=args.voxel_h, universe_path=args.universe,
        config_path=args.config)
    write_npz(args.out, spec, occupancy, args.radius, center_ids)
    print(f"regenerated FP(rho={args.radius:g} A) on a {list(spec.shape)} grid "
          f"at h={spec.h} A ({int(occupancy.sum())} occupied voxels) -> {args.out}")
    return 0


if __name__ == "__main__":                                   # pragma: no cover
    sys.exit(main())
