"""THE frozen, versioned, deterministic footprint API.

This is the single entry point through which a footprint is ever constructed —
Phase C calls it once to produce ``FP_original``, and Phase D calls the *same*
function for the mandatory baseline reproduction and for every perturbation
iteration. That is the whole point: if Phase D used a second implementation, the
comparison would measure the difference between two code paths rather than the
effect of removing centers.

The API is a pure function of its arguments. It performs no I/O, reads no config
file, consults no run context, and takes ``r_hot`` nowhere in its signature — the
mechanical form of F10.

Determinism: identical inputs give byte-identical outputs. The only stochastic
element anywhere in Stages C/D is Phase D's subset sampling, which is seeded by
context string, never by execution order.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..utils.errors import BlockedError, NegativeResult
from ..utils.hashing import code_version, sha256_bytes
from ..utils.multiobjective import SelectionResult
from . import qc as qc_mod
from . import selection as sel
from .domain import FootprintDomain, derive_domain
from .grid import GridSpec, distance_field, make_grid, occupancy_at, seed_indices
from .params import FootprintParams
from .qc import QCVerdict
from .sweep import RadiusState, SweepResult, evaluate_radius, run_sweep

#: Bumped only by a pre-registered methodological revision.
FOOTPRINT_API_VERSION = "1.0.0"


def footprint_code_version() -> str:
    """Content hash of the footprint package alone.

    Scoped deliberately: Phase D must detect drift in *this* code, and would be
    uselessly tripped by an unrelated edit elsewhere in the pipeline.
    """
    return code_version(Path(__file__).resolve().parent)


@dataclass(frozen=True)
class CenterSet:
    """An immutable set of significant hotspot centers.

    Phase D removes *elements* of this set to form ``S_iter``; it never edits one.
    """

    ids: tuple[int, ...]                  # residue indices, ascending
    coords: np.ndarray                    # (n, 3) CA coordinates, same order
    universe_rows: tuple[int, ...]        # row of each center within U_struct

    def __post_init__(self):
        if len(self.ids) != len(self.coords) or len(self.ids) != len(self.universe_rows):
            raise BlockedError("center id / coordinate / universe-row length mismatch")
        if list(self.ids) != sorted(self.ids):
            raise BlockedError("center ids must be in ascending order (determinism)")

    def __len__(self) -> int:
        return len(self.ids)

    def without(self, removed: tuple[int, ...]) -> "CenterSet":
        """``S_iter = S \\ T`` — a new object; ``S`` itself is never mutated."""
        drop = set(int(r) for r in removed)
        keep = [i for i, cid in enumerate(self.ids) if int(cid) not in drop]
        return CenterSet(
            ids=tuple(self.ids[i] for i in keep),
            coords=self.coords[keep],
            universe_rows=tuple(self.universe_rows[i] for i in keep),
        )

    def digest(self) -> str:
        payload = ";".join(f"{i}:{x:.6f},{y:.6f},{z:.6f}"
                           for i, (x, y, z) in zip(self.ids, self.coords))
        return sha256_bytes(payload.encode("utf-8"))


@dataclass(frozen=True)
class Universe:
    """``U_struct`` — every residue with a usable CA. Read-only for this stage."""

    ids: tuple[int, ...]
    coords: np.ndarray

    def __len__(self) -> int:
        return len(self.ids)


@dataclass
class FootprintSolution:
    """The complete Phase C solution, or an honest statement that none exists."""

    domain: FootprintDomain
    sweep: SweepResult
    verdicts: dict[float, QCVerdict]
    selection: SelectionResult
    boundary: dict
    spec: GridSpec
    r_fp: float | None
    selected: RadiusState | None
    excessive_coverage_fraction: float
    coverage_constraint_failure: bool
    n_admissible: int

    # -- derived, decision-bearing views ------------------------------------
    def covered_residue_ids(self, universe: Universe) -> tuple[int, ...]:
        if self.selected is None:
            return ()
        return tuple(int(universe.ids[i])
                     for i in np.nonzero(self.selected.covered)[0])

    def occupancy(self) -> np.ndarray:
        if self.selected is None:
            raise BlockedError("no footprint was selected; occupancy is undefined")
        return occupancy_at(self.sweep.dist, self.selected.rho)

    def digest(self, universe: Universe) -> dict[str, str | float | int | None]:
        """Decision-bearing fingerprint used by the baseline reproduction check."""
        residues = self.covered_residue_ids(universe)
        residue_hash = sha256_bytes(
            ",".join(str(r) for r in residues).encode("utf-8"))
        if self.selected is None:
            occ_hash = "NA"
        else:
            occ = self.occupancy()
            occ_hash = sha256_bytes(np.packbits(occ.ravel()).tobytes())
        return {
            "r_fp_A": self.r_fp,
            "n_covered_residues": len(residues),
            "covered_residues_sha256": residue_hash,
            "occupancy_sha256": occ_hash,
            "volume_A3": self.selected.volume if self.selected else None,
            "surface_area_A2": self.selected.surface_area if self.selected else None,
            "n_components": self.selected.n_components if self.selected else None,
            "grid_shape": list(self.spec.shape),
            "voxel_h_A": self.spec.h,
        }


@dataclass
class FixedReconstruction:
    """Geometry at an externally fixed radius — the II.11 secondary reconstruction.

    It isolates pure center-loss geometry from radius re-selection, so a footprint
    that changes only because ``r_fp^iter`` moved is distinguishable from one that
    changed because the centers that shaped it are gone.
    """

    rho: float
    spec: GridSpec
    dist: np.ndarray
    state: RadiusState
    inside_iteration_domain: bool


def build_footprint(centers: CenterSet, universe: Universe,
                    params: FootprintParams) -> FootprintSolution:
    """Full II.9 -> II.8 -> II.10 construction: domain, sweep, QC, selection.

    Raises :class:`NegativeResult` when the domain itself is empty; returns a
    solution with ``r_fp = None`` when the domain is non-empty but no candidate is
    admissible. Both are scientific outcomes, not exceptions in the ordinary sense.
    """
    if len(centers) == 0:
        raise BlockedError("build_footprint called with an empty center set")

    domain = derive_domain(centers.coords, universe.coords, params)
    sweep = run_sweep(centers.coords, universe.coords, domain, params)

    # QC-C4/C5. Both are mathematically guaranteed — every FP(rho) is a sublevel
    # set of ONE distance field, so the sets are nested and components can only
    # merge. A violation is therefore a numerical fault in this code, and it is
    # raised and investigated, never smoothed away.
    if sweep.monotonicity_violations:
        raise BlockedError(
            f"MONOTONICITY VIOLATION — the sweep is not nested: "
            f"{sweep.monotonicity_violations[:3]}. V(rho) must be non-decreasing and "
            f"n_comp(rho) non-increasing by construction; this is a computational "
            f"fault, not a property of the protein, and it is never smoothed."
        )

    verdicts = {row.rho: qc_mod.evaluate(row, params, domain) for row in sweep.rows}
    constraint_failed, excessive_fraction = qc_mod.coverage_constraint_failure(
        list(verdicts.values()), params)

    candidates = sel.build_candidates(sweep.rows, verdicts)
    result = sel.run_selection(candidates, params, domain.step)
    boundary = sel.run_boundary_diagnostic(domain, result, params)

    selected_state = None
    if result.selected_key is not None:
        selected_state = sweep.by_rho()[result.selected_key]
        qc_mod.assert_qc_f4(selected_state.rho, domain)
        qc_mod.assert_qc_f5(selected_state, np.asarray(centers.universe_rows,
                                                       dtype=np.int64))

    return FootprintSolution(
        domain=domain, sweep=sweep, verdicts=verdicts, selection=result,
        boundary=boundary, spec=sweep.spec,
        r_fp=result.selected_key, selected=selected_state,
        excessive_coverage_fraction=excessive_fraction,
        coverage_constraint_failure=constraint_failed,
        n_admissible=result.n_admissible,
    )


def reconstruct_at(centers: CenterSet, universe: Universe, params: FootprintParams,
                   rho: float) -> FixedReconstruction:
    """Rebuild the footprint at an externally supplied radius (no selection).

    The radius is an input, so QC-F4 is *reported* rather than asserted: the whole
    purpose of the secondary reconstruction is to hold the radius fixed at the
    original ``r_fp`` even when the iteration's own domain would not have offered it.
    """
    if len(centers) == 0:
        raise BlockedError("reconstruct_at called with an empty center set")

    try:
        domain = derive_domain(centers.coords, universe.coords, params)
        inside = bool(domain.rho_min - 1e-9 <= rho <= domain.rho_max + 1e-9)
        pad_radius = max(domain.rho_max, float(rho))
    except NegativeResult:                  # empty domain -> still reconstructable
        inside = False
        pad_radius = float(rho)

    spec = make_grid(centers.coords, pad=pad_radius + params.bbox_padding_A,
                     h=params.voxel_h_A, h_fallback=params.voxel_h_fallback_A,
                     max_voxels=params.voxel_count_fallback_threshold)
    dist = distance_field(spec, centers.coords)
    seed_idx = seed_indices(spec, centers.coords)
    state = evaluate_radius(spec, dist, seed_idx, centers.coords, universe.coords,
                            float(rho), params)
    return FixedReconstruction(rho=float(rho), spec=spec, dist=dist, state=state,
                               inside_iteration_domain=inside)


def occupancy_on_grid(spec: GridSpec, centers: CenterSet, rho: float) -> np.ndarray:
    """``FP(rho)`` for ``centers`` evaluated on an EXTERNAL reference grid.

    Voxel overlap between two footprints is only meaningful in a shared frame, and
    each iteration derives its own bounding box. Both footprints are therefore
    re-evaluated on one common grid before any voxel Jaccard is computed.
    """
    dist = distance_field(spec, centers.coords)
    return occupancy_at(dist, rho)
