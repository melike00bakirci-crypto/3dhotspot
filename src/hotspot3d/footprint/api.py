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

from ..utils.errors import BlockedError
from ..utils.hashing import code_version, sha256_bytes
from ..utils.multiobjective import SelectionResult
from . import qc as qc_mod
from . import selection as sel
from .domain import (FootprintDomain, derive_containment_domain, derive_domain)
from .grid import GridSpec, distance_field, occupancy_at, seed_indices
from .params import FootprintParams
from .qc import QCVerdict
from .sweep import (RadiusState, SweepResult, build_sweep_grid,
                    evaluate_radius, run_sweep)

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
    #: [DECISION-RFP-CONTAINMENT-0001] the containment record; None on the retired
    #: search path, which is now guarded shut.
    containment: dict | None = None

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


#: DECISION-RFP-CONTAINMENT-0001 — reported as a flag, never used to raise
#: the radius.
DETECTION_SCALE_FLOOR_A = 7.0


def build_footprint(centers: CenterSet, universe: Universe,
                    params: FootprintParams) -> FootprintSolution:
    """Full II.9 -> II.8 -> II.10 construction: domain, sweep, QC, selection.

    Raises :class:`NegativeResult` when the domain itself is empty; returns a
    solution with ``r_fp = None`` when the domain is non-empty but no candidate is
    admissible. Both are scientific outcomes, not exceptions in the ordinary sense.
    """
    raise BlockedError(
        "DECISION-RFP-CONTAINMENT-0001 retired the independent r_fp search. "
        "build_footprint() performed domain derivation, a multi-scale sweep, QC "
        "filtering and Pareto/distance-to-utopia selection; r_fp is now DETERMINED "
        "by containment and there is nothing to select. Use "
        "build_footprint_contained(centers, universe, params, evidence). This "
        "function is kept, and kept failing, so that a caller which still expects a "
        "search gets told rather than silently handed a one-candidate result."
    )


def _retired_build_footprint_body(centers: CenterSet, universe: Universe,
                                  params: FootprintParams) -> FootprintSolution:
    """Unreachable. Retained verbatim as the record of the retired procedure."""
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


def build_footprint_contained(centers: CenterSet, universe: Universe,
                              params: FootprintParams,
                              evidence) -> FootprintSolution:
    """DECISION-RFP-CONTAINMENT-0001 — the whole of Stage C's radius determination.

    The domain has exactly one member, so the sweep evaluates one radius, QC judges
    one candidate and the selection is trivial. Every downstream export keeps its
    existing shape: the geometric descriptors are computed and reported, and they
    decide nothing because there is no alternative to prefer.

    ``evidence`` is an ``EvidenceSet``; no radius is passed in and ``r_hot`` is not
    an argument here or in ``footprint.domain``.
    """
    if len(centers) == 0:
        raise BlockedError(
            "build_footprint_contained called with an empty center set")

    domain, record = derive_containment_domain(centers.coords, universe.coords,
                                               evidence, params)
    sweep = run_sweep(centers.coords, universe.coords, domain, params)

    if sweep.monotonicity_violations:
        raise BlockedError(
            f"MONOTONICITY VIOLATION on a one-point domain: "
            f"{sweep.monotonicity_violations[:3]}. With a single radius there is "
            f"nothing to be non-monotone against; this is a computational fault."
        )

    verdicts = {row.rho: qc_mod.evaluate(row, params, domain) for row in sweep.rows}
    constraint_failed, excessive_fraction = qc_mod.coverage_constraint_failure(
        list(verdicts.values()), params)

    candidates = sel.build_candidates(sweep.rows, verdicts)
    result = sel.run_selection(candidates, params, domain.step)

    selected_state = None
    if result.selected_key is not None:
        selected_state = sweep.by_rho()[result.selected_key]
        qc_mod.assert_qc_f4(selected_state.rho, domain)
        qc_mod.assert_qc_f5(selected_state, np.asarray(centers.universe_rows,
                                                       dtype=np.int64))

    # The coverage rule never shrinks the radius: L3 would break. When it fires the
    # run escalates and completes with warnings, carrying the radius it computed.
    record["coverage_escalated"] = bool(constraint_failed)
    record["detection_scale_below_resolution"] = bool(
        record["r_contain_exact_A"] < DETECTION_SCALE_FLOOR_A)

    return FootprintSolution(
        domain=domain, sweep=sweep, verdicts=verdicts, selection=result,
        boundary={}, spec=sweep.spec,
        r_fp=result.selected_key, selected=selected_state,
        excessive_coverage_fraction=excessive_fraction,
        coverage_constraint_failure=constraint_failed,
        n_admissible=result.n_admissible,
        containment=record,
    )


def reconstruct_at(centers: CenterSet, universe: Universe, params: FootprintParams,
                   rho: float) -> FixedReconstruction:
    """Rebuild the footprint at an externally supplied radius (no selection).

    The radius is an input, so QC-F4 is *reported* rather than asserted: the whole
    purpose of the secondary reconstruction is to hold the radius fixed at the
    original ``r_fp``. Under DECISION-RFP-CONTAINMENT-0001 that radius is the
    containment minimum of the UNPERTURBED center set, held fixed across iterations
    exactly as DECISION-STAGE-D-FIXED-RFP-0001 requires.
    """
    if len(centers) == 0:
        raise BlockedError("reconstruct_at called with an empty center set")

    # [DECISION-RFP-CONTAINMENT-0001] The grid is sized from rho alone.
    #
    # This used to call derive_domain() and pad to max(domain.rho_max, rho). Two
    # consequences, both measured on archived runs:
    #
    #   * the retired 25 A ceiling inflated the box whatever rho was — the voxel
    #     count ran 3-10x higher than needed (CACNA1A 11.1M vs 3.7M; PTPN11 2.9M vs
    #     0.3M), and Stage D pays that N_CAP times per gene;
    #   * worse, it made a radius's MEASURED geometry depend on the domain it sat in.
    #     PTPN11's volume came out 7915 against 7947 for the same rho. That is the
    #     defect recorded under "Defect exposed, not caused" in
    #     DECISION-FOOTPRINT-DOMAIN-0001, and sizing from rho closes it: the box now
    #     depends only on the centers and the radius actually being evaluated.
    #
    # There is no domain under containment, so the former inside_iteration_domain
    # flag had no question left to answer. Nothing read it; it is gone.
    spec = build_sweep_grid(centers.coords, float(rho), params)
    dist = distance_field(spec, centers.coords)
    seed_idx = seed_indices(spec, centers.coords)
    state = evaluate_radius(spec, dist, seed_idx, centers.coords, universe.coords,
                            float(rho), params)
    return FixedReconstruction(rho=float(rho), spec=spec, dist=dist, state=state)


def occupancy_on_grid(spec: GridSpec, centers: CenterSet, rho: float) -> np.ndarray:
    """``FP(rho)`` for ``centers`` evaluated on an EXTERNAL reference grid.

    Voxel overlap between two footprints is only meaningful in a shared frame, and
    each iteration derives its own bounding box. Both footprints are therefore
    re-evaluated on one common grid before any voxel Jaccard is computed.
    """
    dist = distance_field(spec, centers.coords)
    return occupancy_at(dist, rho)
