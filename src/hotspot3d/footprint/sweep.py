"""II.8 — the ascending multi-scale sweep with real previous-radius diffs.

¶23 is explicit that a radius "does not represent an isolated solution but rather a
continuation of the previous one". This module therefore walks the domain upward
holding the previous radius' state, and every diff quantity (merge events, neck
widths, dV, expansion rate, dn_comp, ARI) is computed against that retained state
rather than re-derived from scratch. One distance transform underlies all of it.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..utils.geometry import covered_mask
from . import descriptors as desc
from .domain import FootprintDomain
from .grid import GridSpec, distance_field, make_grid, occupancy_at, seed_indices
from .params import FootprintParams


@dataclass
class MergeEvent:
    """One newly formed connection between previously separate regions."""

    rho: float
    component_a_centers: tuple[int, ...]
    component_b_centers: tuple[int, ...]
    neck_half_width_A: float
    neck_width_A: float
    connection_type: str            # natural | artificial | undetermined
    erosion_depth_A: float

    def as_row(self) -> dict:
        return {
            "r_fp_A": self.rho,
            "component_a_center_residues": list(self.component_a_centers),
            "component_b_center_residues": list(self.component_b_centers),
            "neck_half_width_A": self.neck_half_width_A,
            "neck_width_A": self.neck_width_A,
            "connection_type": self.connection_type,
            "erosion_depth_A": self.erosion_depth_A,
        }


@dataclass
class RadiusState:
    """Everything measured at one radius, retained in full (nothing is pruned)."""

    rho: float
    volume: float
    surface_area: float
    n_components: int
    n_covered_residues: int
    coverage: float
    compactness: float
    convexity: float
    connectivity: float
    n_bridges: int
    bridging_index: float
    n_components_eroded: int
    covered: np.ndarray                       # bool over U_struct
    residue_labels: np.ndarray                # component id per residue (0 = none)
    center_labels: np.ndarray                 # component id per center
    mesh_ok: bool
    mesh_failure_reason: str
    mesh: dict | None
    # diff quantities, filled against the retained previous state
    delta_n_components: int = 0
    delta_volume: float = 0.0
    merge_events: list[MergeEvent] = field(default_factory=list)
    # second-pass quantities
    expansion_rate: float = float("nan")
    ari_prev: float = float("nan")
    ari_next: float = float("nan")
    stability: float = float("nan")
    abrupt_transition: bool = False

    @property
    def parsimony(self) -> float:
        return 1.0 - self.coverage

    @property
    def n_merge_events(self) -> int:
        return len(self.merge_events)


@dataclass
class SweepResult:
    spec: GridSpec
    dist: np.ndarray
    seed_idx: np.ndarray
    rows: list[RadiusState]
    merge_events: list[MergeEvent]
    monotonicity_violations: list[dict]

    def by_rho(self) -> dict[float, RadiusState]:
        return {row.rho: row for row in self.rows}


def build_sweep_grid(center_coords: np.ndarray, domain: FootprintDomain,
                     params: FootprintParams) -> GridSpec:
    """Bounding box of ``S`` expanded by ``rho_max + 2 A`` (II.8)."""
    return make_grid(center_coords,
                     pad=domain.rho_max + params.bbox_padding_A,
                     h=params.voxel_h_A,
                     h_fallback=params.voxel_h_fallback_A,
                     max_voxels=params.voxel_count_fallback_threshold)


def evaluate_radius(spec: GridSpec, dist: np.ndarray, seed_idx: np.ndarray,
                    center_coords: np.ndarray, universe_coords: np.ndarray,
                    rho: float, params: FootprintParams) -> RadiusState:
    """All single-radius descriptors. Shared by the sweep and by fixed-radius runs."""
    occ_full = occupancy_at(dist, rho)
    occupied = bool(occ_full.any())
    vol = desc.volume(occ_full, spec.h)

    # Everything below is local to the occupied region; cropping to it is exact.
    sl, offset = desc.crop_slices(occ_full, margin=2)
    occ = occ_full[sl]
    dist_crop = dist[sl]
    crop_origin = tuple(float(spec.origin[d] + offset[d] * spec.h) for d in range(3))

    field_ = desc.components(occ)
    area, mesh, mesh_fail = desc.surface_area(dist_crop, rho, spec.h,
                                              occupied=occupied)
    psi = desc.compactness(vol, area)
    conv = desc.convexity(occ, spec.h, crop_origin, vol)
    conn = desc.connectivity(field_)
    n_bridge, bi, n_eroded = desc.erosion_bridging(
        occ, spec.h, params.erosion_fraction * rho)

    covered = covered_mask(center_coords, universe_coords, rho)
    coverage = float(covered.sum() / len(universe_coords)) if len(universe_coords) else 0.0
    center_labels = desc.center_component_labels(field_, seed_idx - offset)
    residue_labels = desc.residue_component_labels(
        center_coords, universe_coords, covered, rho, center_labels)

    return RadiusState(
        rho=float(rho), volume=vol, surface_area=area,
        n_components=int(field_.n_components),
        n_covered_residues=int(covered.sum()), coverage=coverage,
        compactness=psi, convexity=conv, connectivity=conn,
        n_bridges=n_bridge, bridging_index=bi, n_components_eroded=n_eroded,
        covered=covered, residue_labels=residue_labels, center_labels=center_labels,
        mesh_ok=mesh is not None, mesh_failure_reason=mesh_fail or "NA", mesh=mesh,
    )


def run_sweep(center_coords: np.ndarray, universe_coords: np.ndarray,
              domain: FootprintDomain, params: FootprintParams) -> SweepResult:
    """Ascending sweep over the II.9 grid, tracking evolution radius by radius.

    Meshes are *not* retained across the sweep — only the boolean
    ``mesh_ok`` that QC-F3 depends on. The mesh of the selected radius is
    re-extracted once at export time from the same distance field, which bounds
    memory without losing anything the decision rests on.
    """
    spec = build_sweep_grid(center_coords, domain, params)
    dist = distance_field(spec, center_coords)          # THE single EDT
    seed_idx = seed_indices(spec, center_coords)

    rows: list[RadiusState] = []
    all_events: list[MergeEvent] = []
    violations: list[dict] = []
    previous: RadiusState | None = None

    for rho in domain.grid:
        state = evaluate_radius(spec, dist, seed_idx, center_coords,
                                universe_coords, rho, params)
        if previous is not None:
            state.delta_n_components = state.n_components - previous.n_components
            state.delta_volume = state.volume - previous.volume
            state.merge_events = _detect_merge_events(
                spec, dist, seed_idx, previous, state, rho, params)
            all_events.extend(state.merge_events)
            violations.extend(_check_monotonicity(previous, state))
        state.mesh = None
        rows.append(state)
        previous = state

    _second_pass(rows, domain.step, params)
    return SweepResult(spec=spec, dist=dist, seed_idx=seed_idx, rows=rows,
                       merge_events=all_events, monotonicity_violations=violations)


# --- diffs against the retained previous radius -----------------------------

def _center_groups(labels: np.ndarray) -> dict[int, tuple[int, ...]]:
    groups: dict[int, list[int]] = {}
    for center, label in enumerate(labels):
        groups.setdefault(int(label), []).append(int(center))
    return {k: tuple(sorted(v)) for k, v in sorted(groups.items())}


def _detect_merge_events(spec: GridSpec, dist: np.ndarray, seed_idx: np.ndarray,
                         previous: RadiusState, current: RadiusState, rho: float,
                         params: FootprintParams) -> list[MergeEvent]:
    """Groups of centers that were separate at ``rho - step`` and are now joined.

    Each new connection gets a measured neck and a natural/artificial verdict; a
    connection whose neck half-width is below the erosion depth ``0.15 * rho`` is
    precisely one that the QC-F2 erosion would break, so the two diagnostics agree
    by construction.
    """
    prev_groups = _center_groups(previous.center_labels)
    cur_groups = _center_groups(current.center_labels)

    # map each current component -> the previous groups it absorbed
    absorbed: dict[int, list[tuple[int, ...]]] = {}
    for prev_label, members in prev_groups.items():
        cur_label = int(current.center_labels[members[0]])
        absorbed.setdefault(cur_label, []).append(members)

    events: list[MergeEvent] = []
    depth = params.erosion_fraction * rho
    occ = None
    for cur_label, groups in sorted(absorbed.items()):
        if len(groups) < 2:
            continue
        groups = sorted(groups, key=lambda g: g[0])
        if occ is None:
            occ = occupancy_at(dist, rho)
        for a, b in zip(groups[:-1], groups[1:]):
            a_idx = tuple(int(v) for v in seed_idx[a[0]])
            b_idx = tuple(int(v) for v in seed_idx[b[0]])
            half = desc.bottleneck_radius(occ, spec.h, a_idx, b_idx)
            if np.isnan(half):
                ctype = "undetermined"
            elif half >= depth:
                ctype = "natural"
            else:
                ctype = "artificial"
            events.append(MergeEvent(
                rho=float(rho), component_a_centers=a, component_b_centers=b,
                neck_half_width_A=float(half), neck_width_A=float(2.0 * half),
                connection_type=ctype, erosion_depth_A=float(depth),
            ))
    return events


def _check_monotonicity(previous: RadiusState, current: RadiusState) -> list[dict]:
    """QC-C4/C5. A violation is recorded as a numerical fault — never smoothed."""
    out: list[dict] = []
    if current.n_components > previous.n_components:
        out.append({
            "rule": "QC_C4_n_comp_monotone_non_increasing",
            "r_fp_A": current.rho, "previous_r_fp_A": previous.rho,
            "previous_value": previous.n_components, "value": current.n_components,
        })
    if current.volume < previous.volume - 1e-9:
        out.append({
            "rule": "QC_C5_volume_monotone_non_decreasing",
            "r_fp_A": current.rho, "previous_r_fp_A": previous.rho,
            "previous_value": previous.volume, "value": current.volume,
        })
    return out


def _second_pass(rows: list[RadiusState], step: float,
                 params: FootprintParams) -> None:
    """Expansion rate, ARI stability and abrupt-transition flags.

    ER needs ``rho + step`` and stability needs both neighbours, so these are the
    only quantities that cannot be produced during the ascent itself.
    """
    n = len(rows)
    for i, row in enumerate(rows):
        lo = rows[i - 1] if i > 0 else None
        hi = rows[i + 1] if i + 1 < n else None
        row.expansion_rate = _expansion_rate(lo, row, hi, step)
        if lo is not None:
            shared = lo.covered & row.covered
            row.ari_prev = desc.adjusted_rand(lo.residue_labels, row.residue_labels,
                                              shared)
        if hi is not None:
            shared = hi.covered & row.covered
            row.ari_next = desc.adjusted_rand(hi.residue_labels, row.residue_labels,
                                              shared)
        available = [v for v in (row.ari_prev, row.ari_next) if not np.isnan(v)]
        row.stability = float(np.mean(available)) if available else 0.0

    finite_er = np.array([r.expansion_rate for r in rows
                          if np.isfinite(r.expansion_rate)], dtype=np.float64)
    median_er = float(np.median(finite_er)) if finite_er.size else float("nan")
    for row in rows:
        abrupt = abs(row.delta_n_components) >= params.abrupt_dncomp_threshold
        if np.isfinite(median_er) and median_er > 0 and np.isfinite(row.expansion_rate):
            abrupt = abrupt or (row.expansion_rate
                                >= params.abrupt_er_median_multiple * median_er)
        row.abrupt_transition = bool(abrupt)


def _expansion_rate(lo: RadiusState | None, mid: RadiusState,
                    hi: RadiusState | None, step: float) -> float:
    """``ER = [ln V(rho+step) - ln V(rho-step)] / (2 step)``, one-sided at the ends."""
    def ln(state: RadiusState) -> float:
        return float(np.log(state.volume)) if state.volume > 0 else float("nan")

    if lo is not None and hi is not None:
        a, b = ln(lo), ln(hi)
        return (b - a) / (2.0 * step) if np.isfinite(a) and np.isfinite(b) else float("nan")
    if hi is not None:
        a, b = ln(mid), ln(hi)
        return (b - a) / step if np.isfinite(a) and np.isfinite(b) else float("nan")
    if lo is not None:
        a, b = ln(lo), ln(mid)
        return (b - a) / step if np.isfinite(a) and np.isfinite(b) else float("nan")
    return float("nan")
