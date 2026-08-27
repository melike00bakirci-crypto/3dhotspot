"""Canonical Phase C outputs (Output Contract Part IX).

Two rules shape this module:

* ``r_fp_scan.tsv`` is the **single decision trace**. Every other tabular artifact
  in ``07_FOOTPRINT_RADIUS/`` is a projection of rows that already exist here —
  ``footprint_admissibility.tsv`` is copied verbatim from the scan rows, never
  recomputed (P1). No value is produced twice anywhere in this stage.
* Per-radius voxel grids are **not** stored (P5): they are exactly regenerable
  from the center set and the radius, and the regeneration recipe is written into
  ``footprint_domain.json`` for ``12_REPRODUCIBILITY``. What the decision actually
  rests on — per-radius component membership and every derived metric — is
  retained in full.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..utils.io import write_json, write_text, write_tsv
from .api import FootprintSolution, Universe
from .descriptors import surface_area
from .domain import FootprintDomain
from .params import FootprintParams

R_FP_SCAN_COLUMNS = [
    "r_fp_A", "volume_A3", "surface_area_A2", "n_components", "n_covered_residues",
    "coverage_fraction", "compactness", "convexity", "connectivity",
    "expansion_rate", "n_bridges", "bridging_index", "ari_prev", "ari_next",
    "stability", "parsimony", "n_merge_events", "abrupt_transition",
    "compactness_normalized", "connectivity_normalized", "stability_normalized",
    "parsimony_normalized", "qc_status", "qc_failure_reason",
    "excessive_footprint_coverage", "pareto_member", "distance_to_ideal",
    "selected", "boundary_warning",
]

ADMISSIBILITY_COLUMNS = [
    "r_fp_A", "coverage_fraction", "bridging_index", "volume_A3",
    "qc_status", "qc_failure_reason", "excessive_footprint_coverage",
]

MERGE_EVENT_COLUMNS = [
    "r_fp_A", "component_a_center_residues", "component_b_center_residues",
    "neck_half_width_A", "neck_width_A", "connection_type", "erosion_depth_A",
]

COMPONENT_MEMBERSHIP_COLUMNS = [
    "residue_index", "component_id", "is_significant_center",
    "nearest_center_residue_index", "distance_to_nearest_center_A",
]

FOOTPRINT_RESIDUE_COLUMNS = [
    "residue_index", "component_id", "is_significant_center",
    "nearest_center_residue_index", "distance_to_nearest_center_A", "cohort_class",
]

FOOTPRINT_GEOMETRY_COLUMNS = [
    "r_fp_A", "volume_A3", "surface_area_A2", "n_components",
    "n_covered_residues", "coverage_fraction", "compactness", "convexity",
    "connectivity", "stability", "parsimony", "expansion_rate", "n_bridges",
    "bridging_index", "n_components_eroded", "n_merge_events",
    "abrupt_transition", "voxel_h_A", "voxel_resolution_fallback_applied",
]

FOOTPRINT_COMPONENT_COLUMNS = [
    "component_id", "n_voxels", "volume_A3", "n_residues", "n_centers",
    "center_residues", "centroid_x", "centroid_y", "centroid_z",
]


def _radius_tag(rho: float) -> str:
    """Windows-safe, lexicographically sortable per-radius filename token."""
    return f"{rho:05.2f}"


def build_scan_rows(solution: FootprintSolution) -> list[dict]:
    """THE decision trace. Every candidate, admissible or not, with every metric."""
    result = solution.selection
    band = set(solution.boundary.get("band_members", []))
    rows = []
    for state in solution.sweep.rows:
        verdict = solution.verdicts[state.rho]
        norm = result.normalized.get(state.rho, {})
        rows.append({
            "r_fp_A": state.rho,
            "volume_A3": state.volume,
            "surface_area_A2": state.surface_area,
            "n_components": state.n_components,
            "n_covered_residues": state.n_covered_residues,
            "coverage_fraction": state.coverage,
            "compactness": state.compactness,
            "convexity": state.convexity,
            "connectivity": state.connectivity,
            "expansion_rate": state.expansion_rate,
            "n_bridges": state.n_bridges,
            "bridging_index": state.bridging_index,
            "ari_prev": state.ari_prev,
            "ari_next": state.ari_next,
            "stability": state.stability,
            "parsimony": state.parsimony,
            "n_merge_events": state.n_merge_events,
            "abrupt_transition": state.abrupt_transition,
            # NA for inadmissible rows: normalization is defined over the
            # admissible set only (F15), so a normalized value would be fiction.
            "compactness_normalized": norm.get("compactness"),
            "connectivity_normalized": norm.get("connectivity"),
            "stability_normalized": norm.get("stability"),
            "parsimony_normalized": norm.get("parsimony"),
            "qc_status": verdict.status,
            "qc_failure_reason": verdict.reason_text,
            "excessive_footprint_coverage": verdict.excessive_coverage,
            "pareto_member": state.rho in set(result.pareto_members),
            "distance_to_ideal": result.distances.get(state.rho),
            "selected": (result.selected_key is not None
                         and abs(state.rho - result.selected_key) < 1e-9),
            "boundary_warning": state.rho in band,
        })
    return rows


def write_scan(stage_dir: Path, rows: list[dict]) -> Path:
    return write_tsv(stage_dir / "r_fp_scan.tsv", rows, R_FP_SCAN_COLUMNS)


def write_admissibility(stage_dir: Path, rows: list[dict]) -> Path:
    """A verbatim projection of ``r_fp_scan.tsv`` — nothing is computed here (P1)."""
    projected = [{col: row[col] for col in ADMISSIBILITY_COLUMNS} for row in rows]
    return write_tsv(stage_dir / "footprint_admissibility.tsv", projected,
                     ADMISSIBILITY_COLUMNS)


def write_merge_events(stage_dir: Path, solution: FootprintSolution,
                       center_ids: tuple[int, ...]) -> Path:
    rows = []
    for event in solution.sweep.merge_events:
        row = event.as_row()
        row["component_a_center_residues"] = [int(center_ids[i])
                                              for i in event.component_a_centers]
        row["component_b_center_residues"] = [int(center_ids[i])
                                              for i in event.component_b_centers]
        rows.append(row)
    return write_tsv(stage_dir / "merge_events.tsv", rows, MERGE_EVENT_COLUMNS)


def write_component_memberships(stage_dir: Path, solution: FootprintSolution,
                                universe: Universe,
                                center_ids: tuple[int, ...]) -> list[Path]:
    """Per-radius component membership — retained in full for every radius."""
    out_dir = stage_dir / "footprint_components"
    out_dir.mkdir(parents=True, exist_ok=True)
    center_set = set(int(c) for c in center_ids)
    written = []
    for state in solution.sweep.rows:
        rows = _membership_rows(state, universe, center_set)
        written.append(write_tsv(out_dir / f"r_{_radius_tag(state.rho)}.tsv",
                                 rows, COMPONENT_MEMBERSHIP_COLUMNS))
    return written


def _membership_rows(state, universe: Universe, center_set: set) -> list[dict]:
    from ..utils.geometry import cross_distances

    idx = np.nonzero(state.covered)[0]
    if idx.size == 0:
        return []
    centers_xyz = universe.coords[[i for i, rid in enumerate(universe.ids)
                                   if int(rid) in center_set]]
    center_ids_sorted = [rid for rid in universe.ids if int(rid) in center_set]
    d = cross_distances(universe.coords[idx], centers_xyz)
    nearest = np.argmin(d, axis=1)
    rows = []
    for pos, residue_row in enumerate(idx):
        rid = int(universe.ids[residue_row])
        rows.append({
            "residue_index": rid,
            "component_id": int(state.residue_labels[residue_row]),
            "is_significant_center": rid in center_set,
            "nearest_center_residue_index": int(center_ids_sorted[nearest[pos]]),
            "distance_to_nearest_center_A": float(d[pos, nearest[pos]]),
        })
    return rows


def voxel_h_sensitivity(solution: FootprintSolution, centers_xyz: np.ndarray,
                        params: FootprintParams) -> dict:
    """QC-C15 — how much the reported volume depends on the voxel edge.

    Volume is the one descriptor that scales with discretization, so the decision
    trace has to say by how much. The footprint at ``r_fp`` is rebuilt once at the
    alternative edge and the two volumes are reported side by side; nothing is
    corrected, the number is simply stated.
    """
    from .grid import distance_field, make_grid, occupancy_at

    if solution.selected is None:
        return {"evaluated": False, "reason": "no radius selected"}
    alternative = (params.voxel_h_fallback_A if solution.spec.h == params.voxel_h_A
                   else params.voxel_h_A)
    spec = make_grid(centers_xyz,
                     pad=solution.domain.rho_max + params.bbox_padding_A,
                     h=alternative, h_fallback=alternative,
                     max_voxels=params.voxel_count_fallback_threshold * 1000)
    occ = occupancy_at(distance_field(spec, centers_xyz), solution.r_fp)
    volume_alt = float(int(occ.sum()) * alternative ** 3)
    volume = float(solution.selected.volume)
    return {
        "evaluated": True,
        "voxel_h_used_A": solution.spec.h,
        "voxel_h_alternative_A": alternative,
        "volume_used_A3": volume,
        "volume_alternative_A3": volume_alt,
        "relative_difference": (abs(volume_alt - volume) / volume
                                if volume > 0 else None),
        "fallback_applied": solution.spec.fallback_applied,
        "note": ("reported for transparency; the volume is NOT adjusted and the "
                 "voxel edge is not tuned to a preferred answer"),
    }


def write_domain(stage_dir: Path, solution: FootprintSolution,
                 universe: Universe, center_ids: tuple[int, ...],
                 params: FootprintParams, h_sensitivity: dict | None = None) -> Path:
    domain: FootprintDomain = solution.domain
    payload = domain.as_dict(center_ids=list(center_ids))
    payload["grid_spec"] = solution.spec.as_dict()
    payload["voxel_grid_retention"] = {
        "per_radius_voxel_grids_stored": False,
        "reason": ("P5 — deterministically regenerable from the center set and the "
                   "radius; per-radius component memberships and metrics, which the "
                   "decision rests on, are retained in full."),
        "regeneration_recipe": (
            "from hotspot3d.footprint.grid import make_grid, distance_field, "
            "occupancy_at; spec = make_grid(center_coords, pad=rho_max + "
            f"{params.bbox_padding_A}, h={params.voxel_h_A}, "
            f"h_fallback={params.voxel_h_fallback_A}, "
            f"max_voxels={params.voxel_count_fallback_threshold}); "
            "occ = occupancy_at(distance_field(spec, center_coords), rho)"),
    }
    payload["monotonicity_violations"] = solution.sweep.monotonicity_violations
    payload["n_universe_residues"] = len(universe)
    payload["voxel_h_sensitivity"] = h_sensitivity or {"evaluated": False}
    return write_json(stage_dir / "footprint_domain.json", payload)


def write_objective_matrix(stage_dir: Path, solution: FootprintSolution,
                           params: FootprintParams, redundancy: list[dict]) -> Path:
    result = solution.selection
    payload = {
        "objectives": list(params.objective_names),
        "directions": {name: "maximize" for name in params.objective_names},
        "descriptive_only": ["convexity", "expansion_rate", "n_components",
                             "n_bridges", "merge_events",
                             "per_component_residue_membership"],
        "w_k": params.w_k,
        "weights_source": "config: footprint_selection.w_k (FROZEN, F15)",
        "normalization": "min_max over the admissible set (F15)",
        "utopia_point": result.utopia,
        "degenerate_objectives": result.degenerate_objectives,
        "raw": {f"{k:g}": solution_objectives(solution, k)
                for k in [row.rho for row in solution.sweep.rows]},
        "normalized": {f"{k:g}": v for k, v in result.normalized.items()},
        "correlation_matrix": result.correlation,
        "objective_redundancy_abs_rho": params.objective_redundancy_abs_rho,
        "redundant_objective_pairs": redundancy,
    }
    return write_json(stage_dir / "footprint_objective_matrix.json", payload)


def solution_objectives(solution: FootprintSolution, rho: float) -> dict:
    from .selection import objective_vector
    return objective_vector(solution.sweep.by_rho()[rho])


def write_pareto(stage_dir: Path, solution: FootprintSolution) -> tuple[Path, Path]:
    result = solution.selection
    members = write_json(stage_dir / "footprint_pareto.json", {
        "n_admissible": result.n_admissible,
        "n_pareto_members": len(result.pareto_members),
        "members": [
            {"r_fp_A": k,
             "objectives": solution_objectives(solution, k),
             "normalized": result.normalized.get(k, {}),
             "distance_to_ideal": result.distances.get(k)}
            for k in result.pareto_members
        ],
    })
    dominated = write_json(stage_dir / "footprint_pareto_dominated.json", {
        "n_dominated": len(result.dominators),
        "note": ("Dominated and inadmissible candidates are retained with raw and "
                 "normalized values, distance, QC status and rejection reason "
                 "(agent §9.9)."),
        "dominated": [
            {"r_fp_A": k,
             "dominated_by": v,
             "objectives": solution_objectives(solution, k),
             "normalized": result.normalized.get(k, {}),
             "distance_to_ideal": result.distances.get(k)}
            for k, v in sorted(result.dominators.items())
        ],
        "inadmissible": [
            {"r_fp_A": k, "qc_failure_reason": reason,
             "objectives": solution_objectives(solution, k)}
            for k, reason in sorted(result.rejected.items())
        ],
    })
    return members, dominated


def write_decision(stage_dir: Path, solution: FootprintSolution,
                   params: FootprintParams, escalations: list[dict],
                   diagnostics=None) -> Path:
    result = solution.selection
    selected = result.selected_key
    payload = {
        "selected_r_fp_A": selected,
        "selection_procedure": (
            "admissible set -> Pareto -> min-max normalization over the admissible "
            "set -> utopia (1,1,1,1) -> equally weighted L2 distance -> argmin over "
            "Pareto members -> tie chain"),
        "objectives": list(params.objective_names),
        "w_k": params.w_k,
        "distance_metric": params.distance_metric,
        "utopia_point": result.utopia,
        "normalized_objective_vector": result.normalized.get(selected, {}),
        "distance_to_ideal": result.distances.get(selected),
        "tie_chain_applied": result.tie_chain_applied,
        "near_tie_flag": result.near_tie,
        "near_tie_detail": result.near_tie_detail,
        "n_admissible": result.n_admissible,
        "n_pareto_members": len(result.pareto_members),
        "excessive_coverage_fraction": solution.excessive_coverage_fraction,
        "coverage_constraint_failure": solution.coverage_constraint_failure,
        "r_hot_used_in_selection": False,
        "r_fp_bounded_by_r_hot": False,
        "f10_statement": ("r_fp was selected from geometry alone. r_hot is not an "
                          "argument to the footprint API and r_fp is neither set "
                          "equal to nor bounded below by it."),
        # --- v2 §4/§7 — objective/constraint separation, mirrored from Stage B --
        "degenerate_objectives": result.degenerate_objectives,
        "degenerate_normalized_value": 1.0,
        "effective_objectives_in_distance":
            diagnostics.effective_objectives if diagnostics else list(params.objective_names),
        "degenerate_objectives_excluded_from_distance":
            diagnostics.degenerate_objectives if diagnostics else [],
        "degenerate_exclusion_is_numerically_neutral":
            diagnostics.degenerate_exclusion_is_numerically_neutral if diagnostics else True,
        "degenerate_objective_note": (
            "An r_fp objective identical at every admissible candidate "
            "discriminates nothing. It is excluded from distance-to-utopia and "
            "reported as non-informative (v2 §4/§7). Under the frozen min-max rule "
            "it normalizes to 1.0 everywhere, so exclusion changes no distance — "
            "verified, not assumed."),
        "VACUOUS_PARETO_SELECTION": diagnostics.vacuous_pareto if diagnostics else False,
        "vacuous_pareto_reason": diagnostics.vacuous_reason if diagnostics else "NA",
        "ADMISSIBILITY_DOMINATES_SELECTION":
            diagnostics.admissibility_dominates if diagnostics else False,
        "inadmissible_fraction": diagnostics.inadmissible_fraction if diagnostics else 0.0,
        "inadmissible_by_rule": diagnostics.inadmissible_by_rule if diagnostics else {},
        "significance_may_determine_admissibility": False,
        "admissibility_note": (
            "Footprint admissibility (QC_F1 coverage, QC_F2 bridging, QC_F3 "
            "undefined geometry) is purely geometric. It was never conditioned on "
            "significant-center count in the first place, since r_fp operates on "
            "the already-frozen SIGNIFICANT_HOTSPOT_CENTERS set and never re-tests "
            "significance (v2 §7)."),
        "escalations": escalations,
    }
    return write_json(stage_dir / "footprint_decision.json", payload)


def write_boundary_diagnostic(stage_dir: Path, solution: FootprintSolution) -> Path:
    return write_json(stage_dir / "footprint_domain_boundary_diagnostic.json",
                      solution.boundary)


# --- 08_FINAL_FOOTPRINT -----------------------------------------------------

def write_final_residues(stage_dir: Path, solution: FootprintSolution,
                         universe: Universe, center_ids: tuple[int, ...],
                         plp: frozenset, blb: frozenset) -> Path:
    center_set = set(int(c) for c in center_ids)
    rows = _membership_rows(solution.selected, universe, center_set)
    for row in rows:
        rid = row["residue_index"]
        row["cohort_class"] = "PLP" if rid in plp else ("BLB" if rid in blb else "NA")
    return write_tsv(stage_dir / "footprint_residues.tsv", rows,
                     FOOTPRINT_RESIDUE_COLUMNS)


def write_final_geometry(stage_dir: Path, solution: FootprintSolution) -> Path:
    state = solution.selected
    row = {
        "r_fp_A": state.rho, "volume_A3": state.volume,
        "surface_area_A2": state.surface_area, "n_components": state.n_components,
        "n_covered_residues": state.n_covered_residues,
        "coverage_fraction": state.coverage, "compactness": state.compactness,
        "convexity": state.convexity, "connectivity": state.connectivity,
        "stability": state.stability, "parsimony": state.parsimony,
        "expansion_rate": state.expansion_rate, "n_bridges": state.n_bridges,
        "bridging_index": state.bridging_index,
        "n_components_eroded": state.n_components_eroded,
        "n_merge_events": state.n_merge_events,
        "abrupt_transition": state.abrupt_transition,
        "voxel_h_A": solution.spec.h,
        "voxel_resolution_fallback_applied": solution.spec.fallback_applied,
    }
    return write_tsv(stage_dir / "footprint_geometry.tsv", [row],
                     FOOTPRINT_GEOMETRY_COLUMNS)


def write_final_components(stage_dir: Path, solution: FootprintSolution,
                           universe: Universe,
                           center_ids: tuple[int, ...]) -> Path:
    from .descriptors import components

    state = solution.selected
    occ = solution.occupancy()
    field = components(occ)
    rows = []
    for label in range(1, field.n_components + 1):
        member_residues = np.nonzero(state.residue_labels == label)[0]
        member_centers = [int(center_ids[i])
                          for i, lab in enumerate(state.center_labels)
                          if int(lab) == label]
        idx = np.argwhere(field.labels == label)
        centroid = (np.asarray(solution.spec.origin)
                    + (idx + 0.5) * solution.spec.h).mean(axis=0)
        rows.append({
            "component_id": label,
            "n_voxels": int(field.sizes[label]),
            "volume_A3": float(field.sizes[label] * solution.spec.voxel_volume),
            "n_residues": int(len(member_residues)),
            "n_centers": len(member_centers),
            "center_residues": member_centers,
            "centroid_x": float(centroid[0]),
            "centroid_y": float(centroid[1]),
            "centroid_z": float(centroid[2]),
        })
    return write_tsv(stage_dir / "footprint_components.tsv", rows,
                     FOOTPRINT_COMPONENT_COLUMNS)


def write_occupancy(stage_dir: Path, solution: FootprintSolution) -> Path:
    """Voxel occupancy at the FINAL radius only (P5)."""
    occ = solution.occupancy()
    path = stage_dir / "footprint_occupancy.npz"
    np.savez_compressed(
        path,
        occupancy_packed=np.packbits(occ.ravel()),
        shape=np.asarray(occ.shape, dtype=np.int64),
        origin=np.asarray(solution.spec.origin, dtype=np.float64),
        voxel_h_A=np.asarray([solution.spec.h], dtype=np.float64),
        r_fp_A=np.asarray([solution.r_fp], dtype=np.float64),
    )
    return path


def write_surface_ply(stage_dir: Path, solution: FootprintSolution) -> tuple[Path, str]:
    """Marching-cubes mesh of the final footprint, re-extracted from the same EDT."""
    area, mesh, reason = surface_area(solution.sweep.dist, solution.r_fp,
                                      solution.spec.h, occupied=True)
    path = stage_dir / "footprint_surface.ply"
    if mesh is None:
        write_text(path, "ply\nformat ascii 1.0\ncomment MESH_UNAVAILABLE "
                         f"{reason}\nelement vertex 0\nelement face 0\nend_header")
        return path, reason
    verts = np.asarray(mesh["verts"]) + np.asarray(solution.spec.origin)
    faces = np.asarray(mesh["faces"])
    lines = [
        "ply", "format ascii 1.0",
        f"comment hotspot3d final footprint r_fp={solution.r_fp:g} A",
        f"element vertex {len(verts)}",
        "property float x", "property float y", "property float z",
        f"element face {len(faces)}",
        "property list uchar int vertex_index", "end_header",
    ]
    lines += [f"{v[0]:.4f} {v[1]:.4f} {v[2]:.4f}" for v in verts]
    lines += [f"3 {f[0]} {f[1]} {f[2]}" for f in faces]
    write_text(path, "\n".join(lines))
    return path, "NA"


def write_structures(struct_dir: Path, solution: FootprintSolution,
                     universe: Universe, center_ids: tuple[int, ...],
                     residue_names: dict[int, str], gene: str) -> dict[str, Path]:
    """CA-only footprint model (F1 fixes the residue representation as CA)."""
    struct_dir.mkdir(parents=True, exist_ok=True)
    state = solution.selected
    rows = np.nonzero(state.covered)[0]
    center_set = set(int(c) for c in center_ids)

    pdb_lines = [
        "REMARK   1 hotspot3d final footprint (CA-only, F1 residue representation)",
        f"REMARK   1 GENE {gene}  R_FP {solution.r_fp:.2f} A  "
        f"N_RESIDUES {len(rows)}  N_COMPONENTS {state.n_components}",
        "REMARK   1 B-FACTOR = footprint component id; OCCUPANCY = 1.00 center, "
        "0.50 covered",
    ]
    cif_rows = []
    for serial, row in enumerate(rows, start=1):
        rid = int(universe.ids[row])
        x, y, z = universe.coords[row]
        comp = int(state.residue_labels[row])
        occ = 1.00 if rid in center_set else 0.50
        name = residue_names.get(rid, "UNK")
        pdb_lines.append(
            f"ATOM  {serial:5d}  CA  {name:>3s} A{rid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{float(comp):6.2f}           C"
        )
        cif_rows.append(
            f"ATOM {serial} C CA . {name} A 1 {rid} ? "
            f"{x:.3f} {y:.3f} {z:.3f} {occ:.2f} {float(comp):.2f} 1"
        )
    pdb_lines.append("END")
    pdb = write_text(struct_dir / "final_footprint.pdb", "\n".join(pdb_lines))

    cif_lines = [
        f"data_hotspot3d_footprint_{gene}",
        "#",
        f"_hotspot3d.r_fp_A {solution.r_fp:.2f}",
        f"_hotspot3d.n_components {state.n_components}",
        "#",
        "loop_",
        "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol",
        "_atom_site.label_atom_id", "_atom_site.label_alt_id",
        "_atom_site.label_comp_id", "_atom_site.label_asym_id",
        "_atom_site.label_entity_id", "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code", "_atom_site.Cartn_x", "_atom_site.Cartn_y",
        "_atom_site.Cartn_z", "_atom_site.occupancy", "_atom_site.B_iso_or_equiv",
        "_atom_site.pdbx_PDB_model_num",
    ] + cif_rows + ["#"]
    cif = write_text(struct_dir / "final_footprint.cif", "\n".join(cif_lines))

    residues = sorted(int(universe.ids[r]) for r in rows)
    centers_sorted = sorted(center_set)
    cxc = write_text(struct_dir / "view_footprint.cxc", "\n".join([
        f"# hotspot3d footprint view — {gene}, r_fp = {solution.r_fp:g} A",
        "color lightgray",
        f"select :{','.join(str(r) for r in residues)}",
        "color sel cornflowerblue",
        "show sel atoms",
        f"select :{','.join(str(c) for c in centers_sorted)}",
        "color sel firebrick",
        "style sel sphere",
        "surface :" + ",".join(str(r) for r in residues),
        "~select",
    ]))
    pml = write_text(struct_dir / "view_footprint.pml", "\n".join([
        f"# hotspot3d footprint view — {gene}, r_fp = {solution.r_fp:g} A",
        "hide everything",
        "show cartoon",
        "color grey80",
        f"select footprint, resi {'+'.join(str(r) for r in residues)}",
        "color skyblue, footprint",
        "show surface, footprint",
        f"select centers, resi {'+'.join(str(c) for c in centers_sorted)}",
        "color firebrick, centers",
        "show spheres, centers",
        "deselect",
    ]))
    return {"pdb": pdb, "cif": cif, "cxc": cxc, "pml": pml}
