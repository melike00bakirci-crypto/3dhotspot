"""Synthetic upstream fixtures for Stage C/D tests.

No network, no real gene, no ClinVar payload. The geometry is a deterministic
cubic lattice so that MST merge scales, coverage fractions and component counts
are analytically predictable and a test can assert on them rather than on
whatever the code happened to produce.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from hotspot3d.orchestration.contracts import Handoff
from hotspot3d.utils.config import load_config
from hotspot3d.utils.hashing import manifest_for
from hotspot3d.utils.io import write_json, write_text, write_tsv
from hotspot3d.utils.runctx import RunContext

CONFIG_PATH = Path(__file__).resolve().parents[2] / "config" / "pipeline.yaml"
FIXED_RUN_ID = "20250101T000000Z_testcfg0_testcode"

COORD_COLUMNS = ["residue_index", "x_ca", "y_ca", "z_ca", "plddt", "ca_usable",
                 "residue_aa3"]
UNIVERSE_COLUMNS = ["residue_index", "in_u_struct"]
COHORT_COLUMNS = ["residue_index", "class", "n_records"]
CENTER_COLUMNS = ["center_residue_index", "region_id", "carries_clinvar_variant",
                  "p_emp", "q_value"]
COVERED_COLUMNS = ["residue_index", "center_residue_index"]
HS_VARIANT_COLUMNS = ["residue_index", "center_residue_index", "class"]


def lattice(n_side: int = 6, spacing: float = 5.5) -> tuple[list[int], np.ndarray]:
    """A cubic lattice of CA positions; residue index is the lattice raster order."""
    pts = []
    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                pts.append((i * spacing, j * spacing, k * spacing))
    coords = np.asarray(pts, dtype=np.float64)
    ids = list(range(1, len(coords) + 1))
    return ids, coords


def two_cluster_centers(n_side: int = 6) -> list[int]:
    """Six centers in two triplets that merge *inside* the derived domain.

    Intra-cluster MST edges are one lattice step; the single inter-cluster edge is
    long enough that ``1.25 * rho_all`` binds ``rho_max`` and the merge happens
    partway up the grid, so merge events, neck measurement and the
    natural/artificial classification are all exercised.
    """
    def idx(i, j, k):
        return i * n_side * n_side + j * n_side + k + 1

    return sorted([
        idx(1, 1, 1), idx(1, 1, 2), idx(1, 2, 1),          # cluster A
        idx(3, 3, 3), idx(3, 3, 2), idx(3, 2, 3),          # cluster B
    ])


def make_context(tmp_path: Path, gene: str = "SYNTH",
                 run_id: str = FIXED_RUN_ID) -> RunContext:
    cfg = load_config(CONFIG_PATH)
    return RunContext.create(gene=gene, config=cfg, results_root=tmp_path / "results",
                             synthetic=True, run_id=run_id)


def write_upstream(ctx: RunContext, *, center_ids: list[int],
                   ids: list[int], coords: np.ndarray,
                   r_hot: float = 8.0,
                   n_plp: int = 40, n_blb: int = 40) -> Handoff:
    """Materialize the canonical Stage A / Stage B artifacts this stage consumes."""
    dir_02 = ctx.stage_dir("02_CLINVAR")
    dir_03 = ctx.stage_dir("03_STRUCTURE_QC")
    dir_06 = ctx.stage_dir("06_FINAL_HOTSPOTS")

    write_tsv(dir_03 / "residue_coordinates.tsv", [
        {"residue_index": rid, "x_ca": float(c[0]), "y_ca": float(c[1]),
         "z_ca": float(c[2]), "plddt": 92.0, "ca_usable": True, "residue_aa3": "ALA"}
        for rid, c in zip(ids, coords)
    ], COORD_COLUMNS)
    write_tsv(dir_03 / "positional_universe.tsv",
              [{"residue_index": rid, "in_u_struct": True} for rid in ids],
              UNIVERSE_COLUMNS)
    (dir_03 / "structures").mkdir(exist_ok=True)
    write_text(dir_03 / "structures" / "structure_with_plddt.cif",
               "data_synthetic\n#\n_entry.id SYNTHETIC\n#\n")

    # cohort L: pathogenic residues near the centers, benign residues far away
    center_coords = coords[[ids.index(c) for c in center_ids]]
    d = np.linalg.norm(coords[:, None, :] - center_coords[None, :, :], axis=2).min(axis=1)
    order = np.argsort(d, kind="stable")
    plp = [ids[i] for i in order[:n_plp]]
    blb = [ids[i] for i in order[-n_blb:]]
    write_tsv(dir_02 / "classified_cohort.tsv",
              [{"residue_index": r, "class": "PLP", "n_records": 3} for r in sorted(plp)]
              + [{"residue_index": r, "class": "BLB", "n_records": 2}
                 for r in sorted(blb)],
              COHORT_COLUMNS)

    write_tsv(dir_06 / "significant_hotspot_centers.tsv", [
        {"center_residue_index": c,
         "region_id": f"R{1 if i < len(center_ids) / 2 else 2}",
         "carries_clinvar_variant": c in set(plp) | set(blb),
         "p_emp": 0.0001, "q_value": 0.001}
        for i, c in enumerate(center_ids)
    ], CENTER_COLUMNS)
    covered = []
    for c in center_ids:
        cc = coords[ids.index(c)]
        within = np.nonzero(np.linalg.norm(coords - cc, axis=1) <= r_hot)[0]
        covered += [{"residue_index": ids[i], "center_residue_index": c}
                    for i in within]
    write_tsv(dir_06 / "hotspot_covered_residues.tsv", covered, COVERED_COLUMNS)
    write_tsv(dir_06 / "hotspot_classified_variants.tsv", [
        {"residue_index": r["residue_index"],
         "center_residue_index": r["center_residue_index"],
         "class": "PLP" if r["residue_index"] in set(plp) else "BLB"}
        for r in covered
        if r["residue_index"] in set(plp) | set(blb)
    ], HS_VARIANT_COLUMNS)

    payload = {
        "hotspot_radius": r_hot, "search_domain_source": "pcf",
        "fallback_radius_domain": False, "q": 0.05, "fdr_method": "BH", "B": 10000,
        "n_significant_centers": len(center_ids), "n_centers_without_variant": 0,
        "bh_boundary_p": 0.01, "permutation_resolution_limited": False,
        "b_recommended": 10000, "domain_boundary_warning": False,
        "global_clustering_flag": True, "near_tie_flag": False,
        "loo_sparse_proportion": 0.1, "loo_zero_neighbour_proportion": 0.0,
        "sensitivity_summary": {}, "code_version": "synthetic",
    }
    files = [p for p in (dir_02, dir_03, dir_06)
             for p in p.rglob("*") if p.is_file()]
    handoff = Handoff(name="handoff_02", run_id=ctx.run_id,
                      config_sha256=ctx.config.sha256, qc_status="PASS",
                      manifest=manifest_for(files, ctx.run_root), payload=payload)
    write_json(ctx.handoff_path("handoff_02"), handoff.as_dict())
    return handoff


def synthetic_run(tmp_path: Path, *, n_side: int = 6, spacing: float = 5.5,
                  center_ids: list[int] | None = None, r_hot: float = 8.0,
                  gene: str = "SYNTH", run_id: str = FIXED_RUN_ID,
                  n_plp: int = 40, n_blb: int = 40
                  ) -> tuple[RunContext, Handoff]:
    ctx = make_context(tmp_path, gene=gene, run_id=run_id)
    ids, coords = lattice(n_side, spacing)
    centers = center_ids if center_ids is not None else two_cluster_centers(n_side)
    handoff = write_upstream(ctx, center_ids=centers, ids=ids, coords=coords,
                             r_hot=r_hot, n_plp=n_plp, n_blb=n_blb)
    return ctx, handoff
