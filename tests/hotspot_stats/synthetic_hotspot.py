"""Synthetic Stage A inputs for the Stage B test suite.

`tests/fixtures/synthetic_clinvar.py` (owned by `data-structure`) does not exist yet
and no signature for it is documented anywhere in the repository, so this module
builds the Stage A surface Stage B actually consumes — `residue_coordinates.tsv`,
`positional_universe.tsv`, `classified_cohort.tsv` and `handoff_01.json` with a real
SHA-256 manifest. When the shared fixture lands, :func:`write_stage_a` is the single
call site to redirect.

Everything here is deterministic: a fixed RNG seed, a lattice-based "protein" at a
realistic residue density, and a planted 3D cluster of pathogenic residues so that
the pipeline has something true to find.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml

from hotspot3d.orchestration.contracts import Handoff
from hotspot3d.utils.config import load_config
from hotspot3d.utils.hashing import manifest_for
from hotspot3d.utils.io import write_json, write_tsv
from hotspot3d.utils.runctx import RunContext

COORD_COLUMNS = ["residue_index", "aa", "x_ca", "y_ca", "z_ca", "plddt", "plddt_band",
                 "ca_usable"]
UNIVERSE_COLUMNS = ["residue_index", "in_universe"]
COHORT_COLUMNS = ["residue_index", "class"]
REVIEW_COLUMNS = ["residue_index", "max_star", "n_records_plp", "n_records_blb"]


def synthetic_protein(sphere_radius: float = 18.0, spacing: float = 5.0,
                      jitter: float = 0.6, seed: int = 7) -> np.ndarray:
    """A compact lattice "protein" at a realistic CA density (~1 residue / 125 A^3)."""
    rng = np.random.default_rng(seed)
    steps = int(np.ceil(sphere_radius / spacing))
    grid = np.arange(-steps, steps + 1) * spacing
    points = np.array([[x, y, z] for x in grid for y in grid for z in grid])
    inside = np.linalg.norm(points, axis=1) <= sphere_radius
    coords = points[inside] + rng.normal(0.0, jitter, size=(int(inside.sum()), 3))
    order = np.lexsort((coords[:, 2], coords[:, 1], coords[:, 0]))
    return coords[order]


def uniform_cohort(coords: np.ndarray, n_plp: int = 24, n_blb: int = 40,
                   seed: int = 11) -> tuple[np.ndarray, np.ndarray]:
    """A cohort drawn uniformly at random over the protein — NO planted structure.

    This is the honest negative control for the structure-aware positional null: the
    P/LP set is generated from exactly the distribution the null assumes, so a
    rejection would be a false positive rather than a missed signal. A cohort merely
    pushed to the outer shell would NOT do — non-uniform is not the same as
    unclustered, and the positional null would rightly flag it.
    """
    rng = np.random.default_rng(seed)
    chosen = rng.choice(len(coords), size=n_plp + n_blb, replace=False)
    positions = np.array(sorted(int(p) for p in chosen), dtype=np.int64)
    plp = set(int(p) for p in chosen[:n_plp])
    labels = np.array([1 if int(p) in plp else 0 for p in positions], dtype=np.int64)
    return positions, labels


def planted_cohort(coords: np.ndarray, n_cluster_plp: int = 16, n_scatter_plp: int = 8,
                   n_blb: int = 40, cluster_seed: int = 11
                   ) -> tuple[np.ndarray, np.ndarray]:
    """Positions and labels: one tight P/LP cluster plus scattered P/LP and B/LB."""
    rng = np.random.default_rng(cluster_seed)
    centre = coords[int(np.argmin(np.linalg.norm(coords - coords.mean(axis=0), axis=1)))]
    distance = np.linalg.norm(coords - centre, axis=1)
    cluster = list(np.argsort(distance)[:n_cluster_plp])

    far = [i for i in np.argsort(distance)[::-1] if i not in cluster]
    scatter_plp = list(rng.choice(far[:len(far) // 2], size=n_scatter_plp, replace=False))
    remaining = [i for i in far if i not in scatter_plp]
    blb = list(rng.choice(remaining, size=n_blb, replace=False))

    positions = np.array(sorted(cluster + scatter_plp + blb), dtype=np.int64)
    plp = set(int(i) for i in cluster + scatter_plp)
    labels = np.array([1 if int(p) in plp else 0 for p in positions], dtype=np.int64)
    return positions, labels


def write_stage_a(ctx: RunContext, coords: np.ndarray, positions: np.ndarray,
                  labels: np.ndarray, low_confidence: np.ndarray | None = None,
                  stars: dict[int, int] | None = None) -> dict[str, Path]:
    """Write the Stage A tables Stage B reads, in the documented Stage A schema."""
    qc = ctx.stage_dir("03_STRUCTURE_QC")
    clinvar = ctx.stage_dir("02_CLINVAR")
    indices = np.arange(1, len(coords) + 1, dtype=np.int64)
    plddt = np.full(len(coords), 92.0)
    if low_confidence is not None:
        plddt[low_confidence] = 45.0

    coord_rows = [{
        "residue_index": int(indices[i]), "aa": "ALA",
        "x_ca": float(coords[i, 0]), "y_ca": float(coords[i, 1]),
        "z_ca": float(coords[i, 2]), "plddt": float(plddt[i]),
        "plddt_band": ">90" if plddt[i] > 90 else "<50", "ca_usable": True,
    } for i in range(len(coords))]
    universe_rows = [{"residue_index": int(i), "in_universe": True} for i in indices]
    cohort_rows = [{"residue_index": int(indices[p]),
                    "class": "PLP" if labels[k] == 1 else "BLB"}
                   for k, p in enumerate(positions)]

    paths = {
        "residue_coordinates": write_tsv(qc / "residue_coordinates.tsv", coord_rows,
                                         COORD_COLUMNS),
        "positional_universe": write_tsv(qc / "positional_universe.tsv", universe_rows,
                                         UNIVERSE_COLUMNS),
        "classified_cohort": write_tsv(qc / "classified_cohort.tsv", cohort_rows,
                                       COHORT_COLUMNS),
    }
    if stars is not None:
        review_rows = [{"residue_index": int(indices[p]),
                        "max_star": int(stars.get(int(indices[p]), 1)),
                        "n_records_plp": 1 if labels[k] == 1 else 0,
                        "n_records_blb": 0 if labels[k] == 1 else 1}
                       for k, p in enumerate(positions)]
        paths["variants_residue_level"] = write_tsv(
            clinvar / "variants_residue_level.tsv", review_rows, REVIEW_COLUMNS)
    return paths


def make_context(tmp_path: Path, *, B: int = 200, gene: str = "SYN",
                 run_id: str = "20250101T000000Z_test_test",
                 overrides: dict | None = None) -> RunContext:
    """A synthetic RunContext with a reduced ``B`` supplied through a config overlay.

    ``B`` is injected the way a real deviation would be — through configuration, not
    through a function argument — and the stage refuses a non-default ``B`` unless the
    run is flagged synthetic.
    """
    overlay = {"permutation": {"B_default": B}}
    if overrides:
        overlay.update(overrides)
    tmp_path.mkdir(parents=True, exist_ok=True)
    overlay_path = tmp_path / "overlay.yaml"
    overlay_path.write_text(yaml.safe_dump(overlay), encoding="utf-8")
    cfg = load_config(_pipeline_config(), gene_overlay=overlay_path)
    return RunContext.create(gene=gene, config=cfg, results_root=tmp_path / "results",
                             synthetic=True, run_id=run_id)


def _pipeline_config() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / "config" / "pipeline.yaml"
        if candidate.is_file():
            return candidate
    raise FileNotFoundError("config/pipeline.yaml not found above the test tree")


def make_handoff_01(ctx: RunContext, paths: dict[str, Path], *, qc_status: str = "PASS",
                    n_p: int = 0, n_b: int = 0, review_channel: bool = True) -> Handoff:
    """A handoff_01 with a real SHA-256 manifest, as the consumer contract requires."""
    payload = {
        "gene": ctx.gene, "uniprot_acc": "P00000", "mane_transcript": "NM_000000.1",
        "clinvar_release": "2025-01-01", "alphafold_model_version": "v4",
        "M": 0, "N": n_p + n_b, "N_P": n_p, "N_B": n_b, "n_conflict": 0,
        "stage_b_permitted_columns": ["residue_index", "class", "x_ca", "y_ca", "z_ca",
                                      "plddt", "ca_usable"],
        "forbidden_downstream": ["clinvar_ids", "condition_text", "review_status",
                                 "star_levels", "max_star", "submitter",
                                 "n_records_plp", "n_records_blb"],
    }
    if review_channel and "variants_residue_level" in paths:
        # The shape Stage A actually emits (see data/stage.py handoff_01 payload).
        payload["sensitivity_channel"] = {
            "declared": True,
            "file": str(Path(paths["variants_residue_level"]).relative_to(ctx.run_root)),
            "columns": ["residue_index", "max_star"],
            "may_redefine_primary": False,
        }
    handoff = Handoff(
        name="handoff_01", run_id=ctx.run_id, config_sha256=ctx.config.sha256,
        qc_status=qc_status,
        manifest=manifest_for(list(paths.values()), root=ctx.run_root),
        payload=payload)
    write_json(ctx.handoff_path("handoff_01"), handoff.as_dict())
    return handoff


def full_synthetic_run(tmp_path: Path, *, B: int = 200, sphere_radius: float = 18.0,
                       spacing: float = 5.0, n_cluster_plp: int = 16,
                       n_scatter_plp: int = 8, n_blb: int = 40,
                       low_confidence_fraction: float = 0.0,
                       stars: bool = True, run_id: str = "20250101T000000Z_test_test",
                       qc_status: str = "PASS",
                       dispersed: bool = False, cohort_seed: int = 11,
                       overrides: dict | None = None) -> tuple[RunContext, Handoff, dict]:
    """One call producing (ctx, handoff_01, geometry) ready for ``run_stage_b``.

    ``dispersed=True`` replaces the planted cluster with a uniform draw over the
    protein — the negative control for the positional null.
    """
    coords = synthetic_protein(sphere_radius=sphere_radius, spacing=spacing)
    if dispersed:
        positions, labels = uniform_cohort(
            coords, n_plp=n_cluster_plp + n_scatter_plp, n_blb=n_blb, seed=cohort_seed)
    else:
        positions, labels = planted_cohort(coords, n_cluster_plp, n_scatter_plp, n_blb)
    ctx = make_context(tmp_path, B=B, run_id=run_id, overrides=overrides)

    low = None
    if low_confidence_fraction > 0:
        rng = np.random.default_rng(3)
        low = rng.choice(len(coords), size=int(low_confidence_fraction * len(coords)),
                         replace=False)
    star_map = None
    if stars:
        star_map = {int(positions[k] + 1): (2 if k % 3 == 0 else 1)
                    for k in range(len(positions))}

    paths = write_stage_a(ctx, coords, positions, labels, low, star_map)
    handoff = make_handoff_01(ctx, paths, qc_status=qc_status,
                              n_p=int(labels.sum()), n_b=int((labels == 0).sum()))
    return ctx, handoff, {"coords": coords, "positions": positions, "labels": labels,
                          "paths": paths}
