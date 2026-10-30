"""Upstream input contract and BLOCKING validation — agent §5, §6.

Stage E reads a fixed, declared set of upstream artifacts. :data:`UPSTREAM_PATHS` is
that contract in one place so the Lead can reconcile it against what Stages A–D
actually emit; a rename upstream shows up here as a blocking failure rather than as
a silently missing annotation.

Every failure in this module is BLOCKING and nothing is ever defaulted (agent §6).
The one exception that is *not* a failure is an upstream **scientific negative**:
if discovery ended with no significant hotspot, Stage E does not invent regions —
it reports the negative outcome and stops.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..orchestration.contracts import Handoff, QC_FAIL
from ..utils.errors import BlockedError, NegativeResult
from ..utils.hashing import sha256_file
from ..utils.io import read_json, read_tsv
from ..utils.runctx import RunContext
from .schema import AGENT
from .verbatim import VerbatimTable, load_verbatim_table

#: The declared read contract: logical name -> (stage directory, relative path).
UPSTREAM_PATHS: dict[str, tuple[str, str]] = {
    # Stage B — the statistics Stage E copies verbatim
    "hotspot_regions": ("06_FINAL_HOTSPOTS", "hotspot_regions.tsv"),
    "significant_hotspot_centers": (
        "06_FINAL_HOTSPOTS", "significant_hotspot_centers.tsv"),
    # Internal key keeps the METHOD_SPEC object name; the FILE name is the
    # canonical one fixed by the Output Contract (IX.3/IX.4) and written by Stage B.
    "hotspot_sphere_classified_variants": (
        "06_FINAL_HOTSPOTS", "hotspot_classified_variants.tsv"),
    "hotspot_covered_residues": (
        "06_FINAL_HOTSPOTS", "hotspot_covered_residues.tsv"),
    "all_residue_center_tests": (
        "06_FINAL_HOTSPOTS", "all_residue_center_tests.tsv"),
    # Stage A — structure and cohort
    "plddt_profile": ("03_STRUCTURE_QC", "plddt_profile.tsv"),
    "structure_with_plddt": ("03_STRUCTURE_QC", "structures/structure_with_plddt.cif"),
    "classified_cohort": ("03_STRUCTURE_QC", "classified_cohort.tsv"),
    "variants_residue_level": ("02_CLINVAR", "variants_residue_level.tsv"),
    "review_star_distribution": ("02_CLINVAR", "review_star_distribution.tsv"),
    # Stage C/D — footprint and robustness
    "footprint_residues": ("08_FINAL_FOOTPRINT", "footprint_residues.tsv"),
    "robustness_profile": ("09_ROBUSTNESS", "robustness_profile.json"),
    "center_sensitivity": ("09_ROBUSTNESS", "center_sensitivity.tsv"),
    # Stage B sensitivity — descriptive evidence-quality context only
    "sensitivity_overlap": ("11_SENSITIVITY", "sensitivity_overlap.tsv"),
}


@dataclass
class StageEInputs:
    """Everything Stage E is allowed to see, already validated."""

    gene: str
    uniprot_acc: str
    hotspot_ids: list[str]
    regions: VerbatimTable
    centers: dict[str, list[dict]]
    sphere_variants: dict[str, list[dict]]
    covered_residues: dict[str, list[dict]]
    plddt: dict[int, float]
    coords: dict[int, tuple[float, float, float]]
    structure_path: Path
    robustness_profile: dict
    center_sensitivity: list[dict]
    star_distribution: list[dict]
    sensitivity_overlap: list[dict]
    footprint_residues: list[dict]
    flags: dict[str, Any]
    upstream_hashes: dict[str, str] = field(default_factory=dict)
    paths: dict[str, Path] = field(default_factory=dict)
    #: SHA-256 of every upstream input as read, so §9.11 can re-verify at stage end.
    input_hashes: dict[str, str] = field(default_factory=dict)

    @property
    def all_annotated_residues(self) -> list[int]:
        out: set[int] = set()
        for rows in self.covered_residues.values():
            out.update(int(r["residue_index"]) for r in rows)
        return sorted(out)


# --- validation -------------------------------------------------------------

def validate_handoffs(
    ctx: RunContext, handoff_02: Handoff, handoff_03: Handoff, handoff_04: Handoff
) -> dict[str, str]:
    """Recompute every upstream hash and check QC. Mismatch BLOCKS (agent §6.1–6.2).

    ``Handoff.verify`` recomputes the manifest against what is on disk *now*, so a
    file edited after the handoff was written fails here — trust is never transitive.
    """
    seen: dict[str, str] = {}
    for handoff in (handoff_02, handoff_03, handoff_04):
        if handoff is None:
            raise BlockedError(
                "BLOCKED — Stage E requires handoff_02, handoff_03 and handoff_04. "
                "A missing handoff is never substituted."
            )
        if handoff.qc_status == QC_FAIL:
            raise BlockedError(
                f"BLOCKED — upstream {handoff.name} reports qc_status=FAIL. A failing "
                f"upstream stage never licenses annotation."
            )
        if handoff.run_id != ctx.run_id:
            raise BlockedError(
                f"BLOCKED — {handoff.name} belongs to run {handoff.run_id!r} but this "
                f"run is {ctx.run_id!r}. Handoffs are invalidated only by a new RUN_ID."
            )
        handoff.verify(ctx)
        seen[handoff.name] = _manifest_digest(handoff)

    if not handoff_04.payload.get("robustness_profile"):
        raise BlockedError(
            "BLOCKED — handoff_04 carries no robustness_profile. A missing profile is "
            "a declared blocking failure (agent §11); annotation may not proceed "
            "without the robustness context that must travel with every claim."
        )
    return seen


def _manifest_digest(handoff: Handoff) -> str:
    from ..utils.hashing import sha256_text
    items = ";".join(f"{k}={v}" for k, v in sorted(handoff.manifest.items()))
    return sha256_text(items)


def check_upstream_negative(handoff_02: Handoff, regions: VerbatimTable | None) -> None:
    """Raise :class:`NegativeResult` when discovery legitimately found nothing.

    A valid negative is not a failure and is never repaired by annotating a
    "best non-significant" region (agent §15).
    """
    if handoff_02.negative_result:
        raise NegativeResult(
            condition=handoff_02.negative_result.get("condition", "UPSTREAM_NEGATIVE"),
            detail=handoff_02.negative_result.get(
                "detail", "Stage B reported a scientific negative."
            ),
            upstream="handoff_02",
        )

    n_centers = handoff_02.payload.get("n_significant_centers")
    if n_centers is not None and int(n_centers) == 0:
        raise NegativeResult(
            condition="NO_SIGNIFICANT_HOTSPOTS",
            detail=(
                "Stage B reported n_significant_centers=0. No region is annotated; "
                "Stage E does not invent regions when discovery ends negatively."
            ),
            upstream="handoff_02",
        )

    if regions is not None and not regions.rows:
        raise NegativeResult(
            condition="NO_SIGNIFICANT_HOTSPOTS",
            detail=(
                "hotspot_regions.tsv contains no region. The negative outcome is "
                "reported as-is."
            ),
            upstream="06_FINAL_HOTSPOTS/hotspot_regions.tsv",
        )


def resolve_path(ctx: RunContext, name: str) -> Path:
    """Resolve a declared upstream artifact, enforcing the read barrier."""
    if name not in UPSTREAM_PATHS:
        raise BlockedError(f"BLOCKED — {name!r} is not a declared Stage E input.")
    stage, rel = UPSTREAM_PATHS[name]
    ctx.assert_may_read(AGENT, stage)
    return ctx.full_results / stage / rel


def require(ctx: RunContext, name: str) -> Path:
    path = resolve_path(ctx, name)
    if not path.is_file():
        stage, rel = UPSTREAM_PATHS[name]
        raise BlockedError(
            f"BLOCKED — required Stage E input {name!r} not found at {stage}/{rel}. "
            f"Missing inputs are blocking and are never defaulted (agent §6)."
        )
    return path


def load_inputs(
    ctx: RunContext, handoff_02: Handoff, handoff_03: Handoff, handoff_04: Handoff,
    *, uniprot_acc: str = "NA",
) -> StageEInputs:
    """Load and validate the complete Stage E input set."""
    upstream_hashes = validate_handoffs(ctx, handoff_02, handoff_03, handoff_04)

    regions_path = require(ctx, "hotspot_regions")
    regions = load_verbatim_table(regions_path, "hotspot_id")
    check_upstream_negative(handoff_02, regions)

    paths = {name: require(ctx, name) for name in UPSTREAM_PATHS}

    centers = _group_by_hotspot(read_tsv(paths["significant_hotspot_centers"]),
                                RESIDUE_KEY_BY_OBJECT["centers"])
    sphere_variants = _group_by_hotspot(
        read_tsv(paths["hotspot_sphere_classified_variants"]),
        RESIDUE_KEY_BY_OBJECT["sphere_variants"])
    covered = _group_by_hotspot(read_tsv(paths["hotspot_covered_residues"]),
                                RESIDUE_KEY_BY_OBJECT["covered"])

    plddt = {
        int(r["residue_index"]): float(r["plddt"])
        for r in read_tsv(paths["plddt_profile"])
        if r.get("plddt") is not None
    }

    coords = _load_coords(paths["classified_cohort"])

    profile = read_json(paths["robustness_profile"])

    flags = _collect_flags(handoff_02, handoff_03)

    return StageEInputs(
        gene=ctx.gene,
        uniprot_acc=uniprot_acc,
        hotspot_ids=regions.keys(),
        regions=regions,
        centers=centers,
        sphere_variants=sphere_variants,
        covered_residues=covered,
        plddt=plddt,
        coords=coords,
        structure_path=paths["structure_with_plddt"],
        robustness_profile=profile,
        center_sensitivity=read_tsv(paths["center_sensitivity"]),
        star_distribution=read_tsv(paths["review_star_distribution"]),
        sensitivity_overlap=read_tsv(paths["sensitivity_overlap"]),
        footprint_residues=read_tsv(paths["footprint_residues"]),
        flags=flags,
        upstream_hashes=upstream_hashes,
        paths=paths,
        input_hashes={name: sha256_file(path) for name, path in paths.items()},
    )


def _load_coords(path: Path) -> dict[int, tuple[float, float, float]]:
    """Cα coordinates for the classified cohort.

    Read only to position variants for the gated post hoc overlay. Stage E computes no
    geometry of its own: it never re-derives a footprint, a distance threshold or a
    hotspot boundary from these points.
    """
    coords: dict[int, tuple[float, float, float]] = {}
    for row in read_tsv(path):
        if any(row.get(k) is None for k in ("x_ca", "y_ca", "z_ca")):
            continue
        coords[int(row["residue_index"])] = (
            float(row["x_ca"]), float(row["y_ca"]), float(row["z_ca"])
        )
    return coords


#: The three residue objects deliberately carry DISTINCT key columns (Output
#: Contract IX.4) so that an accidental join across the three concepts fails
#: loudly instead of silently producing a wrong answer. Stage E must therefore
#: name the object it is reading; it may never assume a shared `residue_index`.
RESIDUE_KEY_BY_OBJECT = {
    "centers": "center_residue_index",
    "sphere_variants": "variant_residue_index",
    "covered": "covered_residue_index",
}


def _group_by_hotspot(rows: list[dict], residue_key: str | None = None
                      ) -> dict[str, list[dict]]:
    """Group rows by hotspot_id, exposing the object's key under `residue_index`.

    The canonical distinct key is retained on every row; `residue_index` is added
    as a read-only alias so downstream formatting code has one spelling to use.
    The alias is derived HERE, at the single point where the object's identity is
    known — never guessed at a call site, which is what the distinct key names
    exist to prevent.
    """
    out: dict[str, list[dict]] = {}
    for row in rows:
        if residue_key is not None:
            if residue_key not in row:
                raise BlockedError(
                    f"BLOCKED — expected key column {residue_key!r} is absent "
                    f"(present: {sorted(row)[:8]}). The three residue objects carry "
                    f"distinct key columns by design; Stage E does not guess."
                )
            row = {**row, "residue_index": row[residue_key]}
        # significant_hotspot_centers.tsv carries a single "hotspot_id" (a
        # center belongs to exactly one hotspot by construction);
        # hotspot_classified_variants.tsv and hotspot_covered_residues.tsv
        # carry "hotspot_ids" (plural, ";"-separated) because a variant or a
        # covered residue CAN legitimately fall inside more than one
        # overlapping hotspot's sphere. Both are handled uniformly here: a
        # singular value with no ";" splits into a list of one.
        raw = row.get("hotspot_ids", row.get("hotspot_id"))
        for hid in str(raw).split(";"):
            if hid:
                out.setdefault(hid, []).append(row)
    return out


def _collect_flags(handoff_02: Handoff, handoff_03: Handoff) -> dict[str, Any]:
    """Carry the four upstream flags forward. A missing flag is BLOCKING.

    Flags are never dropped and never defaulted: if Stage B did not state one,
    Stage E cannot invent a value, because "absent" and "FALSE" would then be
    indistinguishable to every downstream reader (agent §8).
    """
    required = {
        "global_clustering_flag": handoff_02,
        "permutation_resolution_limited": handoff_02,
        "fallback_radius_domain": handoff_02,
        "domain_boundary_warning": handoff_02,
    }
    flags: dict[str, Any] = {}
    missing: list[str] = []
    for name, handoff in required.items():
        if name in handoff.payload:
            flags[name] = handoff.payload[name]
        elif name in handoff_03.payload:
            flags[name] = handoff_03.payload[name]
        else:
            missing.append(name)
    if missing:
        raise BlockedError(
            f"BLOCKED — upstream handoffs do not state flag(s) {missing}. Flags are "
            f"never dropped and never defaulted; Stage E cannot distinguish an absent "
            f"flag from a FALSE one."
        )
    return flags
