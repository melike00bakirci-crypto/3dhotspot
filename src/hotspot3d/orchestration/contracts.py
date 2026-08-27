"""Stage API contracts and handoff schemas — the interlock between the four agents.

Lead-owned. Each agent implements the ``run_stage_*`` entry point declared here with
EXACTLY this signature and returns EXACTLY the declared handoff payload. Nothing
else crosses a stage boundary: a consumer reads only the producer's canonical files
plus the handoff JSON, and recomputes every hash it depends on (Part IV).

Handoffs are append-only. Only the Lead invalidates one, by declaring a new RUN_ID.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Protocol

from ..utils.errors import BlockedError
from ..utils.hashing import verify_manifest
from ..utils.io import read_json
from ..utils.runctx import RunContext

QC_PASS = "PASS"
QC_PASS_WITH_WARNINGS = "PASS_WITH_WARNINGS"
QC_FAIL = "FAIL"

# Columns Stage B may load from Stage A in its PRIMARY path (handoff_01 allowlist).
STAGE_B_PERMITTED_COLUMNS = [
    "residue_index", "class", "x_ca", "y_ca", "z_ca", "plddt", "ca_usable",
]
# Columns Stage B's primary path must never see. The >=1*/>=2* sensitivity
# analyses read review metadata through a separate, explicitly non-redefining
# channel declared in handoff_01.
FORBIDDEN_DOWNSTREAM_COLUMNS = [
    "clinvar_ids", "condition_text", "review_status", "star_levels", "max_star",
    "submitter", "n_records_plp", "n_records_blb",
]


@dataclass
class Handoff:
    """Base handoff: identity, QC verdict, manifest. Every consumer verifies."""
    name: str
    run_id: str
    config_sha256: str
    qc_status: str
    manifest: dict[str, str]
    payload: dict[str, Any] = field(default_factory=dict)
    negative_result: dict | None = None

    def as_dict(self) -> dict:
        out = {
            "handoff": self.name, "run_id": self.run_id,
            "config_sha256": self.config_sha256, "qc_status": self.qc_status,
            **self.payload, "manifest_sha256": self.manifest,
        }
        if self.negative_result is not None:
            out["negative_result"] = self.negative_result
        return out

    # -- consumer-side verification -----------------------------------------
    def verify(self, ctx: RunContext, *, require_qc: bool = True) -> None:
        if require_qc and self.qc_status == QC_FAIL:
            raise BlockedError(
                f"BLOCKED — upstream {self.name} reports qc_status=FAIL. "
                f"A failing upstream stage never licenses downstream work."
            )
        if self.config_sha256 != ctx.config.sha256:
            raise BlockedError(
                f"BLOCKED — {self.name} was produced under config_sha256="
                f"{self.config_sha256[:12]} but this run uses {ctx.config.sha256[:12]}."
            )
        bad = verify_manifest(self.manifest, ctx.run_root)
        if bad:
            raise BlockedError(
                f"BLOCKED — {self.name} manifest verification failed "
                f"(trust is never transitive): {bad[:5]}"
            )


def load_handoff(ctx: RunContext, name: str) -> Handoff:
    raw = read_json(ctx.handoff_path(name))
    known = {"schema_version", "handoff", "run_id", "config_sha256", "qc_status",
             "manifest_sha256", "negative_result"}
    return Handoff(
        name=raw.get("handoff", name), run_id=raw["run_id"],
        config_sha256=raw["config_sha256"], qc_status=raw["qc_status"],
        manifest=raw.get("manifest_sha256", {}),
        payload={k: v for k, v in raw.items() if k not in known},
        negative_result=raw.get("negative_result"),
    )


# --- required handoff payload keys, asserted by the Lead integration audit ---
HANDOFF_REQUIRED_KEYS: dict[str, list[str]] = {
    "handoff_01": [
        "gene", "uniprot_acc", "mane_transcript", "clinvar_release",
        "alphafold_model_version", "M", "N", "N_P", "N_B", "n_conflict",
        "stage_b_permitted_columns", "forbidden_downstream",
    ],
    "handoff_02": [
        "hotspot_radius", "search_domain_source", "fallback_radius_domain", "q",
        "fdr_method", "B", "n_significant_centers", "n_centers_without_variant",
        "bh_boundary_p", "permutation_resolution_limited", "b_recommended",
        "domain_boundary_warning", "global_clustering_flag", "near_tie_flag",
        "loo_sparse_proportion", "loo_zero_neighbour_proportion",
        "sensitivity_summary", "code_version",
    ],
    "handoff_03": [
        "hotspot_radius", "footprint_radius", "n_components", "volume_A3",
        "surface_area_A2", "n_footprint_residues", "coverage", "voxel_h",
        "near_tie_flag", "domain_boundary_warning", "code_version",
    ],
    "handoff_04": [
        "footprint_code_version", "robustness_profile", "n_S", "K_MAX",
        "total_subset_space", "total_evaluated", "per_level",
        "robustness_r_fp_invariant",
        "clinical_dataset_unmodified", "geometric_footprint_robustness_only",
    ],
    "handoff_05": [
        "annotation_sources", "n_hotspots_annotated", "mechanism_counts",
        "posthoc_gate_passed", "posthoc_gate_numbers", "robustness_profile",
        "global_clustering_flag", "interpretation_caveats",
    ],
}


def assert_handoff_shape(handoff: Handoff) -> None:
    missing = [k for k in HANDOFF_REQUIRED_KEYS.get(handoff.name, [])
               if k not in handoff.payload]
    if missing:
        raise BlockedError(
            f"BLOCKED — {handoff.name} is missing required payload keys: {missing}"
        )


# --- stage entry-point protocols -------------------------------------------

class StageA(Protocol):
    def __call__(self, ctx: RunContext, *, gene: str,
                 source: Any | None = None) -> Handoff: ...


class StageB(Protocol):
    def __call__(self, ctx: RunContext, *, upstream: Handoff) -> Handoff: ...


class StageCD(Protocol):
    def __call__(self, ctx: RunContext, *,
                 upstream: Handoff) -> tuple[Handoff, Handoff]: ...


class StageE(Protocol):
    def __call__(self, ctx: RunContext, *, handoff_02: Handoff, handoff_03: Handoff,
                 handoff_04: Handoff, source: Any | None = None) -> Handoff: ...
