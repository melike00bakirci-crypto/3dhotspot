"""Frozen footprint parameters, read once from ``config/pipeline.yaml``.

Bundling the constants here keeps the geometric core (:mod:`hotspot3d.footprint.api`)
a *pure function of its arguments*: Phase D re-executes that core tens of thousands
of times and must get bit-identical behaviour without touching the filesystem.

Every value is read with :meth:`FrozenConfig.get`, which raises rather than
defaulting — a missing parameter is BLOCKING (agent §6.4).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ..utils.config import FrozenConfig
from ..utils.errors import BlockedError

#: Config keys the stage validates up front so it fails fast, not mid-sweep.
REQUIRED_CONFIG_KEYS: list[str] = [
    "footprint_domain.rho_min_A",
    "footprint_domain.rho_max_rule",
    "footprint_domain.rho_all_source",
    "footprint_domain.post_merge_factor",
    "footprint_domain.dmax_fraction_cap",
    "footprint_domain.hard_ceiling_A",
    "footprint_domain.step_fp_A",
    "footprint_domain.bound_below_by_r_hot",
    "footprint_geometry.voxel_h_A",
    "footprint_geometry.voxel_h_fallback_A",
    "footprint_geometry.voxel_count_fallback_threshold",
    "footprint_geometry.bbox_padding_A",
    "footprint_geometry.connectivity",
    "footprint_geometry.surface_method",
    "footprint_geometry.erosion_fraction",
    "footprint_geometry.abrupt_dncomp_threshold",
    "footprint_geometry.abrupt_er_median_multiple",
    "footprint_qc.QC_F1_max_coverage",
    "footprint_qc.coverage_denominator",
    "footprint_qc.QC_F2_max_bridging_index",
    "footprint_qc.all_inadmissible_escalation_fraction",
    "footprint_qc.relax_on_failure",
    "footprint_qc.admissibility_dominates_fraction",
    "footprint_selection.objectives",
    "footprint_selection.descriptive_only",
    "footprint_selection.w_k",
    "footprint_selection.normalization",
    "footprint_selection.normalization_domain",
    "footprint_selection.distance_metric",
    "footprint_selection.tie_chain",
    "footprint_selection.near_tie_relative_threshold",
    "footprint_selection.objective_redundancy_abs_rho",
    "footprint_selection.exclude_degenerate_objectives",
    "boundary_diagnostic.band_fraction",
    "boundary_diagnostic.bw2_pareto_fraction",
    "boundary_diagnostic.bw3_top_n",
    "boundary_diagnostic.auto_widen_domain",
]


@dataclass(frozen=True)
class FootprintParams:
    """Immutable snapshot of every constant the footprint core consumes."""

    # II.9 domain
    rho_min_A: float
    post_merge_factor: float
    dmax_fraction_cap: float
    hard_ceiling_A: float
    step_fp_A: float
    # II.8 geometry
    voxel_h_A: float
    voxel_h_fallback_A: float
    voxel_count_fallback_threshold: int
    bbox_padding_A: float
    connectivity: int
    erosion_fraction: float
    abrupt_dncomp_threshold: int
    abrupt_er_median_multiple: float
    # II.10 QC
    qc_f1_max_coverage: float
    qc_f2_max_bridging_index: float
    all_inadmissible_escalation_fraction: float
    admissibility_dominates_fraction: float
    # II.10 selection
    objective_names: tuple[str, ...]
    tie_chain: tuple[str, ...]
    w_k: float
    distance_metric: str
    near_tie_relative_threshold: float
    near_tie_radius_separation_steps: int
    objective_redundancy_abs_rho: float
    exclude_degenerate_objectives: bool
    # boundary diagnostic
    band_fraction: float
    bw2_pareto_fraction: float
    bw3_top_n: int

    @classmethod
    def from_config(cls, cfg: FrozenConfig) -> "FootprintParams":
        cfg.require_all(REQUIRED_CONFIG_KEYS)

        # F10 is structural, not stylistic: r_fp must never be bounded by r_hot.
        if bool(cfg.get("footprint_domain.bound_below_by_r_hot")):
            raise BlockedError(
                "FROZEN methodology violation (F10): "
                "footprint_domain.bound_below_by_r_hot is TRUE. r_fp is an "
                "independent analytical quantity and is neither set equal to nor "
                "bounded below by r_hot."
            )
        if bool(cfg.get("footprint_qc.relax_on_failure")):
            raise BlockedError(
                "FROZEN methodology violation (II.10): footprint_qc.relax_on_failure "
                "is TRUE. The 0.50 coverage rule is never relaxed to rescue a run."
            )
        if bool(cfg.get("boundary_diagnostic.auto_widen_domain")):
            raise BlockedError(
                "FROZEN methodology violation (II.10): automatic domain widening is "
                "prohibited; a boundary optimum is escalated, never self-corrected."
            )
        if cfg.get("footprint_qc.coverage_denominator") != "U_struct":
            raise BlockedError(
                "FROZEN methodology violation (II.10): the QC-F1 coverage denominator "
                "is always |U_struct|."
            )
        if cfg.get("footprint_geometry.surface_method") != "marching_cubes":
            raise BlockedError(
                "FROZEN methodology violation (A18): surface_method must be "
                "marching_cubes."
            )
        if int(cfg.get("footprint_geometry.connectivity")) != 6:
            raise BlockedError(
                "FROZEN methodology violation (A18): voxel connectivity is 6 "
                "(face adjacency)."
            )
        if cfg.get("footprint_selection.normalization") != "min_max":
            raise BlockedError("FROZEN (F15): normalization must be min_max.")
        if cfg.get("footprint_selection.normalization_domain") != "admissible_set":
            raise BlockedError("FROZEN (F15): normalization domain is the admissible set.")

        objectives = _objective_names(cfg.get("footprint_selection.objectives"))
        if objectives != ("compactness", "connectivity", "stability", "parsimony"):
            raise BlockedError(
                "FROZEN methodology violation (II.10): the four r_fp objectives are "
                f"compactness, connectivity, stability, parsimony — config declares "
                f"{objectives}."
            )
        descriptive = tuple(cfg.get("footprint_selection.descriptive_only"))
        overlap = sorted(set(objectives) & set(descriptive))
        if overlap:
            raise BlockedError(
                f"config declares {overlap} as both an objective and descriptive-only."
            )

        return cls(
            rho_min_A=float(cfg.get("footprint_domain.rho_min_A")),
            post_merge_factor=float(cfg.get("footprint_domain.post_merge_factor")),
            dmax_fraction_cap=float(cfg.get("footprint_domain.dmax_fraction_cap")),
            hard_ceiling_A=float(cfg.get("footprint_domain.hard_ceiling_A")),
            step_fp_A=float(cfg.get("footprint_domain.step_fp_A")),
            voxel_h_A=float(cfg.get("footprint_geometry.voxel_h_A")),
            voxel_h_fallback_A=float(cfg.get("footprint_geometry.voxel_h_fallback_A")),
            voxel_count_fallback_threshold=int(
                cfg.get("footprint_geometry.voxel_count_fallback_threshold")),
            bbox_padding_A=float(cfg.get("footprint_geometry.bbox_padding_A")),
            connectivity=int(cfg.get("footprint_geometry.connectivity")),
            erosion_fraction=float(cfg.get("footprint_geometry.erosion_fraction")),
            abrupt_dncomp_threshold=int(
                cfg.get("footprint_geometry.abrupt_dncomp_threshold")),
            abrupt_er_median_multiple=float(
                cfg.get("footprint_geometry.abrupt_er_median_multiple")),
            qc_f1_max_coverage=float(cfg.get("footprint_qc.QC_F1_max_coverage")),
            qc_f2_max_bridging_index=float(cfg.get("footprint_qc.QC_F2_max_bridging_index")),
            all_inadmissible_escalation_fraction=float(
                cfg.get("footprint_qc.all_inadmissible_escalation_fraction")),
            admissibility_dominates_fraction=float(
                cfg.get("footprint_qc.admissibility_dominates_fraction")),
            objective_names=objectives,
            tie_chain=tuple(cfg.get("footprint_selection.tie_chain")),
            w_k=float(cfg.get("footprint_selection.w_k")),
            distance_metric=str(cfg.get("footprint_selection.distance_metric")),
            near_tie_relative_threshold=float(
                cfg.get("footprint_selection.near_tie_relative_threshold")),
            near_tie_radius_separation_steps=int(
                cfg.get_optional("footprint_selection.near_tie_radius_separation_steps", 2)),
            objective_redundancy_abs_rho=float(
                cfg.get("footprint_selection.objective_redundancy_abs_rho")),
            exclude_degenerate_objectives=bool(
                cfg.get("footprint_selection.exclude_degenerate_objectives")),
            band_fraction=float(cfg.get("boundary_diagnostic.band_fraction")),
            bw2_pareto_fraction=float(cfg.get("boundary_diagnostic.bw2_pareto_fraction")),
            bw3_top_n=int(cfg.get("boundary_diagnostic.bw3_top_n")),
        )

    def as_dict(self) -> dict:
        """Provenance-ready parameter record."""
        return {
            "rho_min_A": self.rho_min_A,
            "post_merge_factor": self.post_merge_factor,
            "dmax_fraction_cap": self.dmax_fraction_cap,
            "hard_ceiling_A": self.hard_ceiling_A,
            "step_fp_A": self.step_fp_A,
            "voxel_h_A": self.voxel_h_A,
            "voxel_h_fallback_A": self.voxel_h_fallback_A,
            "voxel_count_fallback_threshold": self.voxel_count_fallback_threshold,
            "bbox_padding_A": self.bbox_padding_A,
            "connectivity": self.connectivity,
            "erosion_fraction": self.erosion_fraction,
            "abrupt_dncomp_threshold": self.abrupt_dncomp_threshold,
            "abrupt_er_median_multiple": self.abrupt_er_median_multiple,
            "QC_F1_max_coverage": self.qc_f1_max_coverage,
            "QC_F2_max_bridging_index": self.qc_f2_max_bridging_index,
            "all_inadmissible_escalation_fraction":
                self.all_inadmissible_escalation_fraction,
            "admissibility_dominates_fraction": self.admissibility_dominates_fraction,
            "objectives": list(self.objective_names),
            "tie_chain": list(self.tie_chain),
            "w_k": self.w_k,
            "distance_metric": self.distance_metric,
            "near_tie_relative_threshold": self.near_tie_relative_threshold,
            "near_tie_radius_separation_steps": self.near_tie_radius_separation_steps,
            "exclude_degenerate_objectives": self.exclude_degenerate_objectives,
            "objective_redundancy_abs_rho": self.objective_redundancy_abs_rho,
            "band_fraction": self.band_fraction,
            "bw2_pareto_fraction": self.bw2_pareto_fraction,
            "bw3_top_n": self.bw3_top_n,
        }


def _objective_names(entries: Sequence[dict]) -> tuple[str, ...]:
    names = []
    for entry in entries:
        if entry.get("direction") != "maximize":
            raise BlockedError(
                f"FROZEN (II.10): every r_fp objective maximizes; {entry} does not."
            )
        names.append(str(entry["name"]))
    return tuple(names)
