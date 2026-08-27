"""Stage B entry point — primary spatial hotspot discovery.

    run_stage_b(ctx, upstream=handoff_01) -> handoff_02

Execution order is METHOD_SPEC's as revised by Workflow v2, without deviation:

    II.1  global K and PCF (structure-aware positional null, both cohorts)
    II.2  candidate radius domain (deterministic, with the four fallback triggers)
    §5.4  POWER CERTIFICATE, pre-flight — before any permutation compute is spent
    II.4  permutation-resolution PRE-FLIGHT
    II.5  radius scan (LOO-MCC, Zg, fold enrichment, stability) — each radius independent
    §4    admissibility — ONLY genuinely undefined quantities; significance never filters
    II.6  Pareto -> min-max normalization -> L2 distance-to-ideal -> r_hot
    II.7  final detection at r_hot under the PRIMARY positional null, BH at q = 0.05,
          the three residue objects, regions; the label null reported alongside
    §5.4  POWER CERTIFICATE, post-hoc — binding
    II.4  permutation-resolution POST-HOC
    II.7  the three post-primary, non-redefining sensitivity analyses

**Then it stops.** No footprint, no r_fp, no robustness, no biology. There is no
post-hoc hotspot LOO validation stage (F9/F13).

**Terminal states (v2 §0).** Exactly one of:

  COMPLETED           significant hotspots found;
  COMPLETED_NEGATIVE  nothing found AND the power certificate PASSED — the design
                      demonstrably could have found something;
  UNDERPOWERED        nothing found and the design could not have found anything.
                      This is NOT a result about the gene, and no claim about the
                      absence of hotspots may be made (v2 Appendix A);
  BLOCKED             a structural, data or QC precondition was not met.

A run may be reported as ``COMPLETED_NEGATIVE`` only after the §5.4 certificate has
been computed and has passed. Reporting an underpowered run as a negative is a
specification violation, and it is the defect that motivated v2.
"""
from __future__ import annotations

import logging
import os
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from ..orchestration.contracts import QC_PASS, QC_PASS_WITH_WARNINGS, Handoff
from ..spatial.domain import derive_domain
from ..spatial.ripley import build_grid, cohort_statistics, intensity, pcf_peaks
from ..utils.errors import BlockedError, NegativeResult
from ..utils.hashing import code_version, manifest_for, sha256_file
from ..utils.io import write_json, write_text, write_tsv
from ..utils.runctx import RunContext
from ..utils.seeds import ctx_final_detection, ctx_global_null, ctx_scan, ctx_sensitivity
from ..utils.status import (
    STAGE_STATUS_SCHEMA,
    OutcomeType,
    Severity,
    StageStatus,
    Status,
    WarningCollector,
)
from . import figures as figs
from . import resolution
from .detection import (
    ALL_CENTER_TEST_COLUMNS,
    BH_TABLE_COLUMNS,
    CLASSIFIED_VARIANT_COLUMNS,
    COVERED_RESIDUE_COLUMNS,
    HOTSPOT_REGION_COLUMNS,
    final_detection,
    significant_subset_matches,
)
from .inputs import (
    ReadGuard,
    assert_no_downstream_dirs_touched,
    assert_owned_dirs_empty,
    load_inputs,
    load_review_channel,
)
from .loo import SPARSE_COLUMNS
from .power import (
    POSITIONAL_NULL,
    UNINFORMATIVE_STATEMENT,
    PowerCertificate,
    UnderpoweredResult,
    certify,
    domain_occupancies,
    preflight as power_preflight_certificates,
)
from .report import render_report
from .scan import (
    ADMISSIBILITY_RULES,
    OBJECTIVE_NAMES,
    QC_RULES,
    R_HOT_SCAN_COLUMNS,
    SCAN_CENTER_COLUMNS,
    SIGNIFICANCE_CONDITIONED_DIAGNOSTICS,
    attach_neighbor_stability,
    scan_one_radius,
)
from .selection import decision_payload, run_selection
from .sensitivity import (
    NON_REDEFINING_TEXT,
    SENSITIVITY_OVERLAP_COLUMNS,
    SensitivityOutcome,
    conclusion_flipped,
    overlap_row,
    run_restricted_detection,
)
from .structures import CaRecord, write_ca_cif, write_ca_pdb, write_chimerax_script, write_pymol_script

AGENT = "hotspot-statistics"

# METHOD_SPEC II.5B, stated once and carried into every output that prints Zg.
# T(r) = sum_c n_P(c)/n_L(c) = sum_i y_i * w_i with w_i = sum_{c covering i} 1/n_L(c),
# so w_i is SMALLER exactly where labelled residues are dense. In the regime this
# pipeline runs in — a small classified cohort L inside a much larger U_struct — a
# genuine 3D cluster therefore drives Zg NEGATIVE. Confirmed independently by the Lead
# (clustered Zg = -1.95 vs dispersed Zg = +1.04 on matched geometry). The statistic is
# FROZEN and is implemented exactly as specified; what must never happen is a reader
# taking a negative Zg as evidence AGAINST a hotspot.
ZG_DIRECTION_NOTE = (
    "DIRECTION OF perm_evidence_zg (II.5B, FROZEN): T(r) = sum_c n_P(c)/n_L(c) is "
    "algebraically sum_i y_i * w_i with the purely geometric weight "
    "w_i = sum_{c covering i} 1/n_L(c), which is SMALLER where labelled residues are "
    "dense. With a small classified cohort inside a much larger U_struct — the regime "
    "this pipeline runs in — genuine 3D clustering of P/LP therefore pushes Zg "
    "NEGATIVE. A negative Zg here is the EXPECTED reading when clustering is present "
    "and is NOT evidence against a hotspot. Zg is one of four equally weighted "
    "selection objectives and is never a significance test; it is derived from the "
    "SECONDARY label-permutation null, while significance comes only from the "
    "per-center structure-aware POSITIONAL null (v2 §5.2) with BH at q."
)

#: v2 §5.5 — the scan and the final detection draw from independent seed contexts, so
#: their center sets can differ. How that difference may be described depends on the
#: §5.4 certificate, and the two readings must never be interchanged.
SCAN_FINAL_NOTE_CERTIFIED = (
    "The radius scan and the final detection draw independent permutations from "
    "different seed contexts (II.12 / v2 §5.5). The FINAL detection is authoritative "
    "and the provisional scan set is never promoted. With the §5.4 power certificate "
    "PASSED, a disagreement means the evidence sits near the BH boundary, where "
    "discrete permutation p-values can move between RUN_IDs."
)
SCAN_FINAL_NOTE_UNCERTIFIED = (
    "The radius scan and the final detection draw independent permutations from "
    "different seed contexts (II.12 / v2 §5.5). The FINAL detection is authoritative. "
    "The §5.4 power certificate FAILED, so this disagreement must NOT be described as "
    "the evidence sitting near the significance boundary: it is evidence that the "
    "test is operating at its floor and could not have rejected regardless of the "
    "data."
)

GLOBAL = "04_GLOBAL_CLUSTERING"
RADIUS = "05_HOTSPOT_RADIUS"
FINAL = "06_FINAL_HOTSPOTS"
SENS = "11_SENSITIVITY"
STAGES = (GLOBAL, RADIUS, FINAL, SENS)

MANIFEST_COLUMNS = ["relative_path", "kind", "status", "reason", "sha256", "size_bytes"]

REQUIRED_CONFIG_KEYS = [
    "global_clustering.grid_step_A", "global_clustering.pcf_shell_width_A",
    "global_clustering.r_max_rule", "global_clustering.null_model",
    "global_clustering.envelope_percentiles", "global_clustering.alpha",
    "global_clustering.global_test_statistic", "global_clustering.compute_for_cohorts",
    "global_clustering.non_significant_terminates_discovery",
    "radius_domain.R_FLOOR", "radius_domain.R_CEIL", "radius_domain.r_cap_rule",
    "radius_domain.step_hot_A", "radius_domain.N_MIN_PCF",
    "radius_domain.pcf_elevation_percentile", "radius_domain.margin_A",
    "radius_domain.min_run_length_A", "radius_domain.min_grid_points",
    "radius_domain.automatic_domain_expansion",
    "permutation.B_default", "permutation.primary_null", "permutation.local_null",
    "permutation.report_secondary_null", "permutation.p_value_estimator",
    "permutation.exclude_uninformative_centers",
    "power_certificate.enabled", "power_certificate.preflight",
    "power_certificate.posthoc", "power_certificate.emit_artifact",
    "fdr.method", "fdr.q", "fdr.secondary_report_method", "fdr.boundary_ties_all_rejected",
    "permutation_resolution_diagnostic.enabled",
    "permutation_resolution_diagnostic.preflight",
    "permutation_resolution_diagnostic.safety_factor_s",
    "permutation_resolution_diagnostic.k_target",
    "permutation_resolution_diagnostic.b_ladder",
    "permutation_resolution_diagnostic.max_recommendable_B",
    "permutation_resolution_diagnostic.auto_change_B",
    "loo_mcc.role", "loo_mcc.kappa", "loo_mcc.zero_neighbour_prediction",
    "loo_mcc.exact_tie_prediction", "loo_mcc.tie_tolerance",
    "loo_mcc.sparse_reporting_cutoff", "loo_mcc.drop_sparse_residues",
    "loo_mcc.undefined_mcc_value",
    "radius_qc.QC_H2_max_coverage", "radius_qc.QC_H3_min_median_n_labeled",
    "radius_qc.QC_H4_max_isolated_center_fraction",
    "radius_qc.significance_may_determine_admissibility",
    "radius_qc.admissibility_constraints",
    "radius_qc.significance_conditioned_diagnostics",
    "radius_qc.admissibility_dominates_fraction",
    "radius_selection.exclude_degenerate_objectives",
    "radius_selection.objectives", "radius_selection.weights_equal",
    "radius_selection.w_k", "radius_selection.normalization",
    "radius_selection.normalization_domain", "radius_selection.utopia_point",
    "radius_selection.distance_metric", "radius_selection.tie_chain",
    "radius_selection.near_tie_relative_threshold",
    "radius_selection.near_tie_radius_separation_steps",
    "radius_selection.objective_redundancy_abs_rho",
    "radius_selection.degenerate_normalized_value",
    "boundary_diagnostic.enabled", "boundary_diagnostic.band_fraction",
    "boundary_diagnostic.bw2_pareto_fraction", "boundary_diagnostic.bw3_top_n",
    "boundary_diagnostic.auto_widen_domain",
    "final_detection.center_universe", "final_detection.region_edge_rule",
    "final_detection.emit_three_residue_objects",
    "sensitivity_analyses.run_at_frozen_r_hot", "sensitivity_analyses.may_redefine_primary",
    "sensitivity_analyses.plddt70", "sensitivity_analyses.review_status_strata",
    "sensitivity_analyses.global_clustering_caveat",
    "plddt.primary_filtering_enabled", "plddt.sensitivity_threshold",
    "plddt.center_universe_min_plddt", "plddt.center_universe_excludes_low_confidence",
    "biological_plausibility.enabled",
    "biological_plausibility.structural_coverage_penalty_threshold",
    "biological_plausibility.structural_coverage_major_threshold",
    "biological_plausibility.cohort_absorption_penalty_threshold",
    "biological_plausibility.r_to_domain_ratio_penalty_threshold",
    "biological_plausibility.enters_pareto_objectives",
    "biological_plausibility.radius_of_gyration_source",
    "seeding.MASTER_SEED", "seeding.generator",
]

FROZEN_STAGE_B = [
    ("fdr.method", "BH"), ("fdr.q", 0.05), ("fdr.boundary_ties_all_rejected", True),
    ("loo_mcc.kappa", 2.0), ("loo_mcc.classification_threshold", None),
    ("loo_mcc.role", "radius_selection_only"), ("loo_mcc.drop_sparse_residues", False),
    ("loo_mcc.zero_neighbour_prediction", "BLB"), ("loo_mcc.exact_tie_prediction", "BLB"),
    ("radius_selection.w_k", 1.0), ("radius_selection.weights_equal", True),
    ("radius_selection.normalization", "min_max"),
    ("radius_selection.normalization_domain", "admissible_set"),
    ("radius_selection.distance_metric", "L2"),
    ("radius_domain.R_FLOOR", 5.0), ("radius_domain.R_CEIL", 25.0),
    ("radius_domain.step_hot_A", 0.5),
    ("radius_domain.automatic_domain_expansion", False),
    ("boundary_diagnostic.auto_widen_domain", False),
    ("permutation_resolution_diagnostic.auto_change_B", False),
    # v2 §5.2 — the PRIMARY per-center null is the same structure-aware positional
    # null as the global analysis; label permutation is the SECONDARY analysis.
    ("permutation.primary_null", "structure_aware_positional"),
    ("permutation.local_null", "label_permutation"),
    ("permutation.report_secondary_null", True),
    # v2 §5.4 — the certificate is never optional and never disabled by a run.
    ("power_certificate.enabled", True),
    ("power_certificate.preflight", True),
    ("power_certificate.posthoc", True),
    ("power_certificate.emit_artifact", True),
    # v2 §4 — significance is an objective, never an admissibility filter.
    ("radius_qc.significance_may_determine_admissibility", False),
    ("radius_selection.exclude_degenerate_objectives", True),
    ("permutation_resolution_diagnostic.k_target", 1),
    ("global_clustering.alpha", 0.05),
    ("global_clustering.null_model", "structure_aware_positional"),
    ("global_clustering.non_significant_terminates_discovery", False),
    ("sensitivity_analyses.may_redefine_primary", False),
    ("plddt.primary_filtering_enabled", False),
    ("final_detection.center_universe", "U_struct"),
    # v2 §2 — the candidate-CENTER universe pLDDT floor is fixed and always applied.
    ("plddt.center_universe_min_plddt", 50.0),
    ("plddt.center_universe_excludes_low_confidence", True),
    # v2 §4 — biological plausibility is computed and reported; F15 keeps exactly
    # four Pareto objectives (Lead ruling in config/pipeline.yaml), so this is
    # asserted rather than trusted.
    ("biological_plausibility.enabled", True),
    ("biological_plausibility.enters_pareto_objectives", False),
    ("biological_plausibility.radius_of_gyration_source",
     "modelled_structure_ca_coordinates"),
]

EXPECTED_OUTPUTS: dict[str, list[str]] = {
    GLOBAL: ["ripleys_k_plp.tsv", "ripleys_k_blb.tsv", "ripleys_k_summary.json",
             "pair_correlation_plp.tsv", "pair_correlation_blb.tsv",
             "pair_correlation_peaks.json", "candidate_radius_domain.json",
             "figures/F1_ripleys_k.png", "figures/F1_ripleys_k.svg",
             "figures/F2_pair_correlation.png", "figures/F2_pair_correlation.svg",
             "stage_status.json", f"warnings_{GLOBAL}.tsv", "stage_manifest.tsv"],
    RADIUS: ["r_hot_scan.tsv", "loo_diagnostics.tsv", "objective_matrix.json",
             "pareto_front.json", "pareto_dominated.json", "radius_decision.json",
             "domain_boundary_diagnostic.json", "power_certificate_preflight.tsv",
             "center_universe.json", "center_universe_exclusions.tsv",
             "figures/F3_radius_scan.png", "figures/F3_radius_scan.svg",
             "figures/F4_pareto_parallel_coordinates.png",
             "figures/F4_pareto_parallel_coordinates.svg",
             "figures/F5_distance_to_ideal.png", "figures/F5_distance_to_ideal.svg",
             "stage_status.json", f"warnings_{RADIUS}.tsv", "stage_manifest.tsv"],
    FINAL: ["all_residue_center_tests.tsv", "significant_hotspot_centers.tsv",
            "hotspot_classified_variants.tsv", "hotspot_covered_residues.tsv",
            "hotspot_regions.tsv", "bh_fdr_table.tsv",
            "permutation_resolution_diagnostic.json",
            "power_certificate.json", "secondary_null_label_permutation.json",
            "structures/final_hotspots.pdb", "structures/final_hotspots.cif",
            "structures/hotspot_centers_bfactor_qvalue.pdb",
            "structures/view_hotspots.cxc", "structures/view_hotspots.pml",
            "figures/F6_hotspot_map.png", "figures/F6_hotspot_map.svg",
            "handoff_02.json", "stage_b_report.md",
            "stage_status.json", f"warnings_{FINAL}.tsv", "stage_manifest.tsv"],
    SENS: ["sensitivity_plddt70.json", "sensitivity_review_status.json",
           "sensitivity_overlap.tsv", "SENSITIVITY_IS_NON_REDEFINING.txt",
           "stage_status.json", f"warnings_{SENS}.tsv", "stage_manifest.tsv"],
}


def _utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


@dataclass
class _StageBook:
    """Per-stage bookkeeping: warnings, timing, created files, status."""

    stage: str
    warnings: WarningCollector
    started: str = field(default_factory=_utc)
    t0: float = field(default_factory=time.perf_counter)
    created: list[Path] = field(default_factory=list)
    status: Status = Status.NOT_RUN
    outcome: OutcomeType = OutcomeType.NOT_APPLICABLE
    reason: str = ""
    recommended: str = ""
    negative: dict | None = None


class StageB:
    """The Stage B run. One instance per invocation; holds no global state."""

    def __init__(self, ctx: RunContext, upstream: Handoff):
        self.ctx = ctx
        self.upstream = upstream
        self.cfg = ctx.config
        self.guard = ReadGuard()
        self.books = {s: _StageBook(s, WarningCollector(s, AGENT)) for s in STAGES}
        self.escalations: list[dict] = []
        self.summary: dict = {}
        self.log = self._logger()
        self.t_start = time.perf_counter()
        self.steps: list[dict] = []

    # -- infrastructure ------------------------------------------------------
    def _logger(self) -> logging.Logger:
        root = Path(os.environ.get("HOTSPOT3D_LOG_ROOT", "logs")) / self.ctx.run_id
        root.mkdir(parents=True, exist_ok=True)
        logger = logging.getLogger(f"hotspot3d.stage_b.{self.ctx.run_id}")
        logger.setLevel(logging.INFO)
        logger.handlers.clear()
        logger.propagate = False
        handler = logging.FileHandler(root / "hotspot_statistics.log", encoding="utf-8")
        handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        logger.addHandler(handler)
        return logger

    def _dir(self, stage: str) -> Path:
        self.ctx.assert_may_write(AGENT, stage)
        return self.ctx.stage_dir(stage)

    def _record(self, stage: str, path: Path) -> Path:
        self.books[stage].created.append(Path(path))
        return Path(path)

    def _step(self, name: str, t0: float) -> None:
        self.steps.append({"step": name, "wall_seconds": round(time.perf_counter() - t0, 3),
                           "ended_utc": _utc()})
        self.log.info("step complete: %s (%.3fs)", name, time.perf_counter() - t0)

    def _escalate(self, ambiguity: str, options: list[str], consequences: list[str],
                  recommendation: str) -> None:
        self.escalations.append({
            "ambiguity": ambiguity, "options": options, "consequences": consequences,
            "recommendation": recommendation, "resolved_by_agent": False,
            "detected_utc": _utc(),
        })
        self.log.warning("ESCALATION: %s", ambiguity)

    # -- II.0 validation -----------------------------------------------------
    def validate(self) -> None:
        t0 = time.perf_counter()
        self.upstream.verify(self.ctx)                     # hashes + qc_status + config
        self.cfg.require_all(REQUIRED_CONFIG_KEYS)
        for key, expected in FROZEN_STAGE_B:
            self.cfg.assert_frozen_value(key, expected)
        assert_owned_dirs_empty(self.ctx, self.guard)

        # v2 §4 — the config DECLARES the rule split and the code APPLIES it. If the
        # two ever drift, a rule could be declared and not enforced (or enforced and
        # not declared), which is exactly the silence the declaration exists to
        # prevent. Asserted rather than assumed.
        declared_constraints = set(self.cfg.get("radius_qc.admissibility_constraints"))
        declared_diagnostics = set(
            self.cfg.get("radius_qc.significance_conditioned_diagnostics"))
        if declared_constraints != set(ADMISSIBILITY_RULES):
            raise BlockedError(
                f"BLOCKED — radius_qc.admissibility_constraints declares "
                f"{sorted(declared_constraints)} but Stage B applies "
                f"{sorted(ADMISSIBILITY_RULES)}. v2 §4 requires each admissibility "
                f"constraint to be declared in configuration with its threshold.")
        if declared_diagnostics != set(SIGNIFICANCE_CONDITIONED_DIAGNOSTICS):
            raise BlockedError(
                f"BLOCKED — radius_qc.significance_conditioned_diagnostics declares "
                f"{sorted(declared_diagnostics)} but Stage B records "
                f"{sorted(SIGNIFICANCE_CONDITIONED_DIAGNOSTICS)}.")
        if declared_constraints & declared_diagnostics:
            raise BlockedError(
                f"BLOCKED — {sorted(declared_constraints & declared_diagnostics)} is "
                f"declared BOTH as an admissibility constraint and as a "
                f"significance-conditioned diagnostic. v2 §4 forbids significance from "
                f"determining admissibility.")

        self.B = int(self.cfg.get("permutation.B_default"))
        if self.B != 100000 and not self.ctx.synthetic:
            raise BlockedError(
                f"BLOCKED — permutation.B_default is {self.B}, not the frozen default "
                f"100000 (DECISION-B-DEFAULT-0001), and this is not a synthetic run. B "
                f"is never changed to chase or avoid significance; a different B "
                f"requires a pre-registered revision under a new RUN_ID.")
        if self.B != 100000:
            self.books[FINAL].warnings.add(
                "SYNTHETIC_INPUT_MODE",
                f"B = {self.B} (frozen default is 100000); permitted only because this "
                f"run is flagged synthetic.",
                potential_consequence="Permutation resolution is coarser than a real run.",
                recommended_action="Never compare these p-values with a production run.",
                affected_output="all_residue_center_tests.tsv", status="ACCEPTED_BY_DESIGN")
        if self.ctx.synthetic:
            self.books[GLOBAL].warnings.add(
                "SYNTHETIC_INPUT_MODE", "Stage B ran on synthetic inputs.",
                status="ACCEPTED_BY_DESIGN")

        self.inputs = load_inputs(self.ctx, self.upstream, self.guard)
        self.log.info("column allowlist asserted: %s", self.inputs.allowlist_used)
        self.log.info("M=%d N=%d N_P=%d N_B=%d", self.inputs.M, self.inputs.N,
                      self.inputs.N_P, self.inputs.N_B)
        self._step("input_validation", t0)

    # -- II.1 ----------------------------------------------------------------
    def global_clustering(self) -> None:
        t0 = time.perf_counter()
        inp = self.inputs
        grid_step = float(self.cfg.get("global_clustering.grid_step_A"))
        shell = float(self.cfg.get("global_clustering.pcf_shell_width_A"))
        lo_p, hi_p = [float(x) for x in self.cfg.get("global_clustering.envelope_percentiles")]
        self.alpha = float(self.cfg.get("global_clustering.alpha"))
        # The configured alpha must remain the two-sided level implied by the frozen
        # envelope percentiles; a silent divergence between the plotted envelope and
        # the decision level would be a methodological drift, so it is asserted.
        derived_alpha = round((lo_p + (100.0 - hi_p)) / 100.0, 10)
        if abs(self.alpha - derived_alpha) > 1e-12:
            raise BlockedError(
                f"BLOCKED — global_clustering.alpha is {self.alpha} but the frozen "
                f"envelope percentiles {[lo_p, hi_p]} imply {derived_alpha}. The plotted "
                f"envelope and the global-test decision level must be the same level.")

        # v2 §4 — radius of gyration of the modelled structure (Cα coordinates of
        # the FULL U_struct, per config biological_plausibility.radius_of_gyration_
        # source), used only for the reported r_to_domain_ratio DIAGNOSTIC — never
        # a Pareto objective, never an admissibility constraint.
        centroid = inp.universe_coords.mean(axis=0)
        self.radius_of_gyration_A = float(
            np.sqrt(np.mean(np.sum((inp.universe_coords - centroid) ** 2, axis=1))))

        self.grid, self.d_max, self.r_max = build_grid(inp.universe_coords, grid_step)
        if len(self.grid) == 0:
            raise BlockedError(
                f"BLOCKED — r_max = floor(D_max/2) = {self.r_max} A leaves no radius on "
                f"the {grid_step} A estimator grid (D_max = {self.d_max:.4g} A).")
        lam, volume, volume_method = intensity(inp.universe_coords)

        self.stats = {}
        for cohort, mask in (("PLP", inp.y == 1), ("BLB", inp.y == 0)):
            context = ctx_global_null(cohort)
            self.stats[cohort] = cohort_statistics(
                cohort, inp.labeled_coords[mask], inp.universe_coords, self.grid, volume,
                shell, self.B, self.ctx.seeds.rng(context), context,
                envelope_percentiles=(lo_p, hi_p), alpha=self.alpha)

        self.peaks = pcf_peaks(self.grid, self.stats["PLP"].g_obs)
        plp = self.stats["PLP"]
        self.global_flag = "significant" if plp.p_global <= self.alpha else "non_significant"

        if self.global_flag == "non_significant":
            self.books[GLOBAL].warnings.add(
                "GLOBAL_CLUSTERING_NON_SIGNIFICANT",
                f"P/LP global test max|Z| = {plp.t_obs:.4g}, p = {plp.p_global:.4g} "
                f"(alpha = {self.alpha:g}) under the structure-aware positional null.",
                potential_consequence=(
                    "Local hotspots may still exist and local discovery CONTINUES (F5); "
                    "every downstream claim must carry this caveat."),
                recommended_action=("Report global_clustering = non_significant alongside "
                                    "every hotspot claim. Do not terminate discovery."),
                affected_output="ripleys_k_summary.json")

        self.volume_detail = {
            "lambda_struct": lam, "V_struct_A3": volume, "volume_method": volume_method,
            "lambda_per_cohort": {c: s.n_points / volume for c, s in self.stats.items()},
            # v2 §4 — radius of gyration of U_struct, the sole input to the reported
            # r_to_domain_ratio DIAGNOSTIC (never a Pareto objective, F15).
            "radius_of_gyration_A": self.radius_of_gyration_A,
            "radius_of_gyration_source": self.cfg.get(
                "biological_plausibility.radius_of_gyration_source"),
            "note": ("K, L and g are normalized by the COHORT's own intensity "
                     "n_C / V_struct, so L(r) - r = 0 and g(r) = 1 are the correct CSR "
                     "reference lines for each cohort. Observed and null cohorts share "
                     "n_C, so the choice cancels in every comparison."),
        }
        self._write_global_outputs()
        self._step("II.1_global_clustering", t0)

    def _curve_rows(self, s) -> tuple[list[dict], list[dict]]:
        k_rows, g_rows = [], []
        for i, r in enumerate(s.grid):
            k_rows.append({
                "radius_A": float(r), "n_points": s.n_points,
                "k_obs": float(s.k_obs[i]), "k_null_mean": float(s.k_null_mean[i]),
                "k_null_sd": float(s.k_null_sd[i]), "k_env_lo": float(s.k_env_lo[i]),
                "k_env_hi": float(s.k_env_hi[i]), "l_obs": float(s.l_obs[i]),
                "l_null_mean": float(s.l_null_mean[i]), "l_env_lo": float(s.l_env_lo[i]),
                "l_env_hi": float(s.l_env_hi[i]),
                "l_minus_r_obs": float(s.l_obs[i] - r),
                "l_minus_r_env_lo": float(s.l_env_lo[i] - r),
                "l_minus_r_env_hi": float(s.l_env_hi[i] - r),
                "z_score": float(s.k_z[i]),
                "outside_envelope": bool(s.k_obs[i] < s.k_env_lo[i] or
                                         s.k_obs[i] > s.k_env_hi[i]),
            })
            g_rows.append({
                "radius_A": float(r), "n_points": s.n_points,
                "g_obs": float(s.g_obs[i]), "g_null_mean": float(s.g_null_mean[i]),
                "g_null_sd": float(s.g_null_sd[i]), "g_env_lo": float(s.g_env_lo[i]),
                "g_env_hi": float(s.g_env_hi[i]), "g_null_p95": float(s.g_null_p95[i]),
                "z_score": float(s.g_z[i]),
                "elevated_above_1": bool(s.g_obs[i] > 1.0),
                "above_null_p95": bool(s.g_obs[i] > s.g_null_p95[i]),
                "in_set_E": (bool(s.g_obs[i] > 1.0 and s.g_obs[i] > s.g_null_p95[i])
                             if s.cohort == "PLP" else None),
            })
        return k_rows, g_rows

    def _write_global_outputs(self) -> None:
        d = self._dir(GLOBAL)
        k_cols = ["radius_A", "n_points", "k_obs", "k_null_mean", "k_null_sd",
                  "k_env_lo", "k_env_hi", "l_obs", "l_null_mean", "l_env_lo", "l_env_hi",
                  "l_minus_r_obs", "l_minus_r_env_lo", "l_minus_r_env_hi", "z_score",
                  "outside_envelope"]
        g_cols = ["radius_A", "n_points", "g_obs", "g_null_mean", "g_null_sd",
                  "g_env_lo", "g_env_hi", "g_null_p95", "z_score", "elevated_above_1",
                  "above_null_p95", "in_set_E"]
        for cohort, suffix in (("PLP", "plp"), ("BLB", "blb")):
            k_rows, g_rows = self._curve_rows(self.stats[cohort])
            self._record(GLOBAL, write_tsv(d / f"ripleys_k_{suffix}.tsv", k_rows, k_cols))
            self._record(GLOBAL, write_tsv(d / f"pair_correlation_{suffix}.tsv",
                                           g_rows, g_cols))

        self._record(GLOBAL, write_json(d / "ripleys_k_summary.json", {
            "run_id": self.ctx.run_id, "gene": self.ctx.gene,
            "D_max_A": self.d_max, "r_max_A": self.r_max,
            "grid_A": [float(r) for r in self.grid],
            "grid_step_A": float(self.cfg.get("global_clustering.grid_step_A")),
            "pcf_shell_width_A": float(self.cfg.get("global_clustering.pcf_shell_width_A")),
            "null_model": "structure_aware_positional",
            "null_draws_without_replacement": True,
            "B": self.B, "envelope_percentiles": self.cfg.get(
                "global_clustering.envelope_percentiles"),
            "alpha": self.alpha,
            "alpha_source": ("config:global_clustering.alpha, asserted equal to the "
                             "two-sided level implied by the frozen envelope percentiles"),
            "intensity": self.volume_detail,
            "global_test": {
                cohort: {
                    "statistic": "T = max_r |Z_K(r)|", "T_obs": s.t_obs,
                    "p_global": s.p_global, "significant": s.p_global <= self.alpha,
                    "n_points": s.n_points, "seed_context": s.seed_context,
                    "seed": self.ctx.seeds.seed(s.seed_context),
                } for cohort, s in self.stats.items()},
            "global_clustering_flag": self.global_flag,
            "non_significant_terminates_discovery": False,
            "note": ("Global assessment only. Ripley's K never selects r_hot, and a "
                     "non-significant global result does not terminate local discovery "
                     "(F5); it sets the caveat flag carried by every downstream claim."),
        }))
        self._record(GLOBAL, write_json(d / "pair_correlation_peaks.json", self.peaks))
        for p in figs.figure_f1_ripley(d, self.stats):
            self._record(GLOBAL, p)

    # -- II.2 ----------------------------------------------------------------
    def candidate_domain(self) -> None:
        t0 = time.perf_counter()
        self.domain = derive_domain(
            self.grid, self.stats["PLP"].g_obs, self.stats["PLP"].g_null_p95,
            n_plp=self.inputs.N_P, d_max=self.d_max,
            r_floor=float(self.cfg.get("radius_domain.R_FLOOR")),
            r_ceil=float(self.cfg.get("radius_domain.R_CEIL")),
            step=float(self.cfg.get("radius_domain.step_hot_A")),
            margin=float(self.cfg.get("radius_domain.margin_A")),
            n_min_pcf=int(self.cfg.get("radius_domain.N_MIN_PCF")),
            min_run_length=float(self.cfg.get("radius_domain.min_run_length_A")),
            min_grid_points=int(self.cfg.get("radius_domain.min_grid_points")))

        d = self._dir(GLOBAL)
        self._record(GLOBAL, write_json(d / "candidate_radius_domain.json",
                                        self.domain.as_json()))
        for p in figs.figure_f2_pcf(d, self.stats, self.peaks, self.domain):
            self._record(GLOBAL, p)

        if self.domain.fallback:
            self.books[GLOBAL].warnings.add(
                "FALLBACK_RADIUS_DOMAIN",
                f"Candidate radius domain fell back to [{self.domain.lo:g}, "
                f"{self.domain.hi:g}] A; trigger {self.domain.trigger_id} "
                f"(all fired: {','.join(self.domain.triggers_fired)}).",
                potential_consequence=(
                    "The domain was not derived from the observed PCF, so the scanned "
                    "range is the pre-registered envelope rather than data-driven."),
                recommended_action=("Report FALLBACK_RADIUS_DOMAIN = TRUE with the "
                                    "trigger ID alongside r_hot."),
                affected_output="candidate_radius_domain.json")
            self._escalate(
                f"FALLBACK_RADIUS_DOMAIN fired (trigger {self.domain.trigger_id}).",
                ["Accept the fallback-envelope result as reported.",
                 "Re-run under a new RUN_ID after the cohort or the PCF conditions change."],
                ["The scanned domain is the frozen envelope, not PCF-derived; r_hot is "
                 "still selected by the full II.6 procedure inside it.",
                 "A new RUN_ID is a Lead decision; the agent never re-derives the domain."],
                "Accept and report the fallback prominently; do not re-derive the domain.")

        if not self.domain.admissible:
            raise NegativeResult(
                "RADIUS_DOMAIN_INADMISSIBLE",
                self.domain.detail.get("inadmissible_reason", ""),
                stage=RADIUS, domain=self.domain.as_json())

        self.log.info("domain: %s [%g, %g] with %d grid points", self.domain.source,
                      self.domain.lo, self.domain.hi, len(self.domain.grid))
        self._step("II.2_candidate_radius_domain", t0)

    # -- v2 §5.4 pre-flight power certificate --------------------------------
    def power_preflight(self) -> None:
        """The certificate that decides whether this run may make a claim at all.

        Computed BEFORE any permutation compute is spent, per candidate radius, from
        sphere occupancies alone. The run is terminated pre-flight only when **no**
        radius in the domain could have rejected — the representative certificate is
        therefore the most favourable radius, never the worst.
        """
        t0 = time.perf_counter()
        inp = self.inputs
        q = float(self.cfg.get("fdr.q"))
        self.primary_null = str(self.cfg.get("permutation.primary_null"))

        self.occupancies = domain_occupancies(
            inp.universe_coords, inp.labeled_positions,
            [float(r) for r in self.domain.grid], center_coords=inp.center_coords)
        self.occupancy_by_radius = {round(o.radius_A, 6): o for o in self.occupancies}

        self.preflight_certificate, self.preflight_rows = power_preflight_certificates(
            self.primary_null, B=self.B, q=q, N_P=inp.N_P, N_B=inp.N_B,
            n_universe=inp.M, n_center_universe=inp.n_center_universe,
            occupancies=self.occupancies)
        # Per-radius certificates, keyed for the scan trace.
        self.certificate_by_radius = {round(row["radius_A"], 6): row
                                      for row in self.preflight_rows}

        d = self._dir(RADIUS)
        self._record(RADIUS, write_tsv(
            d / "power_certificate_preflight.tsv", self.preflight_rows,
            ["schema_version", "radius_A", "m", "p_res", "p_comb_best", "p_floor",
             "c_1", "p_floor_over_c_1", "binding_floor", "passes"]))
        self._write_center_universe_report()

        cert = self.preflight_certificate
        self.log.info("pre-flight power certificate: %s", cert.headline())
        if not cert.passes:
            self._raise_underpowered(cert, RADIUS)
        self._step("v2_5.4_power_certificate_preflight", t0)

    def _write_center_universe_report(self) -> None:
        """v2 §2 — U_center vs U_struct, and which residues the pLDDT floor removed.

        The number of residues removed from U_center, and whether they form
        contiguous stretches of the sequence, is an explicit v2 §2 requirement.
        "Contiguous" is defined on consecutive integer ``residue_index`` values —
        exactly how a low-confidence AlphaFold region manifests in practice.
        """
        inp = self.inputs
        excluded_mask = ~inp.center_universe_mask
        excluded_index = sorted(int(i) for i in inp.universe_index[excluded_mask])
        plddt_by_index = {int(i): float(p) for i, p in
                          zip(inp.universe_index, inp.universe_plddt)}

        regions: list[list[int]] = []
        for idx in excluded_index:
            if regions and idx == regions[-1][-1] + 1:
                regions[-1].append(idx)
            else:
                regions.append([idx])
        region_rows = [{
            "region_id": f"R{n}", "residue_index_min": run[0], "residue_index_max": run[-1],
            "n_residues": len(run),
        } for n, run in enumerate(regions, start=1)]
        region_of_index = {idx: row["region_id"] for row in region_rows
                           for idx in range(row["residue_index_min"],
                                           row["residue_index_max"] + 1)}

        d = self._dir(RADIUS)
        self._record(RADIUS, write_json(d / "center_universe.json", {
            "run_id": self.ctx.run_id, "gene": self.ctx.gene,
            "specification": "Workflow v2 §2",
            "min_plddt_threshold": inp.center_universe_min_plddt,
            "n_universe_U_struct": inp.M,
            "n_center_universe_U_center": inp.n_center_universe,
            "n_excluded_from_center_universe": inp.M - inp.n_center_universe,
            "n_contiguous_excluded_regions": len(region_rows),
            "excluded_regions": region_rows,
            "note": ("Residues below the pLDDT floor are excluded ONLY as candidate "
                     "sphere CENTERS. They remain in U_struct (the positional "
                     "reference, F9/§3 unchanged) and in the classified cohort L if "
                     "they carry a variant (F2 unaffected)."),
        }))
        self._record(RADIUS, write_tsv(
            d / "center_universe_exclusions.tsv",
            [{"schema_version": "1.0.0", "residue_index": idx,
              "plddt": plddt_by_index.get(idx), "region_id": region_of_index.get(idx, "NA")}
             for idx in excluded_index],
            ["schema_version", "residue_index", "plddt", "region_id"]))

    def _raise_underpowered(self, cert: PowerCertificate, stage: str) -> None:
        """Emit the BLOCKING warning and escalate, then terminate UNDERPOWERED."""
        self.books[stage].warnings.add(
            "TEST_CANNOT_REJECT",
            f"{cert.headline()} No center can be declared significant under any "
            f"realization of the data.",
            potential_consequence=UNINFORMATIVE_STATEMENT,
            recommended_action=(
                "; ".join(cert.remedies) if cert.remedies else
                "Lead decision — the design must change before this gene is re-run."),
            affected_output="power_certificate.json")
        options = list(cert.remedies) or ["Lead decision on the design."]
        self._escalate(
            f"POWER CERTIFICATE FAILED ({cert.phase}): p_floor = {cert.p_floor:.6g} > "
            f"c_1 = q/m = {cert.c_1:.6g} (ratio {cert.ratio:.6g}); binding floor "
            f"{cert.binding_floor}; N_P = {cert.N_P}, N_B = {cert.N_B}, m = {cert.m}, "
            f"B = {cert.B}, null = {cert.null_model}.",
            options,
            ["The run is UNDERPOWERED and uninformative about this gene; it is NOT a "
             "negative result and no absence claim may be made (v2 Appendix A).",
             "Increasing B is futile when the combinatorial floor binds; only a "
             "larger benign cohort, a different null or a smaller test family helps."],
            (f"Re-run under a NEW RUN_ID at B = {cert.B_to_clear_c_1} "
             if not cert.increasing_B_is_futile else
             "Do NOT raise B — the floor is combinatorial. ") +
            "Report the run as UNDERPOWERED; the agent changes neither B nor the "
            "FDR method.")
        raise UnderpoweredResult(cert, stage)

    # -- II.4 pre-flight -----------------------------------------------------
    def preflight_resolution(self) -> None:
        t0 = time.perf_counter()
        self.preflight = resolution.evaluate(
            "preflight", self.B, float(self.cfg.get("fdr.q")), self.inputs.M, None, None,
            int(self.cfg.get("permutation_resolution_diagnostic.safety_factor_s")),
            list(self.cfg.get("permutation_resolution_diagnostic.b_ladder")),
            int(self.cfg.get("permutation_resolution_diagnostic.max_recommendable_B")),
            k_target=int(self.cfg.get("permutation_resolution_diagnostic.k_target")))
        self.log.info("pre-flight resolution: C1=%s R=%.4g B_req=%d B_rec=%s",
                      self.preflight.conditions["C1"], self.preflight.R,
                      self.preflight.B_req, self.preflight.B_rec)
        self._step("II.4_preflight_resolution", t0)

    def _note_resolution_limited(self, diag) -> None:
        """II.4 — emit PERMUTATION_RESOLUTION_LIMITED once, whichever branch reaches it.

        Called from the normal path with the post-hoc diagnostic, and from the
        UNDERPOWERED terminal handler with the pre-flight one, so a run that stops
        before the scan still surfaces the flag it already knows about.
        """
        if self.books[FINAL].warnings.has("PERMUTATION_RESOLUTION_LIMITED"):
            return
        self.books[FINAL].warnings.add(
            "PERMUTATION_RESOLUTION_LIMITED",
            f"B = {self.B}, p_res = {diag.p_res:.6g}, m = {diag.m}, "
            f"n_floor = {diag.n_floor}, R = {diag.R:.4g}; conditions fired: "
            f"pre-flight {self.preflight.conditions}, {diag.phase} {diag.conditions}. "
            f"B_req = {diag.B_req}, B_rec = {diag.B_rec}.",
            potential_consequence=("The permutation resolution, not the data, may be "
                                   "bounding the rejection set. An empty result is "
                                   "NOT evidence of absence."),
            recommended_action=(f"Lead decision: re-run under a NEW RUN_ID with "
                                f"B = {diag.B_rec}. B is never changed automatically "
                                f"and never to chase significance."),
            affected_output="all_residue_center_tests.tsv")
        self._escalate(
            f"PERMUTATION_RESOLUTION_LIMITED = TRUE (B = {self.B}, m = {diag.m}, "
            f"p_res = {diag.p_res:.6g}, R = {diag.R:.4g}, n_floor = {diag.n_floor}).",
            [f"Re-run under a NEW RUN_ID with B = {diag.B_rec}.",
             "Accept the result with the limitation stated prominently."],
            ["A larger B costs compute but restores the ability to reject at the "
             "achieved family size.",
             "With the flag set, an empty center set must not be read as evidence "
             "of absence."],
            f"Re-run with B = {diag.B_rec} under a new RUN_ID; the agent has changed "
            f"neither B nor the FDR method.")

    def _write_resolution_diagnostic(self) -> None:
        """II.4 artefact — written on the normal path and on the UNDERPOWERED path."""
        post = getattr(self, "posthoc", None)
        payload = {
            "preflight": self.preflight.as_json(),
            "posthoc": (post.as_json() if post is not None else
                        {"computed": False,
                         "reason": "the run terminated before the final detection"}),
            "PERMUTATION_RESOLUTION_LIMITED": bool(
                getattr(self, "resolution_limited", self.preflight.limited)),
            "B_used": self.B,
            "B_recommended": (post.B_rec if post is not None else self.preflight.B_rec),
            "B_was_changed_by_the_agent": False,
            "fdr_method_was_changed_by_the_agent": False,
            "bh_margin": getattr(self, "bh_margin", None),
        }
        self._record(FINAL, write_json(
            self._dir(FINAL) / "permutation_resolution_diagnostic.json", payload))

    # -- II.5 ----------------------------------------------------------------
    def scan(self) -> None:
        t0 = time.perf_counter()
        inp = self.inputs
        q = float(self.cfg.get("fdr.q"))
        self.results = []
        for r in self.domain.grid:
            context = ctx_scan(float(r))
            # v2 §5.2/§5.5 — the primary (positional) and secondary (label) nulls
            # never share a random stream; each has its own derived seed context.
            label_context = f"{context}|secondary_label_null"
            self.results.append(scan_one_radius(
                float(r), inp.universe_coords, inp.universe_index, inp.labeled_positions,
                inp.labeled_coords, inp.y.astype(np.float64),
                self.ctx.seeds.rng(context),
                label_rng=self.ctx.seeds.rng(label_context), B=self.B, q=q,
                kappa=float(self.cfg.get("loo_mcc.kappa")),
                tie_tolerance=float(self.cfg.get("loo_mcc.tie_tolerance")),
                sparse_cutoff=int(self.cfg.get("loo_mcc.sparse_reporting_cutoff")),
                undefined_mcc=float(self.cfg.get("loo_mcc.undefined_mcc_value")),
                qc_h2_max_coverage=float(self.cfg.get("radius_qc.QC_H2_max_coverage")),
                qc_h3_min_median=float(self.cfg.get("radius_qc.QC_H3_min_median_n_labeled")),
                qc_h4_max_isolated=float(
                    self.cfg.get("radius_qc.QC_H4_max_isolated_center_fraction")),
                r_floor=self.domain.r_floor, r_cap=self.domain.r_cap,
                seed_context=context,
                center_coords=inp.center_coords, center_index=inp.center_index,
                radius_of_gyration_A=self.radius_of_gyration_A))
            self.results[-1].certificate = self.certificate_by_radius.get(
                round(float(r), 6))
            self.log.info("radius %.1f: |S|=%d mcc=%.4f zg=%.4f admissible=%s "
                          "diagnostics=%s", r,
                          len(self.results[-1].centers), self.results[-1].loo.mcc,
                          self.results[-1].perm.zg, self.results[-1].admissible,
                          self.results[-1].diagnostics_reason)

        attach_neighbor_stability(self.results, float(self.cfg.get("radius_domain.step_hot_A")))

        for r in self.results:
            if not r.boundary_neutral_verified:
                raise BlockedError(
                    f"BLOCKED — kappa boundary-neutrality failed at r = {r.radius_A} A. "
                    f"The LOO smoothing must never shift the decision boundary (F13).")

        if not any(r.admissible for r in self.results):
            self._write_scan_outputs(selection=None)
            raise NegativeResult(
                "NO_ADMISSIBLE_RADIUS",
                f"All {len(self.results)} scanned radii failed II.5E admissibility; "
                f"reasons: " + "; ".join(
                    f"{r.radius_A:g}:{r.qc_failure_reason}" for r in self.results),
                stage=RADIUS)
        self._step("II.5_radius_scan", t0)

    # -- II.6 ----------------------------------------------------------------
    def select_radius(self) -> None:
        t0 = time.perf_counter()
        self.selection = run_selection(
            self.results, [float(r) for r in self.domain.grid],
            tie_chain=list(self.cfg.get("radius_selection.tie_chain")),
            w_k=float(self.cfg.get("radius_selection.w_k")),
            metric=str(self.cfg.get("radius_selection.distance_metric")),
            near_tie_rel=float(self.cfg.get("radius_selection.near_tie_relative_threshold")),
            near_tie_sep_steps=int(
                self.cfg.get("radius_selection.near_tie_radius_separation_steps")),
            step=float(self.cfg.get("radius_domain.step_hot_A")),
            band_fraction=float(self.cfg.get("boundary_diagnostic.band_fraction")),
            bw2_fraction=float(self.cfg.get("boundary_diagnostic.bw2_pareto_fraction")),
            bw3_top_n=int(self.cfg.get("boundary_diagnostic.bw3_top_n")),
            redundancy_rho=float(self.cfg.get("radius_selection.objective_redundancy_abs_rho")),
            admissibility_dominates_fraction=float(
                self.cfg.get("radius_qc.admissibility_dominates_fraction")),
            exclude_degenerate_objectives=bool(
                self.cfg.get("radius_selection.exclude_degenerate_objectives")))

        if not self.selection.degenerate_exclusion_is_numerically_neutral:
            raise BlockedError(
                "BLOCKED — excluding the degenerate objectives from distance-to-utopia "
                "changes a distance. Under the frozen min-max rule a degenerate "
                "objective normalizes to 1.0 everywhere and must contribute exactly 0; "
                "a discrepancy means the selection machinery has drifted.")

        if not self.selection.independent_pareto_check:
            raise BlockedError(
                "BLOCKED — the independent dominance re-check disagrees with the Pareto "
                "front returned by utils.multiobjective. Selection cannot be trusted.")

        self.r_hot = float(self.selection.selected)
        if not any(abs(self.r_hot - float(g)) < 1e-9 for g in self.domain.grid):
            raise BlockedError(
                f"BLOCKED — selected r_hot = {self.r_hot} is not on the scanned grid.")

        sel = self.selection

        # --- v2 §4 guards on the selection itself ------------------------------
        if sel.degenerate_objectives:
            self.books[RADIUS].warnings.add(
                "DEGENERATE_OBJECTIVE",
                f"Objective(s) {sel.degenerate_objectives} take the identical value at "
                f"every admissible radius and discriminate nothing. Excluded from "
                f"distance-to-utopia; the effective objective set is "
                f"{sel.effective_objectives}.",
                potential_consequence=(
                    "Fewer objectives are actually deciding than the four the design "
                    "declares, so 'multi-objective' overstates the decision."),
                recommended_action=("Report the effective objective set alongside "
                                    "r_hot. The objective set is frozen for this run."),
                affected_output="radius_decision.json")
        if sel.vacuous_pareto:
            self.books[RADIUS].warnings.add(
                "VACUOUS_PARETO_SELECTION",
                f"VACUOUS_PARETO_SELECTION: {sel.vacuous_reason}.",
                potential_consequence=(
                    "r_hot was not chosen by trading objectives off against each "
                    "other. Combined with a boundary optimum this is a symptom of a "
                    "power failure (v2 §4/§5.4), not a property of the radius domain."),
                recommended_action=("Report the Pareto set size and the effective "
                                    "objective set prominently beside r_hot."),
                affected_output="pareto_front.json")
            self._escalate(
                f"VACUOUS_PARETO_SELECTION at r_hot = {self.r_hot:g} A: "
                f"{sel.vacuous_reason}.",
                ["Accept r_hot and report the vacuity prominently (current behaviour).",
                 "Lead review of whether the objective set discriminates at all for "
                 "this gene, under a NEW RUN_ID."],
                ["The reported radius carries no trade-off information.",
                 "Changing the objective set after seeing results would be "
                 "outcome-dependent and is prohibited."],
                "Accept and report; the objective set is frozen for this run.")
        if sel.admissibility_dominates:
            self.books[RADIUS].warnings.add(
                "ADMISSIBILITY_DOMINATES_SELECTION",
                f"{sel.inadmissible_fraction:.1%} of scanned radii were ruled "
                f"inadmissible (threshold "
                f"{self.cfg.get('radius_qc.admissibility_dominates_fraction'):.0%}); "
                f"by rule: {sel.inadmissible_by_rule}.",
                potential_consequence=("Admissibility, not the objectives, decided the "
                                       "outcome."),
                recommended_action=("Report which constraint did the eliminating "
                                    "(v2 §4). Significance is never one of them."),
                affected_output="r_hot_scan.tsv")

        if sel.result.near_tie:
            det = sel.result.near_tie_detail
            self.books[RADIUS].warnings.add(
                "NEAR_TIE_RADIUS_SELECTION",
                f"r = {det['first']:g} A and r = {det['second']:g} A differ by "
                f"{det['relative_gap']:.4g} in distance-to-ideal (< "
                f"{self.cfg.get('radius_selection.near_tie_relative_threshold')}) while "
                f"being more than 2 steps apart.",
                potential_consequence="A materially different radius was nearly selected.",
                recommended_action="Non-blocking. Report both radii and their objectives.",
                affected_output="radius_decision.json")
            self._escalate(
                f"Near-tie between r_hot = {det['first']:g} A and {det['second']:g} A.",
                ["Accept the tie-chain result (current behaviour).",
                 "Lead review of both candidates' objective vectors."],
                ["The tie chain is deterministic and pre-registered, so the result is "
                 "reproducible either way.",
                 "Changing the selection post hoc would be outcome-dependent and is "
                 "prohibited."],
                "Accept the tie-chain result and report the runner-up alongside it.")

        if sel.redundant_pairs:
            pairs = ", ".join(f"{p['objective_a']}~{p['objective_b']}(rho={p['rho']:.4f})"
                              for p in sel.redundant_pairs)
            self.books[RADIUS].warnings.add(
                "OBJECTIVE_REDUNDANCY",
                f"Objective correlation exceeds "
                f"{self.cfg.get('radius_selection.objective_redundancy_abs_rho')}: {pairs}",
                potential_consequence=("Two objectives carry nearly the same information, "
                                       "so the effective weighting is not equal."),
                recommended_action="Lead decision; the objective set is frozen for this run.",
                affected_output="objective_matrix.json")
            self._escalate(
                f"Objective redundancy |rho| > "
                f"{self.cfg.get('radius_selection.objective_redundancy_abs_rho')}: {pairs}",
                ["Report and proceed with the frozen objective set.",
                 "Pre-register a revised objective set under a new RUN_ID."],
                ["Equal weights are effectively applied to fewer than four independent "
                 "objectives.",
                 "Changing objectives now would be outcome-dependent."],
                "Report and proceed; the objective set is frozen for this run.")

        if sel.boundary.get("DOMAIN_BOUNDARY_WARNING"):
            self.books[RADIUS].warnings.add(
                "BOUNDARY_OPTIMUM_WARNING",
                f"DOMAIN_BOUNDARY_WARNING: BW1={sel.boundary['BW1_selected_in_band']} "
                f"BW2={sel.boundary['BW2_pareto_majority_in_band']} "
                f"BW3={sel.boundary['BW3_top3_same_band']}; band "
                f"{sel.boundary['band_members']}, end {sel.boundary['which_end']}.",
                potential_consequence=("The optimum sits against a domain edge; a better "
                                       "radius may lie outside the pre-registered domain."),
                recommended_action=("Lead-authorized re-run under an explicitly widened, "
                                    "pre-registered domain and a NEW RUN_ID. The agent "
                                    "never widens the domain itself."),
                affected_output="domain_boundary_diagnostic.json")
            binding = ("R_FLOOR" if "lower" in sel.boundary["which_end"]
                       else "r_cap = min(R_CEIL, 0.25*D_max)")
            self._escalate(
                f"DOMAIN_BOUNDARY_WARNING at the {sel.boundary['which_end']} end of "
                f"[{self.domain.lo:g}, {self.domain.hi:g}] A; binding constraint {binding}.",
                ["Accept r_hot as selected inside the pre-registered domain.",
                 "Lead-authorized re-run under a widened, pre-registered domain "
                 "(NEW RUN_ID)."],
                ["The reported r_hot may be an edge artefact of the domain rather than a "
                 "true interior optimum.",
                 "Automatic widening is prohibited — it would be domain shopping."],
                "Report the warning; the Lead decides whether to authorize a widened re-run.")

        # --- v2 §4: biological-plausibility DIAGNOSTICS at the SELECTED radius ----
        # Computed at every candidate radius (r_hot_scan.tsv) but only WARNED on
        # here, at r_hot — these never remove a radius from the Pareto set and
        # never enter the objective vector (F15 stays at exactly four objectives;
        # see config/pipeline.yaml `biological_plausibility` for the Lead ruling).
        bio = next(r for r in self.results if r.radius_A == self.r_hot)
        cov_major = float(
            self.cfg.get("biological_plausibility.structural_coverage_major_threshold"))
        cov_soft = float(
            self.cfg.get("biological_plausibility.structural_coverage_penalty_threshold"))
        absorb_soft = float(
            self.cfg.get("biological_plausibility.cohort_absorption_penalty_threshold"))
        ratio_soft = float(
            self.cfg.get("biological_plausibility.r_to_domain_ratio_penalty_threshold"))
        if bio.structural_coverage > cov_major:
            self.books[RADIUS].warnings.add(
                "EXCESSIVE_STRUCTURAL_COVERAGE",
                f"structural_coverage = {bio.structural_coverage:.4g} at the selected "
                f"r_hot = {self.r_hot:g} A exceeds the MAJOR bound {cov_major:g} "
                f"(soft penalty threshold {cov_soft:g}).",
                potential_consequence=("A sphere covering most of U_struct localizes "
                                       "nothing; the finding may not be a spatially "
                                       "restricted hotspot."),
                recommended_action=("Report prominently. Never selected against — "
                                    "coverage is a diagnostic, not an admissibility "
                                    "constraint or a Pareto objective (v2 §4)."),
                affected_output="r_hot_scan.tsv")
        if bio.cohort_absorption > absorb_soft:
            self.books[RADIUS].warnings.add(
                "HIGH_COHORT_ABSORPTION",
                f"cohort_absorption = {bio.cohort_absorption:.4g} at the selected "
                f"r_hot = {self.r_hot:g} A exceeds {absorb_soft:g}: the single largest "
                f"significant sphere holds most of the classified cohort L.",
                potential_consequence=("A radius whose sphere absorbs most of the "
                                       "cohort does not localize anything."),
                recommended_action="Report as a caveat; never used to reselect r_hot.",
                affected_output="r_hot_scan.tsv")
        if bio.r_to_domain_ratio > ratio_soft:
            self.books[RADIUS].warnings.add(
                "HIGH_R_TO_DOMAIN_RATIO",
                f"r_to_domain_ratio = {bio.r_to_domain_ratio:.4g} at the selected "
                f"r_hot = {self.r_hot:g} A exceeds {ratio_soft:g} "
                f"(radius of gyration = {self.radius_of_gyration_A:.4g} A).",
                potential_consequence=("The radius is large relative to the modelled "
                                       "structure's own scale."),
                recommended_action="Report as a caveat; never used to reselect r_hot.",
                affected_output="r_hot_scan.tsv")

        self._write_scan_outputs(selection=self.selection)
        self.log.info("selected r_hot = %.2f A (distance %.6g, Pareto members %s)",
                      self.r_hot, sel.result.distances[self.r_hot],
                      sel.result.pareto_members)
        self._step("II.6_radius_selection", t0)

    def _scan_row(self, r, selection) -> dict:
        norm = (selection.result.normalized.get(r.radius_A)
                if (selection and r.admissible) else None)
        boundary = bool(selection.boundary.get("DOMAIN_BOUNDARY_WARNING")) if selection else None
        m = r.loo.metrics
        cert = r.certificate or {}
        dv = r.diagnostic_values
        return {
            "radius_A": r.radius_A,
            "n_labeled_in_universe": r.n_labeled_in,
            "loo_mcc": r.loo.mcc, "loo_sens": m["sens"], "loo_spec": m["spec"],
            "loo_ppv": m["ppv"], "loo_npv": m["npv"], "loo_f1": m["f1"],
            "loo_balacc": m["balacc"],
            "loo_n_zero_neighbour": r.loo.n_zero_neighbour,
            "loo_prop_zero_neighbour": r.loo.prop_zero_neighbour,
            "loo_n_sparse": r.loo.n_sparse, "loo_n_tie": r.loo.n_tie,
            "perm_evidence_raw": r.perm.zg,
            "perm_evidence_normalized": norm["perm_evidence_zg"] if norm else None,
            "fold_enrichment": r.fold_enrichment,
            "fold_enrichment_normalized": norm["fold_enrichment"] if norm else None,
            "neighbor_stability": r.neighbor_stability,
            "neighbor_stability_normalized": norm["neighbor_stability"] if norm else None,
            "loo_mcc_normalized": norm["loo_mcc"] if norm else None,
            "n_significant_centers": len(r.centers),
            "coverage_fraction": r.coverage_fraction,
            "n_singleton_centers": r.n_singleton_centers,
            "qc_status": "PASS" if r.admissible else "FAIL",
            "qc_failure_reason": r.qc_failure_reason,
            "pareto_member": ((r.radius_A in selection.result.pareto_members)
                              if (selection and r.admissible) else None),
            "distance_to_ideal": (selection.result.distances.get(r.radius_A)
                                  if selection else None),
            "selected": bool(selection and r.radius_A == selection.selected),
            "search_domain_source": self.domain.source,
            "boundary_warning": boundary,
            # --- v2 §4: the admissibility verdict, its rule and its reason -----
            "admissible": r.admissible,
            "admissibility_rule": r.qc_failure_reason,
            "admissibility_reason": r.admissibility_reason,
            "significance_diagnostics_fired": r.diagnostics_reason,
            "qc_h1_zero_significant_centers": dv["QC_H1_zero_significant_centers"],
            "qc_h2_coverage_exceeds_max": dv["QC_H2_exceeds_max"],
            "qc_h3_median_n_labeled": dv["QC_H3_median_n_labeled"],
            "qc_h3_below_threshold": dv["QC_H3_below_threshold"],
            "qc_h4_isolated_center_fraction": dv["QC_H4_isolated_center_fraction"],
            "qc_h4_exceeds_max": dv["QC_H4_exceeds_max"],
            "fold_enrichment_defined": r.fold_enrichment_defined,
            "n_in_test_family": r.n_in_family,
            # --- v2 §5.2 / §5.4 ------------------------------------------------
            "primary_null": r.extra["primary_null"],
            "secondary_null": r.extra["secondary_null"],
            "power_p_res": cert.get("p_res"),
            "power_p_comb_best": cert.get("p_comb_best"),
            "power_p_floor": cert.get("p_floor"),
            "power_c_1": cert.get("c_1"),
            "power_certificate_passes": cert.get("passes"),
            "power_binding_floor": cert.get("binding_floor"),
            # --- v2 §4: biological-plausibility DIAGNOSTICS, every candidate radius --
            "structural_coverage": r.structural_coverage,
            "cohort_absorption": r.cohort_absorption,
            "r_to_domain_ratio": r.r_to_domain_ratio,
        }

    def _write_scan_outputs(self, selection) -> None:
        d = self._dir(RADIUS)
        rows = [self._scan_row(r, selection) for r in self.results]
        self._record(RADIUS, write_tsv(d / "r_hot_scan.tsv", rows, R_HOT_SCAN_COLUMNS))

        for r in self.results:
            path = d / "radius_scan_centers" / f"r_{r.radius_A:.1f}.tsv"
            self._record(RADIUS, write_tsv(path, r.center_rows, SCAN_CENTER_COLUMNS))

        loo_rows = [row for r in self.results for row in r.loo.diagnostics]
        self._record(RADIUS, write_tsv(d / "loo_diagnostics.tsv", loo_rows, SPARSE_COLUMNS))

        correlation = selection.result.correlation if selection else {}
        self._record(RADIUS, write_json(d / "objective_matrix.json", {
            "objectives": OBJECTIVE_NAMES,
            "directions": {n: "maximize" for n in OBJECTIVE_NAMES},
            "weights": {n: 1.0 for n in OBJECTIVE_NAMES},
            "weights_frozen_equal_w_k_1": True,
            "normalization": "min_max", "normalization_domain": "admissible_set",
            "correlation_matrix": correlation,
            "redundancy_threshold_abs_rho": float(
                self.cfg.get("radius_selection.objective_redundancy_abs_rho")),
            "correlation_informative": (selection.correlation_informative
                                        if selection else False),
            "correlation_informative_rule": (
                "A Pearson correlation over fewer than three admissible candidates is "
                "+/-1 by arithmetic; the redundancy diagnostic is then reported as "
                "uninformative rather than raised as a finding."),
            "redundant_pairs": selection.redundant_pairs if selection else [],
            "per_objective_argmax": selection.per_objective_argmax if selection else {},
            "raw_matrix": {f"{r.radius_A:g}": r.objectives() for r in self.results},
            "admissible": {f"{r.radius_A:g}": r.admissible for r in self.results},
            "qc_rules": QC_RULES,
            "admissibility_constraints": list(
                self.cfg.get("radius_qc.admissibility_constraints")),
            "significance_conditioned_diagnostics": list(
                self.cfg.get("radius_qc.significance_conditioned_diagnostics")),
            "significance_may_determine_admissibility": False,
            "effective_objectives_in_distance": (selection.effective_objectives
                                                 if selection else []),
            "degenerate_objectives_excluded": (selection.degenerate_objectives
                                               if selection else []),
            "perm_evidence_zg_direction": ZG_DIRECTION_NOTE,
            "perm_evidence_zg_null_model": "label_permutation (SECONDARY, v2 §5.2)",
            "note": ("Objectives are never scalarized beyond the frozen min-max "
                     "normalization; LOO-MCC is a SELECTION objective only and is never "
                     "hotspot validation (F13). v2 §4: significance is an objective, "
                     "never an admissibility filter."),
        }))
        self._record(RADIUS, write_json(d / "pareto_front.json", {
            "pareto_members": selection.result.pareto_members if selection else [],
            "n_admissible": selection.result.n_admissible if selection else 0,
            "independent_dominance_check_passed": (
                selection.independent_pareto_check if selection else None),
            "normalized": ({f"{k:g}": v for k, v in selection.result.normalized.items()}
                           if selection else {}),
            "distances": ({f"{k:g}": v for k, v in selection.result.distances.items()}
                          if selection else {}),
            "degenerate_objectives": selection.result.degenerate_objectives if selection else [],
            "effective_objectives_in_distance": (selection.effective_objectives
                                                 if selection else []),
            "VACUOUS_PARETO_SELECTION": bool(selection.vacuous_pareto) if selection else None,
            "vacuous_pareto_reason": selection.vacuous_reason if selection else "NA",
            "utopia": selection.result.utopia if selection else {},
        }))
        self._record(RADIUS, write_json(d / "pareto_dominated.json", {
            "dominated": ({f"{k:g}": v for k, v in selection.result.dominators.items()}
                          if selection else {}),
            "inadmissible": {f"{r.radius_A:g}": r.qc_failure_reason
                             for r in self.results if not r.admissible},
            "significance_conditioned_diagnostics_fired": {
                f"{r.radius_A:g}": r.diagnostics_reason for r in self.results
                if r.diagnostics_fired},
            "ADMISSIBILITY_DOMINATES_SELECTION": (bool(selection.admissibility_dominates)
                                                  if selection else None),
            "inadmissible_fraction": (selection.inadmissible_fraction
                                      if selection else None),
            "inadmissible_by_rule": (selection.inadmissible_by_rule if selection else {}),
            "note": ("Dominated and inadmissible candidates are retained in full with "
                     "their dominating vectors or failing rules — a scan that drops "
                     "radii cannot be audited. v2 §4: a radius is inadmissible ONLY "
                     "when a declared quantity is genuinely undefined; zero significant "
                     "centers is a diagnostic (QC_H1), never an admissibility failure."),
        }))
        self._record(RADIUS, write_json(d / "radius_decision.json", decision_payload(
            selection, self.results, self.domain.source, self.domain.fallback)
            if selection else {
                "selected_r_hot_A": None,
                "search_domain_source": self.domain.source,
                "FALLBACK_RADIUS_DOMAIN": self.domain.fallback,
                "n_candidates": len(self.results), "n_admissible": 0,
                "negative_result": "NO_ADMISSIBLE_RADIUS",
                "candidates": [{"radius_A": r.radius_A, "admissible": False,
                                "qc_failure_reason": r.qc_failure_reason,
                                "objectives_raw": r.objectives()} for r in self.results]}))
        self._record(RADIUS, write_json(d / "domain_boundary_diagnostic.json", {
            **(selection.boundary if selection else {"DOMAIN_BOUNDARY_WARNING": False,
                                                     "reason": "no admissible candidate"}),
            "evaluated": True,
            "grid_A": [float(g) for g in self.domain.grid],
            "domain_lo_A": self.domain.lo, "domain_hi_A": self.domain.hi,
            "search_domain_source": self.domain.source,
            "FALLBACK_RADIUS_DOMAIN": self.domain.fallback,
            "automatic_domain_expansion": False,
        }))

        step = float(self.cfg.get("radius_domain.step_hot_A"))
        selected = selection.selected if selection else None
        for p in figs.figure_f3_scan(d, rows, selected, step):
            self._record(RADIUS, p)
        for p in figs.figure_f4_pareto(
                d, selection.result.normalized if selection else {},
                selection.result.pareto_members if selection else [],
                OBJECTIVE_NAMES, selected):
            self._record(RADIUS, p)
        for p in figs.figure_f5_distance(
                d, selection.result.distances if selection else {},
                selection.result.pareto_members if selection else [], selected,
                selection.boundary if selection else {}, step):
            self._record(RADIUS, p)

    # -- II.7 ----------------------------------------------------------------
    def detect(self) -> None:
        t0 = time.perf_counter()
        inp = self.inputs
        q = float(self.cfg.get("fdr.q"))
        context = ctx_final_detection()
        label_context = f"{context}|secondary_label_null"
        self.detection = final_detection(
            self.r_hot, inp.universe_coords, inp.universe_index, inp.labeled_positions,
            inp.labeled_coords, inp.y.astype(np.float64), self.ctx.seeds.rng(context),
            B=self.B, q=q, seed_context=context,
            label_rng=(self.ctx.seeds.rng(label_context)
                       if self.cfg.get("permutation.report_secondary_null") else None),
            label_seed_context=label_context,
            center_coords=inp.center_coords, center_index=inp.center_index)
        # final_detection() runs the B-draw positional-null permutation test then
        # Benjamini-Hochberg FDR internally (detection.py: positional_pass() then
        # benjamini_hochberg(), sequential, fused into one call) -- timed here as
        # one combined step; splitting further would mean instrumenting inside
        # detection.py itself, which this timing pass does not touch.
        self._step("II.7_final_detection_permutation_and_fdr", t0)
        t0_cert = time.perf_counter()

        det = self.detection
        if not significant_subset_matches(det.center_rows, det.significant_rows):
            raise BlockedError(
                "BLOCKED — significant_hotspot_centers is not the exact significant=TRUE "
                "subset of all_residue_center_tests.")

        # How close the evidence came to the BH bar, stated so a reader can see the
        # margin without recomputing it. The rank-1 critical value q/m is the hardest
        # bar in the procedure: a single isolated center must clear it alone.
        m_family = len(det.family_positions)
        rank1 = (q / m_family) if m_family else None
        smallest_p = min((r["p_emp"] for r in det.bh_rows), default=None)

        # --- v2 §5.4 POST-HOC power certificate, on the realized family --------
        self.posthoc_certificate = certify(
            "posthoc", null_model=self.primary_null, B=self.B, q=q, m=m_family,
            N_P=inp.N_P, N_B=inp.N_B, n_universe=inp.M,
            n_center_universe=inp.n_center_universe,
            occupancies_labeled=det.perm.n_labeled[det.family_positions],
            occupancies_universe=det.perm.n_universe[det.family_positions],
            radius_A=self.r_hot,
            extra={"smallest_observed_p_emp": smallest_p,
                   "smallest_observed_p_exact": (
                       float(det.perm.p_exact[det.family_positions].min())
                       if m_family else None),
                   "n_significant_bh": det.n_significant})
        self.log.info("post-hoc power certificate: %s",
                      self.posthoc_certificate.headline())
        self.certificate_passes = bool(self.posthoc_certificate.passes)

        self.bh_margin = {
            "m_family": m_family,
            "rank1_critical_value": rank1,
            "smallest_observed_p": smallest_p,
            "p_resolution_floor": 1.0 / (self.B + 1),
            "margin_ratio": ((smallest_p / rank1)
                             if (rank1 and smallest_p is not None and rank1 > 0) else None),
            "power_certificate_passed": self.certificate_passes,
            "margin_ratio_note": (
                "smallest observed p divided by the rank-1 BH critical value q/m. "
                "Below 1 the smallest p clears the hardest bar on its own; above 1 a "
                "rejection requires several centers to clear their joint bar together. "
                "A ratio near 1 means the result is close to the boundary and can move "
                "between RUN_IDs, because permutation p-values are discrete."
                if self.certificate_passes else
                "smallest observed p divided by the rank-1 BH critical value q/m. The "
                "v2 §5.4 power certificate FAILED for this configuration, so this "
                "ratio must NOT be read as proximity to the significance boundary: "
                "the test was operating at its floor and could not have rejected "
                "regardless of the data."),
            "top_ranks": [
                {"rank": r["rank"], "p_emp": r["p_emp"],
                 "bh_critical_value": r["bh_critical_value"], "reject_bh": r["reject_bh"]}
                for r in det.bh_rows[:10]],
        }

        self.posthoc = resolution.evaluate(
            "posthoc", self.B, q, len(det.family_positions),
            det.perm.p_emp[det.family_positions], det.fdr.reject_bh,
            int(self.cfg.get("permutation_resolution_diagnostic.safety_factor_s")),
            list(self.cfg.get("permutation_resolution_diagnostic.b_ladder")),
            int(self.cfg.get("permutation_resolution_diagnostic.max_recommendable_B")),
            k_target=int(self.cfg.get("permutation_resolution_diagnostic.k_target")))
        self.resolution_limited = bool(self.preflight.limited or self.posthoc.limited)

        if self.resolution_limited:
            self._note_resolution_limited(self.posthoc)
        if self.posthoc.escalate_above_ladder:
            self._escalate(
                f"B_req = {self.posthoc.B_req} exceeds the maximum recommendable B "
                f"({self.cfg.get('permutation_resolution_diagnostic.max_recommendable_B')}).",
                ["Reduce the test family by a pre-registered restriction of U_struct.",
                 "Accept that single-center detection is not attainable at this family size."],
                ["Any restriction of U_struct is a methodological change requiring "
                 "pre-registration.",
                 "The reported result remains valid but is resolution-bounded."],
                "Escalate to the Lead; no ladder value is recommended above the cap.")

        if det.fdr.n_reject_bh != det.fdr.n_reject_by:
            self.books[FINAL].warnings.add(
                "BH_BY_DISAGREEMENT",
                f"BH rejects {det.fdr.n_reject_bh} center(s); BY rejects "
                f"{det.fdr.n_reject_by} (c(m) = {det.fdr.by_constant:.4g}).",
                potential_consequence=("Under arbitrary dependence the discovery set would "
                                       "be smaller. BH remains the primary and only "
                                       "decision procedure (A10)."),
                recommended_action="Report both counts; never substitute BY for BH.",
                affected_output="bh_fdr_table.tsv")

        # The provisional S(r_hot) computed during the scan and the final S at the same
        # radius come from DIFFERENT draws (II.12 / v2 §5.5 give the scan and the final
        # detection different seed contexts), so they can disagree. How that
        # disagreement may be DESCRIBED depends on the §5.4 certificate.
        scan_at_r_hot = next(r for r in self.results if r.radius_A == self.r_hot)
        self.n_provisional_at_r_hot = len(scan_at_r_hot.centers)
        self.scan_final_disagreement = (self.n_provisional_at_r_hot != det.n_significant)
        self.scan_final_note = (SCAN_FINAL_NOTE_CERTIFIED if self.certificate_passes
                                else SCAN_FINAL_NOTE_UNCERTIFIED)

        # --- v2 §5.4: the certificate is BINDING at the post-hoc instance too ---
        if not self.certificate_passes:
            self._write_final_outputs()
            if det.n_significant:
                # p_floor > c_1 rules out a RANK-1 rejection; BH can still reject
                # several centers jointly at rank k > 1. v2 §5.4 nevertheless states
                # the run "must terminate as UNDERPOWERED", so the rejection is
                # recorded and escalated rather than published as a finding.
                self._escalate(
                    f"The post-hoc §5.4 certificate failed (p_floor = "
                    f"{self.posthoc_certificate.p_floor:.6g} > c_1 = "
                    f"{self.posthoc_certificate.c_1:.6g}) yet BH rejected "
                    f"{det.n_significant} center(s) jointly at rank > 1.",
                    ["Terminate UNDERPOWERED as v2 §5.4 directs (current behaviour); "
                     "the rejection is recorded in all_residue_center_tests.tsv but "
                     "not published as a finding.",
                     "Lead ruling on whether §5.4 was intended to bind only the rank-1 "
                     "case, under a NEW RUN_ID."],
                    ["A potentially real finding is withheld from downstream stages.",
                     "Publishing it would contradict the pre-registered v2 §5.4 rule, "
                     "which is stated without a rank-1 exception."],
                    "Terminate UNDERPOWERED and escalate the tension to the Lead; the "
                    "agent never relaxes a pre-registered rule on its own.")
            self._raise_underpowered(self.posthoc_certificate, FINAL)

        if det.n_significant == 0:
            # Reached ONLY with a PASSED certificate, so this is a real negative:
            # the design demonstrably could have rejected and did not.
            self.books[FINAL].warnings.add(
                "NO_SIGNIFICANT_HOTSPOTS",
                f"No center of U_struct survives BH at q = {q} (family size "
                f"{len(det.family_positions)}, B = {self.B}, null = "
                f"{self.primary_null}). The v2 §5.4 power certificate PASSED "
                f"(p_floor = {self.posthoc_certificate.p_floor:.6g} <= c_1 = "
                f"{self.posthoc_certificate.c_1:.6g}), so the design could have "
                f"rejected. The scan found {self.n_provisional_at_r_hot} provisional "
                f"center(s) at this radius under its own draw (seed context "
                f"{ctx_scan(self.r_hot)!r} vs {context!r}).",
                potential_consequence=(
                    "No hotspot exists at r_hot under the pre-registered procedure."
                    if not self.resolution_limited else
                    "The empty set may be a resolution artefact, NOT evidence of absence."),
                recommended_action=("Report as a valid negative result, with the power "
                                    "certificate attached as evidence (v2 §5.4)."
                                    if not self.resolution_limited else
                                    "Report as RESOLUTION-LIMITED and escalate."),
                affected_output="significant_hotspot_centers.tsv")

        self._write_final_outputs()
        self.log.info("final detection: |S| = %d, %d without a ClinVar variant, "
                      "boundary p = %s", det.n_significant, det.n_centers_without_variant,
                      det.boundary_p)
        self._step("v2_5.4_power_certificate_posthoc_and_bookkeeping", t0_cert)

    def _write_final_outputs(self) -> None:
        d = self._dir(FINAL)
        det = self.detection
        self._record(FINAL, write_tsv(d / "all_residue_center_tests.tsv",
                                      det.center_rows, ALL_CENTER_TEST_COLUMNS))
        self._record(FINAL, write_tsv(d / "significant_hotspot_centers.tsv",
                                      det.significant_rows, ALL_CENTER_TEST_COLUMNS))
        self._record(FINAL, write_tsv(d / "hotspot_classified_variants.tsv",
                                      det.classified_rows, CLASSIFIED_VARIANT_COLUMNS))
        self._record(FINAL, write_tsv(d / "hotspot_covered_residues.tsv",
                                      det.covered_rows, COVERED_RESIDUE_COLUMNS))
        self._record(FINAL, write_tsv(d / "hotspot_regions.tsv",
                                      det.region_rows, HOTSPOT_REGION_COLUMNS))
        self._record(FINAL, write_tsv(d / "bh_fdr_table.tsv", det.bh_rows, BH_TABLE_COLUMNS))
        self._write_resolution_diagnostic()
        self._write_power_certificate()
        self._record(FINAL, write_json(d / "secondary_null_label_permutation.json", {
            "run_id": self.ctx.run_id, "gene": self.ctx.gene,
            "radius_A": self.r_hot,
            "primary_null": self.primary_null,
            "secondary_null": "label_permutation",
            "primary_n_significant_bh": self.detection.n_significant,
            "primary_significant_center_residue_indices": sorted(
                r["center_residue_index"] for r in self.detection.significant_rows),
            **self.detection.secondary,
            "v2_5_2": ("The two nulls answer different questions and must not be mixed "
                       "across the stages of one run. This block is reported so a "
                       "reader can see BOTH, and is never aggregated with, substituted "
                       "for, or promoted over the primary result."),
        }))
        self._write_structures(d)
        for p in figs.figure_f6_hotspot_map(d, det.center_rows,
                                            float(self.cfg.get("fdr.q")),
                                            det.boundary_p, self.r_hot):
            self._record(FINAL, p)

    def _write_power_certificate(self) -> None:
        """v2 §11 — the certificate is a first-class artefact in EVERY terminal state.

        Written from whatever exists: a run blocked before the pre-flight still gets
        a file saying so, because "the certificate is missing" and "the certificate
        passed" must never look the same on disk.
        """
        pre = getattr(self, "preflight_certificate", None)
        post = getattr(self, "posthoc_certificate", None)
        passed = None
        if post is not None:
            passed = bool(post.passes)
        elif pre is not None:
            passed = bool(pre.passes)
        payload = {
            "run_id": self.ctx.run_id, "gene": self.ctx.gene,
            "specification": "Workflow v2 §5.4",
            "enabled": bool(self.cfg.get("power_certificate.enabled")),
            "primary_null": getattr(self, "primary_null", None),
            "preflight": pre.as_json() if pre is not None else {
                "computed": False,
                "reason": "the run terminated before the pre-flight certificate"},
            "posthoc": post.as_json() if post is not None else {
                "computed": False,
                "reason": "the run terminated before the final detection"},
            "PASSES": passed,
            "TEST_CANNOT_REJECT": (None if passed is None else not passed),
            "terminal_state_rule": (
                "A run may be reported as COMPLETED_NEGATIVE only after this "
                "certificate has been computed and has passed (v2 §0/§5.4). A failed "
                "certificate terminates the run as UNDERPOWERED, which is NOT a result "
                "about the gene."),
            "uninformative_statement_if_failed": UNINFORMATIVE_STATEMENT,
            "B_was_changed_by_the_agent": False,
            "fdr_method_was_changed_by_the_agent": False,
        }
        path = write_json(self._dir(FINAL) / "power_certificate.json", payload)
        if path not in self.books[FINAL].created:
            self._record(FINAL, path)

    def _write_structures(self, d: Path) -> None:
        det = self.detection
        inp = self.inputs
        hotspots: dict[str, list[int]] = {}
        for row in det.region_rows:
            hotspots[row["hotspot_id"]] = [int(x) for x in
                                           row["center_residue_indices"].split(";")]

        centers = [CaRecord(residue_index=r["center_residue_index"], x=r["x_ca"],
                            y=r["y_ca"], z=r["z_ca"],
                            bfactor=float(int(r["hotspot_id"][1:]) if
                                          r["hotspot_id"] != "NA" else 0))
                   for r in det.significant_rows]
        omitted = write_ca_pdb(
            d / "structures" / "final_hotspots.pdb", centers,
            f"STAGE B HOTSPOT CENTERS r_hot={self.r_hot:g}A",
            remarks=[f"run_id {self.ctx.run_id}",
                     f"{len(centers)} significant centers, BH q<={self.cfg.get('fdr.q')}",
                     "B-factor column carries the hotspot region number",
                     "resname UNK: Stage B never loads the amino-acid identity column"])
        self._record(FINAL, d / "structures" / "final_hotspots.pdb")
        self._record(FINAL, write_ca_cif(
            d / "structures" / "final_hotspots.cif", centers,
            comments=[f"hotspot3d Stage B significant centers, r_hot = {self.r_hot:g} A",
                      f"run_id {self.ctx.run_id}"]))

        # v2 §2 — det.center_rows now covers ONLY U_center; a residue of U_struct
        # outside U_center was never tested and gets b=0 (indistinguishable on
        # this map from a tested-but-non-significant residue, which is correct:
        # neither carries a q-value).
        q_map = {r["center_residue_index"]: r for r in det.center_rows}
        bfactor_records = []
        for i, idx in enumerate(inp.universe_index):
            row = q_map.get(int(idx))
            if row is None:
                b = 0.0
            else:
                qv = row["q_bh"]
                b = 0.0 if (not row["significant"] or qv is None or qv <= 0) else float(
                    min(-np.log10(max(qv, 1e-300)), 999.99))
            bfactor_records.append(CaRecord(
                residue_index=int(idx), x=float(inp.universe_coords[i, 0]),
                y=float(inp.universe_coords[i, 1]), z=float(inp.universe_coords[i, 2]),
                bfactor=b))
        omitted += write_ca_pdb(
            d / "structures" / "hotspot_centers_bfactor_qvalue.pdb", bfactor_records,
            "B-FACTOR = -LOG10(Q_BH), 0 FOR NON-SIGNIFICANT OR UNTESTED",
            remarks=["every residue of U_center is a candidate sphere center (F9, "
                     "restricted from U_struct by v2 sec.2 pLDDT >= "
                     f"{self.inputs.center_universe_min_plddt:g})",
                     "residues of U_struct outside U_center were never tested and "
                     "carry b=0, same as a tested non-significant residue",
                     f"r_hot = {self.r_hot:g} A, q = {self.cfg.get('fdr.q')}"])
        self._record(FINAL, d / "structures" / "hotspot_centers_bfactor_qvalue.pdb")
        self._record(FINAL, write_chimerax_script(
            d / "structures" / "view_hotspots.cxc", hotspots, self.r_hot))
        self._record(FINAL, write_pymol_script(
            d / "structures" / "view_hotspots.pml", hotspots, self.r_hot))

        if omitted:
            self.books[FINAL].warnings.add(
                "RENDERING_UNAVAILABLE",
                f"{omitted} residue(s) exceed the PDB 4-digit residue-number field and "
                f"were omitted from the .pdb exports; the .cif file is complete.",
                potential_consequence="PDB viewers see a truncated model.",
                recommended_action="Use final_hotspots.cif for the complete model.",
                affected_output="structures/final_hotspots.pdb")

    # -- II.7 sensitivity ----------------------------------------------------
    def sensitivity(self) -> None:
        t0 = time.perf_counter()
        inp = self.inputs
        q = float(self.cfg.get("fdr.q"))
        primary = {r["center_residue_index"] for r in self.detection.significant_rows}
        self.sens_outcomes: list[SensitivityOutcome] = []

        if not primary:
            reason = ("primary SIGNIFICANT_HOTSPOT_CENTERS is empty; a non-redefining "
                      "sensitivity analysis cannot establish a positive result and is "
                      "therefore not applicable")
            self.sens_outcomes = [SensitivityOutcome(n, "NOT_APPLICABLE", reason)
                                  for n in ("plddt70", "star1", "star2")]
            self._write_sensitivity_outputs()
            self._step("II.7_sensitivity", t0)
            return

        threshold = float(self.cfg.get("plddt.sensitivity_threshold"))
        keep_u = np.nan_to_num(inp.universe_plddt, nan=-1.0) >= threshold
        keep_l = keep_u[inp.labeled_positions]
        context = ctx_sensitivity("plddt70")
        plddt = run_restricted_detection(
            "plddt70", self.r_hot, inp.universe_coords, inp.universe_index, keep_u,
            inp.labeled_positions, inp.labeled_coords, inp.y.astype(np.float64), keep_l,
            self.ctx.seeds.rng(context), B=self.B, q=q, seed_context=context,
            primary_centers=primary)
        plddt.detail.setdefault("plddt_threshold", threshold)
        self.sens_outcomes.append(plddt)

        channel = load_review_channel(self.ctx, self.upstream, self.guard)
        self.review_channel = channel
        for stratum in [int(s) for s in self.cfg.get("sensitivity_analyses.review_status_strata")]:
            name = f"star{stratum}"
            if not channel["available"]:
                self.sens_outcomes.append(
                    SensitivityOutcome(name, "NOT_EVALUABLE", channel["reason"]))
                continue
            stars = channel["stars"]
            keep_l_star = np.array([stars.get(int(i), -1.0) >= stratum
                                    for i in inp.labeled_index], dtype=bool)
            context = ctx_sensitivity(name)
            outcome = run_restricted_detection(
                name, self.r_hot, inp.universe_coords, inp.universe_index,
                np.ones(inp.M, dtype=bool), inp.labeled_positions, inp.labeled_coords,
                inp.y.astype(np.float64), keep_l_star, self.ctx.seeds.rng(context),
                B=self.B, q=q, seed_context=context, primary_centers=primary)
            outcome.detail.setdefault("min_review_stars", stratum)
            self.sens_outcomes.append(outcome)

        for o in self.sens_outcomes:
            if not conclusion_flipped(o):
                continue
            code = ("STRUCTURAL_CONFIDENCE_SENSITIVE" if o.name == "plddt70"
                    else "REVIEW_STATUS_SENSITIVE")
            self.books[SENS].warnings.add(
                code,
                f"Sensitivity analysis {o.name} at the frozen r_hot = {self.r_hot:g} A "
                f"yields NO significant center, while the primary analysis yields "
                f"{len(primary)}.",
                potential_consequence=("The primary result depends on residues that the "
                                       "restriction removes."),
                recommended_action=("Report as a caveat on the primary result. The "
                                    "sensitivity analysis NEVER redefines it."),
                affected_output=f"sensitivity_{o.name}.json")
            self._escalate(
                f"Sensitivity analysis {o.name} flips the qualitative conclusion.",
                ["Report the primary result with the caveat attached (current behaviour).",
                 "Lead review of whether the primary cohort is adequate."],
                ["The primary result stands; the caveat travels with every claim.",
                 "Adopting the sensitivity result would violate the non-redefining rule."],
                "Report the caveat; never adopt the sensitivity result as primary.")

        self._write_sensitivity_outputs()
        self._step("II.7_sensitivity", t0)

    def _write_sensitivity_outputs(self) -> None:
        d = self._dir(SENS)
        by_name = {o.name: o for o in self.sens_outcomes}
        common = {
            "run_at_frozen_r_hot": getattr(self, "r_hot", None),
            "IS_NON_REDEFINING": True,
            "may_redefine_primary": False,
            "primary_centers": sorted(
                r["center_residue_index"] for r in self.detection.significant_rows),
        }
        plddt = by_name.get("plddt70")
        self._record(SENS, write_json(d / "sensitivity_plddt70.json", {
            **common, "analysis": "plddt70", "status": plddt.status,
            "reason": plddt.reason, "centers": plddt.centers, **plddt.detail,
            "note": ("pLDDT is recorded and is NEVER a primary filter (F2). This analysis "
                     "shows what the result would have been under a >= 70 pLDDT filter."),
        }))
        self._record(SENS, write_json(d / "sensitivity_review_status.json", {
            **common, "analysis": "review_status",
            "channel": getattr(self, "review_channel", {"available": False}),
            "strata": {name: {"status": by_name[name].status,
                              "reason": by_name[name].reason,
                              "centers": by_name[name].centers, **by_name[name].detail}
                       for name in by_name if name.startswith("star")},
            "note": ("ClinVar review stars are NEVER an inclusion criterion (F12). Star "
                     "metadata is read only through the separate non-redefining channel "
                     "declared in handoff_01, never through the Stage B primary path."),
        }))
        rows = [overlap_row(o) for o in self.sens_outcomes]
        rows.append({
            "analysis": "global_clustering_caveat",
            "status": "COMPLETED", "radius_A": getattr(self, "r_hot", None),
            "n_labeled_used": self.inputs.N, "n_plp_used": self.inputs.N_P,
            "n_blb_used": self.inputs.N_B,
            "n_centers_universe": self.inputs.n_center_universe,
            "n_centers_primary": len(self.detection.significant_rows),
            "n_centers_sensitivity": None, "n_overlap": None,
            "overlap_fraction_of_primary": None, "jaccard": None,
            "conclusion_preserved": None, "is_non_redefining": True,
            "reason": (f"global_clustering = {self.global_flag}; carried as a caveat on "
                       f"every downstream claim (F5). Local discovery is never terminated "
                       f"by a non-significant global result."),
        })
        self._record(SENS, write_tsv(d / "sensitivity_overlap.tsv", rows,
                                     SENSITIVITY_OVERLAP_COLUMNS))
        self._record(SENS, write_text(d / "SENSITIVITY_IS_NON_REDEFINING.txt",
                                      NON_REDEFINING_TEXT))

    # -- closing -------------------------------------------------------------
    def finish(self) -> Handoff:
        det = self.detection
        q = float(self.cfg.get("fdr.q"))
        sel = self.selection
        selected_norm = sel.result.normalized[self.r_hot]

        loo_at_r = next(r for r in self.results if r.radius_A == self.r_hot).loo
        sens_summary = {o.name: {"status": o.status, "reason": o.reason,
                                 "n_centers": len(o.centers),
                                 "overlap_fraction_of_primary":
                                     o.detail.get("overlap_fraction_of_primary"),
                                 "jaccard": o.detail.get("jaccard"),
                                 "is_non_redefining": True}
                        for o in self.sens_outcomes}
        sens_summary["global_clustering_caveat"] = {
            "status": "COMPLETED", "flag": self.global_flag, "is_non_redefining": True}

        negative = None
        if det.n_significant == 0:
            # Only reachable with a PASSED §5.4 certificate — detect() terminates
            # UNDERPOWERED otherwise, so a negative here is always a certified one.
            negative = {
                "condition": "NO_SIGNIFICANT_HOTSPOT_CENTERS",
                "detail": (f"No center of U_struct survives BH at q = {q} at "
                           f"r_hot = {self.r_hot:g} A."),
                "settled_negative": not self.resolution_limited,
                "interpretation": (
                    "Valid, complete scientific negative: the v2 §5.4 power "
                    "certificate PASSED, so the design demonstrably could have "
                    "rejected, and no 3D hotspot is detectable under the "
                    "pre-registered procedure."
                    if not self.resolution_limited else
                    "RESOLUTION-LIMITED: the permutation resolution may bound the "
                    "rejection set. This is explicitly NOT evidence of absence."),
                "power_certificate_passed": True,
                "power_certificate": self.posthoc_certificate.as_json(),
                "family_size": len(det.family_positions),
                "B": self.B, "B_recommended": self.posthoc.B_rec,
                "n_provisional_centers_at_r_hot": self.n_provisional_at_r_hot,
                "scan_final_disagreement": self.scan_final_disagreement,
                "scan_final_disagreement_note": self.scan_final_note,
                "licenses_no_method_modification": True,
            }

        self.summary = self._build_summary(selected_norm, loo_at_r, sens_summary, negative)

        payload = {
            "gene": self.ctx.gene, "code_version": self.ctx.code_version,
            "terminal_state": ("COMPLETED_NEGATIVE" if det.n_significant == 0
                               else "COMPLETED"),
            "hotspot_radius": self.r_hot,
            "primary_null": self.primary_null,
            "secondary_null": "label_permutation",
            "secondary_null_summary": det.secondary,
            "power_certificate_passed": self.certificate_passes,
            "power_certificate": self.posthoc_certificate.as_json(),
            "power_certificate_preflight": self.preflight_certificate.as_json(),
            "search_domain_source": self.domain.source,
            "fallback_radius_domain": self.domain.fallback,
            "fallback_trigger_id": self.domain.trigger_id or "NA",
            "q": q, "fdr_method": "BH", "B": self.B,
            "n_significant_centers": det.n_significant,
            "n_centers_without_variant": det.n_centers_without_variant,
            "bh_boundary_p": det.boundary_p,
            "permutation_resolution_limited": self.resolution_limited,
            "b_recommended": self.posthoc.B_rec,
            "domain_boundary_warning": bool(sel.boundary.get("DOMAIN_BOUNDARY_WARNING")),
            "global_clustering_flag": self.global_flag,
            "near_tie_flag": bool(sel.result.near_tie),
            "loo_sparse_proportion": loo_at_r.prop_sparse,
            "loo_zero_neighbour_proportion": loo_at_r.prop_zero_neighbour,
            "perm_evidence_zg_at_r_hot": next(
                r.perm.zg for r in self.results if r.radius_A == self.r_hot),
            "perm_evidence_zg_direction": ZG_DIRECTION_NOTE,
            "bh_rank1_critical_value": self.bh_margin["rank1_critical_value"],
            "bh_smallest_observed_p": self.bh_margin["smallest_observed_p"],
            "bh_margin_ratio": self.bh_margin["margin_ratio"],
            "sensitivity_summary": sens_summary,
            "derived_seeds": self.ctx.seeds.issued,
            "n_hotspot_regions": len(det.region_rows),
            "n_covered_residues": len(det.covered_rows),
            "n_classified_in_hotspots": len(det.classified_rows),
            "n_provisional_centers_at_r_hot": self.n_provisional_at_r_hot,
            "scan_final_disagreement": self.scan_final_disagreement,
            "scan_final_disagreement_note": self.scan_final_note,
            "vacuous_pareto_selection": bool(sel.vacuous_pareto),
            "admissibility_dominates_selection": bool(sel.admissibility_dominates),
            "degenerate_objectives": sel.degenerate_objectives,
            "effective_objectives_in_distance": sel.effective_objectives,
            # v2 §2 — |U_struct| and |U_center| reported distinctly wherever the
            # test-family construction is documented.
            "n_universe_U_struct": self.inputs.M,
            "n_center_universe_U_center": self.inputs.n_center_universe,
            "escalations": self.escalations,
        }

        for stage in (GLOBAL, RADIUS, FINAL, SENS):
            book = self.books[stage]
            if book.status == Status.NOT_RUN:
                book.status = Status.COMPLETED
                book.outcome = OutcomeType.COMPLETED
        if det.n_significant == 0:
            book = self.books[FINAL]
            if self.resolution_limited:
                book.status, book.outcome = Status.BLOCKED, OutcomeType.TECHNICAL_FAILURE
                book.reason = negative["interpretation"]
                book.recommended = (f"Lead decision: re-run under a NEW RUN_ID with "
                                    f"B = {self.posthoc.B_rec}. Do not report this as a "
                                    f"settled negative.")
            else:
                book.status, book.outcome = (Status.COMPLETED_NEGATIVE,
                                             OutcomeType.SCIENTIFIC_NEGATIVE)
                book.reason = negative["detail"]
                book.recommended = ("Report as a valid negative result. The chain stops "
                                    "here; no footprint is constructed.")
            book.negative = negative
            self.books[SENS].status = Status.NOT_RUN
            self.books[SENS].outcome = OutcomeType.NOT_APPLICABLE
            self.books[SENS].reason = (
                "The primary center set is empty, so a non-redefining sensitivity "
                "analysis has nothing to be sensitive about and cannot establish a "
                "positive result.")

        report = render_report(self.summary)
        self._record(FINAL, write_text(self._dir(FINAL) / "stage_b_report.md", report))

        handoff = Handoff(
            name="handoff_02", run_id=self.ctx.run_id,
            config_sha256=self.cfg.sha256,
            qc_status=self._qc_status(), manifest={}, payload=payload,
            negative_result=negative)
        self._finalize_outputs(handoff)
        return handoff

    def _qc_status(self) -> str:
        for book in self.books.values():
            if book.warnings.counts()[Severity.BLOCKING.value]:
                return "FAIL"
        major = sum(b.warnings.counts()[Severity.MAJOR.value] for b in self.books.values())
        return QC_PASS_WITH_WARNINGS if major else QC_PASS

    def _build_summary(self, selected_norm, loo_at_r, sens_summary, negative) -> dict:
        sel = self.selection
        det = self.detection
        rejected = sorted(
            ((c.key, _reason(c, sel)) for c in sel.candidates if c.key != self.r_hot),
            key=lambda t: sel.result.distances.get(t[0], float("inf")))
        return {
            "run_id": self.ctx.run_id, "gene": self.ctx.gene,
            "agent": AGENT, "code_version": self.ctx.code_version,
            "config_sha256": self.cfg.sha256, "synthetic": self.ctx.synthetic,
            "started_utc": self.ctx.started_utc, "ended_utc": _utc(),
            "wall_seconds": round(time.perf_counter() - self.t_start, 3),
            "steps": self.steps,
            "inputs": {
                "M": self.inputs.M, "N": self.inputs.N, "N_P": self.inputs.N_P,
                "N_B": self.inputs.N_B, "D_max_A": self.d_max, "r_max_A": self.r_max,
                "n_center_universe": self.inputs.n_center_universe,
                "allowlist": self.inputs.allowlist_used,
                "load_report": self.inputs.load_report,
                "upstream_handoff": self.upstream.name,
                "upstream_qc": self.upstream.qc_status,
            },
            "global": {
                "flag": self.global_flag, "alpha": self.alpha,
                "PLP": {"T_obs": self.stats["PLP"].t_obs,
                        "p_global": self.stats["PLP"].p_global},
                "BLB": {"T_obs": self.stats["BLB"].t_obs,
                        "p_global": self.stats["BLB"].p_global},
                "principal_pcf_peak": self.peaks.get("principal_peak"),
            },
            "domain": self.domain.as_json(),
            "power_certificate": {
                "preflight": self.preflight_certificate.as_json(),
                "posthoc": self.posthoc_certificate.as_json(),
                "passes": self.certificate_passes,
            },
            "primary_null": self.primary_null,
            "secondary_null": det.secondary,
            "selection": {
                "r_hot": self.r_hot,
                "normalized": selected_norm,
                "distance_to_ideal": sel.result.distances[self.r_hot],
                "pareto_members": sel.result.pareto_members,
                "vacuous_pareto": sel.vacuous_pareto,
                "vacuous_reason": sel.vacuous_reason,
                "admissibility_dominates": sel.admissibility_dominates,
                "inadmissible_fraction": sel.inadmissible_fraction,
                "inadmissible_by_rule": sel.inadmissible_by_rule,
                "effective_objectives": sel.effective_objectives,
                "n_admissible": sel.result.n_admissible,
                "n_candidates": len(sel.candidates),
                "near_tie": sel.result.near_tie,
                "near_tie_detail": sel.result.near_tie_detail,
                "tie_chain": sel.result.tie_chain_applied,
                "per_objective_argmax": sel.per_objective_argmax,
                "boundary": sel.boundary,
                "degenerate_objectives": sel.result.degenerate_objectives,
                "correlation": sel.result.correlation,
                "top_rejected": rejected[:5],
                "raw_objectives": {f"{r.radius_A:g}": r.objectives() for r in self.results},
                "inadmissible": {f"{r.radius_A:g}": r.qc_failure_reason
                                 for r in self.results if not r.admissible},
            },
            "detection": {
                **det.counts,
                "bh_boundary_p": det.boundary_p,
                "n_provisional_centers_at_r_hot": self.n_provisional_at_r_hot,
                "scan_final_disagreement": self.scan_final_disagreement,
                "scan_final_disagreement_note": self.scan_final_note,
                "coverage_at_r_hot": next(
                    r.coverage_fraction for r in self.results if r.radius_A == self.r_hot),
                "loo_prop_sparse": loo_at_r.prop_sparse,
                "loo_prop_zero_neighbour": loo_at_r.prop_zero_neighbour,
                "loo_n_tie": loo_at_r.n_tie,
                "by_constant": det.fdr.by_constant,
                "n_significant_by": det.fdr.n_reject_by,
                "bh_margin": self.bh_margin,
                "perm_evidence_zg_at_r_hot": next(
                    r.perm.zg for r in self.results if r.radius_A == self.r_hot),
            },
            "resolution": {"preflight": self.preflight.as_json(),
                           "posthoc": self.posthoc.as_json(),
                           "limited": self.resolution_limited},
            "sensitivity": sens_summary,
            "warnings": [w.as_row() for b in self.books.values() for w in b.warnings.items],
            "escalations": self.escalations,
            "negative_result": negative,
            "B": self.B, "q": float(self.cfg.get("fdr.q")),
            "kappa": float(self.cfg.get("loo_mcc.kappa")),
            "perm_evidence_zg_direction": ZG_DIRECTION_NOTE,
            "seeds": self.ctx.seeds.issued,
        }

    def _finalize_outputs(self, handoff: Handoff) -> None:
        """Manifests, statuses, provenance and handoff_02 — written for every stage."""
        barrier = assert_no_downstream_dirs_touched(self.ctx, self.guard)
        for stage in STAGES:
            self._write_stage_manifest(stage)
            book = self.books[stage]
            book.warnings.write(self.ctx.full_results / stage)
            status = StageStatus(
                stage=stage, agent_owner=AGENT, status=book.status,
                outcome_type=book.outcome, reason=book.reason or "NA",
                started_utc=book.started, ended_utc=_utc(),
                wall_seconds=time.perf_counter() - book.t0,
                upstream_handoff=self.upstream.name,
                upstream_qc_status=self.upstream.qc_status,
                n_outputs_expected=len(EXPECTED_OUTPUTS[stage]),
                n_outputs_created=sum(
                    1 for f in EXPECTED_OUTPUTS[stage]
                    if (self.ctx.full_results / stage / f).is_file()),
                warnings=book.warnings.counts(),
                recommended_action=book.recommended or "NA",
                negative_result=book.negative,
            )
            _write_stage_status(status, self.ctx.full_results / stage)

        files = [p for stage in STAGES
                 for p in (self.ctx.full_results / stage).rglob("*") if p.is_file()]
        handoff.manifest = manifest_for(files, root=self.ctx.run_root)
        write_json(self.ctx.handoff_path("handoff_02"), handoff.as_dict())

        for stage in STAGES:
            self._write_provenance(stage, barrier)

    def _write_stage_manifest(self, stage: str) -> None:
        directory = self.ctx.full_results / stage
        rows = []
        created = {str(Path(p).relative_to(directory)) for p in self.books[stage].created
                   if str(p).startswith(str(directory))}
        for rel in sorted(set(EXPECTED_OUTPUTS[stage]) | created):
            path = directory / rel
            if path.is_file():
                rows.append({"relative_path": rel, "kind": path.suffix.lstrip(".") or "txt",
                             "status": "CREATED", "reason": "NA",
                             "sha256": sha256_file(path),
                             "size_bytes": path.stat().st_size})
            else:
                rows.append({"relative_path": rel, "kind": Path(rel).suffix.lstrip(".") or "txt",
                             "status": "NOT_CREATED",
                             "reason": self.books[stage].reason or
                             "stage did not reach this output",
                             "sha256": None, "size_bytes": None})
        write_tsv(directory / "stage_manifest.tsv", rows, MANIFEST_COLUMNS)

    def _write_provenance(self, stage: str, barrier: dict) -> None:
        book = self.books[stage]
        self.ctx.record_provenance(
            stage, AGENT,
            inputs={**getattr(self.inputs, "paths", {}),
                    "upstream_handoff": self.upstream.name,
                    "upstream_manifest_sha256": self.upstream.manifest,
                    "column_allowlist": getattr(self.inputs, "allowlist_used", []),
                    "load_report": getattr(self.inputs, "load_report", {})},
            outputs={"created": sorted(str(Path(p).relative_to(self.ctx.run_root))
                                       for p in book.created)},
            parameters=self._parameters(),
            commands=[f"run_stage_b(ctx, upstream={self.upstream.name})"],
            warnings=[w.as_row() for w in book.warnings.items],
            extra={"information_barrier": barrier,
                   "steps": self.steps,
                   "escalations": self.escalations,
                   "hotspot_code_version": code_version(
                       Path(__file__).resolve().parent),
                   "spatial_code_version": code_version(
                       Path(__file__).resolve().parents[1] / "spatial"),
                   "seed_contexts": self.ctx.seeds.issued,
                   "stage_status": book.status.value})

    def _parameters(self) -> dict:
        keys = ["global_clustering", "radius_domain", "permutation", "fdr",
                "permutation_resolution_diagnostic", "loo_mcc", "radius_qc",
                "radius_selection", "boundary_diagnostic", "final_detection",
                "sensitivity_analyses", "seeding", "plddt", "biological_plausibility"]
        return {k: self.cfg.get(k) for k in keys} | {"B_used": self.B}


def _reason(candidate, sel) -> str:
    from .selection import _rejection_reason
    return _rejection_reason(candidate, sel)


#: Workflow v2 introduced ``OutcomeType.UNINFORMATIVE``; the Lead-owned
#: ``utils.status._write_not_run`` paragraph writer does not yet carry a branch for
#: it and would raise ``KeyError``. Stage B therefore writes the UNINFORMATIVE
#: paragraph itself and defers to the shared writer for every other outcome, so no
#: Lead-owned file is edited and no stage is ever left without a status file.
_UNINFORMATIVE_PARAGRAPH = (
    "This run is UNINFORMATIVE. It is NOT a scientific negative and NOT a technical "
    "failure. The per-center test could not have rejected any center regardless of "
    "the data, because the smallest p-value it can return (p_floor) exceeds the "
    "rank-1 Benjamini-Hochberg critical value c_1 = q/m. See power_certificate.json "
    "for p_floor, c_1, their ratio, N_P, N_B, m, B and which floor binds.\n\n"
    "NOTHING may be concluded about the presence or absence of hotspots in this gene "
    "from this run (Workflow v2 Appendix A). In particular the statements 'no 3D "
    "hotspot is detectable in this gene', 'this is a valid, complete scientific "
    "negative' and 'this result licenses no modification of the method' are all "
    "prohibited here."
)


def _write_stage_status(status: StageStatus, stage_dir: Path) -> Path:
    """Write ``stage_status.json`` (+ ``NOT_RUN.txt``) for any outcome type."""
    if status.outcome_type is not OutcomeType.UNINFORMATIVE:
        return status.write(stage_dir)
    stage_dir = Path(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)
    path = write_json(stage_dir / "stage_status.json", status.as_dict(),
                      schema_version=STAGE_STATUS_SCHEMA)
    (stage_dir / "NOT_RUN.txt").write_text(
        f"STAGE NOT COMPLETED: {status.stage}\n"
        f"{'=' * (21 + len(status.stage))}\n\n"
        f"Status       : {status.status.value}\n"
        f"Outcome type : {status.outcome_type.value}\n"
        f"Owner        : {status.agent_owner}\n"
        f"Detected     : {status.ended_utc or status.started_utc or 'NA'}\n\n"
        f"What happened\n-------------\n{status.reason or 'No reason recorded.'}\n\n"
        f"How to read this\n----------------\n{_UNINFORMATIVE_PARAGRAPH}\n\n"
        f"Recommended action\n------------------\n"
        f"{status.recommended_action or 'None recorded.'}\n",
        encoding="utf-8", newline="\n")
    return path


# --- entry point -------------------------------------------------------------

def run_stage_b(ctx: RunContext, *, upstream: Handoff) -> Handoff:
    """Stage B: primary spatial hotspot discovery (METHOD_SPEC II.1-II.7)."""
    stage = StageB(ctx, upstream)
    try:
        stage.validate()
        stage.global_clustering()
        stage.candidate_domain()
        # II.4 first, so the permutation-resolution diagnostic exists even when the
        # §5.4 certificate terminates the run at the pre-flight instance. The two
        # overlap (C1 is p_res > q/m) and the reader needs both, not the first one
        # that happened to fire.
        stage.preflight_resolution()
        stage.power_preflight()                    # v2 §5.4 — binding, before compute
        stage.scan()
        stage.select_radius()
        stage.detect()                             # v2 §5.4 post-hoc is inside detect()
        stage.sensitivity()
        return stage.finish()
    except UnderpoweredResult as up:
        return _finish_underpowered(stage, up)
    except NegativeResult as neg:
        return _finish_negative(stage, neg)
    except Exception as exc:                       # BlockedError / LeakageError / defects
        _finish_blocked(stage, exc)
        raise


def is_underpowered(handoff: Handoff) -> bool:
    """The signal an orchestrator reads to map Stage B onto ``RunStatus.UNDERPOWERED``.

    Stage B returns a normal ``handoff_02`` in this terminal state rather than
    raising, so the run's artefacts are complete and auditable. The handoff carries
    ``qc_status = FAIL`` (which blocks every consumer), ``negative_result = None``
    (because this is NOT a negative) and ``terminal_state = "UNDERPOWERED"``.
    """
    return handoff.payload.get("terminal_state") == Status.UNDERPOWERED.value


def _finish_underpowered(stage: StageB, up: UnderpoweredResult) -> Handoff:
    """Terminal state UNDERPOWERED (v2 §0/§5.4) — uninformative, never a negative.

    Deliberately shaped so that nothing downstream, and no reader, can mistake this
    for a result about the gene:

      * ``negative_result`` is ``None`` — there is no negative to report;
      * ``qc_status = FAIL`` — every consumer contract blocks on it, so the chain
        stops without a stage having to interpret the science;
      * ``terminal_state = UNDERPOWERED`` — the signal the orchestrator maps onto
        ``RunStatus.UNDERPOWERED`` (exit code 3);
      * the power certificate is written first, because v2 §11 requires it in every
        terminal state and it is the whole evidence for the verdict.
    """
    cert = up.certificate
    stage.log.error("UNDERPOWERED: %s", cert.headline())
    stage._write_power_certificate()
    if getattr(stage, "preflight", None) is not None:
        if stage.preflight.limited:
            stage._note_resolution_limited(stage.preflight)
        stage._write_resolution_diagnostic()

    order = list(STAGES)
    reached = order.index(up.stage)
    for i, s in enumerate(order):
        book = stage.books[s]
        if i < reached:
            book.status, book.outcome = Status.COMPLETED, OutcomeType.COMPLETED
        elif i == reached:
            book.status, book.outcome = Status.UNDERPOWERED, OutcomeType.UNINFORMATIVE
            book.reason = cert.headline()
            book.recommended = ("; ".join(cert.remedies) if cert.remedies else
                                "Lead decision — the design must change before this "
                                "gene is re-run.")
        else:
            book.status, book.outcome = Status.NOT_RUN, OutcomeType.NOT_APPLICABLE
            book.reason = (f"upstream stage {up.stage} terminated UNDERPOWERED: the "
                           f"test could not have rejected any center")
    # The FINAL book carries the BLOCKING warning even when the certificate failed
    # pre-flight, so a reader who opens only 06_FINAL_HOTSPOTS still sees it.
    if up.stage != FINAL and not stage.books[FINAL].warnings.has("TEST_CANNOT_REJECT"):
        stage.books[FINAL].warnings.add(
            "TEST_CANNOT_REJECT", cert.headline(),
            potential_consequence=UNINFORMATIVE_STATEMENT,
            recommended_action="See power_certificate.json.",
            affected_output="power_certificate.json")

    payload = {
        "gene": stage.ctx.gene, "code_version": stage.ctx.code_version,
        "terminal_state": Status.UNDERPOWERED.value,
        "outcome_type": OutcomeType.UNINFORMATIVE.value,
        "power_certificate_passed": False,
        "power_certificate": cert.as_json(),
        "power_certificate_preflight": (
            stage.preflight_certificate.as_json()
            if getattr(stage, "preflight_certificate", None) is not None else None),
        "uninformative_statement": UNINFORMATIVE_STATEMENT,
        "prohibited_inferences": [
            "No 3D hotspot is detectable in this gene.",
            "This is a valid, complete scientific negative.",
            "This result licenses no modification of the method.",
        ],
        "hotspot_radius": getattr(stage, "r_hot", None),
        "primary_null": getattr(stage, "primary_null", None),
        "search_domain_source": (stage.domain.source if hasattr(stage, "domain")
                                 else "NA"),
        "fallback_radius_domain": getattr(getattr(stage, "domain", None), "fallback",
                                          False),
        "fallback_trigger_id": (getattr(stage.domain, "trigger_id", None) or "NA"
                                if hasattr(stage, "domain") else "NA"),
        "q": float(stage.cfg.get("fdr.q")), "fdr_method": "BH",
        "B": getattr(stage, "B", None),
        "B_to_clear_c_1": cert.B_to_clear_c_1,
        "increasing_B_is_futile": cert.increasing_B_is_futile,
        # Nothing is CLAIMED, so the published count is zero and every consumer
        # contract stops here. What BH mechanically rejected is still recorded — the
        # evidence survives on disk in all_residue_center_tests.tsv, it is simply not
        # a finding. Reporting 0 while the file held rows, with no reconciling number,
        # would be the same quiet inconsistency v2 exists to remove.
        "n_significant_centers": 0,
        "n_bh_rejections_recorded_not_published": int(
            getattr(getattr(stage, "detection", None), "n_significant", 0) or 0),
        "withheld_rejection_note": (
            "p_floor > c_1 rules out a rank-1 rejection; BH can still reject several "
            "centers jointly at rank > 1. v2 §5.4 directs termination as UNDERPOWERED "
            "without a rank-1 exception, so any such rejection is recorded and "
            "escalated, never published as a finding."),
        "n_centers_without_variant": 0, "bh_boundary_p": None,
        "permutation_resolution_limited": getattr(
            getattr(stage, "preflight", None), "limited", False),
        "b_recommended": getattr(getattr(stage, "preflight", None), "B_rec", None),
        "domain_boundary_warning": False,
        "global_clustering_flag": getattr(stage, "global_flag", "NA"),
        "near_tie_flag": False, "loo_sparse_proportion": None,
        "loo_zero_neighbour_proportion": None,
        "sensitivity_summary": {
            "status": "NOT_APPLICABLE",
            "reason": ("the run is UNDERPOWERED; a non-redefining sensitivity analysis "
                       "cannot rescue a test that could not reject")},
        "n_universe_U_struct": stage.inputs.M,
        "n_center_universe_U_center": stage.inputs.n_center_universe,
        "escalations": stage.escalations,
        "derived_seeds": stage.ctx.seeds.issued,
    }
    handoff = Handoff(name="handoff_02", run_id=stage.ctx.run_id,
                      config_sha256=stage.cfg.sha256, qc_status="FAIL", manifest={},
                      payload=payload, negative_result=None)

    stage.summary = {
        "run_id": stage.ctx.run_id, "gene": stage.ctx.gene, "agent": AGENT,
        "terminal_state": Status.UNDERPOWERED.value,
        "underpowered": True,
        "power_certificate": {
            "preflight": (stage.preflight_certificate.as_json()
                          if getattr(stage, "preflight_certificate", None) is not None
                          else None),
            "posthoc": (stage.posthoc_certificate.as_json()
                        if getattr(stage, "posthoc_certificate", None) is not None
                        else None),
            "passes": False, "binding": cert.as_json(),
        },
        "primary_null": getattr(stage, "primary_null", None),
        "escalations": stage.escalations,
        "warnings": [w.as_row() for b in stage.books.values() for w in b.warnings.items],
        "domain": (stage.domain.as_json() if hasattr(stage, "domain") else {}),
        "global": {"flag": getattr(stage, "global_flag", "NA")},
        "inputs": {"M": stage.inputs.M, "N": stage.inputs.N, "N_P": stage.inputs.N_P,
                   "N_B": stage.inputs.N_B,
                   "n_center_universe": stage.inputs.n_center_universe,
                   "D_max_A": getattr(stage, "d_max", None),
                   "r_max_A": getattr(stage, "r_max", None),
                   "allowlist": stage.inputs.allowlist_used,
                   "upstream_handoff": stage.upstream.name,
                   "upstream_qc": stage.upstream.qc_status},
        "code_version": stage.ctx.code_version, "config_sha256": stage.cfg.sha256,
        "synthetic": stage.ctx.synthetic,
        "stopped_at_stage": up.stage,
        "B": getattr(stage, "B", None), "q": float(stage.cfg.get("fdr.q")),
        "kappa": float(stage.cfg.get("loo_mcc.kappa")),
        "steps": stage.steps, "seeds": stage.ctx.seeds.issued, "ended_utc": _utc(),
    }
    write_text(stage._dir(FINAL) / "stage_b_report.md", render_report(stage.summary))
    stage.books[FINAL].created.append(stage.ctx.full_results / FINAL / "stage_b_report.md")
    stage._finalize_outputs(handoff)
    return handoff


def _finish_negative(stage: StageB, neg: NegativeResult) -> Handoff:
    """A valid, complete scientific negative — never a technical failure (P4)."""
    stopped_at = neg.context.get("stage", RADIUS)
    negative = {
        "condition": neg.condition, "detail": neg.detail,
        "settled_negative": True,
        "interpretation": ("The analysis ran correctly and the pre-registered procedure "
                           "yields no result at this stage. This licenses no method "
                           "modification: B is not increased, the correction is not "
                           "loosened, the domain is not widened and no residue is dropped."),
        "licenses_no_method_modification": True,
        **{k: v for k, v in neg.context.items() if k != "stage"},
    }
    stage.log.warning("NEGATIVE RESULT [%s]: %s", neg.condition, neg.detail)
    # v2 §11 — the certificate is a first-class artefact in EVERY terminal state.
    stage._write_power_certificate()
    stage._escalate(
        f"Negative result: {neg.condition}.",
        ["Report the negative result and stop the chain (current behaviour).",
         "Lead review of cohort size or structure suitability under a NEW RUN_ID."],
        ["Stage C is not entered; no footprint is constructed.",
         "Re-running with altered parameters to obtain a positive result is prohibited."],
        "Report the negative result as final for this RUN_ID.")

    order = list(STAGES)
    reached = order.index(stopped_at)
    for i, s in enumerate(order):
        book = stage.books[s]
        if i < reached:
            book.status, book.outcome = Status.COMPLETED, OutcomeType.COMPLETED
        elif i == reached:
            book.status = Status.COMPLETED_NEGATIVE
            book.outcome = OutcomeType.SCIENTIFIC_NEGATIVE
            book.reason = f"{neg.condition}: {neg.detail}"
            book.negative = negative
            book.recommended = "Report as a valid negative result; the chain stops here."
        else:
            book.status, book.outcome = Status.NOT_RUN, OutcomeType.NOT_APPLICABLE
            book.reason = f"upstream stage {stopped_at} ended in {neg.condition}"

    payload = {
        "gene": stage.ctx.gene, "code_version": stage.ctx.code_version,
        "terminal_state": Status.COMPLETED_NEGATIVE.value,
        "power_certificate_passed": (
            bool(stage.preflight_certificate.passes)
            if getattr(stage, "preflight_certificate", None) is not None else None),
        "power_certificate": (
            stage.preflight_certificate.as_json()
            if getattr(stage, "preflight_certificate", None) is not None else None),
        "primary_null": getattr(stage, "primary_null", None),
        "hotspot_radius": None,
        "search_domain_source": getattr(stage, "domain", None).source
        if hasattr(stage, "domain") else "NA",
        "fallback_radius_domain": getattr(getattr(stage, "domain", None), "fallback", False),
        "fallback_trigger_id": (getattr(stage.domain, "trigger_id", None) or "NA"
                                if hasattr(stage, "domain") else "NA"),
        "q": float(stage.cfg.get("fdr.q")), "fdr_method": "BH",
        "B": getattr(stage, "B", None), "n_significant_centers": 0,
        "n_centers_without_variant": 0, "bh_boundary_p": None,
        "permutation_resolution_limited": getattr(
            getattr(stage, "preflight", None), "limited", False),
        "b_recommended": getattr(getattr(stage, "preflight", None), "B_rec", None),
        "domain_boundary_warning": False,
        "global_clustering_flag": getattr(stage, "global_flag", "NA"),
        "near_tie_flag": False, "loo_sparse_proportion": None,
        "loo_zero_neighbour_proportion": None,
        "sensitivity_summary": {"status": "NOT_APPLICABLE",
                                "reason": "negative result upstream of the frozen r_hot"},
        "n_universe_U_struct": stage.inputs.M,
        "n_center_universe_U_center": stage.inputs.n_center_universe,
        "escalations": stage.escalations,
        "derived_seeds": stage.ctx.seeds.issued,
    }
    handoff = Handoff(name="handoff_02", run_id=stage.ctx.run_id,
                      config_sha256=stage.cfg.sha256, qc_status=QC_PASS, manifest={},
                      payload=payload, negative_result=negative)

    preflight_cert = getattr(stage, "preflight_certificate", None)
    stage.summary = {"run_id": stage.ctx.run_id, "gene": stage.ctx.gene, "agent": AGENT,
                     "negative_result": negative, "escalations": stage.escalations,
                     "primary_null": getattr(stage, "primary_null", None),
                     "power_certificate": (
                         {"preflight": preflight_cert.as_json(),
                          "passes": bool(preflight_cert.passes)}
                         if preflight_cert is not None else {}),
                     "warnings": [w.as_row() for b in stage.books.values()
                                  for w in b.warnings.items],
                     "domain": (stage.domain.as_json() if hasattr(stage, "domain") else {}),
                     "global": {"flag": getattr(stage, "global_flag", "NA")},
                     "inputs": {"M": stage.inputs.M, "N": stage.inputs.N,
                                "N_P": stage.inputs.N_P, "N_B": stage.inputs.N_B,
                                "n_center_universe": stage.inputs.n_center_universe},
                     "B": getattr(stage, "B", None),
                     "q": float(stage.cfg.get("fdr.q")),
                     "steps": stage.steps, "seeds": stage.ctx.seeds.issued,
                     "ended_utc": _utc()}
    write_text(stage._dir(FINAL) / "stage_b_report.md", render_report(stage.summary))
    stage.books[FINAL].created.append(stage.ctx.full_results / FINAL / "stage_b_report.md")
    stage._finalize_outputs(handoff)
    return handoff


def _finish_blocked(stage: StageB, exc: Exception) -> None:
    """Technical failure: record it everywhere, then let the exception propagate."""
    stage.log.error("BLOCKED/FAILED: %s: %s", type(exc).__name__, exc)
    try:
        # v2 §11 — even a blocked run states, on disk, what the certificate was or
        # that it was never reached. Silence is not an acceptable third state.
        stage._write_power_certificate()
    except Exception:                                           # pragma: no cover
        pass
    for s in STAGES:
        book = stage.books[s]
        if book.status in (Status.NOT_RUN,):
            book.status = Status.BLOCKED
            book.outcome = OutcomeType.TECHNICAL_FAILURE
            book.reason = f"{type(exc).__name__}: {exc}"
            book.recommended = ("Resolve the blocking condition and re-run. Nothing can "
                                "be concluded about the biology from this stage.")
        try:
            book.warnings.write(stage.ctx.full_results / s)
            StageStatus(stage=s, agent_owner=AGENT, status=book.status,
                        outcome_type=book.outcome, reason=book.reason,
                        started_utc=book.started, ended_utc=_utc(),
                        wall_seconds=time.perf_counter() - book.t0,
                        upstream_handoff=stage.upstream.name,
                        upstream_qc_status=stage.upstream.qc_status,
                        n_outputs_expected=len(EXPECTED_OUTPUTS[s]),
                        n_outputs_created=0,
                        warnings=book.warnings.counts(),
                        recommended_action=book.recommended,
                        ).write(stage.ctx.full_results / s)
        except Exception:                                       # pragma: no cover
            pass
