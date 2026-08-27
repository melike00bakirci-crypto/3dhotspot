"""Phase D — geometric footprint robustness (METHOD_SPEC II.11).

One question is being answered: **is ``FP_original`` excessively dependent on one
or a small subset of significant hotspot centers, or is the same spatial footprint
preserved when centers are removed?**

The product is a continuous robustness profile. Low preservation is a complete,
publishable result and is reported as limited robustness, never repaired. No
categorical ROBUST / NOT_ROBUST verdict is emitted as a primary output, and every
statement produced here is *geometric footprint robustness* — it is not validation
of hotspot discovery.
"""
from __future__ import annotations

import gzip
import io
import tarfile
import time
from pathlib import Path

import numpy as np

from ..footprint.api import (FOOTPRINT_API_VERSION, build_footprint,
                             footprint_code_version, reconstruct_at)
from ..footprint.freeze import (load_freeze, verify_freeze, verify_upstream_unchanged)
from ..footprint.staging import stage_files, write_stage_manifest
from ..orchestration.contracts import (QC_PASS, QC_PASS_WITH_WARNINGS, Handoff,
                                       assert_handoff_shape)
from ..utils.errors import BlockedError, EscalationRequired
from ..utils.hashing import manifest_for, sha256_text
from ..utils.io import write_json, write_tsv
from ..utils.runctx import RunContext
from ..utils.seeds import ctx_bootstrap, ctx_iteration
from ..utils.status import (OutcomeType, Severity, StageStatus, Status,
                            WarningCollector, utc_now)
from . import figures as fig_mod
from . import report as report_mod
from .aggregate import (AGGREGATE_COLUMNS, distribution_summary,
                        preservation_frequency, summarize)
from .classify import CLASSIFICATION_FIELDS
from .compare import GEOMETRY_FIELDS
from .design import build_design, k_max
from .params import RobustnessParams
from .parallel import run_iterations as run_iterations_parallel
from .perturb import FAILURE_COLUMNS
from .sensitivity import (CENTER_SENSITIVITY_COLUMNS, RECURRENCE_DEFINITION,
                          center_influence)

AGENT = "footprint-robustness"
STAGE = "09_ROBUSTNESS"

PERTURBATION_BASE_COLUMNS = [
    "iteration_index", "k", "removed_center_residues", "n_removed", "n_remaining",
    "enumeration_mode", "seed_context", "status", "outcome_type", "failure_reason",
    "robustness_r_fp",
]
# Per-iteration wall time is recorded in provenance (agent §13), NOT here: a
# timing column would make an otherwise byte-reproducible decision table differ
# between identical runs, and determinism is asserted by re-run.
#
# DECISION-STAGE-D-FIXED-RFP-0001: every iteration reconstructs at the frozen
# Stage C r_fp (docs/decisions/DECISION-STAGE-D-FIXED-RFP-0001.md), so there is
# exactly one geometry block per iteration, not a primary/fixed pair.
PERTURBATION_COLUMNS = (
    PERTURBATION_BASE_COLUMNS
    + GEOMETRY_FIELDS
    + CLASSIFICATION_FIELDS
    + [f"{c}_hs" for c in CLASSIFICATION_FIELDS]
)

PERTURBATION_UNIVERSE_COLUMNS = [
    "center_residue_index", "region_id", "carries_clinvar_variant",
    "x_ca", "y_ca", "z_ca", "in_final_footprint",
]
DESIGN_COLUMNS = ["k", "total_subsets_T_k", "n_evaluated", "mode",
                  "sampling_fraction", "seed_context", "seed"]

EXPECTED_09 = [
    ("baseline_reproduction.json", "mandatory byte-identical baseline reproduction"),
    ("perturbation_universe.tsv", "S with region membership and carries_clinvar_variant"),
    ("perturbation_design.tsv", "per-level design, recorded before iteration 1"),
    ("subset_design.json", "run-level design scalars (copied from the TSV)"),
    ("perturbation_results.tsv", "every iteration: the fixed-r_fp geometry block "
                                 "and both classification families"),
    ("iteration_failures.tsv", "failed iterations with status codes"),
    ("center_sensitivity.tsv", "per-center influence ranking"),
    ("robustness_summary.tsv", "aggregates with n, mean, SD, CI, CV, median, IQR"),
    ("robustness_profile.json", "the continuous robustness profile"),
    ("perturbation_footprints.tar.gz", "per-iteration residue lists, bundled"),
    ("handoff_04.json", "handoff to biological-annotation"),
    ("stage_d_report.md", "narrative Phase D report"),
    ("stage_status.json", "stage status (always present)"),
    (f"warnings_{STAGE}.tsv", "stage warnings"),
    ("stage_manifest.tsv", "this manifest"),
    ("figures/F9_robustness_similarity_distributions.png", "F9 (300 dpi)"),
    ("figures/F9_robustness_similarity_distributions.svg", "F9 (vector)"),
    ("figures/F10_center_influence_ranking.png", "F10 (300 dpi)"),
    ("figures/F10_center_influence_ranking.svg", "F10 (vector)"),
]

FORBIDDEN_IMPORTS = ("hotspot3d.hotspot_stats", "hotspot3d.hotspot",
                     "hotspot3d.spatial")

#: v2 §8 [NEW] investigated and found NOT APPLICABLE under this run's frozen
#: design — see hotspot3d.robustness.perturb's module docstring for the full
#: reasoning. Always present in robustness_profile.json (and, through it,
#: handoff_04.json), unconditionally, so a reader never has to infer its absence.
V2_SECTION8_UNDERPOWERED_APPLICABLE = False
V2_SECTION8_UNDERPOWERED_NOTE = (
    "v2 Sec8's per-iteration UNDERPOWERED clause ('re-running the analysis inside "
    "a validation iteration must re-run Sec5.4... recorded as UNDERPOWERED') is "
    "NOT APPLICABLE under this run's frozen design: "
    "robustness.rerun_hotspot_pipeline = false (FROZEN, F11/A7) means Phase D "
    "never re-runs the per-center significance test, so Sec5.4's power "
    "certificate is never computed inside an iteration and there is nothing an "
    "UNDERPOWERED verdict could report. See hotspot3d.robustness.perturb's module "
    "docstring for the full investigation.")


# --------------------------------------------------------------------------

def run_phase_d(ctx: RunContext, *, phase_c, handoff_03: Handoff) -> Handoff:
    started = utc_now()
    t0 = time.perf_counter()
    stage_dir = ctx.stage_dir(STAGE)
    warn = WarningCollector(STAGE, AGENT)

    up = phase_c.inputs
    fp_params = phase_c.params
    params = RobustnessParams.from_config(ctx.config)

    # -- §6.6/§9.20 structural assertions ---------------------------------
    assert_no_hotspot_statistics_import()
    assert_no_radius_optimization_in_perturb()
    cohort_fingerprint = _cohort_fingerprint(up.cohort_plp, up.cohort_blb)

    # -- §6.8/§6.9 freeze and code-version integrity ----------------------
    freeze_payload = load_freeze(ctx)
    verify_freeze(ctx, freeze_payload, when="Phase D entry")
    installed_version = footprint_code_version()
    if freeze_payload.get("footprint_code_version") != installed_version:
        raise BlockedError(
            f"BLOCKED — footprint code_version drift: the freeze records "
            f"{str(freeze_payload.get('footprint_code_version'))[:12]} but the "
            f"installed package hashes to {installed_version[:12]}. The baseline "
            f"would not be comparable."
        )

    # -- §11.1 consumer-side verification of handoff_03 --------------------
    _verify_handoff_03_claims(phase_c, up, handoff_03)

    # -- §6.10 mandatory baseline reproduction ----------------------------
    baseline = _baseline_reproduction(ctx, stage_dir, phase_c, up, fp_params)

    # -- §6.11 evaluability -----------------------------------------------
    n_s = len(up.centers)
    if n_s < params.min_n_s:
        return _not_evaluable(ctx, stage_dir, warn, phase_c, handoff_03, params,
                              n_s, baseline, started, t0)

    _write_universe(stage_dir, phase_c, up)

    # -- §7.9 design, recorded IN FULL before iteration 1 ------------------
    design = build_design(n_s, params.n_cap, ctx.seeds, up.centers.ids)
    write_tsv(stage_dir / "perturbation_design.tsv",
              [level.as_row() for level in design.levels], DESIGN_COLUMNS)
    # DECISION-STAGE-D-FIXED-RFP-0001: an iteration now performs one
    # reconstruct_at() call, not build_footprint()'s full sweep, so the
    # per-iteration cost is probed with the SAME operation an iteration
    # actually runs -- reconstructing the full, unperturbed center set at the
    # frozen r_fp -- rather than reusing the baseline-reproduction sweep time,
    # which is no longer representative of iteration cost.
    t_probe = time.perf_counter()
    reconstruct_at(up.centers, up.universe, fp_params, phase_c.r_fp)
    seconds_each = time.perf_counter() - t_probe
    estimate = seconds_each * design.total_evaluated
    write_json(stage_dir / "subset_design.json", {
        **design.as_dict(),
        "planned_iterations": design.total_evaluated,
        "estimated_wall_seconds": estimate,
        "estimated_seconds_per_iteration": seconds_each,
        "levels_note": ("per-level rows live in perturbation_design.tsv; this file "
                        "holds only run-level scalars copied from it (P1)"),
        "master_seed": int(ctx.config.get("seeding.MASTER_SEED")),
        "derived_seeds": {level.seed_context: level.seed
                          for level in design.levels if level.seed is not None},
    })
    if design.K_MAX == 1 and design.total_evaluated <= 2:
        warn.add("ROBUSTNESS_EVIDENCE_MINIMAL",
                 f"n_S = {n_s}: K_MAX = 1 and only {design.total_evaluated} "
                 f"iteration(s) are possible. The profile is valid but rests on "
                 f"that many observations and is reported honestly as minimal.",
                 potential_consequence="the robustness estimate has very few "
                                       "degrees of freedom",
                 recommended_action="read the profile as minimal evidence, not as "
                                    "a precise estimate",
                 affected_output="robustness_profile.json",
                 status="ACCEPTED_BY_DESIGN")
    if params.max_wall_seconds is not None and estimate > params.max_wall_seconds:
        raise EscalationRequired(
            ambiguity=(f"the recorded design plans {design.total_evaluated:,} "
                       f"iterations (subset space {design.total_subset_space:,}, "
                       f"N_CAP {params.n_cap:,}); estimated wall time "
                       f"{estimate / 3600:.1f} h exceeds the configured budget of "
                       f"{params.max_wall_seconds / 3600:.1f} h"),
            options=["Lead lowers N_CAP and the design is rebuilt under a new RUN_ID",
                     "Lead raises the wall-time budget and this run proceeds"],
            consequences=["a smaller N_CAP samples more sparsely at high k",
                          "a longer run occupies the budget"],
            recommendation="the Lead sets N_CAP. The design is never silently "
                           "truncated and the budget is never quietly exceeded.",
        )

    # -- §7.10 iterations ----------------------------------------------------
    # Independent by construction (II.12): distributed across a bounded,
    # niced worker pool (execution.n_jobs / execution.worker_nice) that never
    # auto-scales to the node's CPU count. ProcessPoolExecutor.map() returns
    # results in TASK order regardless of completion order or worker count,
    # so the result list is bit-identical to a serial run — see
    # hotspot3d.robustness.parallel's module docstring.
    original_state = phase_c.solution.selected
    original_residues = frozenset(phase_c.solution.covered_residue_ids(up.universe))
    mode_by_k = {level.k: level.mode for level in design.levels}

    tasks = [(index, removed, mode_by_k[len(removed)], ctx_iteration(index))
             for index, removed in enumerate(design.subsets)]
    results = run_iterations_parallel(
        tasks, n_workers=params.n_workers, worker_nice=params.worker_nice,
        centers=up.centers, universe=up.universe, params=fp_params,
        original_state=original_state, original_residues=original_residues,
        r_fp_original=phase_c.r_fp, plp=up.cohort_plp, blb=up.cohort_blb,
        r_hot_value=up.r_hot.value_for_coverage_query,
        undefined_mcc=params.undefined_mcc)

    planned = design.total_evaluated
    completed = sum(1 for r in results if not r.failed)
    failed = sum(1 for r in results if r.failed)
    if planned != completed + failed:
        raise BlockedError(
            f"iteration accounting failure: planned {planned} != completed "
            f"{completed} + failed {failed}. No iteration is ever removed from a "
            f"denominator."
        )

    # DECISION-STAGE-D-FIXED-RFP-0001's invariant, checked at the wiring level
    # (run_iteration already asserts it per-call): every completed iteration
    # reconstructed at exactly the frozen Stage C r_fp, never a re-selected one.
    off_invariant = [r.index for r in results
                     if not r.failed and r.robustness_r_fp != phase_c.r_fp]
    if off_invariant:
        raise BlockedError(
            f"BLOCKED — DECISION-STAGE-D-FIXED-RFP-0001 invariant violated: "
            f"iterations {off_invariant[:10]} reconstructed at a radius other "
            f"than the frozen r_fp={phase_c.r_fp}."
        )

    # -- outputs -----------------------------------------------------------
    write_tsv(stage_dir / "perturbation_results.tsv",
              [_ordered_row(r.row) for r in results], PERTURBATION_COLUMNS)
    write_tsv(stage_dir / "iteration_failures.tsv", [
        {"iteration_index": r.index, "k": r.k,
         "removed_center_residues": list(r.removed), "status": r.status,
         "outcome_type": r.outcome_type, "failure_reason": r.failure_reason}
        for r in results if r.failed], FAILURE_COLUMNS)

    sensitivity = center_influence(results, up.centers.ids, up.center_regions,
                                   up.center_carries_variant)
    write_tsv(stage_dir / "center_sensitivity.tsv", sensitivity,
              CENTER_SENSITIVITY_COLUMNS)

    aggregates = _aggregate_all(ctx, results, planned, params)
    write_tsv(stage_dir / "robustness_summary.tsv", aggregates, AGGREGATE_COLUMNS)

    profile = _build_profile(results, aggregates, sensitivity, design, params,
                             phase_c, planned, completed, failed)
    write_json(stage_dir / "robustness_profile.json", profile)

    _bundle_footprints(stage_dir / "perturbation_footprints.tar.gz", results)

    jaccards = [r.row.get("jaccard_residues") for r in results]
    dices = [r.row.get("dice_residues") for r in results]
    fig_mod.figure_f9(stage_dir / "figures", jaccards, dices,
                      (params.preservation_primary, params.preservation_secondary))
    fig_mod.figure_f10(stage_dir / "figures", sensitivity)

    _emit_profile_warnings(warn, profile, results, params)

    # -- §9.20 end-of-stage integrity -------------------------------------
    verify_freeze(ctx, freeze_payload, when="Phase D stage end")
    verify_upstream_unchanged(up.hashes, up.paths, when="Phase D stage end")
    if _cohort_fingerprint(up.cohort_plp, up.cohort_blb) != cohort_fingerprint:
        raise BlockedError(
            "BLOCKED — the classified cohort changed during Phase D. No ClinVar "
            "record or class label may be altered at any point in this stage."
        )

    report_mod.write_stage_d_report(
        stage_dir / "stage_d_report.md", profile, design, params, results,
        sensitivity, phase_c, warn.items, baseline, ctx.run_id, ctx.gene)

    handoff = _build_handoff_04(ctx, stage_dir, profile, design, params, phase_c,
                                results, warn)
    _write_status(stage_dir, warn, planned, completed, failed, started, t0,
                  handoff_03)
    warn.write(stage_dir)
    write_stage_manifest(stage_dir, EXPECTED_09, extra_created=stage_files(stage_dir))
    _record_provenance(ctx, up, fp_params, params, design, phase_c, results, warn)
    return handoff


# --------------------------------------------------------------------------
# structural assertions
# --------------------------------------------------------------------------

def assert_no_hotspot_statistics_import() -> None:
    """§9.20 — the robustness procedure must not import the hotspot-statistics code.

    Checked against the *source* of the two packages rather than ``sys.modules``,
    which would false-positive whenever the orchestrator has already run Stage B in
    the same process. This is the mechanical form of "no hotspot statistic is
    recomputed here".
    """
    roots = [Path(__file__).resolve().parent,
             Path(__file__).resolve().parents[1] / "footprint"]
    offenders = []
    for root in roots:
        for path in sorted(root.rglob("*.py")):
            text = path.read_text(encoding="utf-8")
            for token in FORBIDDEN_IMPORTS:
                package = token.split(".")[-1]
                if f"import {token}" in text or f"from {token}" in text \
                        or f"from ..{package} import" in text \
                        or f"from .{package} import" in text:
                    offenders.append(f"{path.name}: {token}")
    if offenders:
        raise BlockedError(
            f"BLOCKED — the footprint/robustness packages import hotspot-statistics "
            f"code: {offenders}. Phase D is geometric footprint robustness; it never "
            f"re-runs Ripley's K, pair correlation, the radius scan, r_hot selection, "
            f"permutation testing or FDR."
        )


def assert_no_radius_optimization_in_perturb() -> None:
    """DECISION-STAGE-D-FIXED-RFP-0001 — a Phase D iteration never re-selects r_fp.

    Checked against the *parsed AST* of ``perturb.py`` rather than by tracing a
    single run (catches the defect on every run, not only when a particular
    execution happens to exercise it) and rather than a raw text search (which
    would false-positive on this module's own docstring naming
    ``build_footprint`` in prose). ``build_footprint`` (the domain-derivation ->
    multi-scale-sweep -> QC -> Pareto -> distance-to-ideal selection entry point)
    must be neither imported nor called anywhere in ``perturb.py``'s actual code;
    only ``reconstruct_at`` (no selection) may be.
    """
    import ast

    path = Path(__file__).resolve().parent / "perturb.py"
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

    imported_names = set()
    called_names = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            imported_names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            called_names.add(node.func.id)

    if "build_footprint" in imported_names or "build_footprint" in called_names:
        raise BlockedError(
            f"BLOCKED — {path.name} imports or calls build_footprint: a Phase D "
            f"iteration must reconstruct only at the frozen Stage C r_fp via "
            f"reconstruct_at, never re-select r_fp via the full sweep/Pareto "
            f"selection path (DECISION-STAGE-D-FIXED-RFP-0001)."
        )
    if "reconstruct_at" not in imported_names:
        raise BlockedError(
            f"BLOCKED — {path.name} no longer imports reconstruct_at; the "
            f"fixed-radius reconstruction this decision requires is missing."
        )


def _verify_handoff_03_claims(phase_c, up, handoff_03: Handoff) -> None:
    """Trust is never transitive: the consumer re-checks what it depends on.

    Phase C asserted these; Phase D asserts them again from the artifacts, because
    a robustness comparison against a footprint that does not contain its own
    centers, or against an off-grid radius, would be meaningless.
    """
    footprint = set(phase_c.solution.covered_residue_ids(up.universe))
    missing = sorted(set(up.centers.ids) - footprint)
    if missing:
        raise BlockedError(
            f"BLOCKED — FP_original does not contain significant centers {missing[:10]}. "
            f"A final footprint that excludes a significant center is a computational "
            f"fault and is escalated to the Lead, never worked around."
        )
    domain = phase_c.solution.domain
    if phase_c.r_fp not in domain.grid:
        raise BlockedError(
            f"BLOCKED — r_fp = {phase_c.r_fp} is not on the II.9 grid "
            f"[{domain.rho_min}, {domain.rho_max}] step {domain.step}."
        )
    declared = handoff_03.payload.get("footprint_radius")
    if declared is None or abs(float(declared) - phase_c.r_fp) > 1e-9:
        raise BlockedError(
            f"BLOCKED — handoff_03 declares footprint_radius={declared} but the "
            f"frozen solution selected {phase_c.r_fp}."
        )
    for stage in ("10_ANNOTATION",):
        assert stage not in {p for path in up.paths.values() for p in path.parts}, (
            "10_ANNOTATION must never be opened in either phase")


def _cohort_fingerprint(plp: frozenset, blb: frozenset) -> str:
    return sha256_text("PLP:" + ",".join(map(str, sorted(plp)))
                       + "|BLB:" + ",".join(map(str, sorted(blb))))


# --------------------------------------------------------------------------
# baseline reproduction
# --------------------------------------------------------------------------

def _baseline_reproduction(ctx: RunContext, stage_dir: Path, phase_c, up,
                           fp_params) -> dict:
    """Re-run the frozen API on the FULL center set; require byte-identical output."""
    t0 = time.perf_counter()
    reproduced = build_footprint(up.centers, up.universe, fp_params)
    elapsed = time.perf_counter() - t0

    original_digest = phase_c.solution.digest(up.universe)
    reproduced_digest = reproduced.digest(up.universe)
    mismatches = sorted(k for k in original_digest
                        if original_digest[k] != reproduced_digest[k])

    payload = {
        "performed": True,
        "footprint_code_version": footprint_code_version(),
        "footprint_api_version": FOOTPRINT_API_VERSION,
        "center_set_sha256": up.centers.digest(),
        "original": original_digest,
        "reproduced": reproduced_digest,
        "byte_identical": not mismatches,
        "mismatched_fields": mismatches,
        "wall_seconds": elapsed,
        "statement": ("The frozen footprint API reproduced FP_original from the "
                      "original inputs before any iteration ran. Without this, a "
                      "robustness comparison would measure code drift."),
    }
    write_json(stage_dir / "baseline_reproduction.json", payload)
    if mismatches:
        raise BlockedError(
            f"BLOCKED — baseline reproduction failed on {mismatches}. The frozen "
            f"footprint API did not reproduce FP_original byte-identically from the "
            f"original inputs; this is a reproducibility fault and Phase D cannot "
            f"produce a meaningful comparison."
        )
    return payload


# --------------------------------------------------------------------------
# outputs
# --------------------------------------------------------------------------

def _write_universe(stage_dir: Path, phase_c, up) -> None:
    footprint_residues = set(phase_c.solution.covered_residue_ids(up.universe))
    rows = []
    for i, cid in enumerate(up.centers.ids):
        x, y, z = up.centers.coords[i]
        rows.append({
            "center_residue_index": int(cid),
            "region_id": up.center_regions.get(int(cid), "NA"),
            "carries_clinvar_variant": bool(
                up.center_carries_variant.get(int(cid), False)),
            "x_ca": float(x), "y_ca": float(y), "z_ca": float(z),
            "in_final_footprint": int(cid) in footprint_residues,
        })
    write_tsv(stage_dir / "perturbation_universe.tsv", rows,
              PERTURBATION_UNIVERSE_COLUMNS)


def _ordered_row(row: dict) -> dict:
    return {c: row.get(c) for c in PERTURBATION_COLUMNS}


def _aggregate_all(ctx: RunContext, results, planned: int,
                   params: RobustnessParams) -> list[dict]:
    """Every ¶47 and ¶49 metric is aggregated — including the unflattering ones."""
    rows = []
    families = [
        ("fixed_radius", GEOMETRY_FIELDS, ""),
        ("classification_footprint", CLASSIFICATION_FIELDS, ""),
        ("classification_hotspot_coverage", CLASSIFICATION_FIELDS, "_hs"),
    ]
    for reconstruction, fields, suffix in families:
        for field_name in fields:
            column = f"{field_name}{suffix}"
            values = [r.row.get(column) for r in results]
            context = ctx_bootstrap(f"{reconstruction}|{field_name}")
            rows.append(summarize(field_name, reconstruction, values, planned,
                                  ctx.seeds.rng(context),
                                  params.bootstrap_resamples,
                                  params.bootstrap_confidence, context))
    # DECISION-STAGE-D-FIXED-RFP-0001's invariant, in the same aggregate table a
    # reader already checks: degenerate_constant / SD=0 / mean=r_fp IS the proof.
    context = ctx_bootstrap("diagnostic|robustness_r_fp")
    rows.append(summarize("robustness_r_fp", "diagnostic",
                          [r.robustness_r_fp for r in results], planned,
                          ctx.seeds.rng(context), params.bootstrap_resamples,
                          params.bootstrap_confidence, context))
    return rows


def _build_profile(results, aggregates, sensitivity, design, params,
                   phase_c, planned: int, completed: int, failed: int) -> dict:
    jaccards = [r.row.get("jaccard_residues") for r in results]
    dices = [r.row.get("dice_residues") for r in results]

    p50, hits50, n50 = preservation_frequency(jaccards, params.preservation_primary)
    p70, hits70, n70 = preservation_frequency(jaccards, params.preservation_secondary)

    robustness_r_fp = [r.robustness_r_fp for r in results
                       if r.robustness_r_fp is not None]
    n_at_original = sum(1 for v in robustness_r_fp if v == phase_c.r_fp)

    return {
        "result_type": "geometric_footprint_robustness_profile",
        "primary_result": ("the continuous robustness profile below; the "
                           "preservation frequencies exist solely to compute the "
                           "para-51 preservation proportion"),
        "categorical_verdict_emitted": False,
        "jaccard_residues": distribution_summary(
            [0.0 if v is None else v for v in jaccards]),
        "dice_residues": distribution_summary(
            [0.0 if v is None else v for v in dices]),
        # Top-level aliases read verbatim by REVIEW_PACK and gene_summary. They are
        # the same numbers as the "preservation" block below, surfaced under the
        # contract's names rather than recomputed.
        "preservation_freq_050": p50,
        "preservation_freq_070": p70,
        "preservation": {
            "threshold_primary": params.preservation_primary,
            "threshold_secondary": params.preservation_secondary,
            "P_jaccard_ge_primary": p50,
            "n_iterations_ge_primary": hits50,
            "P_jaccard_ge_secondary": p70,
            "n_iterations_ge_secondary": hits70,
            "denominator": n50,
            "note": ("the denominator is EVERY evaluated iteration; failed "
                     "iterations contribute J = 0 and are not removed"),
        },
        "mcc_footprint_family": distribution_summary(
            [r.row.get("mcc") for r in results]),
        "mcc_hotspot_coverage_family": distribution_summary(
            [r.row.get("mcc_hs") for r in results]),
        "component_preservation_rate": _rate(
            [r.row.get("component_preserved") for r in results], planned),
        "geometric_continuity_rate": _rate(
            [r.row.get("continuity") for r in results], planned),
        "centroid_shift_A": distribution_summary(
            [r.row.get("centroid_shift_A") for r in results]),
        "hausdorff95_A": distribution_summary(
            [r.row.get("hausdorff95_A") for r in results]),
        "volume_ratio": distribution_summary(
            [r.row.get("volume_ratio") for r in results]),
        "center_influence_ranking": sensitivity,
        "center_recurrence_definition": RECURRENCE_DEFINITION,
        "robustness_r_fp_invariant": {
            "r_fp_original": phase_c.r_fp,
            "n_iterations_completed": len(robustness_r_fp),
            "n_equal_to_original": n_at_original,
            "invariant_holds": bool(n_at_original == len(robustness_r_fp)),
            "distribution": distribution_summary(robustness_r_fp),
            "interpretation": ("DECISION-STAGE-D-FIXED-RFP-0001: every completed "
                               "iteration reconstructs at exactly the frozen "
                               "Stage C r_fp; radius selection is never re-run "
                               "inside Phase D. Asserted at run_iteration's "
                               "origin and re-checked here over all results."),
        },
        "iterations": {"planned": planned, "completed": completed, "failed": failed},
        "evidence_flags": _evidence_flags(design, planned, failed),
        "design": design.as_dict(),
        "clinical_dataset_unmodified": True,
        "perturbation_acted_on": "geometric center positions (SIGNIFICANT_HOTSPOT_CENTERS)",
        "geometric_footprint_robustness_only": True,
        "not_validation_of_hotspot_discovery": (
            "This result is geometric footprint robustness only. It is not "
            "validation of hotspot discovery and must never be presented as "
            "hotspot leave-one-out validation."),
        "v2_section8_underpowered_iteration_applicable":
            V2_SECTION8_UNDERPOWERED_APPLICABLE,
        "v2_section8_underpowered_iteration_note": V2_SECTION8_UNDERPOWERED_NOTE,
    }


def _evidence_flags(design, planned: int, failed: int) -> dict:
    """How much evidence the profile actually rests on.

    Both conditions below change how the numbers should be read, and neither has a
    code in the FROZEN warning enum, so they are reported structurally here and in
    the QC section of ``stage_d_report.md``. Reusing an unrelated Stage B code
    would send a reader to the wrong stage; inventing one is not mine to do.
    """
    fraction = float(failed / planned) if planned else 0.0
    minimal = bool(design.K_MAX == 1 and design.total_evaluated <= 2)
    return {
        "minimal_design": minimal,
        "minimal_design_detail": (
            f"n_S = {design.n_S}, K_MAX = {design.K_MAX}, "
            f"{design.total_evaluated} iteration(s) evaluated; the profile rests on "
            f"that many observations."
            + (" This is the minimum possible design and the profile is reported "
               "honestly as minimal." if minimal else "")),
        "n_iterations_failed": int(failed),
        "failed_fraction": fraction,
        "iterations_failed_en_masse": bool(fraction >= 0.50),
        "failure_detail": (
            f"{failed} of {planned} iterations failed ({fraction:.1%}); all remain "
            f"in every denominator and contribute J = 0. See "
            f"iteration_failures.tsv."
            + (" Failures dominate: the profile reflects reconstruction failure "
               "more than measured geometric change and must be read as such."
               if fraction >= 0.50 else "")),
    }


def _rate(values, planned: int) -> dict:
    hits = sum(1 for v in values if v is True)
    return {"n_true": hits, "denominator": planned,
            "rate": float(hits / planned) if planned else None}


def _bundle_footprints(path: Path, results) -> None:
    """One archive, not fifty thousand directories (inode control, IX.3).

    Written with fixed mtime/uid/gid and sorted member order so the archive is
    byte-identical across runs.
    """
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w", format=tarfile.GNU_FORMAT) as tar:
        for result in sorted(results, key=lambda r: r.index):
            body = ("\n".join(str(r) for r in sorted(result.residues)) + "\n") \
                if result.residues else ""
            data = body.encode("utf-8")
            info = tarfile.TarInfo(name=f"iter_{result.index:06d}.txt")
            info.size = len(data)
            info.mtime = 0
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            info.mode = 0o644
            tar.addfile(info, io.BytesIO(data))
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        with gzip.GzipFile(filename="", mode="wb", fileobj=fh, mtime=0) as gz:
            gz.write(buffer.getvalue())


def _emit_profile_warnings(warn: WarningCollector, profile: dict, results,
                           params: RobustnessParams) -> None:
    failed = profile["iterations"]["failed"]
    planned = profile["iterations"]["planned"]
    if failed and planned and failed / planned >= 0.50:
        warn.add("ROBUSTNESS_ITERATIONS_FAILED_EN_MASSE",
                 f"{failed} of {planned} iterations failed "
                 f"({failed / planned:.1%}). The profile reflects reconstruction "
                 f"failure more than measured geometric change and must be read as "
                 f"such. All failures remain in every denominator and contribute "
                 f"J = 0.",
                 potential_consequence="the robustness profile is dominated by "
                                       "failures rather than by geometry",
                 recommended_action="Lead review: inspect iteration_failures.tsv "
                                    "for a common cause before interpreting the "
                                    "profile",
                 affected_output="perturbation_results.tsv")

    preservation = profile["preservation"]["P_jaccard_ge_primary"]
    if preservation is not None and np.isfinite(preservation) \
            and preservation < params.preservation_primary:
        warn.add("LIMITED_GEOMETRIC_ROBUSTNESS",
                 f"the footprint was preserved at Jaccard >= "
                 f"{params.preservation_primary:g} in only {preservation:.1%} of "
                 f"iterations. The footprint depends substantially on which "
                 f"significant centers are present. Reported as measured; not "
                 f"repaired.",
                 potential_consequence="every downstream statement about this "
                                       "footprint must carry the robustness profile",
                 recommended_action="carry the profile with every hotspot claim",
                 affected_output="robustness_profile.json")


def _not_evaluable(ctx: RunContext, stage_dir: Path, warn: WarningCollector,
                   phase_c, handoff_03: Handoff, params: RobustnessParams,
                   n_s: int, baseline: dict, started: str, t0: float) -> Handoff:
    """``n_S < 2`` — reported honestly, never faked (§3.15)."""
    reason = (f"ROBUSTNESS_NOT_EVALUABLE: n_S = {n_s} significant hotspot center(s). "
              f"Perturbation removes geometric centers from S; with fewer than "
              f"{params.min_n_s} centers no non-empty proper subset exists, so no "
              f"iteration is possible. No surrogate analysis is substituted.")
    warn.add("ROBUSTNESS_NOT_EVALUABLE", reason,
             potential_consequence="no geometric robustness estimate exists for this "
                                   "footprint",
             recommended_action="report the footprint with an explicit "
                                "not-evaluable robustness statement",
             affected_output="robustness_profile.json")
    profile = {
        "result_type": "geometric_footprint_robustness_profile",
        "evaluable": False, "reason": reason, "n_S": n_s,
        "categorical_verdict_emitted": False,
        # The contract keys are present but null: a consumer must be able to read
        # them without a KeyError and must not be able to mistake "not evaluable"
        # for a low score.
        "jaccard_residues": distribution_summary([]),
        "dice_residues": distribution_summary([]),
        "preservation_freq_050": None,
        "preservation_freq_070": None,
        "iterations": {"planned": 0, "completed": 0, "failed": 0},
        "clinical_dataset_unmodified": True,
        "geometric_footprint_robustness_only": True,
        "v2_section8_underpowered_iteration_applicable":
            V2_SECTION8_UNDERPOWERED_APPLICABLE,
        "v2_section8_underpowered_iteration_note": V2_SECTION8_UNDERPOWERED_NOTE,
    }
    write_json(stage_dir / "robustness_profile.json", profile)
    _write_universe(stage_dir, phase_c, phase_c.inputs)
    write_tsv(stage_dir / "perturbation_design.tsv", [], DESIGN_COLUMNS)
    write_json(stage_dir / "subset_design.json", {
        "n_S": n_s, "K_MAX": k_max(n_s) if n_s >= 2 else 0, "N_CAP": params.n_cap,
        "total_subset_space_sum_T_k": 0, "total_evaluated": 0,
        "recorded_before_iteration_1": True,
        "reason_not_evaluable": reason,
    })
    write_tsv(stage_dir / "perturbation_results.tsv", [], PERTURBATION_COLUMNS)
    write_tsv(stage_dir / "iteration_failures.tsv", [], FAILURE_COLUMNS)
    write_tsv(stage_dir / "center_sensitivity.tsv", [], CENTER_SENSITIVITY_COLUMNS)
    write_tsv(stage_dir / "robustness_summary.tsv", [], AGGREGATE_COLUMNS)
    _bundle_footprints(stage_dir / "perturbation_footprints.tar.gz", [])
    fig_mod.figure_f9(stage_dir / "figures", [], [],
                      (params.preservation_primary, params.preservation_secondary))
    fig_mod.figure_f10(stage_dir / "figures", [])
    report_mod.write_not_evaluable_report(stage_dir / "stage_d_report.md", reason,
                                          phase_c, baseline, ctx.run_id, ctx.gene)

    payload = _handoff_payload(profile, None, params, phase_c, [])
    handoff = Handoff(name="handoff_04", run_id=ctx.run_id,
                      config_sha256=ctx.config.sha256,
                      qc_status=QC_PASS_WITH_WARNINGS,
                      manifest=manifest_for(stage_files(stage_dir), ctx.run_root),
                      payload=payload)
    assert_handoff_shape(handoff)
    write_json(ctx.handoff_path("handoff_04"), handoff.as_dict())
    StageStatus(
        stage=STAGE, agent_owner=AGENT, status=Status.COMPLETED_NEGATIVE,
        outcome_type=OutcomeType.NOT_APPLICABLE, reason=reason,
        started_utc=started, ended_utc=utc_now(),
        wall_seconds=time.perf_counter() - t0, upstream_handoff="handoff_03",
        upstream_qc_status=handoff_03.qc_status,
        n_outputs_expected=len(EXPECTED_09),
        n_outputs_created=len(stage_files(stage_dir)), warnings=warn.counts(),
        recommended_action="Report the footprint with an explicit not-evaluable "
                           "robustness statement.",
        negative_result={"condition": "ROBUSTNESS_NOT_EVALUABLE", "detail": reason},
    ).write(stage_dir)
    warn.write(stage_dir)
    write_stage_manifest(stage_dir, EXPECTED_09, extra_created=stage_files(stage_dir))
    return handoff


# --------------------------------------------------------------------------
# handoff / status / provenance
# --------------------------------------------------------------------------

def _handoff_payload(profile: dict, design, params: RobustnessParams, phase_c,
                     results) -> dict:
    return {
        "footprint_code_version": footprint_code_version(),
        "footprint_api_version": FOOTPRINT_API_VERSION,
        "robustness_profile": profile,
        "n_S": profile.get("n_S", design.n_S if design else None),
        "K_MAX": design.K_MAX if design else None,
        "total_subset_space": design.total_subset_space if design else 0,
        "total_evaluated": design.total_evaluated if design else 0,
        "per_level": [level.as_row() for level in design.levels] if design else [],
        "robustness_r_fp_invariant": profile.get("robustness_r_fp_invariant", {}),
        "r_fp": phase_c.r_fp,
        "clinical_dataset_unmodified": True,
        "geometric_footprint_robustness_only": True,
        "statement": ("Perturbation acted on geometric hotspot center positions. No "
                      "ClinVar record was removed and no class label was altered. "
                      "This result is geometric footprint robustness only and is NOT "
                      "validation of hotspot discovery."),
    }


def _build_handoff_04(ctx: RunContext, stage_dir: Path, profile: dict, design,
                      params: RobustnessParams, phase_c, results,
                      warn: WarningCollector) -> Handoff:
    payload = _handoff_payload(profile, design, params, phase_c, results)
    payload["n_S"] = design.n_S
    counts = warn.counts()
    qc = QC_PASS_WITH_WARNINGS if (counts[Severity.MAJOR.value]
                                   or counts[Severity.ADVISORY.value]) else QC_PASS
    handoff = Handoff(name="handoff_04", run_id=ctx.run_id,
                      config_sha256=ctx.config.sha256, qc_status=qc,
                      manifest=manifest_for(stage_files(stage_dir), ctx.run_root),
                      payload=payload)
    assert_handoff_shape(handoff)
    write_json(ctx.handoff_path("handoff_04"), handoff.as_dict())
    return handoff


def _write_status(stage_dir: Path, warn: WarningCollector, planned: int,
                  completed: int, failed: int, started: str, t0: float,
                  handoff_03: Handoff) -> None:
    StageStatus(
        stage=STAGE, agent_owner=AGENT, status=Status.COMPLETED,
        outcome_type=OutcomeType.COMPLETED,
        reason=(f"geometric footprint robustness profile computed over {planned} "
                f"iterations (completed {completed}, failed {failed})"),
        started_utc=started, ended_utc=utc_now(),
        wall_seconds=time.perf_counter() - t0, upstream_handoff="handoff_03",
        upstream_qc_status=handoff_03.qc_status,
        n_outputs_expected=len(EXPECTED_09),
        n_outputs_created=len(stage_files(stage_dir)), warnings=warn.counts(),
    ).write(stage_dir)


def _record_provenance(ctx: RunContext, up, fp_params, params: RobustnessParams,
                       design, phase_c, results, warn: WarningCollector) -> None:
    from ..utils.runctx import software_versions

    ctx.record_provenance(
        STAGE, AGENT,
        inputs={"upstream_files_sha256": up.hashes,
                "center_set_sha256": up.centers.digest(),
                "freeze_files_sha256": phase_c.freeze_payload.get("files_sha256", {}),
                "r_fp_original": phase_c.r_fp,
                "r_hot_recorded_readonly": up.r_hot.value_for_coverage_query},
        outputs={"planned": design.total_evaluated,
                 "completed": sum(1 for r in results if not r.failed),
                 "failed": sum(1 for r in results if r.failed)},
        parameters={**params.as_dict(),
                    "footprint_parameters": fp_params.as_dict(),
                    "footprint_code_version": phase_c.code_version,
                    "n_S": design.n_S, "K_MAX": design.K_MAX,
                    "per_level": [level.as_row() for level in design.levels],
                    "total_subset_space": design.total_subset_space,
                    "total_evaluated": design.total_evaluated,
                    "hotspot_statistics_recomputed": False,
                    "clinvar_records_modified": False,
                    "software_versions": software_versions()},
        commands=[f"hotspot3d.robustness.stage.run_phase_d(run_id={ctx.run_id})"]
        + [f"iteration {r.index:06d}: remove {list(r.removed)} "
           f"[seed_context={r.seed_context}] -> {r.status} "
           f"({r.wall_seconds:.3f}s)" for r in results],
        warnings=[w.as_row() for w in warn.items],
    )
