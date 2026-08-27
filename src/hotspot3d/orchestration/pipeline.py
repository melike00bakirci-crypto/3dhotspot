"""Pipeline orchestration — the Lead's own code.

Runs A -> B -> C/D -> E, enforces the handoff contracts, and assembles the
standardized output package. The orchestrator NEVER authors a scientific number
(P2): every value it writes into run-level files is copied from an agent-owned
artifact.

Negative-result design (IX.8): a stage that terminates the chain scientifically
still yields a complete, valid package. Downstream stages are still created as
directories with their own stage_status.json referencing the upstream cause, so
no mysteriously empty folder ever exists.
"""
from __future__ import annotations

import shutil
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from ..reporting.archives import build_archives, write_readme
from ..reporting.manifest import (ALL_STAGES, aggregate_warnings, build_manifest,
                                  dir_size, write_gene_summary, write_run_manifest,
                                  write_run_warnings)
from ..reporting.review_pack import build_review_pack
from ..reporting.terminal import render_terminal_summary
from ..utils.config import FrozenConfig, assert_frozen_methodology, load_config
from ..utils.errors import (BlockedError, EscalationRequired, LeakageError,
                            NegativeResult)
from ..utils.io import read_json, write_json, write_text, write_tsv
from ..utils.runctx import STAGE_OWNERS, RunContext, software_versions, uv_pip_freeze
from ..utils.status import (BANNERS, EXIT_CODES, OutcomeType, RunStatus, Status,
                            StageStatus, utc_now)
from .contracts import Handoff, assert_handoff_shape

# Which stage directories belong to which pipeline step, for NOT_RUN propagation.
STEP_STAGES = {
    "A": ["01_INPUT_RAW", "02_CLINVAR", "03_STRUCTURE_QC"],
    "B": ["04_GLOBAL_CLUSTERING", "05_HOTSPOT_RADIUS", "06_FINAL_HOTSPOTS", "11_SENSITIVITY"],
    "CD": ["07_FOOTPRINT_RADIUS", "08_FINAL_FOOTPRINT", "09_ROBUSTNESS"],
    "E": ["10_ANNOTATION"],
}
STEP_AGENT = {"A": "data-structure", "B": "hotspot-statistics",
              "CD": "footprint-robustness", "E": "biological-annotation"}


@dataclass
class RunResult:
    run_root: Path
    run_status: RunStatus
    outcome_type: OutcomeType
    exit_code: int
    terminal_summary: str
    metrics: dict = field(default_factory=dict)
    negative_result: dict | None = None
    escalations: list[dict] = field(default_factory=list)
    error: str | None = None


class _StageAAdapter:
    """Presents a variant source and a structure source as ONE Stage A provider.

    Stage A's contract takes a single provider that satisfies both protocols. When
    the caller supplies one object for both slots (the usual case) it is passed
    through untouched; when it supplies two, this adapter delegates each method to
    whichever object implements it.
    """

    def __init__(self, variant_source: Any, structure_source: Any):
        self._variant = variant_source
        self._structure = structure_source

    def __getattr__(self, name: str) -> Any:
        for provider in (self._variant, self._structure):
            if provider is not None and hasattr(provider, name):
                return getattr(provider, name)
        raise AttributeError(
            f"neither the variant source nor the structure source provides {name!r}")


@dataclass
class SourceBundle:
    """Providers injected by the caller. Synthetic runs pass mocks; a live run
    would pass the real clients. Keeping them injectable is what lets the whole
    pipeline be tested without touching the network."""
    variant_source: Any = None
    structure_source: Any = None
    annotation_source: Any = None

    @property
    def stage_a_source(self) -> Any:
        """The single provider Stage A's contract expects."""
        if self.structure_source is None or self.structure_source is self.variant_source:
            return self.variant_source
        if self.variant_source is None:
            return self.structure_source
        return _StageAAdapter(self.variant_source, self.structure_source)


def run_pipeline(gene: str, *, config: FrozenConfig | None = None,
                 config_path: str | Path = "config/pipeline.yaml",
                 results_root: str | Path = "results",
                 sources: SourceBundle | None = None,
                 synthetic: bool = False,
                 run_id: str | None = None,
                 stage_overrides: dict[str, Callable] | None = None) -> RunResult:
    cfg = config or load_config(config_path)
    assert_frozen_methodology(cfg)          # drift cannot enter silently

    ctx = RunContext.create(gene=gene, config=cfg, results_root=results_root,
                            synthetic=synthetic, run_id=run_id)
    sources = sources or SourceBundle()
    stages = _resolve_stages(stage_overrides)

    for stage in ALL_STAGES:               # P3 — every stage directory exists
        ctx.stage_dir(stage)

    _snapshot_reproducibility(ctx)

    negative: dict | None = None
    uninformative: dict | None = None
    escalations: list[dict] = []
    error: str | None = None
    not_created: dict[str, str] = {}
    handoffs: dict[str, Handoff] = {}
    completed_steps: list[str] = []

    order = [
        ("A", lambda: stages["A"](ctx, gene=gene, source=sources.stage_a_source)),
        ("B", lambda: stages["B"](ctx, upstream=handoffs["handoff_01"])),
        ("CD", lambda: stages["CD"](ctx, upstream=handoffs["handoff_02"])),
        ("E", lambda: stages["E"](ctx, handoff_02=handoffs["handoff_02"],
                                  handoff_03=handoffs["handoff_03"],
                                  handoff_04=handoffs["handoff_04"],
                                  source=sources.annotation_source)),
    ]

    for step, call in order:
        started = time.time()
        try:
            result = call()
            produced = result if isinstance(result, tuple) else (result,)
            for h in produced:
                assert_handoff_shape(h)
                handoffs[h.name] = h
            completed_steps.append(step)

            # Workflow v2 §0/§5.4 — UNDERPOWERED is checked FIRST, and is never a
            # negative. The test rejected nothing AND could not have rejected
            # anything, so the run says nothing about the gene. Stage B returns a
            # complete handoff rather than raising, so its artefacts (including the
            # power certificate) stay on disk; without this branch the orchestrator
            # would call Stage CD, which refuses on qc_status=FAIL, and the run
            # would be typed BLOCKED / TECHNICAL_FAILURE — misreporting an
            # uninformative run as a pipeline fault, the mirror of the misreport
            # that motivated v2.
            underpowered = next(
                (h for h in produced
                 if h.payload.get("terminal_state") == "UNDERPOWERED"), None)
            if underpowered is not None:
                uninformative = {
                    "stage": step, "agent": STEP_AGENT[step],
                    "terminal_state": "UNDERPOWERED",
                    **{k: underpowered.payload.get(k) for k in
                       ("power_certificate", "power_certificate_passed")
                       if underpowered.payload.get(k) is not None},
                }
                remaining = _next_step(step)
                if remaining is not None:
                    _mark_not_run(
                        ctx, step, remaining_after=remaining,
                        reason=("TEST_CANNOT_REJECT: the per-center test could not "
                                "have rejected any center regardless of the data "
                                "(p_floor > q/m). See power_certificate.json."),
                        outcome=OutcomeType.UNINFORMATIVE,
                        not_created=not_created)
                break

            # A stage may DECLARE a settled scientific negative in its handoff and
            # return normally, rather than raising NegativeResult. Both routes end
            # the chain for the same scientific reason and must be typed the same.
            #
            # Only the raised route was handled, so a stage that completed
            # COMPLETED_NEGATIVE still had its successor invoked; that successor
            # correctly refused to build on an empty result, but it refused with a
            # BlockedError, which types as TECHNICAL_FAILURE. The first real KCNA2
            # run came back BLOCKED / TECHNICAL_FAILURE / negative_result: null even
            # though Stage B had recorded a settled negative — a legitimate
            # scientific finding reported as a pipeline fault. (Under v2 that same
            # run is UNDERPOWERED and is caught by the branch above instead.)
            settled = next((h.negative_result for h in produced
                            if (h.negative_result or {}).get("settled_negative")), None)
            if settled:
                negative = {"stage": step, "agent": STEP_AGENT[step], **settled}
                remaining = _next_step(step)
                if remaining is not None:
                    _mark_not_run(
                        ctx, step, remaining_after=remaining,
                        reason=f"{settled.get('condition')}: {settled.get('detail')}",
                        outcome=OutcomeType.SCIENTIFIC_NEGATIVE,
                        not_created=not_created)
                break

        except NegativeResult as nr:
            negative = {"condition": nr.condition, "detail": nr.detail,
                        "stage": step, "agent": STEP_AGENT[step], **nr.context}
            _mark_not_run(ctx, step, remaining_after=step,
                          reason=f"{nr.condition}: {nr.detail}",
                          outcome=OutcomeType.SCIENTIFIC_NEGATIVE,
                          not_created=not_created)
            break

        except EscalationRequired as esc:
            escalations.append({
                "stage": step, "agent": STEP_AGENT[step], "ambiguity": esc.ambiguity,
                "options": esc.options, "consequences": esc.consequences,
                "recommendation": esc.recommendation,
            })
            _mark_not_run(ctx, step, remaining_after=step,
                          reason=f"ESCALATION: {esc.ambiguity}",
                          outcome=OutcomeType.NOT_APPLICABLE, not_created=not_created)
            break

        except (BlockedError, LeakageError) as exc:
            error = f"{type(exc).__name__}: {exc}"
            _mark_not_run(ctx, step, remaining_after=step, reason=error,
                          outcome=OutcomeType.TECHNICAL_FAILURE, not_created=not_created)
            break

        except Exception as exc:                              # unexpected fault
            error = f"{type(exc).__name__}: {exc}\n{traceback.format_exc(limit=6)}"
            _mark_not_run(ctx, step, remaining_after=step, reason=error,
                          outcome=OutcomeType.TECHNICAL_FAILURE, not_created=not_created)
            break

        finally:
            _ensure_step_status(ctx, step, wall=time.time() - started)

    return _finalize(ctx, handoffs=handoffs, negative=negative, error=error,
                     escalations=escalations, completed_steps=completed_steps,
                     not_created=not_created, uninformative=uninformative)


# --------------------------------------------------------------------------- #
# stage resolution                                                            #
# --------------------------------------------------------------------------- #

def _resolve_stages(overrides: dict[str, Callable] | None) -> dict[str, Callable]:
    """Late-bind the agent entry points so a partially implemented tree still runs."""
    overrides = overrides or {}
    resolved: dict[str, Callable] = {}
    spec = {
        "A": ("hotspot3d.data.stage", "run_stage_a"),
        "B": ("hotspot3d.hotspot.stage", "run_stage_b"),
        "CD": ("hotspot3d.footprint.stage", "run_stage_cd"),
        "E": ("hotspot3d.annotation.stage", "run_stage_e"),
    }
    for key, (module, func) in spec.items():
        if key in overrides:
            resolved[key] = overrides[key]
            continue
        try:
            mod = __import__(module, fromlist=[func])
            resolved[key] = getattr(mod, func)
        except (ImportError, AttributeError) as exc:
            resolved[key] = _missing_stage(key, module, func, exc)
    return resolved


def _missing_stage(key: str, module: str, func: str, exc: Exception):
    def _raise(*_args, **_kwargs):
        raise BlockedError(
            f"Stage {key} entry point {module}.{func} is not available ({exc}). "
            f"The owning agent has not delivered its implementation."
        )
    return _raise


# --------------------------------------------------------------------------- #
# status propagation                                                          #
# --------------------------------------------------------------------------- #

def _ensure_step_status(ctx: RunContext, step: str, wall: float) -> None:
    """Guarantee stage_status.json exists for every stage of a step (P3)."""
    for stage in STEP_STAGES[step]:
        path = ctx.stage_dir(stage) / "stage_status.json"
        if path.exists():
            continue
        StageStatus(
            stage=stage, agent_owner=STEP_AGENT[step], status=Status.NOT_RUN,
            outcome_type=OutcomeType.NOT_APPLICABLE,
            reason="Stage did not report a status; the owning step did not reach it.",
            started_utc=utc_now(), ended_utc=utc_now(), wall_seconds=wall,
        ).write(ctx.stage_dir(stage))


#: Chain order, shared by the not-run marker and the settled-negative handler.
STEP_ORDER = ["A", "B", "CD", "E"]


def _next_step(step: str) -> str | None:
    """The step after ``step``, or None if it is the last.

    Used when a stage COMPLETED but declared a settled negative: it has already
    written its own COMPLETED_NEGATIVE status, so only what follows is NOT_RUN.
    """
    index = STEP_ORDER.index(step)
    return STEP_ORDER[index + 1] if index + 1 < len(STEP_ORDER) else None


def _mark_not_run(ctx: RunContext, failed_step: str, *, remaining_after: str,
                  reason: str, outcome: OutcomeType,
                  not_created: dict[str, str]) -> None:
    """Every stage that could not run gets an explicit status referencing the cause."""
    steps = STEP_ORDER
    start = steps.index(remaining_after)
    for step in steps[start:]:
        for stage in STEP_STAGES[step]:
            status_path = ctx.stage_dir(stage) / "stage_status.json"
            if status_path.exists() and step == failed_step:
                continue
            is_origin = step == failed_step
            StageStatus(
                stage=stage, agent_owner=STEP_AGENT[step],
                status=Status.COMPLETED_NEGATIVE if (
                    is_origin and outcome is OutcomeType.SCIENTIFIC_NEGATIVE)
                    else Status.NOT_RUN,
                outcome_type=outcome if is_origin else OutcomeType.NOT_APPLICABLE,
                reason=reason if is_origin else
                    f"Upstream step {failed_step} terminated the chain: {reason}",
                started_utc=utc_now(), ended_utc=utc_now(),
                recommended_action=(
                    "No action — this is a valid scientific outcome."
                    if outcome is OutcomeType.SCIENTIFIC_NEGATIVE else
                    "Review the originating stage's NOT_RUN.txt and the Lead escalation."),
            ).write(ctx.stage_dir(stage))
            not_created[stage] = reason


# --------------------------------------------------------------------------- #
# reproducibility snapshot                                                    #
# --------------------------------------------------------------------------- #

def _snapshot_reproducibility(ctx: RunContext) -> None:
    repro = ctx.stage_dir("12_REPRODUCIBILITY")
    cfg = ctx.config

    # config.yaml is the EFFECTIVE merged configuration — every parameter this run
    # actually read, overlay included. It used to be a copy of the base file, which
    # for an overlaid run displayed values the run never used while the manifest
    # still called it "frozen config snapshot".
    #
    # These are the same bytes the config's sha256 digests, so the snapshot is the
    # configuration that digest identifies, not a re-rendering of it.
    write_text(repro / "config.yaml", cfg.effective_yaml())

    # ... and the sources are kept verbatim beside it, because the effective dump
    # cannot carry YAML comments and this repository states each parameter's
    # FROZEN provenance in exactly those comments.
    shutil.copy2(cfg.path, repro / "config_base.yaml")
    if cfg.overlay is not None:
        shutil.copy2(cfg.overlay, repro / "config_overlay.yaml")

    write_json(repro / "run_metadata.json", {
        "run_id": ctx.run_id, "gene": ctx.gene,
        "config_sha256": cfg.sha256, "code_version": ctx.code_version,
        # Byte-level provenance of each source, so an edit confined to comments
        # stays detectable even though it leaves config_sha256 unmoved.
        "config_effective_sha256": cfg.sha256,
        "config_base_path": str(cfg.path),
        "config_base_sha256": cfg.source_sha256,
        "config_overlay_path": str(cfg.overlay) if cfg.overlay else None,
        "config_overlay_sha256": cfg.overlay_sha256,
        "started_utc": ctx.started_utc, "synthetic_input_mode": ctx.synthetic,
        "run_root": str(ctx.run_root),
    })
    write_tsv(repro / "software_versions.tsv",
              [{"package": k, "version": v} for k, v in software_versions().items()],
              ["package", "version"])
    write_text(repro / "uv_pip_freeze.txt", uv_pip_freeze())
    write_text(repro / "INTERMEDIATE_RETENTION.md", _retention_doc())
    # P3 — absence is explicit: EVERY stage directory always carries a status,
    # including the two the Lead owns. Without this the reproducibility and summary
    # directories were the only ones a reader could find with no status at all.
    StageStatus(
        stage="12_REPRODUCIBILITY", agent_owner="lead", status=Status.COMPLETED,
        outcome_type=OutcomeType.COMPLETED,
        reason="Config snapshot, environment, seeds and retention policy recorded.",
        started_utc=ctx.started_utc, ended_utc=utc_now(),
    ).write(repro)


def _retention_doc() -> str:
    return """# Intermediate retention policy

Everything that is **not** deterministically regenerable is retained in full (P5).

## Retained
- Per-radius provisional center sets — `05_HOTSPOT_RADIUS/radius_scan_centers/r_<R>.tsv`
- Per-radius footprint component memberships — `07_FOOTPRINT_RADIUS/footprint_components/r_<R>.tsv`
- Every scanned candidate with raw values, normalized values, QC status, Pareto
  membership, distance-to-ideal and rejection reason
- Per-iteration footprint residue lists — `09_ROBUSTNESS/perturbation_footprints.tar.gz`

## Omitted, with a regeneration recipe

**Per-radius voxel occupancy grids.** Deterministically regenerable from the
significant center set and the radius. Only the final radius grid is stored
(`08_FINAL_FOOTPRINT/footprint_occupancy.npz`).

Regenerate any radius with:

```bash
python -m hotspot3d.footprint.regenerate_grid \\
    --centers FULL_RESULTS/06_FINAL_HOTSPOTS/significant_hotspot_centers.tsv \\
    --coords  FULL_RESULTS/03_STRUCTURE_QC/residue_coordinates.tsv \\
    --radius  <RHO> --voxel-h <H> --out grid_<RHO>.npz
```

`<H>` is recorded in `07_FOOTPRINT_RADIUS/footprint_domain.json`. The grid is a
pure function of (centers, radius, h, bounding box), so the regenerated grid is
byte-identical to the one used during the run.

## Rationale
Up to 10,000 robustness iterations x a directory each would consume ~50,000
inodes and break HPC inode quotas. Bundling costs one inode instead. Nothing
decision-bearing is lost.
"""


# --------------------------------------------------------------------------- #
# finalization                                                                #
# --------------------------------------------------------------------------- #

def _finalize(ctx: RunContext, *, handoffs: dict[str, Handoff], negative: dict | None,
              error: str | None, escalations: list[dict], completed_steps: list[str],
              not_created: dict[str, str],
              uninformative: dict | None = None) -> RunResult:
    warnings_rows = aggregate_warnings(ctx)
    n_blocking = sum(1 for w in warnings_rows if w.get("severity") == "BLOCKING")
    n_major = sum(1 for w in warnings_rows if w.get("severity") == "MAJOR")

    if error is not None:
        run_status = RunStatus.BLOCKED if "Blocked" in error or "BLOCKED" in error \
            else RunStatus.FAILED
        outcome = OutcomeType.TECHNICAL_FAILURE
    elif uninformative is not None:
        # Ranked ABOVE escalations and negatives on purpose. An underpowered run
        # carries a BLOCKING TEST_CANNOT_REJECT warning and an escalation, and it
        # has an empty center set; every one of those would otherwise capture it
        # and report it as PARTIAL or, worse, as COMPLETED_NEGATIVE. Workflow v2
        # Appendix A makes the negative reading a specification violation.
        run_status, outcome = RunStatus.UNDERPOWERED, OutcomeType.UNINFORMATIVE
    elif escalations:
        run_status, outcome = RunStatus.PARTIAL, OutcomeType.NOT_APPLICABLE
    elif negative is not None:
        run_status, outcome = RunStatus.COMPLETED_NEGATIVE, OutcomeType.SCIENTIFIC_NEGATIVE
    elif len(completed_steps) < 4:
        run_status, outcome = RunStatus.PARTIAL, OutcomeType.NOT_APPLICABLE
    elif n_major or n_blocking:
        run_status, outcome = RunStatus.COMPLETED_WITH_WARNINGS, OutcomeType.COMPLETED
    else:
        run_status, outcome = RunStatus.COMPLETED, OutcomeType.COMPLETED

    identity = _identity(ctx, handoffs)
    metrics = _metrics(ctx, handoffs, n_major=n_major, n_blocking=n_blocking)

    write_run_warnings(ctx, warnings_rows)
    manifest_rows = build_manifest(ctx, not_created)

    pack_info = build_review_pack(
        ctx, run_status=run_status, outcome_type=outcome.value, metrics=metrics,
        identity=identity, warnings_rows=warnings_rows,
        negative_result=negative, manifest_rows=manifest_rows)

    archives = build_archives(ctx)
    write_readme(ctx, run_status.value, BANNERS[run_status].format(gene=ctx.gene))

    paths = {
        "index_html": str((ctx.review_pack / "index.html").resolve()),
        "full_results": str(ctx.full_results.resolve()),
        "review_pack_zip": archives.get("review_pack_zip", "NA"),
        "full_results_tar": archives.get("full_results_tar", "NA"),
    }
    terminal = render_terminal_summary(
        gene=ctx.gene, run_id=ctx.run_id, run_status=run_status,
        outcome_type=outcome.value, identity=identity, metrics=metrics,
        warnings_rows=warnings_rows, paths=paths)

    summary_dir = ctx.stage_dir("00_RUN_SUMMARY")
    StageStatus(
        stage="00_RUN_SUMMARY", agent_owner="lead",
        status=Status.COMPLETED if run_status not in (RunStatus.FAILED,
                                                      RunStatus.BLOCKED)
               else Status.COMPLETED,
        outcome_type=OutcomeType.COMPLETED,
        reason=f"Run-level summary assembled; run_status={run_status.value}.",
        started_utc=ctx.started_utc, ended_utc=utc_now(),
    ).write(summary_dir)
    write_text(summary_dir / "terminal_summary.txt", terminal)
    write_text(summary_dir / "pipeline_flow.md",
               _pipeline_flow(ctx, completed_steps, negative, error, escalations))

    flat_metrics = {k: v.get("value") for k, v in metrics.items()}
    flat_metrics.update({"run_status": run_status.value, "outcome_type": outcome.value,
                         "uniprot_acc": identity.get("UniProt"),
                         "mane_transcript": identity.get("MANE transcript"),
                         "clinvar_release": identity.get("ClinVar release"),
                         "alphafold_model_version": identity.get("Structure"),
                         "review_pack_bytes": pack_info["bytes"],
                         "full_results_bytes": dir_size(ctx.full_results)})
    write_gene_summary(ctx, flat_metrics)

    run_summary = {
        "run_id": ctx.run_id, "gene": ctx.gene, "run_status": run_status.value,
        "outcome_type": outcome.value, "identity": identity,
        "metrics": metrics, "flags": _flag_summary(warnings_rows),
        "warning_counts": {"MAJOR": n_major, "BLOCKING": n_blocking,
                           "total": len(warnings_rows)},
        "stages": _stage_statuses(ctx), "escalations": escalations,
        "negative_result": negative, "error": error,
        "review_pack_bytes": pack_info["bytes"],
        "paths": paths, "archives": archives,
    }
    write_json(ctx.run_root / "run_summary.json", run_summary)
    write_json(summary_dir / "run_summary.json", run_summary)

    # manifest last so it sees every file it can, then re-write with pack rows
    write_run_manifest(ctx, build_manifest(ctx, not_created))

    if pack_info["over_budget"]:
        warnings_rows.append({
            "severity": "ADVISORY", "warning_code": "REVIEW_PACK_SIZE_EXCEEDED",
            "stage": "REVIEW_PACK", "agent_owner": "lead",
            "message": f"REVIEW_PACK is {pack_info['bytes']} bytes, over the "
                       f"{pack_info['budget']} byte budget. Nothing was dropped.",
            "potential_consequence": "Slower transfer", "recommended_action":
                "Review key_figures/ sizes", "affected_output": "REVIEW_PACK",
            "status": "ACCEPTED_BY_DESIGN", "detected_utc": utc_now(),
        })
        write_run_warnings(ctx, warnings_rows)

    final_root = ctx.finalize(ok=run_status not in (RunStatus.FAILED,))
    terminal = terminal.replace(str(ctx.run_root), str(final_root))

    return RunResult(run_root=final_root, run_status=run_status, outcome_type=outcome,
                     exit_code=EXIT_CODES[run_status], terminal_summary=terminal,
                     metrics=flat_metrics, negative_result=negative,
                     escalations=escalations, error=error)


def _identity(ctx: RunContext, handoffs: dict[str, Handoff]) -> dict:
    h1 = handoffs.get("handoff_01")
    p = h1.payload if h1 else {}
    return {
        "Gene": ctx.gene, "UniProt": p.get("uniprot_acc", "NA"),
        "MANE transcript": p.get("mane_transcript", "NA"),
        "RUN_ID": ctx.run_id,
        "Structure": p.get("alphafold_model_version", "NA"),
        "ClinVar release": p.get("clinvar_release", "NA"),
        "config_sha256": ctx.config.sha256[:16],
        "code_version": ctx.code_version[:16],
        "Input mode": "SYNTHETIC" if ctx.synthetic else "LIVE",
    }


def _metrics(ctx: RunContext, handoffs: dict[str, Handoff], *, n_major: int,
             n_blocking: int) -> dict:
    """Copy headline numbers VERBATIM from handoffs — the Lead authors nothing."""
    def g(name, key, default=None):
        h = handoffs.get(name)
        return h.payload.get(key, default) if h else default

    def entry(value, unit="", stage="", source=""):
        return {"value": value, "unit": unit, "stage": stage, "source": source}

    h1, h2 = "handoff_01", "handoff_02"
    h3, h4 = "handoff_03", "handoff_04"
    prof = g(h4, "robustness_profile", {}) or {}
    jac = prof.get("jaccard_residues", {}) if isinstance(prof, dict) else {}

    return {
        "n_missense_retrieved": entry(g(h1, "n_missense_retrieved"), "records", "02_CLINVAR"),
        "n_residues_plp": entry(g(h1, "N_P"), "residues", "02_CLINVAR"),
        "n_residues_blb": entry(g(h1, "N_B"), "residues", "02_CLINVAR"),
        "n_residue_class_conflict": entry(g(h1, "n_conflict"), "residues", "02_CLINVAR"),
        "n_positional_universe": entry(g(h1, "M"), "residues", "03_STRUCTURE_QC"),
        "r_hot": entry(g(h2, "hotspot_radius"), "A", "05_HOTSPOT_RADIUS"),
        "r_hot_domain_source": entry(g(h2, "search_domain_source"), "", "04_GLOBAL_CLUSTERING"),
        "fallback_radius_domain": entry(g(h2, "fallback_radius_domain"), "", "04_GLOBAL_CLUSTERING"),
        "permutation_count": entry(g(h2, "B"), "permutations", "06_FINAL_HOTSPOTS"),
        "permutation_resolution_limited": entry(g(h2, "permutation_resolution_limited"),
                                                "", "06_FINAL_HOTSPOTS"),
        "b_recommended": entry(g(h2, "b_recommended"), "", "06_FINAL_HOTSPOTS"),
        "fdr_method": entry(g(h2, "fdr_method"), "", "06_FINAL_HOTSPOTS"),
        "fdr_q": entry(g(h2, "q"), "", "06_FINAL_HOTSPOTS"),
        "n_significant_centers": entry(g(h2, "n_significant_centers"), "centers",
                                       "06_FINAL_HOTSPOTS"),
        "n_hotspot_regions": entry(g(h2, "n_hotspot_regions"), "regions", "06_FINAL_HOTSPOTS"),
        "n_centers_without_variant": entry(g(h2, "n_centers_without_variant"), "centers",
                                           "06_FINAL_HOTSPOTS"),
        "global_clustering_plp_p": entry(g(h2, "global_clustering_plp_p"), "",
                                         "04_GLOBAL_CLUSTERING"),
        "global_clustering_blb_p": entry(g(h2, "global_clustering_blb_p"), "",
                                         "04_GLOBAL_CLUSTERING"),
        "r_fp": entry(g(h3, "footprint_radius"), "A", "07_FOOTPRINT_RADIUS"),
        "footprint_volume_A3": entry(g(h3, "volume_A3"), "A^3", "08_FINAL_FOOTPRINT"),
        "footprint_surface_A2": entry(g(h3, "surface_area_A2"), "A^2", "08_FINAL_FOOTPRINT"),
        "footprint_n_residues": entry(g(h3, "n_footprint_residues"), "residues",
                                      "08_FINAL_FOOTPRINT"),
        "footprint_coverage": entry(g(h3, "coverage"), "fraction", "08_FINAL_FOOTPRINT"),
        "footprint_n_components": entry(g(h3, "n_components"), "components",
                                        "08_FINAL_FOOTPRINT"),
        "n_perturbations_possible": entry(g(h4, "total_subset_space"), "subsets",
                                          "09_ROBUSTNESS"),
        "n_perturbations_evaluated": entry(g(h4, "total_evaluated"), "iterations",
                                           "09_ROBUSTNESS"),
        "jaccard_median": entry(jac.get("median"), "", "09_ROBUSTNESS"),
        "jaccard_iqr_lo": entry(jac.get("iqr_lo"), "", "09_ROBUSTNESS"),
        "jaccard_iqr_hi": entry(jac.get("iqr_hi"), "", "09_ROBUSTNESS"),
        "preservation_freq_050": entry(prof.get("preservation_freq_050"), "",
                                       "09_ROBUSTNESS"),
        "preservation_freq_070": entry(prof.get("preservation_freq_070"), "",
                                       "09_ROBUSTNESS"),
        "n_warnings_major": entry(n_major, "warnings", "run"),
        "n_warnings_blocking": entry(n_blocking, "warnings", "run"),
    }


def _flag_summary(warnings_rows: list[dict]) -> list[dict]:
    seen, out = set(), []
    for w in warnings_rows:
        code = str(w.get("warning_code"))
        if code not in seen:
            seen.add(code)
            out.append({"warning_code": code, "severity": w.get("severity")})
    return out


def _stage_statuses(ctx: RunContext) -> list[dict]:
    out = []
    for stage in STAGE_OWNERS:
        path = ctx.full_results / stage / "stage_status.json"
        if path.is_file():
            data = read_json(path)
            out.append({"stage": stage, "status": data.get("status"),
                        "outcome_type": data.get("outcome_type"),
                        "agent_owner": data.get("agent_owner"),
                        "reason": data.get("reason", "")})
    return out


def _pipeline_flow(ctx: RunContext, completed: list[str], negative, error,
                   escalations) -> str:
    names = {"A": "data-structure", "B": "hotspot-statistics",
             "CD": "footprint-robustness (Phase C -> FREEZE GATE -> Phase D)",
             "E": "biological-annotation"}
    lines = [f"# Pipeline flow — {ctx.gene} / {ctx.run_id}", ""]
    for step in ["A", "B", "CD", "E"]:
        mark = "x" if step in completed else " "
        lines.append(f"- [{mark}] **{step}** — {names[step]}")
    lines += ["", "## Termination", ""]
    if error:
        lines.append(f"Technical failure: `{error.splitlines()[0]}`")
    elif escalations:
        lines.append(f"Escalated to the Lead: {escalations[0]['ambiguity']}")
    elif negative:
        lines.append(f"Scientific negative: **{negative['condition']}** — "
                     f"{negative['detail']}")
    else:
        lines.append("All four steps completed.")
    return "\n".join(lines) + "\n"
