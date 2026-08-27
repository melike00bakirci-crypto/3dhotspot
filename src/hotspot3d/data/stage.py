"""Stage A entry point — ``run_stage_a``.

Composes the whole of Stage A: retrieval, the four ClinVar strata, structural QC,
pLDDT recording, the CA positional universe ``U_struct``, the classified cohort
``L``, and ``handoff_01.json``.

Stage A is deterministic. No randomness enters it, so every derived-seed field is
written explicitly as ``null`` with the reason ``not_applicable_deterministic_stage``.

Stage A is also blind downstream: nothing here reads, requests or reasons about
any Stage B/C/D/E artifact. Knowing where hotspots landed would bias curation,
mapping and exclusion decisions, so the gate is one-way by construction.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..orchestration.contracts import (
    FORBIDDEN_DOWNSTREAM_COLUMNS,
    QC_FAIL,
    QC_PASS,
    QC_PASS_WITH_WARNINGS,
    STAGE_B_PERMITTED_COLUMNS,
    Handoff,
)
from ..structure.cif import write_structure_with_plddt
from ..structure.mapping import (
    CLASSIFIED_COHORT_COLUMNS,
    POSITIONAL_UNIVERSE_COLUMNS,
    RESIDUE_COORDINATE_COLUMNS,
    map_cohort,
    positional_universe,
    residue_coordinates,
)
from ..structure.plddt import (
    PLDDT_PROFILE_COLUMNS,
    PLDDT_REGION_COLUMNS,
    band_labels,
    low_confidence_regions,
    plddt_profile,
    summarize,
)
from ..structure.qc import assert_numbering_usable, structure_qc
from ..structure.sources import assert_model_admissible, assert_oligomeric_assembly_permitted
from ..utils.errors import BlockedError, EscalationRequired, NegativeResult
from ..utils.hashing import manifest_for, sha256_bytes, sha256_file
from ..utils.io import write_json, write_text, write_tsv
from ..utils.runctx import RunContext, uv_pip_freeze
from ..utils.status import OutcomeType, Severity, StageStatus, Status, WarningCollector, utc_now
from . import figures as figs
from .audit import assert_no_star_based_inclusion
from .classify import (
    CLASS_BLB,
    CLASS_CONFLICT,
    CLASS_PLP,
    REASON_NON_MISSENSE_VARIANT_TYPE,
    REASON_UNPARSEABLE_HGVS,
    SignificancePolicy,
    compound_policy_cost,
)
from .clinvar import normalize_records, raw_tsv_text
from .cohort import (
    CONFLICT_COLUMNS,
    STRATUM1_COLUMNS,
    STRATUM2_COLUMNS,
    STRATUM3_COLUMNS,
    STRATUM4_COLUMNS,
    apply_structural_exclusions,
    assert_cohort_viable,
    cohort_counts,
    collapse_to_residues,
    evaluate_records,
    prove_conservation,
    review_star_distribution,
    star_summary,
)
from .sources import (
    OLIGOMERIC_STATE_OLIGOMERIC,
    assert_release_metadata,
    assert_source_complete,
    raw_payload_of,
)

AGENT = "data-structure"
STAGE_RAW = "01_INPUT_RAW"
STAGE_CLINVAR = "02_CLINVAR"
STAGE_STRUCT = "03_STRUCTURE_QC"
STAGES = (STAGE_RAW, STAGE_CLINVAR, STAGE_STRUCT)

#: warnings_<stage>.tsv suffixes
WARN_SUFFIX = {STAGE_RAW: "01_input_raw", STAGE_CLINVAR: "02_clinvar",
               STAGE_STRUCT: "03_structure_qc"}

#: ``qc_status`` is decided from the FROZEN severity vocabulary (IX.7) alone, with
#: no per-code exception list. BLOCKING means "a stage or the run stopped; result
#: unusable as-is", and that is exactly what qc_status = FAIL communicates. A list
#: of codes exempted from the verdict would be a methodology change in
#: implementation clothing: the next exemption extends it, and the severity
#: vocabulary stops meaning what IX.7 says it means.

MANIFEST_COLUMNS = ("path", "stage", "description", "status", "sha256",
                    "size_bytes", "reason")

#: Every file Stage A is expected to produce (Output Contract IX.3 / P3: anything
#: not produced still appears, with status NOT_CREATED and a reason).
EXPECTED_FILES: dict[str, tuple[tuple[str, str], ...]] = {
    STAGE_RAW: (
        ("clinvar_raw.tsv", "immutable ClinVar payload as retrieved"),
        ("clinvar_raw.tsv.sha256", "SHA-256 of the ClinVar payload"),
        ("uniprot_canonical.fasta", "UniProt canonical isoform sequence"),
        ("uniprot_canonical.fasta.sha256", "SHA-256 of the UniProt sequence"),
        ("alphafold_model.cif", "AlphaFold model as served"),
        ("alphafold_model.cif.sha256", "SHA-256 of the AlphaFold model"),
        ("retrieval_log.json", "URLs, queries, release dates, status, retries"),
        ("stage_status.json", "stage outcome (always written)"),
        ("warnings_01_input_raw.tsv", "stage warnings"),
        ("stage_manifest.tsv", "expected vs created files with hashes"),
    ),
    STAGE_CLINVAR: (
        ("variants_missense_all.tsv", "stratum 1 — every retrieved record, unfiltered"),
        ("variants_residue_level.tsv", "stratum 2 — residue-level collapse"),
        ("variants_excluded_from_primary.tsv", "stratum 3 — exclusions with reasons"),
        ("review_star_distribution.tsv", "stratum 4 — stars, described never used"),
        ("residue_class_conflicts.tsv", "F3 RESIDUE_CLASS_CONFLICT residues"),
        ("cohort_summary.json", "cohort headline counts and audits"),
        ("figures/significance_composition.png", "records by significance class"),
        ("figures/review_star_distribution.png", "review-star composition"),
        ("figures/cohort_along_sequence.png", "cohort position along the sequence"),
        ("stage_status.json", "stage outcome (always written)"),
        ("warnings_02_clinvar.tsv", "stage warnings"),
        ("stage_manifest.tsv", "expected vs created files with hashes"),
    ),
    STAGE_STRUCT: (
        ("structure_qc.json", "structural QC checks and verdict"),
        ("residue_coordinates.tsv", "CA/CB/side-chain coordinates and pLDDT"),
        ("positional_universe.tsv", "U_struct — every residue with a usable CA"),
        ("classified_cohort.tsv", "L — the binary P/LP vs B/LB cohort"),
        ("plddt_profile.tsv", "per-residue pLDDT and band"),
        ("plddt_regions.tsv", "contiguous low-confidence regions (never excluded)"),
        ("mapping_report.json", "ClinVar-to-structure mapping and unmapped residues"),
        ("structures/structure_with_plddt.cif", "model with B-factor = pLDDT"),
        ("figures/plddt_profile.png", "pLDDT profile with band edges"),
        ("figures/plddt_bands.png", "pLDDT band distribution"),
        ("figures/universe_projection.png", "U_struct projection with L overlaid"),
        ("stage_a_report.md", "Stage A narrative report"),
        ("handoff_01.json", "handoff to hotspot-statistics"),
        ("stage_status.json", "stage outcome (always written)"),
        ("warnings_03_structure_qc.tsv", "stage warnings"),
        ("stage_manifest.tsv", "expected vs created files with hashes"),
    ),
}

#: FROZEN constants Stage A depends on. Drift here is a methodological revision.
FROZEN_STAGE_A_ASSERTIONS: tuple[tuple[str, Any], ...] = (
    ("representation.coordinate_atom", "CA"),
    ("representation.distance_metric", "euclidean"),
    ("clinvar.review_star_filter", None),
    ("clinvar.use_review_stars_for_inclusion", False),
    ("clinvar.non_binary_disposition", "preserve_outside_primary"),
    ("clinvar.residue_class_conflict.policy", "exclude_from_primary"),
    ("clinvar.residue_class_conflict.preserve_all_source_records", True),
    ("transcript.reference", "MANE_SELECT"),
    ("transcript.require_ref_aa_match", True),
    ("structure.allow_homology_models", False),
    ("structure.allow_ortholog_structures", False),
    ("structure.multi_fragment_policy", "block_and_escalate"),
    ("plddt.primary_filtering_enabled", False),
    ("plddt.pae_threshold", None),
)

REQUIRED_CONFIG_KEYS = [
    "clinvar.pathogenic_labels", "clinvar.benign_labels",
    "clinvar.exclusion_reason_enum", "clinvar.source", "clinvar.variant_type",
    "clinvar.compound_significance_policy", "clinvar.min_classified_residues",
    "transcript.fallback", "structure.source", "structure.model_template",
    "structure.missing_entry_policy", "structure.allow_monomer_for_obligate_oligomer",
    "structure.oligomeric_state_source", "plddt.bands", "plddt.sensitivity_threshold",
    "representation.emit_unused_atoms", "output.full_results.figure_dpi",
]


# --- bookkeeping ------------------------------------------------------------

@dataclass
class _Emitted:
    stage: str
    name: str
    path: Path
    description: str
    status: str = "CREATED"
    reason: str = "NA"


@dataclass
class _Run:
    """Mutable Stage A bookkeeping: what was written, warned and decided."""

    ctx: RunContext
    gene: str
    started_utc: str = field(default_factory=utc_now)
    t0: float = field(default_factory=time.perf_counter)
    warnings: dict = field(default_factory=dict)
    emitted: list[_Emitted] = field(default_factory=list)
    commands: list[str] = field(default_factory=list)
    step_seconds: dict[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.warnings = {s: WarningCollector(WARN_SUFFIX[s], AGENT) for s in STAGES}

    # -- paths ---------------------------------------------------------------
    def dir(self, stage: str) -> Path:
        self.ctx.assert_may_write(AGENT, stage)
        return self.ctx.stage_dir(stage)

    def rel(self, path: Path) -> str:
        return str(Path(path).relative_to(self.ctx.run_root))

    def note(self, stage: str, name: str, path: Path, description: str) -> Path:
        self.emitted.append(_Emitted(stage, name, Path(path), description))
        return Path(path)

    def skip(self, stage: str, name: str, description: str, reason: str) -> None:
        self.emitted.append(_Emitted(stage, name, self.dir(stage) / name, description,
                                     status="NOT_CREATED", reason=reason))

    def created_names(self, stage: str) -> set[str]:
        return {e.name for e in self.emitted if e.stage == stage and e.status == "CREATED"}

    def paths(self) -> list[Path]:
        return [e.path for e in self.emitted if e.status == "CREATED" and e.path.is_file()]

    def warn(self, stage: str, code: str, message: str, **kwargs):
        return self.warnings[stage].add(code, message, **kwargs)

    def all_warnings(self) -> list:
        """Every Stage A warning, most severe first (IX.7 reading order)."""
        order = {Severity.BLOCKING: 0, Severity.MAJOR: 1,
                 Severity.ADVISORY: 2, Severity.INFO: 3}
        items = [w for stage in STAGES for w in self.warnings[stage].items]
        return sorted(items, key=lambda w: (order[w.severity], w.stage, w.warning_code))

    def warning_counts(self) -> dict[str, int]:
        out = {s.value: 0 for s in Severity}
        for w in self.all_warnings():
            out[w.severity.value] += 1
        return out

    def timed(self, label: str, start: float) -> None:
        self.step_seconds[label] = round(time.perf_counter() - start, 6)


# --- entry point ------------------------------------------------------------

def run_stage_a(ctx: RunContext, *, gene: str, source: Any | None = None) -> Handoff:
    """Execute Stage A and return ``handoff_01``.

    ``source`` is a provider satisfying :class:`hotspot3d.data.sources.StageASource`
    (variants + canonical sequence + structure). Synthetic runs always pass one.
    With ``source=None`` and ``ctx.synthetic=False`` the live clients would be
    used; they are not enabled in the implementation phase and refuse.
    """
    run = _Run(ctx=ctx, gene=gene)
    try:
        return _execute(run, source)
    except NegativeResult:
        raise
    except (BlockedError, EscalationRequired) as exc:
        _write_terminal_status(run, Status.BLOCKED, OutcomeType.TECHNICAL_FAILURE, str(exc))
        raise
    except Exception as exc:                                    # pragma: no cover
        _write_terminal_status(run, Status.FAILED, OutcomeType.TECHNICAL_FAILURE,
                               f"{type(exc).__name__}: {exc}")
        raise


def _write_terminal_status(run: _Run, status: Status, outcome: OutcomeType,
                           reason: str) -> None:
    """P3 — stage_status.json is written for every stage, always, even on failure."""
    for stage in STAGES:
        try:
            stage_dir = run.ctx.stage_dir(stage)
            run.warnings[stage].write(stage_dir)
            StageStatus(
                stage=stage, agent_owner=AGENT, status=status, outcome_type=outcome,
                reason=reason, started_utc=run.started_utc, ended_utc=utc_now(),
                wall_seconds=time.perf_counter() - run.t0,
                n_outputs_expected=len(EXPECTED_FILES[stage]),
                n_outputs_created=len(run.created_names(stage)),
                warnings=run.warning_counts(),
                recommended_action=(
                    "Stage A stopped before producing a usable cohort. Read this "
                    "reason, resolve it with the Lead, and re-run under a new RUN_ID; "
                    "no downstream stage may run against this directory."),
            ).write(stage_dir)
        except Exception:                                       # pragma: no cover
            continue


# --- the stage itself -------------------------------------------------------

def _execute(run: _Run, source: Any | None) -> Handoff:
    ctx, gene = run.ctx, run.gene
    cfg = ctx.config

    # -- 1. preflight --------------------------------------------------------
    step = time.perf_counter()
    _validate_inputs(run)
    star_audit = assert_no_star_based_inclusion()
    policy = SignificancePolicy.from_config(cfg)
    bands = [float(b) for b in cfg.get("plddt.bands")]
    sensitivity_threshold = float(cfg.get("plddt.sensitivity_threshold"))
    source = _resolve_source(run, source)
    run.timed("preflight", step)

    if ctx.synthetic:
        run.warn(STAGE_RAW, "SYNTHETIC_INPUT_MODE",
                 "Stage A ran against a synthetic source; no network retrieval occurred.",
                 potential_consequence="Results describe fixture data, not ClinVar.",
                 recommended_action="Do not interpret this run biologically.",
                 affected_output="FULL_RESULTS/01_INPUT_RAW/")

    # -- 2. retrieval --------------------------------------------------------
    step = time.perf_counter()
    resolution = source.resolve_gene(gene)
    sequence_meta = assert_release_metadata(dict(source.sequence_metadata()), "UniProt")
    raw_rows = list(source.fetch_missense(gene))
    clinvar_meta = assert_release_metadata(dict(source.release_metadata()), "ClinVar")
    model = source.fetch_model(resolution.uniprot_acc)
    model_meta = dict(source.model_metadata())
    run.commands.append(
        f"run_stage_a(gene={gene!r}, source={type(source).__name__}, "
        f"synthetic={ctx.synthetic})")
    run.timed("retrieval", step)

    if str(model.model_version).strip() in ("", "NA", "None"):
        raise BlockedError(
            f"BLOCKED — the AlphaFold model for {resolution.uniprot_acc} carries no "
            f"served version. An unversioned payload is never analysed (§6.5).")

    try:
        assert_model_admissible(model)
    except BlockedError as exc:
        run.warn(STAGE_STRUCT, "MULTI_FRAGMENT_AFDB_ENTRY", str(exc),
                 potential_consequence="Coordinates would not span the analysed sequence.",
                 recommended_action="Escalate to the Lead; fragment stitching is prohibited.",
                 affected_output="FULL_RESULTS/03_STRUCTURE_QC/")
        raise

    # -- 2b. oligomeric assembly (Workflow v2 §2) ----------------------------
    # UniProt's own subunit annotation, resolved alongside the gene (§6.2). A
    # multimer model is never fetched (no such client exists); the two
    # remaining paths are BLOCKED (allow_monomer=FALSE) or the declared,
    # advisory monomer exception (allow_monomer=TRUE, the current default).
    allow_monomer_oligomer = bool(cfg.get("structure.allow_monomer_for_obligate_oligomer"))
    try:
        assert_oligomeric_assembly_permitted(resolution.oligomeric_state,
                                             allow_monomer=allow_monomer_oligomer)
    except BlockedError as exc:
        run.warn(STAGE_STRUCT, "OLIGOMERIC_ASSEMBLY_UNAVAILABLE", str(exc),
                 potential_consequence="A monomer model omits every residue "
                                       "contributed by a neighbouring subunit; "
                                       "any spatial neighbourhood spanning a subunit "
                                       "interface is invisible to this run.",
                 recommended_action="Escalate to the Lead. This pipeline has no "
                                    "multimer-fetching capability: either authorize "
                                    "the monomer exception "
                                    "(structure.allow_monomer_for_obligate_oligomer) "
                                    "or supply a multimer model out of band.",
                 affected_output="FULL_RESULTS/03_STRUCTURE_QC/structure_qc.json")
        raise

    intra_subunit_only = resolution.oligomeric_state == OLIGOMERIC_STATE_OLIGOMERIC
    if intra_subunit_only:
        run.warn(STAGE_STRUCT, "OLIGOMERIC_ASSEMBLY_MONOMER_DECLARED",
                 f"UniProt's subunit annotation for {resolution.uniprot_acc} reads as an "
                 f"obligate homo-/hetero-oligomer (matched phrase: "
                 f"{resolution.oligomeric_state_matched_phrase!r}). Proceeding on the "
                 f"monomeric AlphaFold model per "
                 f"structure.allow_monomer_for_obligate_oligomer=TRUE. Every hotspot and "
                 f"footprint downstream of this run must be labelled intra-subunit only.",
                 potential_consequence="Coordinates omit residues contributed by "
                                       "neighbouring subunits.",
                 recommended_action="Read structure_qc.json > oligomeric_assembly; treat "
                                    "every downstream hotspot/footprint as intra-subunit "
                                    "only, per the declared limitation.",
                 affected_output="FULL_RESULTS/03_STRUCTURE_QC/structure_qc.json",
                 status="ACCEPTED_BY_DESIGN")

    # -- 3. immutable raw payloads ------------------------------------------
    step = time.perf_counter()
    _write_raw_inputs(run, source, resolution, raw_rows, model, clinvar_meta,
                      sequence_meta, model_meta)
    run.timed("write_01_input_raw", step)

    # -- 4. ClinVar evaluation and residue collapse -------------------------
    step = time.perf_counter()
    records = normalize_records(raw_rows, gene=gene)
    evals = evaluate_records(
        records, policy=policy, canonical_sequence=resolution.sequence,
        mane_transcript=resolution.mane_transcript,
        require_ref_aa_match=bool(cfg.get("transcript.require_ref_aa_match")),
    )
    residue_rows = collapse_to_residues(evals, canonical_sequence=resolution.sequence)
    compound_cost = compound_policy_cost(
        records, pathogenic_labels=cfg.get("clinvar.pathogenic_labels"),
        benign_labels=cfg.get("clinvar.benign_labels"),
        configured_policy=cfg.get("clinvar.compound_significance_policy"))
    run.timed("clinvar_evaluation", step)

    # -- 5. structure QC, pLDDT, coordinates --------------------------------
    step = time.perf_counter()
    qc = structure_qc(model, canonical_sequence=resolution.sequence)
    assert_numbering_usable(qc)
    coordinate_rows = residue_coordinates(model, bands)
    universe_rows = positional_universe(coordinate_rows)
    profile = plddt_profile(model, bands)
    plddt_summary = summarize(profile, bands)
    run.timed("structure_qc", step)

    # -- 6. structural exclusions, then U_struct / L ------------------------
    usable = {row["residue_index"] for row in universe_rows}
    excluded_no_ca = apply_structural_exclusions(
        residue_rows, [r.residue_index for r in residue_rows
                       if r.residue_index not in usable])
    mapping = map_cohort(residue_rows, coordinate_rows)
    counts = cohort_counts(residue_rows)
    conservation = prove_conservation(evals)
    regions = low_confidence_regions(
        profile, sorted({bands[0], sensitivity_threshold}),
        cohort_residues=[row["residue_index"] for row in mapping.cohort])

    _raise_on_structural_failure(run, qc, mapping, conservation, counts)

    # -- 7. 02_CLINVAR ------------------------------------------------------
    step = time.perf_counter()
    star_rows = review_star_distribution(evals)
    _write_clinvar_outputs(run, evals, residue_rows, star_rows, counts, conservation,
                           policy, resolution, clinvar_meta, star_audit,
                           excluded_no_ca, plddt_summary, sensitivity_threshold, profile,
                           mapping, compound_cost)
    run.timed("write_02_clinvar", step)

    # -- 8. 03_STRUCTURE_QC -------------------------------------------------
    step = time.perf_counter()
    _write_structure_outputs(run, model, qc, coordinate_rows, universe_rows,
                             mapping, profile, regions, plddt_summary, bands,
                             model_meta, resolution, intra_subunit_only,
                             allow_monomer_oligomer)
    run.timed("write_03_structure_qc", step)

    # -- 9. cohort viability — BLOCKED, but the evidence stays complete ------
    negative: dict | None = None
    try:
        assert_cohort_viable(counts)
    except NegativeResult as exc:
        negative = {
            "condition": exc.condition, "detail": exc.detail, "context": exc.context,
            "outcome_type": OutcomeType.SCIENTIFIC_NEGATIVE.value,
            "negative_kind": "input_cannot_support_the_analysis",
            "nothing_was_tested": True,
            "interpretation": (
                "The evidence base for this gene cannot support the analysis: a binary "
                "class is empty or |L| < 2, so no spatial test was ever run. This is "
                "NOT a finding of 'no clustering' — nothing was looked for and nothing "
                "was found. The cause is ClinVar's content rather than a fault in the "
                "code, which is why outcome_type is SCIENTIFIC_NEGATIVE and not "
                "TECHNICAL_FAILURE, but the run is BLOCKED and licenses no method "
                "modification: no VUS or conflicting record was admitted, no conflict "
                "residue was reassigned, no non-canonical transcript was accepted and "
                "no review-star filter was introduced to change it."),
            "downstream": ("hotspot-statistics must NOT run; the A->B handoff is "
                           "blocked (qc_status = FAIL)."),
        }
        run.warn(STAGE_CLINVAR, "INSUFFICIENT_CLASSIFIED_RESIDUES", exc.detail,
                 potential_consequence=("Nothing was tested; no statement about "
                                        "clustering in this gene is available."),
                 recommended_action=("Report the cohort as unanalysable; do not relax "
                                     "any criterion to manufacture one."),
                 affected_output="FULL_RESULTS/03_STRUCTURE_QC/classified_cohort.tsv")

    _post_cohort_warnings(run, counts, mapping, plddt_summary, profile,
                          sensitivity_threshold, negative)

    # -- 10. status, handoff, manifests, provenance, report -----------------
    qc_status = _qc_status(run)
    status, outcome = _stage_outcome(negative)
    for stage in STAGES:
        stage_dir = run.dir(stage)
        run.note(stage, f"warnings_{WARN_SUFFIX[stage]}.tsv",
                 run.warnings[stage].write(stage_dir), "stage warnings")
        run.note(stage, "stage_status.json", StageStatus(
            stage=stage, agent_owner=AGENT, status=status, outcome_type=outcome,
            reason=_stage_reason(stage, counts, mapping, negative),
            started_utc=run.started_utc, ended_utc=utc_now(),
            wall_seconds=time.perf_counter() - run.t0,
            n_outputs_expected=len(EXPECTED_FILES[stage]),
            n_outputs_created=len(run.created_names(stage)),
            warnings=run.warnings[stage].counts(),
            recommended_action=_recommended_action(negative),
            negative_result=negative,
        ).write(stage_dir), "stage outcome (always written)")

    _write_report(run, resolution, clinvar_meta, model, qc, counts, mapping,
                  conservation, plddt_summary, star_audit, negative, qc_status,
                  excluded_no_ca, status, star_summary(evals), compound_cost)
    handoff = _write_handoff(run, resolution, clinvar_meta, model, counts, mapping,
                             qc_status, negative, intra_subunit_only,
                             allow_monomer_oligomer)
    _write_manifests(run)
    _write_provenance(run, resolution, clinvar_meta, sequence_meta, model_meta,
                      counts, mapping, star_audit, conservation, len(evals))
    return handoff


# --- preflight --------------------------------------------------------------

def _validate_inputs(run: _Run) -> None:
    """§6 — every failure BLOCKS and nothing is ever defaulted."""
    ctx, cfg = run.ctx, run.ctx.config
    if run.gene != ctx.gene:
        raise BlockedError(
            f"BLOCKED — Stage A was asked to run {run.gene!r} inside a RunContext for "
            f"{ctx.gene!r}. One run root holds exactly one gene.")

    cfg.require_all(REQUIRED_CONFIG_KEYS)
    for dotted, expected in FROZEN_STAGE_A_ASSERTIONS:
        cfg.assert_frozen_value(dotted, expected)

    # Two distinct freeze checks, because they catch different edits.
    #
    # 1. The EFFECTIVE parameter tree. Re-read every source, re-merge, re-digest.
    #    Asking the config to do this is the only correct way: the digest is a
    #    function of the merged tree, and a bare hash of cfg.path is neither the
    #    same value nor even the same KIND of value once an overlay is in play.
    # 2. The source BYTES. This repository keeps its FROZEN annotations in the
    #    comments of config/pipeline.yaml, and a comment-only edit leaves the
    #    effective tree — and so check 1 — completely unmoved.
    recomputed = cfg.recompute_sha256()
    if recomputed is not None and recomputed != cfg.sha256:
        raise BlockedError(
            f"BLOCKED — config_sha256 mismatch: the Lead froze {cfg.sha256} but "
            f"{cfg.describe_source()} now yields {recomputed}. "
            f"The effective configuration changed after it was frozen.")

    base_now, overlay_now = cfg.recompute_source_digests()
    for label, source, frozen, current in (
            ("config", cfg.path, cfg.source_sha256, base_now),
            ("gene overlay", cfg.overlay, cfg.overlay_sha256, overlay_now)):
        if frozen is None or current is None:
            continue                      # no overlay, or the source is gone
        if current != frozen:
            raise BlockedError(
                f"BLOCKED — {label} file changed after it was read: {source} was "
                f"{frozen} and is now {current}. The effective parameters may be "
                f"unaffected (a comment or formatting edit), but this repository "
                f"keeps its FROZEN annotations in those comments, so the edit is "
                f"reported rather than ignored.")

    rerun = bool(getattr(ctx, "rerun", False))
    for stage in STAGES:
        stage_dir = ctx.full_results / stage
        existing = [p for p in stage_dir.rglob("*") if p.is_file()] if stage_dir.exists() else []
        if existing and not rerun:
            raise BlockedError(
                f"BLOCKED — {stage} already contains {len(existing)} file(s) for "
                f"RUN_ID {ctx.run_id}. 01_INPUT_RAW is write-once and a stage directory "
                f"is never silently overwritten; the Lead must declare a re-run or a "
                f"new RUN_ID.")
        ctx.assert_may_write(AGENT, stage)


def _resolve_source(run: _Run, source: Any | None):
    if source is None:
        if run.ctx.synthetic:
            raise BlockedError(
                "BLOCKED — a synthetic run must be given an explicit source; Stage A "
                "never fabricates input data.")
        from .sources import LiveStageASource

        # The gate and the model version both come from the frozen config: the
        # bundle picks neither. allow_network=FALSE reaches every client, so each
        # external retrieval refuses with its own BLOCKED message instead of
        # substituting a cached or hand-made payload.
        source = LiveStageASource(
            allow_network=bool(run.ctx.config.get_optional("execution.allow_network", False)),
            model_template=run.ctx.config.get("structure.model_template"))
    assert_source_complete(source)
    return source


# --- 01_INPUT_RAW -----------------------------------------------------------

def _write_raw_inputs(run: _Run, source, resolution, raw_rows, model,
                      clinvar_meta: dict, sequence_meta: dict, model_meta: dict) -> None:
    """Write the immutable payloads plus their SHA-256 sidecars and the log."""
    stage_dir = run.dir(STAGE_RAW)
    log: list[dict] = []
    failed: list[dict] = []

    def emit(name: str, description: str, payload: bytes, meta: dict, verbatim: bool):
        path = stage_dir / name
        path.write_bytes(payload)
        digest = sha256_bytes(payload)
        (stage_dir / f"{name}.sha256").write_text(f"{digest}  {name}\n",
                                                  encoding="utf-8", newline="\n")
        run.note(STAGE_RAW, name, path, description)
        run.note(STAGE_RAW, f"{name}.sha256", stage_dir / f"{name}.sha256",
                 f"SHA-256 of {name}")
        log.append({
            "artifact": name, "sha256": digest, "n_bytes": len(payload),
            "payload_is_source_verbatim": verbatim,
            "url": meta.get("url", "NA"), "query": meta.get("query", "NA"),
            "release_date": meta.get("release_date", "NA"),
            "http_status": meta.get("http_status", "NA"),
            "retries": meta.get("retries", 0),
            "retrieved_utc": utc_now(),
            "synthetic": bool(run.ctx.synthetic),
            "metadata": meta,
        })
        # Every attempt that failed before the payload arrived, as the source
        # recorded it. A retrieval that needed retries says so here; this list
        # was previously hardcoded empty, which reported a clean transfer
        # whatever actually happened (§13).
        failed.extend({"artifact": name, **dict(record)}
                      for record in (meta.get("failed_attempts") or []))

    clinvar_bytes = raw_payload_of(source, "clinvar")
    emit("clinvar_raw.tsv", "immutable ClinVar payload as retrieved",
         clinvar_bytes if clinvar_bytes is not None
         else raw_tsv_text(raw_rows).encode("utf-8"),
         clinvar_meta, clinvar_bytes is not None)

    uniprot_bytes = raw_payload_of(source, "uniprot")
    emit("uniprot_canonical.fasta", "UniProt canonical isoform sequence",
         uniprot_bytes if uniprot_bytes is not None
         else resolution.fasta().encode("utf-8"),
         sequence_meta, uniprot_bytes is not None)

    structure_bytes = raw_payload_of(source, "structure")
    if structure_bytes is None and model.raw_cif_text is not None:
        structure_bytes = model.raw_cif_text.encode("utf-8")
    if structure_bytes is None:
        from ..structure.cif import render_mmcif

        structure_bytes = render_mmcif(model).encode("utf-8")
    emit("alphafold_model.cif", "AlphaFold model as served", structure_bytes,
         model_meta, model.raw_cif_text is not None)

    run.note(STAGE_RAW, "retrieval_log.json",
             write_json(stage_dir / "retrieval_log.json", {
                 "run_id": run.ctx.run_id, "gene": run.gene, "agent_owner": AGENT,
                 "network_used": False if run.ctx.synthetic else None,
                 "offline_mode": bool(run.ctx.config.get_optional("execution.offline_mode",
                                                                 False)),
                 "gene_resolution": resolution.as_dict(),
                 "n_records_retrieved": len(raw_rows),
                 "retrievals": log,
                 "failed_retrievals": failed,
                 "substitution_policy": (
                     "A failed retrieval is never replaced by a cached, hand-made or "
                     "previous-run file."),
             }, run_root=run.ctx.run_root),
             "URLs, queries, release dates, status, retries")


# --- 02_CLINVAR -------------------------------------------------------------

def _non_missense_true_cause(detail: str) -> str:
    """Workflow v2 §1 — the fine-grained true cause behind a
    ``non_missense_variant_type`` or ``unparseable_hgvs`` exclusion, for
    REPORTING only (``exclusion_reason_detail_counts``, Lead ruling 2026-08-18).

    This groups the SAME ``rec.change.reason`` value the record was already
    excluded with (``hotspot3d.data.hgvs.parse_protein_hgvs`` /
    ``classify.non_missense_exclusion_reason``); it never reclassifies a record
    and never changes ``exclusion_reason``. It exists because the top-level
    reason alone still groups several distinct causes together — e.g.
    ``non_missense_variant_type`` covers synonymous, nonsense, frameshift AND
    other in-frame indel/dup/ext changes — and the per-reason breakdown v2 §1
    requires needs that granularity visible.
    """
    text = str(detail or "")
    if text == "synonymous":
        return "synonymous"
    if text == "nonsense":
        return "nonsense"
    if text.startswith("not_a_substitution:fs"):
        return "frameshift"
    if text.startswith("not_a_substitution:"):
        return "indel_or_other_non_substitution"
    if text == "reference_is_termination":
        return "stop_loss_or_extension"
    if text.startswith("non_standard_reference_aa:") or text.startswith(
            "non_standard_alternate_aa:"):
        return "non_standard_amino_acid"
    if text == "non_positive_position":
        return "invalid_position"
    if text in ("hgvs_p_pattern_unmatched", "not_protein_level_hgvs",
                "empty_protein_hgvs", "no_protein_hgvs"):
        return "true_parse_failure"
    return "other_non_missense"                                  # pragma: no cover


def _write_clinvar_outputs(run: _Run, evals, residue_rows, star_rows, counts,
                           conservation, policy, resolution, clinvar_meta,
                           star_audit, excluded_no_ca, plddt_summary,
                           sensitivity_threshold, profile, mapping,
                           compound_cost: dict) -> None:
    stage_dir = run.dir(STAGE_CLINVAR)
    cfg = run.ctx.config

    run.note(STAGE_CLINVAR, "variants_missense_all.tsv",
             write_tsv(stage_dir / "variants_missense_all.tsv",
                       [e.stratum1_row() for e in evals], STRATUM1_COLUMNS,
                       run_root=run.ctx.run_root),
             "stratum 1 — every retrieved record, unfiltered")

    run.note(STAGE_CLINVAR, "variants_residue_level.tsv",
             write_tsv(stage_dir / "variants_residue_level.tsv",
                       [r.stratum2_row() for r in residue_rows], STRATUM2_COLUMNS,
                       run_root=run.ctx.run_root),
             "stratum 2 — residue-level collapse")

    run.note(STAGE_CLINVAR, "variants_excluded_from_primary.tsv",
             write_tsv(stage_dir / "variants_excluded_from_primary.tsv",
                       [e.stratum3_row() for e in evals if not e.eligible_primary],
                       STRATUM3_COLUMNS, run_root=run.ctx.run_root),
             "stratum 3 — exclusions with reasons")

    run.note(STAGE_CLINVAR, "review_star_distribution.tsv",
             write_tsv(stage_dir / "review_star_distribution.tsv", star_rows,
                       STRATUM4_COLUMNS, run_root=run.ctx.run_root),
             "stratum 4 — stars, described never used")

    conflicts = [r for r in residue_rows if r.residue_class == CLASS_CONFLICT]
    run.note(STAGE_CLINVAR, "residue_class_conflicts.tsv",
             write_tsv(stage_dir / "residue_class_conflicts.tsv",
                       [r.conflict_row() for r in conflicts], CONFLICT_COLUMNS,
                       run_root=run.ctx.run_root),
             "F3 RESIDUE_CLASS_CONFLICT residues")

    reason_counts: dict[str, int] = {}
    for ev in evals:
        if not ev.eligible_primary:
            reason_counts[ev.exclusion_reason] = reason_counts.get(ev.exclusion_reason, 0) + 1
    # [NEW — Workflow v2 §1, Lead ruling 2026-08-18] the per-reason breakdown at
    # the granularity of the true cause, nested under the two reasons that each
    # still cover more than one HGVS.parse_protein_hgvs() outcome.
    exclusion_reason_detail_counts: dict[str, dict[str, int]] = {}
    for ev in evals:
        if ev.exclusion_reason in (REASON_UNPARSEABLE_HGVS, REASON_NON_MISSENSE_VARIANT_TYPE):
            bucket = _non_missense_true_cause(ev.exclusion_detail)
            per_reason = exclusion_reason_detail_counts.setdefault(ev.exclusion_reason, {})
            per_reason[bucket] = per_reason.get(bucket, 0) + 1
    sig_counts: dict[str, int] = {}
    for ev in evals:
        key = ev.significance.significance_class
        sig_counts[key] = sig_counts.get(key, 0) + 1
    inexact = [e.record_id for e in evals
               if e.significance.eligible_primary and not e.significance.exact_match]

    # Read with get_optional: its ABSENCE is itself meaningful (no configured
    # minimum). REQUIRED_CONFIG_KEYS still demands the key be present, so an
    # explicit null is a recorded choice rather than an accidental omission.
    minimum = cfg.get_optional("clinvar.min_classified_residues")
    low_confidence_cohort = _cohort_below_threshold(profile, mapping.cohort,
                                                    sensitivity_threshold)

    run.note(STAGE_CLINVAR, "cohort_summary.json",
             write_json(stage_dir / "cohort_summary.json", {
                 "run_id": run.ctx.run_id, "gene": run.gene, "agent_owner": AGENT,
                 "config_sha256": cfg.sha256,
                 "uniprot_acc": resolution.uniprot_acc,
                 "mane_transcript": resolution.mane_transcript,
                 "mane_source": resolution.mane_source,
                 "mane_fallback_used": resolution.mane_fallback_used,
                 "clinvar_release": clinvar_meta.get("release_date", "NA"),
                 "clinvar_query": clinvar_meta.get("query", "NA"),
                 # [NEW — Workflow v2 §1] N_B reported at the top level, not only
                 # nested in residue_counts below: it sets the resolution of any
                 # P/LP-vs-B/LB contrast and must be evaluated by Stage B's power
                 # certificate (§5.4) before any negative result is issued.
                 "N_B": counts.N_B, "N_P": counts.N_P, "N": counts.N,
                 "benign_cohort_adequacy_note": (
                     "N_B (benign/likely-benign residues) determines the resolution of "
                     "any test contrasting pathogenic against benign residues. "
                     "hotspot-statistics must compute the §5.4 power certificate against "
                     "this N_B before reporting any negative result as COMPLETED_NEGATIVE."
                 ),
                 "strata": {
                     "stratum_1_all_retrieved": len(evals),
                     "stratum_2_residues": len(residue_rows),
                     "stratum_3_excluded_records": conservation.n_excluded,
                     "stratum_4_star_distribution_rows": len(star_rows),
                 },
                 "record_conservation": conservation.as_dict(),
                 "exclusion_reason_counts": dict(sorted(reason_counts.items())),
                 "exclusion_reason_enum": list(cfg.get("clinvar.exclusion_reason_enum")),
                 "reason_precedence": (
                     "A record is reported under the first fact that made it unusable: "
                     "unparseable_hgvs/non_missense_variant_type, non_canonical_transcript, "
                     "ref_aa_mismatch, significance, residue_class_conflict, "
                     "no_ca_coordinate."),
                 # [NEW — Workflow v2 §1, Lead ruling 2026-08-18] the true-cause
                 # breakdown WITHIN each of the two reasons that still cover more
                 # than one parse_protein_hgvs() outcome: non_missense_variant_type
                 # (synonymous / nonsense / frameshift / other in-frame indel-dup-ext
                 # / stop-loss) and unparseable_hgvs (true_parse_failure /
                 # non_standard_amino_acid / invalid_position). No record's
                 # exclusion_reason or eligibility is changed to produce this.
                 "exclusion_reason_detail_counts": {
                     reason: dict(sorted(detail.items()))
                     for reason, detail in sorted(exclusion_reason_detail_counts.items())
                 },
                 "significance_class_counts": dict(sorted(sig_counts.items())),
                 "significance_policy": policy.as_dict(),
                 "n_binary_records_matched_component_wise": len(inexact),
                 # [NEW — Workflow v2 §1] counterfactual cost of the frozen
                 # compound-significance policy — never a second inclusion path.
                 "compound_significance_policy_cost": compound_cost,
                 "residue_counts": {
                     "n_residues_total": counts.n_residues_total,
                     "N": counts.N, "N_P": counts.N_P, "N_B": counts.N_B,
                     "n_conflict": counts.n_conflict,
                     "n_excluded_no_ca_coordinate": len(excluded_no_ca),
                     "residues_excluded_no_ca_coordinate": excluded_no_ca,
                 },
                 "review_star_distribution_overall": star_summary(evals),
                 "star_filter_audit": star_audit,
                 "review_stars_used_for_inclusion": False,
                 "vus_or_conflicting_admitted_to_binary_classes": False,
                 "plddt_used_as_filter": False,
                 "cohort_plddt": {
                     "sensitivity_threshold": sensitivity_threshold,
                     "n_cohort_below_threshold": low_confidence_cohort,
                     "note": "Reported only. F2 — pLDDT never removes a residue.",
                 },
                 "min_classified_residues_configured": minimum,
                 "underpowered": (None if minimum is None else counts.N < int(minimum)),
                 "structure_plddt_summary": plddt_summary.as_dict(),
             }, run_root=run.ctx.run_root),
             "cohort headline counts and audits")

    _emit_figures(run, STAGE_CLINVAR, {
        "figures/significance_composition.png": lambda p: figs.figure_significance_composition(
            [e.stratum1_row() for e in evals], p, run.gene),
        "figures/review_star_distribution.png": lambda p: figs.figure_star_distribution(
            star_rows, p, run.gene),
        "figures/cohort_along_sequence.png": lambda p: figs.figure_cohort_along_sequence(
            [r.stratum2_row() for r in residue_rows], p, run.gene,
            resolution.sequence_length),
    })


def _cohort_below_threshold(profile, cohort, threshold: float) -> int:
    lookup = {row["residue_index"]: row["plddt"] for row in profile}
    return sum(1 for row in cohort
               if lookup.get(row["residue_index"]) is not None
               and lookup[row["residue_index"]] < threshold)


# --- 03_STRUCTURE_QC --------------------------------------------------------

def _write_structure_outputs(run: _Run, model, qc, coordinate_rows, universe_rows,
                             mapping, profile, regions, plddt_summary, bands,
                             model_meta, resolution, intra_subunit_only: bool,
                             allow_monomer_oligomer: bool) -> None:
    stage_dir = run.dir(STAGE_STRUCT)

    cif_record = write_structure_with_plddt(model, stage_dir / "structures" /
                                            "structure_with_plddt.cif")
    run.note(STAGE_STRUCT, "structures/structure_with_plddt.cif",
             stage_dir / "structures" / "structure_with_plddt.cif",
             "model with B-factor = pLDDT")

    run.note(STAGE_STRUCT, "structure_qc.json",
             write_json(stage_dir / "structure_qc.json", {
                 "run_id": run.ctx.run_id, "gene": run.gene, "agent_owner": AGENT,
                 **qc.as_dict(),
                 "model_metadata": model_meta,
                 "uniprot_acc": resolution.uniprot_acc,
                 "uniprot_entry_version": resolution.uniprot_entry_version,
                 "canonical_isoform": resolution.canonical_isoform,
                 "structure_export": cif_record,
                 "plddt_summary": plddt_summary.as_dict(),
                 "plddt_bands": bands,
                 "plddt_band_labels": band_labels(bands),
                 "representation": {
                     "coordinate_atom": "CA", "distance_metric": "euclidean",
                     "emitted_unused": ["CB", "SC_CENTROID"],
                     "frozen_decision": "F1",
                 },
                 # Workflow v2 §2 — determined from UniProt's own subunit comment
                 # only (structure.oligomeric_state_source); never a ComplexPortal
                 # lookup, which is explicitly out of scope. "unknown" means no
                 # SUBUNIT comment was present/parseable at all — a distinct,
                 # honest state from a confirmed "monomer" reading.
                 "oligomeric_assembly": {
                     "oligomeric_state": resolution.oligomeric_state,
                     "oligomeric_state_source": resolution.oligomeric_state_source,
                     "subunit_comment_text": resolution.oligomeric_state_evidence,
                     "matched_phrase": resolution.oligomeric_state_matched_phrase,
                     "allow_monomer_for_obligate_oligomer": allow_monomer_oligomer,
                     "multimer_model_used": False,
                     "intra_subunit_only": intra_subunit_only,
                     "note": (
                         "No multimer-fetching capability exists in this pipeline; "
                         "every structural model is the monomeric AlphaFold prediction. "
                         "When oligomeric_state == 'oligomeric', proceeding is the "
                         "declared exception of config/pipeline.yaml "
                         "structure.allow_monomer_for_obligate_oligomer, and every "
                         "downstream hotspot/footprint must be read as intra-subunit "
                         "only."
                     ),
                 },
             }, run_root=run.ctx.run_root),
             "structural QC checks and verdict")

    run.note(STAGE_STRUCT, "residue_coordinates.tsv",
             write_tsv(stage_dir / "residue_coordinates.tsv", coordinate_rows,
                       RESIDUE_COORDINATE_COLUMNS, run_root=run.ctx.run_root),
             "CA/CB/side-chain coordinates and pLDDT")

    run.note(STAGE_STRUCT, "positional_universe.tsv",
             write_tsv(stage_dir / "positional_universe.tsv", universe_rows,
                       POSITIONAL_UNIVERSE_COLUMNS, run_root=run.ctx.run_root),
             "U_struct — every residue with a usable CA")

    run.note(STAGE_STRUCT, "classified_cohort.tsv",
             write_tsv(stage_dir / "classified_cohort.tsv", mapping.cohort,
                       CLASSIFIED_COHORT_COLUMNS, run_root=run.ctx.run_root),
             "L — the binary P/LP vs B/LB cohort")

    run.note(STAGE_STRUCT, "plddt_profile.tsv",
             write_tsv(stage_dir / "plddt_profile.tsv", profile, PLDDT_PROFILE_COLUMNS,
                       run_root=run.ctx.run_root),
             "per-residue pLDDT and band")

    run.note(STAGE_STRUCT, "plddt_regions.tsv",
             write_tsv(stage_dir / "plddt_regions.tsv", regions, PLDDT_REGION_COLUMNS,
                       run_root=run.ctx.run_root),
             "contiguous low-confidence regions (never excluded)")

    run.note(STAGE_STRUCT, "mapping_report.json",
             write_json(stage_dir / "mapping_report.json", {
                 "run_id": run.ctx.run_id, "gene": run.gene, "agent_owner": AGENT,
                 "mane_transcript": resolution.mane_transcript,
                 "uniprot_canonical_length": resolution.sequence_length,
                 **mapping.report,
                 "alignment": qc.alignment.as_dict() if qc.alignment else None,
                 "numbering_offset_applied": False,
                 "silent_numbering_repair": False,
                 "coordinate_imputation": False,
                 "class_imputation": False,
                 "cohort_balancing": False,
             }, run_root=run.ctx.run_root),
             "ClinVar-to-structure mapping and unmapped residues")

    _emit_figures(run, STAGE_STRUCT, {
        "figures/plddt_profile.png": lambda p: figs.figure_plddt_profile(
            profile, p, run.gene, bands,
            [row["residue_index"] for row in mapping.cohort]),
        "figures/plddt_bands.png": lambda p: figs.figure_plddt_bands(
            plddt_summary.as_dict(), p, run.gene),
        "figures/universe_projection.png": lambda p: figs.figure_universe_projection(
            universe_rows, mapping.cohort, p, run.gene),
    })


def _emit_figures(run: _Run, stage: str, builders: dict) -> None:
    stage_dir = run.dir(stage)
    for name, build in builders.items():
        path = stage_dir / name
        description = dict(EXPECTED_FILES[stage]).get(name, "figure")
        try:
            build(path)
            run.note(stage, name, path, description)
        except Exception as exc:
            # A figure never blocks a stage: absence is recorded explicitly (P3)
            # and the numeric outputs it illustrates are already written.
            run.skip(stage, name, description, f"{type(exc).__name__}: {exc}")
            if not run.warnings[stage].has("RENDERING_UNAVAILABLE"):
                run.warn(stage, "RENDERING_UNAVAILABLE",
                         f"Figure rendering unavailable: {type(exc).__name__}: {exc}",
                         potential_consequence="Numeric outputs are unaffected.",
                         recommended_action="Install/repair matplotlib and regenerate.",
                         affected_output=f"FULL_RESULTS/{stage}/figures/")


# --- QC gates and warnings --------------------------------------------------

def _raise_on_structural_failure(run: _Run, qc, mapping, conservation, counts) -> None:
    blocking = [c for c in qc.failures if c.severity == "BLOCKING"]
    if blocking:
        detail = "; ".join(f"{c.check_id} {c.description} ({c.detail})" for c in blocking)
        run.warn(STAGE_STRUCT, "STRUCTURE_MAPPING_FAILURE",
                 f"Blocking structural QC failure: {detail}",
                 potential_consequence="Coordinates cannot be trusted to correspond "
                                       "to the analysed sequence.",
                 recommended_action="Escalate to the Lead; no silent repair is permitted.",
                 affected_output="FULL_RESULTS/03_STRUCTURE_QC/structure_qc.json")
        raise BlockedError(f"STRUCTURE_MAPPING_FAILURE — {detail}")

    non_blocking = [c for c in qc.failures if c.severity != "BLOCKING"]
    if non_blocking:
        run.warn(STAGE_STRUCT, "STRUCTURAL_CONFIDENCE_SENSITIVE",
                 "Non-blocking structural QC findings: "
                 + "; ".join(f"{c.check_id} ({c.detail})" for c in non_blocking),
                 potential_consequence="Some residues cannot enter U_struct.",
                 recommended_action="Read structure_qc.json; nothing was repaired "
                                    "or trimmed to make these pass.",
                 affected_output="FULL_RESULTS/03_STRUCTURE_QC/structure_qc.json")

    if not conservation.conserved:
        raise BlockedError(
            f"BLOCKED — record-count conservation failed: {conservation.as_dict()}. "
            f"Every retrieved record must land in exactly one stratum.")

    if not mapping.report["cohort_subset_of_universe"]:
        run.warn(STAGE_STRUCT, "STRUCTURE_MAPPING_FAILURE",
                 "A classified residue is absent from U_struct.",
                 affected_output="FULL_RESULTS/03_STRUCTURE_QC/classified_cohort.tsv")
        raise BlockedError("STRUCTURE_MAPPING_FAILURE — L is not a subset of U_struct.")

    if mapping.M == 0:
        raise BlockedError("STRUCTURE_MAPPING_FAILURE — U_struct is empty.")


def _post_cohort_warnings(run: _Run, counts, mapping, plddt_summary, profile,
                          threshold: float, negative: dict | None) -> None:
    if counts.n_conflict:
        run.warn(STAGE_CLINVAR, "RESIDUE_CLASS_CONFLICT",
                 f"{counts.n_conflict} residue(s) carry both P/LP and B/LB evidence and "
                 f"are excluded from the primary comparison with all source records kept.",
                 potential_consequence="The binary cohort is smaller than the record count "
                                       "suggests.",
                 recommended_action="Read residue_class_conflicts.tsv; F3 forbids assigning "
                                    "these residues to either class.",
                 affected_output="FULL_RESULTS/02_CLINVAR/residue_class_conflicts.tsv",
                 status="ACCEPTED_BY_DESIGN")

    if mapping.unmapped:
        run.warn(STAGE_STRUCT, "STRUCTURAL_CONFIDENCE_SENSITIVE",
                 f"{len(mapping.unmapped)} cohort residue(s) have no usable CA and are "
                 f"recorded in mapping_report.json with a reason.",
                 potential_consequence="Those residues cannot enter U_struct or L.",
                 recommended_action="Read mapping_report.json; they were not silently dropped.",
                 affected_output="FULL_RESULTS/03_STRUCTURE_QC/mapping_report.json")

    # Two majority rules, both on the same bright line so no new constant is
    # invented: most of the cohort, or most of the model, below the configured
    # sensitivity threshold. Reporting only — F2 forbids filtering on pLDDT.
    below_cohort = _cohort_below_threshold(profile, mapping.cohort, threshold)
    below_model = sum(1 for row in profile
                      if row["plddt"] is not None and row["plddt"] < threshold)
    if mapping.N and below_cohort * 2 >= mapping.N:
        run.warn(STAGE_STRUCT, "STRUCTURAL_CONFIDENCE_SENSITIVE",
                 f"{below_cohort} of {mapping.N} classified residues sit below "
                 f"pLDDT {threshold:g}.",
                 potential_consequence="Spatial conclusions may rest on low-confidence "
                                       "coordinates.",
                 recommended_action="Read the pLDDT >= 70 sensitivity analysis downstream; "
                                    "F2 forbids filtering on pLDDT here.",
                 affected_output="FULL_RESULTS/03_STRUCTURE_QC/plddt_regions.tsv")
    elif profile and below_model * 2 >= len(profile):
        run.warn(STAGE_STRUCT, "STRUCTURAL_CONFIDENCE_SENSITIVE",
                 f"{below_model} of {len(profile)} modelled residues sit below "
                 f"pLDDT {threshold:g}.",
                 potential_consequence="The model as a whole is low confidence, even "
                                       "where the cohort is not.",
                 recommended_action="Read plddt_regions.tsv; no residue was excluded.",
                 affected_output="FULL_RESULTS/03_STRUCTURE_QC/plddt_regions.tsv")

    minimum = run.ctx.config.get_optional("clinvar.min_classified_residues")
    if minimum is not None and counts.N < int(minimum) and negative is None:
        run.warn(STAGE_CLINVAR, "UNDERPOWERED_COHORT",
                 f"N={counts.N} is below the configured minimum {int(minimum)}.",
                 potential_consequence="Downstream tests may have little power.",
                 recommended_action="Escalate to the Lead; this is not an error to "
                                    "engineer away.",
                 affected_output="FULL_RESULTS/02_CLINVAR/cohort_summary.json")


def _qc_status(run: _Run) -> str:
    """Decide ``qc_status`` from the frozen severity vocabulary alone (IX.7).

    Any BLOCKING warning -> FAIL. Any other ADVISORY-or-above warning ->
    PASS_WITH_WARNINGS. INFO-level notes do not degrade the verdict.

    Stage A's insufficient-cohort branch lands in FAIL by that rule, and that is
    correct. Two different negatives must not be conflated:

      (A) the analysis ran and found nothing — every center tested, none survived
          BH-FDR. The stage did its job; the answer is "no effect". PASS.
      (B) the input cannot support the analysis — a class is empty or |L| < 2, so
          nothing was ever tested. Reporting this as PASS would assert "we looked
          and found nothing" when in fact we could not look.

    Stage A's own negative is always (B): an empty class BLOCKS the A->B handoff.
    ``outcome_type`` still says SCIENTIFIC_NEGATIVE, because the cause is the
    evidence base rather than a fault in the code, but the run is BLOCKED and
    downstream must not proceed.
    """
    warnings = run.all_warnings()
    if any(w.severity == Severity.BLOCKING for w in warnings):
        return QC_FAIL
    notable = [w for w in warnings if w.severity != Severity.INFO]
    return QC_PASS_WITH_WARNINGS if notable else QC_PASS


def _stage_outcome(negative: dict | None) -> tuple[Status, OutcomeType]:
    """BLOCKED, not COMPLETED_NEGATIVE, when the input cannot support the analysis.

    Stage A never reaches a case-(A) negative — it runs no test whose answer could
    be "no effect". Its only negative is "the evidence base cannot be analysed",
    which blocks the A->B handoff. ``outcome_type`` stays SCIENTIFIC_NEGATIVE
    because the cause is ClinVar's content, not a fault in the code, so this is
    still not a TECHNICAL_FAILURE.
    """
    if negative is None:
        return Status.COMPLETED, OutcomeType.COMPLETED
    return Status.BLOCKED, OutcomeType.SCIENTIFIC_NEGATIVE


def _stage_reason(stage: str, counts, mapping, negative: dict | None) -> str:
    if negative is not None:
        # NOT_RUN.txt prints this verbatim, so it must name the cause plainly and
        # must not imply that anything was tested.
        return (f"{negative['detail']} No spatial analysis was attempted: the "
                f"classified cohort could not be assembled, so nothing was tested "
                f"and no statement about clustering in this gene is available "
                f"from this run.")
    return (f"Stage A completed: M={mapping.M}, N={counts.N}, N_P={counts.N_P}, "
            f"N_B={counts.N_B}, n_conflict={counts.n_conflict}.")


def _recommended_action(negative: dict | None) -> str:
    if negative is None:
        return ("Verify every manifest hash, then proceed to hotspot-statistics via "
                "handoff_01.json.")
    return ("Downstream stages must NOT run against this handoff. Report the cohort "
            "as unanalysable and say so as such — this is not a finding of 'no "
            "clustering'. Do not admit VUS or conflicting records, do not reassign "
            "conflict residues, do not extend to a non-canonical transcript and do "
            "not introduce a review-star filter in order to reach a usable cohort.")


# --- handoff, manifests, provenance, report ---------------------------------

def _write_handoff(run: _Run, resolution, clinvar_meta, model, counts, mapping,
                   qc_status: str, negative: dict | None, intra_subunit_only: bool,
                   allow_monomer_oligomer: bool) -> Handoff:
    ctx = run.ctx
    manifest = manifest_for(run.paths(), root=ctx.run_root)
    handoff = Handoff(
        name="handoff_01", run_id=ctx.run_id, config_sha256=ctx.config.sha256,
        qc_status=qc_status, manifest=manifest, negative_result=negative,
        payload={
            "gene": run.gene,
            "uniprot_acc": resolution.uniprot_acc,
            "uniprot_entry_version": resolution.uniprot_entry_version,
            "mane_transcript": resolution.mane_transcript,
            "mane_fallback_used": resolution.mane_fallback_used,
            "clinvar_release": clinvar_meta.get("release_date", "NA"),
            "alphafold_model_version": model.model_version,
            "alphafold_model_id": model.model_id(),
            "alphafold_db_version": model.db_version,
            "M": mapping.M, "N": counts.N, "N_P": counts.N_P, "N_B": counts.N_B,
            "n_conflict": counts.n_conflict,
            "D_max": mapping.report["D_max"],
            "n_unmapped": len(mapping.unmapped),
            "code_version": ctx.code_version,
            "synthetic_input_mode": bool(ctx.synthetic),
            "derived_seeds": None,
            "derived_seeds_reason": "not_applicable_deterministic_stage",
            "stage_b_permitted_columns": list(STAGE_B_PERMITTED_COLUMNS),
            "forbidden_downstream": list(FORBIDDEN_DOWNSTREAM_COLUMNS),
            "stage_b_universe_columns": list(POSITIONAL_UNIVERSE_COLUMNS),
            "primary_inputs": {
                "positional_universe": "FULL_RESULTS/03_STRUCTURE_QC/positional_universe.tsv",
                "classified_cohort": "FULL_RESULTS/03_STRUCTURE_QC/classified_cohort.tsv",
                "residue_coordinates": "FULL_RESULTS/03_STRUCTURE_QC/residue_coordinates.tsv",
            },
            "sensitivity_channel": {
                "declared": True,
                "purpose": ("METHOD_SPEC II.7 >=1* / >=2* sensitivity analyses only — "
                            "post-primary and explicitly non-redefining."),
                "file": "FULL_RESULTS/02_CLINVAR/variants_residue_level.tsv",
                "columns": ["residue_index", "star_levels", "max_star"],
                "may_redefine_primary": False,
                "note": ("These columns are forbidden to the primary path. Reading them "
                         "through this channel never changes L, N_P or N_B."),
            },
            "consumer_must_verify": [
                "every manifest SHA-256",
                "both classes non-empty",
                "every L residue present in U_struct with a usable CA",
                "qc_status != FAIL",
            ],
            # Workflow v2 §2 — declared from UniProt's own subunit annotation
            # (structure.oligomeric_state_source), never a ComplexPortal lookup.
            # No multimer model was fetched or used; a monomer AlphaFold model
            # backs every coordinate in this handoff regardless of the state below.
            "oligomeric_assembly": {
                "oligomeric_state": resolution.oligomeric_state,
                "oligomeric_state_source": resolution.oligomeric_state_source,
                "allow_monomer_for_obligate_oligomer": allow_monomer_oligomer,
                "intra_subunit_only": intra_subunit_only,
                "note": ("hotspot-statistics and footprint-robustness must propagate "
                         "intra_subunit_only=TRUE into every hotspot/footprint label "
                         "when set; it means the biological assembly is a declared "
                         "obligate oligomer analysed on its monomer, so no inter-subunit "
                         "neighbourhood is represented anywhere in U_struct or L."),
            },
        },
    )
    path = write_json(ctx.handoff_path("handoff_01"), handoff.as_dict(),
                      run_root=ctx.run_root)
    run.note(STAGE_STRUCT, "handoff_01.json", path, "handoff to hotspot-statistics")
    return handoff


def _write_manifests(run: _Run) -> None:
    """stage_manifest.tsv lists EXPECTED files, so absence is explicit (P3)."""
    created = {(e.stage, e.name): e for e in run.emitted if e.status == "CREATED"}
    skipped = {(e.stage, e.name): e for e in run.emitted if e.status == "NOT_CREATED"}

    for stage in STAGES:
        stage_dir = run.dir(stage)
        rows = []
        for name, description in EXPECTED_FILES[stage]:
            emitted = created.get((stage, name))
            if emitted is not None and emitted.path.is_file():
                rows.append({
                    "path": run.rel(emitted.path), "stage": stage,
                    "description": description, "status": "CREATED",
                    "sha256": sha256_file(emitted.path),
                    "size_bytes": emitted.path.stat().st_size, "reason": "NA",
                })
            elif name == "stage_manifest.tsv":
                rows.append({
                    "path": run.rel(stage_dir / name), "stage": stage,
                    "description": description, "status": "DEFERRED",
                    "sha256": None, "size_bytes": None,
                    "reason": "self-referential; hashed in the run-level manifest",
                })
            else:
                reason = (skipped[(stage, name)].reason if (stage, name) in skipped
                          else "not produced by this run")
                rows.append({
                    "path": run.rel(stage_dir / name), "stage": stage,
                    "description": description, "status": "NOT_CREATED",
                    "sha256": None, "size_bytes": None, "reason": reason,
                })
        write_tsv(stage_dir / "stage_manifest.tsv", rows, MANIFEST_COLUMNS,
                  run_root=run.ctx.run_root)


def _write_provenance(run: _Run, resolution, clinvar_meta, sequence_meta, model_meta,
                      counts, mapping, star_audit, conservation,
                      n_records_retrieved: int) -> None:
    ctx = run.ctx
    freeze = uv_pip_freeze()
    shared = {
        "derived_seeds": None,
        "derived_seeds_reason": "not_applicable_deterministic_stage",
        "deterministic_stage": True,
        "uv_pip_freeze": freeze,
        "step_wall_seconds": run.step_seconds,
        "config_path": str(ctx.config.path),
        "information_barrier": {
            "stage": "A",
            "downstream_directories_read": [],
            "biology_used_in_any_decision": False,
        },
    }
    warnings_rows = [w.as_row() for w in run.all_warnings()]

    stage_inputs = {
        STAGE_RAW: {"clinvar": clinvar_meta, "uniprot": sequence_meta,
                    "alphafold": model_meta},
        STAGE_CLINVAR: {"clinvar_release": clinvar_meta.get("release_date", "NA"),
                        "clinvar_query": clinvar_meta.get("query", "NA"),
                        "gene_resolution": resolution.as_dict()},
        STAGE_STRUCT: {"alphafold": model_meta,
                       "uniprot_canonical_length": resolution.sequence_length},
    }
    stage_outputs = {
        stage: {e.name: run.rel(e.path) for e in run.emitted
                if e.stage == stage and e.status == "CREATED"}
        for stage in STAGES
    }
    stage_parameters = {
        STAGE_RAW: {"structure.model_template": ctx.config.get("structure.model_template"),
                    "clinvar.source": ctx.config.get("clinvar.source")},
        STAGE_CLINVAR: {
            "clinvar.pathogenic_labels": ctx.config.get("clinvar.pathogenic_labels"),
            "clinvar.benign_labels": ctx.config.get("clinvar.benign_labels"),
            "clinvar.review_star_filter": None,
            "clinvar.use_review_stars_for_inclusion": False,
            "transcript.reference": ctx.config.get("transcript.reference"),
            "transcript.require_ref_aa_match": ctx.config.get(
                "transcript.require_ref_aa_match"),
        },
        STAGE_STRUCT: {
            "plddt.bands": ctx.config.get("plddt.bands"),
            "plddt.primary_filtering_enabled": False,
            "plddt.pae_threshold": None,
            "representation.coordinate_atom": "CA",
            "representation.distance_metric": "euclidean",
        },
    }
    stage_extra = {
        STAGE_RAW: {"n_records_retrieved": n_records_retrieved,
                    "payloads_immutable": True},
        STAGE_CLINVAR: {"star_filter_audit": star_audit,
                        "record_conservation": conservation.as_dict(),
                        "cohort": {"N": counts.N, "N_P": counts.N_P, "N_B": counts.N_B,
                                   "n_conflict": counts.n_conflict}},
        STAGE_STRUCT: {"M": mapping.M, "N": mapping.N,
                       "n_unmapped": len(mapping.unmapped),
                       "D_max": mapping.report["D_max"]},
    }

    for stage in STAGES:
        ctx.record_provenance(
            stage, AGENT,
            inputs=stage_inputs[stage], outputs=stage_outputs[stage],
            parameters=stage_parameters[stage],
            commands=run.commands,
            warnings=[w for w in warnings_rows if w["stage"] == WARN_SUFFIX[stage]],
            extra={**shared, **stage_extra[stage]},
        )


def _write_report(run: _Run, resolution, clinvar_meta, model, qc, counts, mapping,
                  conservation, plddt_summary, star_audit, negative, qc_status,
                  excluded_no_ca, status, star_tally: dict, compound_cost: dict) -> None:
    stage_dir = run.dir(STAGE_STRUCT)
    bands = plddt_summary.band_counts
    stars = ", ".join(f"{k}*: {v}" for k, v in star_tally.items()) or "none"
    failed = [f"{c.check_id} {c.description} — {c.detail}" for c in qc.failures] or ["none"]
    warnings_lines = [
        f"- {w.severity.value} {w.warning_code}: {w.message}" for w in run.all_warnings()
    ] or ["- none"]

    text = f"""# Stage A report — {run.gene}

RUN_ID `{run.ctx.run_id}` · agent `{AGENT}` · config_sha256 `{run.ctx.config.sha256[:16]}`

STATUS: {status.value} (qc_status = {qc_status})

COHORT ADEQUACY [NEW — Workflow v2 §1]: N_P = {counts.N_P} · N_B = {counts.N_B} · N = {counts.N}
- N_B is reported here, prominently, because it sets the resolution of any P/LP-vs-B/LB
  contrast; hotspot-statistics must evaluate its §5.4 power certificate against this N_B
  before any negative result may be reported as COMPLETED_NEGATIVE.

OLIGOMERIC ASSEMBLY [NEW — Workflow v2 §2]: state = {resolution.oligomeric_state}
  (source: {resolution.oligomeric_state_source}; evidence: {resolution.oligomeric_state_evidence!r})
- No multimer model was fetched or used (no such client exists); every coordinate in this
  run is from the monomeric AlphaFold model regardless of this state.
- If state == oligomeric, this is the declared exception of
  structure.allow_monomer_for_obligate_oligomer, and every downstream hotspot/footprint
  must be labelled intra-subunit only. See structure_qc.json > oligomeric_assembly.

INPUTS USED:
- ClinVar {clinvar_meta.get('source', 'NA')} release {clinvar_meta.get('release_date', 'NA')}
  (query: `{clinvar_meta.get('query', 'NA')}`)
- UniProt {resolution.uniprot_acc} entry version {resolution.uniprot_entry_version},
  canonical isoform {resolution.canonical_isoform}, {resolution.sequence_length} aa
- Transcript {resolution.mane_transcript} ({resolution.mane_source},
  fallback_used={resolution.mane_fallback_used})
- AlphaFold {model.model_id()} model version {model.model_version},
  database version {model.db_version}, fragments {model.n_fragments}
- Synthetic input mode: {bool(run.ctx.synthetic)}

METHODS EXECUTED:
1. Gene resolved to exactly one UniProt accession, canonical isoform and transcript.
2. Every ClinVar missense record retrieved at every review star level (F12). No star
   filter exists in any Stage A code path; this is asserted mechanically over
   {star_audit['n_modules_scanned']} modules by an AST audit, not by prose.
3. Significance normalized to P/LP and B/LB only; every other category kept in strata
   1 and 3 with a reason from the frozen enum.
4. Protein HGVS parsed, reference AA verified against the canonical sequence, non-MANE
   records remapped only when position and reference AA both agree.
5. Records collapsed to residues; F3 applied to residues carrying both classes.
6. AlphaFold model retrieved and QC'd; numbering offset established by explicit
   full-sequence alignment and reported, never silently applied.
7. pLDDT extracted, banded and profiled for contiguous low-confidence regions —
   recorded, never used to exclude (F2).
8. U_struct and L emitted on CA coordinates (F1); CB and side-chain centroid columns
   emitted and marked unused.

OUTPUTS GENERATED:
{chr(10).join('- ' + run.rel(e.path) for e in run.emitted if e.status == 'CREATED')}
- FULL_RESULTS/03_STRUCTURE_QC/stage_a_report.md (this file)
- FULL_RESULTS/03_STRUCTURE_QC/handoff_01.json
- FULL_RESULTS/{{01_INPUT_RAW,02_CLINVAR,03_STRUCTURE_QC}}/stage_manifest.tsv
- FULL_RESULTS/12_REPRODUCIBILITY/provenance/provenance_0{{1,2,3}}_*.json
  (the four lines above are written after this report; each stage_manifest.tsv
   carries the hashes of its own directory, and handoff_01.json carries the
   run-level manifest)

QC RESULTS:
- Stratum 1 (all retrieved records): {conservation.n_stratum1}
- Stratum 2 (residues): {counts.n_residues_total}
- Stratum 3 (excluded records): {conservation.n_excluded}
- Stratum 4 (star distribution): emitted, non-empty
- Record conservation `stratum1 = eligible(stratum2) union stratum3`: {conservation.conserved}
- M = |U_struct| = {mapping.M}
- N = |L| = {counts.N} · N_P = {counts.N_P} · N_B = {counts.N_B}
- n_conflict (RESIDUE_CLASS_CONFLICT) = {counts.n_conflict}
- Unmapped cohort residues = {len(mapping.unmapped)}
  (excluded for no usable CA: {len(excluded_no_ca)})
- pLDDT bands: {bands}
- Review stars (all retrieved records): {stars}
- Structural QC: {qc.as_dict()['n_checks'] - qc.as_dict()['n_failed']}/{qc.as_dict()['n_checks']} checks passed
- Numbering offset: {qc.alignment.offset if qc.alignment else 'NA'} (applied: False)
- L subset of U_struct: {mapping.report['cohort_subset_of_universe']}

SCIENTIFIC DECISIONS:
- Review stars were recorded in all four strata and used for nothing else (F12).
- VUS, conflicting and every other non-binary category stayed out of the binary
  classes; none was relabelled to enlarge a cohort.
- RESIDUE_CLASS_CONFLICT residues keep all source records and were assigned to
  neither class (F3).
- pLDDT was recorded and banded; no residue was excluded by it and no PAE threshold
  exists (F2).
- Residue representation is CA and the distance is CA-CA Euclidean (F1).
- No numbering repair, sequence trimming, coordinate or class imputation, or cohort
  balancing was performed.
- Compound-significance policy [NEW — Workflow v2 §1]: {compound_cost['configured_policy']}
  (frozen by Lead ruling). Cost of this ruling versus {compound_cost['counterfactual_policy']}
  (counterfactual comparator only, never applied): {compound_cost['n_excluded_by_configured_but_kept_by_counterfactual']}
  record(s) excluded here that the counterfactual policy would have kept. See
  cohort_summary.json > compound_significance_policy_cost.

WARNINGS:
{chr(10).join(warnings_lines)}

FAILED OR REJECTED ANALYSES:
{chr(10).join('- ' + f for f in failed)}

UNRESOLVED ISSUES:
{_unresolved_block(negative)}

HANDOFF:
- `FULL_RESULTS/03_STRUCTURE_QC/handoff_01.json` -> hotspot-statistics
- Consumer must verify every manifest hash, both classes non-empty, every L residue
  present in U_struct with a usable CA, and qc_status != FAIL.
"""
    run.note(STAGE_STRUCT, "stage_a_report.md",
             write_text(stage_dir / "stage_a_report.md", text,
                        run_root=run.ctx.run_root),
             "Stage A narrative report")


def _unresolved_block(negative: dict | None) -> str:
    if negative is None:
        return ("- none. The significance-modifier reading policy "
                "(component-wise matching of compound ClinicalSignificance strings) is "
                "recorded in cohort_summary.json and is fully reversible from stratum 1 "
                "via significance_exact_match.")
    return (f"- {negative['condition']}: {negative['detail']}\n"
            f"- The INPUT could not support the analysis, so nothing was tested. "
            f"This is NOT a finding of 'no clustering' — read it as 'this gene's "
            f"ClinVar content cannot be analysed', never as 'this gene has no "
            f"hotspots'.\n"
            f"- The cause is the evidence base rather than a fault in the code "
            f"(outcome_type SCIENTIFIC_NEGATIVE), but the run is BLOCKED and "
            f"downstream stages must not proceed.")
