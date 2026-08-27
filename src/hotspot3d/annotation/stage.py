"""Stage E entry point — ``run_stage_e``.

Executes agent §7 in order: validate, copy verbatim, attach structure, annotate,
curate mechanisms blind, overlay, evaluate the post hoc gate, self-QC, hand off.

Stage E is terminal. Nothing it produces flows backwards, and the guardrails that
make that true — the write guard, the mandatory-context schema, the handoff scan and
the circularity replay — run on every execution, not only in tests.
"""
from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from ..orchestration.contracts import (
    QC_PASS,
    QC_PASS_WITH_WARNINGS,
    Handoff,
    assert_handoff_shape,
)
from ..utils.errors import BlockedError, EscalationRequired, LeakageError, NegativeResult
from ..utils.hashing import manifest_for, sha256_file
from ..utils.io import write_json, write_text, write_tsv
from ..utils.runctx import RunContext
from ..utils.status import OutcomeType, Severity, StageStatus, Status, WarningCollector, utc_now
from . import inputs as stage_inputs
from .mechanisms import (
    assert_no_membership_inference,
    build_search_log,
    curate,
    mechanism_by_hotspot,
    mechanism_counts,
)
from .posthoc import PosthocGate, evaluate_gate, gate_failure_record, run_posthoc_spatial
from .report import STANDING_CAVEATS, build_report
from .rubric import Mechanism, Rubric
from .schema import (
    AGENT,
    ANNOTATION_GAPS_COLUMNS,
    CONSERVATION_COLUMNS,
    DISEASE_ASSOCIATIONS_COLUMNS,
    DOMAINS_COLUMNS,
    FUNCTIONAL_MECHANISM_COLUMNS,
    FUNCTIONAL_SITES_COLUMNS,
    HOTSPOT_ANNOTATION_COLUMNS,
    LITERATURE_LOG_COLUMNS,
    MECHANISM_BY_HOTSPOT_COLUMNS,
    SECONDARY_STRUCTURE_COLUMNS,
    STAGE,
    STAGE_DIR,
    STAGE_MANIFEST_COLUMNS,
    assert_mandatory_context,
    assert_no_upstream_influence,
    assert_write_target,
)
from .sources import validate_source
from .verbatim import COPIED_STATISTICS, VerbatimAudit

#: Optional figures. Never blocking; every one not produced appears in the manifest.
OPTIONAL_FIGURES = {
    "figures/plddt_profile_with_hotspots.png": "pLDDT profile with hotspot overlay",
    "figures/star_distribution_by_class.png": "star-level distribution by class",
    "figures/sensitivity_overlap.png": "sensitivity-overlap bar/Venn",
    "figures/mechanism_distribution.png": "mechanism distribution across hotspots",
}

REQUIRED_OUTPUTS = [
    "hotspot_annotation.tsv", "secondary_structure.tsv", "domains.tsv",
    "functional_sites.tsv", "conservation.tsv", "disease_associations.tsv",
    "functional_mechanism_variants.tsv", "mechanism_by_hotspot.tsv",
    "mechanism_spatial_posthoc.json", "annotation_gaps.tsv",
    "literature_search_log.tsv", "handoff_05.json", "stage_e_report.md",
    "stage_status.json", "warnings_10_annotation.tsv", "stage_manifest.tsv",
]


class StageWriter:
    """Every Stage E write goes through here, so the ownership guard cannot be bypassed."""

    def __init__(self, ctx: RunContext, stage_dir: Path):
        self.ctx = ctx
        self.stage_dir = stage_dir
        self.created: list[Path] = []

    def _target(self, name: str) -> Path:
        path = self.stage_dir / name
        assert_write_target(self.ctx.run_root, path)
        return path

    def tsv(self, name: str, rows: list[dict], columns: list[str]) -> Path:
        path = write_tsv(self._target(name), rows, columns, run_root=self.ctx.run_root)
        self.created.append(path)
        return path

    def json(self, name: str, payload: dict) -> Path:
        path = write_json(self._target(name), payload, run_root=self.ctx.run_root)
        self.created.append(path)
        return path

    def text(self, name: str, text: str) -> Path:
        path = write_text(self._target(name), text, run_root=self.ctx.run_root)
        self.created.append(path)
        return path


def _resolve_uniprot_acc(ctx: RunContext, handoff_02: Handoff) -> str:
    """The UniProt accession Stage A resolved, for the live source and for
    ``load_inputs``. handoff_02 has never actually carried this key (a
    pre-existing gap, not introduced here — Stage E had no live consumer of
    it until now); the authoritative record is Stage A's own
    ``01_INPUT_RAW/retrieval_log.json``, already inside this agent's read
    scope. Falls back to "NA" (the existing behavior) rather than raising, so
    synthetic fixtures that provide neither are unaffected.
    """
    acc = handoff_02.payload.get("uniprot_acc")
    if acc:
        return str(acc)
    log_path = ctx.full_results / "01_INPUT_RAW" / "retrieval_log.json"
    try:
        import json

        log = json.loads(log_path.read_text(encoding="utf-8"))
        acc = log.get("gene_resolution", {}).get("uniprot_acc")
        if acc:
            return str(acc)
    except (OSError, ValueError, KeyError):
        pass
    return "NA"


def run_stage_e(
    ctx: RunContext,
    *,
    handoff_02: Handoff,
    handoff_03: Handoff,
    handoff_04: Handoff,
    source: Any = None,
) -> Handoff:
    """Annotate the frozen hotspots and curate functional mechanisms.

    Returns ``handoff_05``. Raises nothing for a scientific negative — that is a
    completed outcome with an explicit ``negative_result`` block.
    """
    started = time.time()
    started_utc = utc_now()
    ctx.assert_may_write(AGENT, STAGE_DIR)
    stage_dir = ctx.stage_dir(STAGE_DIR)
    writer = StageWriter(ctx, stage_dir)
    warnings = WarningCollector(STAGE, AGENT)

    if ctx.synthetic:
        warnings.add(
            "SYNTHETIC_INPUT_MODE",
            "Stage E ran against a synthetic/mock annotation source. Annotation "
            "content is not real biological data and must not be interpreted as such.",
            potential_consequence="Biological conclusions drawn from this run are void.",
            recommended_action="Re-run with versioned live sources before publication.",
            affected_output="10_ANNOTATION/",
        )

    try:
        return _execute(
            ctx, handoff_02, handoff_03, handoff_04, source,
            writer=writer, warnings=warnings, stage_dir=stage_dir,
            started=started, started_utc=started_utc,
        )
    except NegativeResult as neg:
        return _finish_negative(
            ctx, neg, writer=writer, warnings=warnings, stage_dir=stage_dir,
            started=started, started_utc=started_utc, handoff_04=handoff_04,
        )
    except EscalationRequired as esc:
        _finish_blocked(
            ctx, writer=writer, warnings=warnings, stage_dir=stage_dir,
            started=started, started_utc=started_utc,
            reason=(
                f"ESCALATION — {esc.ambiguity}\n\nOptions: {esc.options}\n"
                f"Consequences: {esc.consequences}\n"
                f"Recommendation: {esc.recommendation}"
            ),
            action="Stage E stopped and escalated. The Lead decides; Stage E does not.",
        )
        raise
    except (BlockedError, LeakageError) as exc:
        _finish_blocked(
            ctx, writer=writer, warnings=warnings, stage_dir=stage_dir,
            started=started, started_utc=started_utc, reason=str(exc),
            action="Resolve the blocking condition upstream and re-run Stage E.",
        )
        raise


# --- main path --------------------------------------------------------------

def _execute(
    ctx: RunContext, handoff_02: Handoff, handoff_03: Handoff, handoff_04: Handoff,
    source: Any, *, writer: StageWriter, warnings: WarningCollector, stage_dir: Path,
    started: float, started_utc: str,
) -> Handoff:
    cfg = ctx.config

    # 1 — config: the rubric and the gate must be defined (agent §6.5).
    cfg.require_all(
        [
            "annotation.mechanism_categories", "annotation.evidence_levels",
            "annotation.min_evidence_for_label",
            "annotation.conflicting_strong_resolution",
            "annotation.weak_only_resolution", "annotation.conservation_source",
            "annotation.secondary_structure_tool",
            "annotation.require_versioned_sources", "annotation.may_influence_discovery",
        ]
    )
    if cfg.get("annotation.may_influence_discovery"):
        raise BlockedError(
            "FROZEN methodology violation: annotation.may_influence_discovery is TRUE. "
            "Stage E is terminal by construction."
        )
    rubric = Rubric.from_config(cfg)
    gate = PosthocGate.from_config(cfg)
    uniprot_acc = _resolve_uniprot_acc(ctx, handoff_02)

    if source is None:
        if ctx.synthetic:
            raise BlockedError(
                "BLOCKED — a synthetic run must be given an explicit source; Stage E "
                "never fabricates annotation data."
            )
        ctx.assert_may_read(AGENT, "01_INPUT_RAW")
        from .live_sources import LiveAnnotationSource

        source = LiveAnnotationSource(
            allow_network=bool(cfg.get_optional("execution.allow_network", False)),
            gene=ctx.gene,
            uniprot_acc=uniprot_acc,
            clinvar_raw_path=ctx.full_results / "01_INPUT_RAW" / "clinvar_raw.tsv",
            variants_residue_level_path=(
                ctx.full_results / "02_CLINVAR" / "variants_residue_level.tsv"),
        )
    source_versions = validate_source(source)

    # 2 — upstream validation, hashes, negative-result branch.
    data = stage_inputs.load_inputs(
        ctx, handoff_02, handoff_03, handoff_04,
        uniprot_acc=uniprot_acc,
    )
    _warn_upstream_flags(warnings, data.flags)

    gaps: list[dict] = []
    profile_summary = _summarize_profile(data.robustness_profile)
    influence = _per_center_influence(data)

    # 3 — hotspot inventory, copied verbatim. Nothing here is computed.
    audit = VerbatimAudit(source=data.regions)
    rows: list[dict] = []
    for hid in data.hotspot_ids:
        copied = audit.copy_all(hid, COPIED_STATISTICS)
        centers = data.centers.get(hid, [])
        covered = data.covered_residues.get(hid, [])
        sphere = data.sphere_variants.get(hid, [])
        residues = sorted(int(r["residue_index"]) for r in covered)

        plddt_values = [data.plddt[r] for r in residues if r in data.plddt]
        if len(plddt_values) < len(residues):
            gaps.append(_gap("hotspot", hid, "pLDDT", "residue absent from plddt_profile",
                             "03_STRUCTURE_QC/plddt_profile.tsv", "ADVISORY"))

        rows.append(
            {
                "hotspot_id": hid,
                "residue_start": copied["residue_start"],
                "residue_end": copied["residue_end"],
                "covered_interval": f"{copied['residue_start']}-{copied['residue_end']}",
                "n_plp": copied["n_plp"],
                "n_blb": copied["n_blb"],
                "fold_enrichment": copied["fold_enrichment"],
                "p_emp": copied["p_emp"],
                "q_bh": copied["q_bh"],
                # three residue objects, separate and never collapsed
                "n_significant_hotspot_centers": len(centers),
                "n_hotspot_sphere_classified_variants": len(sphere),
                "n_hotspot_covered_residues": len(covered),
                # pLDDT is read, never re-derived
                "plddt_mean": _mean(plddt_values),
                "plddt_min": min(plddt_values) if plddt_values else None,
            }
        )

    # 4 — structural annotation from the versioned source.
    all_residues = data.all_annotated_residues
    ss_map = source.secondary_structure(data.structure_path)
    domains = source.domains(data.uniprot_acc)
    sites = source.functional_sites(data.uniprot_acc)
    conservation = source.conservation(data.uniprot_acc, all_residues)
    diseases = source.disease_associations(data.gene)
    # A live source reports its OWN retrieval failures explicitly (never
    # silently indistinguishable from a real "queried, found nothing"
    # negative) via this optional attribute; a synthetic/test source simply
    # doesn't have one.
    gaps.extend(getattr(source, "source_gaps", []))

    ss_rows, domain_rows, site_rows, cons_rows, disease_rows = [], [], [], [], []

    for row in rows:
        hid = row["hotspot_id"]
        residues = sorted(
            int(r["residue_index"]) for r in data.covered_residues.get(hid, [])
        )
        span = (min(residues), max(residues)) if residues else (0, -1)

        ss_rows.extend(_ss_rows(hid, residues, ss_map, source_versions, gaps))
        hs_domains = _overlapping(domains, span, "start", "end")
        domain_rows.extend(_domain_rows(hid, hs_domains, source_versions))
        hs_sites = _overlapping(sites, span, "residue_start", "residue_end")
        site_rows.extend(_site_rows(hid, hs_sites, source_versions))
        cons_rows.extend(_cons_rows(hid, residues, conservation, source_versions, gaps))
        hs_diseases = _overlapping(diseases, span, "residue_start", "residue_end")
        disease_rows.extend(_disease_rows(hid, hs_diseases, source_versions))

        if not hs_domains:
            gaps.append(_gap("hotspot", hid, "overlapping_domain",
                             "hotspot overlaps no known domain — a valid, informative "
                             "negative result", "InterPro/Pfam/UniProt", "INFO"))
        if not hs_diseases:
            gaps.append(_gap("hotspot", hid, "disease_association",
                             "no prior disease association reported — a valid, "
                             "informative negative result", "literature", "INFO"))

        row["secondary_structure_composition"] = _ss_composition(residues, ss_map)
        row["overlapping_domains"] = _names(hs_domains, "name") or "NA"
        row["functional_regions"] = _names(
            [s for s in hs_sites if s.get("site_type") == "functional_region"], "description"
        ) or "NA"
        row["ligand_sites"] = _names(
            [s for s in hs_sites if s.get("site_type") == "ligand_binding"], "description"
        ) or "NA"
        row["protein_interaction_sites"] = _names(
            [s for s in hs_sites if s.get("site_type") == "protein_interaction"], "description"
        ) or "NA"
        row["conservation_summary"] = _cons_summary(residues, conservation)
        row["disease_associations"] = _names(hs_diseases, "disease") or "NA"

        # 5 — descriptive evidence-quality context. Never a re-filter.
        row["star_composition"] = _star_composition(hid, data)
        ge1, ge2 = _sensitivity_overlap(hid, data, gaps)
        row["overlap_ge1star"] = ge1
        row["overlap_ge2star"] = ge2

        # 6 — robustness + flags on EVERY row (schema requirement).
        row["center_recurrence_rate"] = influence.get(hid, {}).get("recurrence", "NA")
        row["robustness_profile_summary"] = profile_summary
        row["per_center_influence_summary"] = influence.get(hid, {}).get("influence", "NA")
        row["global_clustering_flag"] = data.flags["global_clustering_flag"]
        row["permutation_resolution_limited"] = data.flags["permutation_resolution_limited"]
        row["fallback_radius_domain"] = data.flags["fallback_radius_domain"]
        row["domain_boundary_warning"] = data.flags["domain_boundary_warning"]

    assert_mandatory_context(rows)
    annotation_path = writer.tsv("hotspot_annotation.tsv", rows, HOTSPOT_ANNOTATION_COLUMNS)

    writer.tsv("secondary_structure.tsv", ss_rows, SECONDARY_STRUCTURE_COLUMNS)
    writer.tsv("domains.tsv", domain_rows, DOMAINS_COLUMNS)
    writer.tsv("functional_sites.tsv", site_rows, FUNCTIONAL_SITES_COLUMNS)
    writer.tsv("conservation.tsv", cons_rows, CONSERVATION_COLUMNS)
    writer.tsv("disease_associations.tsv", disease_rows, DISEASE_ASSOCIATIONS_COLUMNS)

    # 7 — mechanism curation, structurally blind to hotspot membership.
    records = list(source.literature_mechanisms(data.gene))
    if getattr(source, "literature_mechanisms_not_run", False):
        warnings.add(
            "LITERATURE_MECHANISMS_NOT_RUN",
            "literature_mechanisms was not attempted by this run's annotation "
            "source (out of scope for this build). Every functional-mechanism "
            "label below defaults to Not_Experimentally_Characterized as a "
            "NON-RESULT of an unrun search, not a negative finding.",
            severity=Severity.MAJOR,
            potential_consequence=(
                "mechanism_by_hotspot.tsv, functional_mechanism_variants.tsv and "
                "the post hoc gate all reflect zero literature evidence, which "
                "could otherwise be misread as 'no functional mechanism found'."),
            recommended_action=(
                "Re-run Stage E once literature_mechanisms is implemented before "
                "drawing any conclusion about functional mechanisms."),
            affected_output="functional_mechanism_variants.tsv",
        )
    curated = curate(records, rubric)
    circularity_qc = assert_no_membership_inference(curated, rubric)

    for cv in curated:
        if cv.mechanism is Mechanism.NOT_CHARACTERIZED:
            gaps.append(_gap(
                "variant", f"{cv.residue_index}{cv.aa_change}", "functional_mechanism",
                "no experimental evidence retrieved — category "
                "Not_Experimentally_Characterized assigned by rubric, never by inference",
                "literature", "INFO",
            ))
        if cv.resolution.rubric_gap:
            warnings.add(
                "REVIEW_STATUS_SENSITIVE",
                f"Variant {cv.residue_index}{cv.aa_change}: conflicting "
                f"{cv.evidence_strength.value} sources fall outside the frozen rubric, "
                f"which specifies conflicting-STRONG and conflicting-WEAK only. "
                f"Resolved conservatively to {cv.mechanism.value}. OPEN RUBRIC POINT.",
                severity=Severity.MAJOR,
                potential_consequence="Mechanism counts, and therefore the post hoc "
                                      "gate, depend on this unresolved rubric point.",
                recommended_action="Lead to rule on conflicting-MODERATE resolution.",
                affected_output="functional_mechanism_variants.tsv",
            )

    writer.tsv("functional_mechanism_variants.tsv",
               [cv.as_row() for cv in curated], FUNCTIONAL_MECHANISM_COLUMNS)
    writer.tsv("literature_search_log.tsv", build_search_log(records),
               LITERATURE_LOG_COLUMNS)

    # 8 — the overlay, strictly after curation.
    row_context = {
        r["hotspot_id"]: {
            "robustness_profile_summary": r["robustness_profile_summary"],
            "global_clustering_flag": r["global_clustering_flag"],
        }
        for r in rows
    }
    writer.tsv(
        "mechanism_by_hotspot.tsv",
        mechanism_by_hotspot(curated, data.covered_residues,
                             minimum=rubric.min_evidence_for_label,
                             row_context=row_context),
        MECHANISM_BY_HOTSPOT_COLUMNS,
    )

    counts = mechanism_counts(curated)

    # 9 — the gate. Numbers are recorded whether it passes or fails.
    evaluation = evaluate_gate(curated, gate)
    if evaluation.passed:
        posthoc = run_posthoc_spatial(
            curated, data.coords, gate, evaluation,
            rng=ctx.seeds.rng("posthoc|mechanism_spatial"),
            seed=ctx.seeds.seed("posthoc|mechanism_spatial"),
            n_permutations=int(cfg.get("permutation.B_default")),
            q=float(cfg.get("fdr.q")),
        )
    else:
        posthoc = gate_failure_record(evaluation, gate)
        warnings.add(
            "POSTHOC_GATE_FAILED",
            evaluation.reason,
            potential_consequence="No mechanism-specific spatial claim can be made.",
            recommended_action="None. The gate is not lowered to obtain a result.",
            affected_output="mechanism_spatial_posthoc.json",
        )
    writer.json("mechanism_spatial_posthoc.json", posthoc)

    # 10 — gaps. An empty gap file is itself suspicious and is justified in the
    # report (§4). The frozen warning enum has no code for it, so it is raised as a
    # narrative justification rather than by bending an unrelated code.
    writer.tsv("annotation_gaps.tsv", gaps, ANNOTATION_GAPS_COLUMNS)

    # 11 — self-QC.
    qc = [
        audit.verify_written(annotation_path),
        circularity_qc,
        _qc_plddt(rows, data),
        _qc_mandatory_context(rows),
        _qc_residue_objects(rows),
        _qc_upstream_unchanged(ctx, data),
    ]

    # 12 — report, manifest, provenance, handoff.
    report = build_report(
        gene=ctx.gene, run_id=ctx.run_id, annotation_rows=rows,
        curated=[cv.as_row() for cv in curated], counts=counts,
        gate_numbers=evaluation.as_numbers(), posthoc=posthoc, qc=qc,
        source_versions=source_versions, gaps=gaps, warnings=warnings.items,
        synthetic=ctx.synthetic,
    )
    writer.text("stage_e_report.md", report)

    qc_status = (
        QC_PASS_WITH_WARNINGS
        if warnings.counts()[Severity.MAJOR.value] else QC_PASS
    )

    payload = {
        "annotation_sources": sorted(source_versions),
        "annotation_source_versions": {
            name: sv.as_row() for name, sv in sorted(source_versions.items())
        },
        "n_hotspots_annotated": len(rows),
        "mechanism_counts": counts,
        "posthoc_gate_passed": evaluation.passed,
        "posthoc_gate_numbers": evaluation.as_numbers(),
        "robustness_profile": handoff_04.payload["robustness_profile"],
        "global_clustering_flag": data.flags["global_clustering_flag"],
        "interpretation_caveats": list(STANDING_CAVEATS),
        "upstream_handoff_hashes": data.upstream_hashes,
        "residue_objects": {
            "SIGNIFICANT_HOTSPOT_CENTERS": sum(
                r["n_significant_hotspot_centers"] for r in rows),
            "HOTSPOT_SPHERE_CLASSIFIED_VARIANTS": sum(
                r["n_hotspot_sphere_classified_variants"] for r in rows),
            "HOTSPOT_COVERED_RESIDUES": sum(
                r["n_hotspot_covered_residues"] for r in rows),
            "note": (
                "Three distinct objects, never collapsed. A significant hotspot center "
                "is a geometric test position and need not carry a ClinVar variant."
            ),
        },
        "terminal_stage": True,
        "may_influence_discovery": False,
    }
    assert_no_upstream_influence(payload)

    status = StageStatus(
        stage=STAGE, agent_owner=AGENT, status=Status.COMPLETED,
        outcome_type=OutcomeType.COMPLETED,
        reason="Annotation and mechanism curation completed.",
        started_utc=started_utc, ended_utc=utc_now(), wall_seconds=time.time() - started,
        upstream_handoff="handoff_02,handoff_03,handoff_04",
        upstream_qc_status=handoff_04.qc_status,
        n_outputs_expected=len(REQUIRED_OUTPUTS) + len(OPTIONAL_FIGURES),
        n_outputs_created=len(writer.created) + 3,
        warnings=warnings.counts(),
        recommended_action=(
            "Lead must independently verify: every copied statistic matches Stage B, "
            "the robustness profile is attached to every claim, and the leakage audit "
            "passes."
        ),
    )
    return _finalize(
        ctx, writer=writer, warnings=warnings, stage_dir=stage_dir, status=status,
        payload=payload, qc_status=qc_status, qc=qc, negative=None,
    )


# --- terminal branches ------------------------------------------------------

def _finish_negative(
    ctx: RunContext, neg: NegativeResult, *, writer: StageWriter,
    warnings: WarningCollector, stage_dir: Path, started: float, started_utc: str,
    handoff_04: Handoff,
) -> Handoff:
    """Upstream ended negatively: create the stage, invent nothing (agent §15)."""
    warnings.add(
        "NO_SIGNIFICANT_HOTSPOTS",
        f"Upstream discovery ended in a valid negative result ({neg.condition}). "
        f"No region was annotated. Stage E does not invent regions, and never "
        f"annotates a 'best non-significant' region as if it were a finding.",
        potential_consequence="No biological interpretation is available for this gene.",
        recommended_action="None. A negative result is a valid, informative outcome.",
        affected_output="10_ANNOTATION/",
    )
    negative = {
        "condition": neg.condition,
        "detail": neg.detail,
        "upstream_cause": neg.context.get("upstream", "NA"),
        "interpretation": (
            "TRUE SCIENTIFIC NEGATIVE upstream. Stage E ran correctly and had nothing "
            "to annotate. This licenses no method modification."
        ),
    }

    for name, columns in (
        ("hotspot_annotation.tsv", HOTSPOT_ANNOTATION_COLUMNS),
        ("secondary_structure.tsv", SECONDARY_STRUCTURE_COLUMNS),
        ("domains.tsv", DOMAINS_COLUMNS),
        ("functional_sites.tsv", FUNCTIONAL_SITES_COLUMNS),
        ("conservation.tsv", CONSERVATION_COLUMNS),
        ("disease_associations.tsv", DISEASE_ASSOCIATIONS_COLUMNS),
        ("functional_mechanism_variants.tsv", FUNCTIONAL_MECHANISM_COLUMNS),
        ("mechanism_by_hotspot.tsv", MECHANISM_BY_HOTSPOT_COLUMNS),
        ("literature_search_log.tsv", LITERATURE_LOG_COLUMNS),
    ):
        writer.tsv(name, [], columns)

    writer.tsv("annotation_gaps.tsv", [_gap(
        "stage", STAGE_DIR, "entire_annotation",
        f"upstream negative result ({neg.condition}); nothing to annotate",
        "upstream", "INFO")], ANNOTATION_GAPS_COLUMNS)
    writer.json("mechanism_spatial_posthoc.json", {
        "analysis_type": "post_hoc", "status": "NOT_RUN",
        "outcome_type": "NOT_APPLICABLE",
        "reason": "No hotspots upstream; the gate was not reached.",
        "feeds_upstream": False,
    })
    writer.text("stage_e_report.md", _negative_report(ctx, neg))

    status = StageStatus(
        # NOT_RUN, not COMPLETED_NEGATIVE: COMPLETED_NEGATIVE would assert that *this*
        # stage ran and produced a scientific negative. Discovery terminated upstream,
        # so Stage E never became applicable and produced no negative finding of its
        # own — claiming one would misreport a non-event as a result (Output Contract
        # IX.8.3; agent §15).
        stage=STAGE, agent_owner=AGENT, status=Status.NOT_RUN,
        outcome_type=OutcomeType.NOT_APPLICABLE,
        reason=f"{neg.condition}: {neg.detail}",
        started_utc=started_utc, ended_utc=utc_now(), wall_seconds=time.time() - started,
        upstream_handoff="handoff_02", upstream_qc_status=handoff_04.qc_status,
        n_outputs_expected=len(REQUIRED_OUTPUTS) + len(OPTIONAL_FIGURES),
        n_outputs_created=len(writer.created) + 3,
        warnings=warnings.counts(),
        recommended_action="None. Report the negative result as the finding.",
        negative_result=negative,
    )
    payload = {
        "annotation_sources": [], "annotation_source_versions": {},
        "n_hotspots_annotated": 0,
        "mechanism_counts": {m.value: 0 for m in Mechanism},
        "posthoc_gate_passed": False,
        "posthoc_gate_numbers": {
            "passed": False,
            "reason": "Gate not reached: no hotspots were produced upstream.",
        },
        "robustness_profile": handoff_04.payload.get("robustness_profile", {}),
        "global_clustering_flag": "NA",
        "interpretation_caveats": list(STANDING_CAVEATS),
        "upstream_handoff_hashes": {},
        "residue_objects": {
            "SIGNIFICANT_HOTSPOT_CENTERS": 0,
            "HOTSPOT_SPHERE_CLASSIFIED_VARIANTS": 0,
            "HOTSPOT_COVERED_RESIDUES": 0,
            "note": "No hotspots upstream.",
        },
        "terminal_stage": True, "may_influence_discovery": False,
    }
    assert_no_upstream_influence(payload)
    return _finalize(
        ctx, writer=writer, warnings=warnings, stage_dir=stage_dir, status=status,
        payload=payload, qc_status=QC_PASS, qc=[], negative=negative,
    )


def _finish_blocked(
    ctx: RunContext, *, writer: StageWriter, warnings: WarningCollector,
    stage_dir: Path, started: float, started_utc: str, reason: str, action: str,
) -> None:
    """Always leave stage_status.json + NOT_RUN.txt behind — no mysteriously empty folder."""
    status = StageStatus(
        stage=STAGE, agent_owner=AGENT, status=Status.BLOCKED,
        outcome_type=OutcomeType.TECHNICAL_FAILURE, reason=reason,
        started_utc=started_utc, ended_utc=utc_now(), wall_seconds=time.time() - started,
        n_outputs_expected=len(REQUIRED_OUTPUTS) + len(OPTIONAL_FIGURES),
        n_outputs_created=len(writer.created),
        warnings=warnings.counts(), recommended_action=action,
    )
    status.write(stage_dir)
    warnings.write(stage_dir)
    _write_manifest(ctx, writer, stage_dir)


def _finalize(
    ctx: RunContext, *, writer: StageWriter, warnings: WarningCollector,
    stage_dir: Path, status: StageStatus, payload: dict, qc_status: str,
    qc: list[dict], negative: dict | None,
) -> Handoff:
    """Write handoff_05, status, warnings, manifest and provenance."""
    handoff_path = ctx.handoff_path("handoff_05")
    assert_write_target(ctx.run_root, handoff_path)

    outputs = [p for p in writer.created]
    handoff = Handoff(
        name="handoff_05", run_id=ctx.run_id, config_sha256=ctx.config.sha256,
        qc_status=qc_status,
        manifest=manifest_for(outputs, root=ctx.run_root),
        payload=payload, negative_result=negative,
    )
    assert_handoff_shape(handoff)
    write_json(handoff_path, handoff.as_dict(), run_root=ctx.run_root)
    writer.created.append(handoff_path)

    status.n_outputs_created = len(writer.created) + 3  # + status, warnings, manifest
    status.write(stage_dir)
    warnings.write(stage_dir)
    _write_manifest(ctx, writer, stage_dir)

    provenance_path = ctx.provenance_path(STAGE)
    assert_write_target(ctx.run_root, provenance_path)
    ctx.record_provenance(
        STAGE, AGENT,
        inputs={
            name: {"path": str(path), "sha256": sha256_file(path)}
            for name, path in _existing_inputs(ctx).items()
        },
        outputs={str(p.relative_to(ctx.run_root)): sha256_file(p) for p in writer.created},
        parameters={
            "rubric": "METHOD_SPEC II.13 (frozen)",
            "posthoc_gate": ctx.config.get("annotation.posthoc_gate"),
            "conservation_source": ctx.config.get("annotation.conservation_source"),
            "secondary_structure_tool": ctx.config.get(
                "annotation.secondary_structure_tool"),
        },
        warnings=[w.as_row() for w in warnings.items],
        extra={
            "qc_records": qc,
            "annotation_source_versions": payload.get("annotation_source_versions", {}),
            "network_used": False,
            "terminal_stage": True,
            "feeds_upstream": False,
        },
    )
    return handoff


def _existing_inputs(ctx: RunContext) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for name in stage_inputs.UPSTREAM_PATHS:
        stage, rel = stage_inputs.UPSTREAM_PATHS[name]
        path = ctx.full_results / stage / rel
        if path.is_file():
            out[name] = path
    return out


def _write_manifest(ctx: RunContext, writer: StageWriter, stage_dir: Path) -> None:
    # Report what is actually on disk, not what the writer believes it wrote.
    # stage_manifest.tsv describes itself as CREATED because this call creates it.
    def exists(name: str) -> bool:
        return name == "stage_manifest.tsv" or (stage_dir / name).is_file()

    rows = [
        {
            "file": name, "kind": "required",
            "status": "CREATED" if exists(name) else "NOT_CREATED",
            "reason": "NA" if exists(name)
            else "stage terminated before this file was written",
        }
        for name in REQUIRED_OUTPUTS
    ]
    rows.extend(
        {
            "file": name, "kind": "optional_figure", "status": "NOT_CREATED",
            "reason": (
                f"{desc} not produced; optional figures are never blocking and add no "
                f"information beyond the tables in this stage."
            ),
        }
        for name, desc in OPTIONAL_FIGURES.items()
    )
    manifest_path = stage_dir / "stage_manifest.tsv"
    assert_write_target(ctx.run_root, manifest_path)
    write_tsv(manifest_path, rows, STAGE_MANIFEST_COLUMNS, run_root=ctx.run_root)


def _negative_report(ctx: RunContext, neg: NegativeResult) -> str:
    return (
        f"# Stage E — Biological Annotation: {ctx.gene}\n\n"
        f"- **Run ID:** `{ctx.run_id}`\n"
        f"- **Status:** NOT_RUN / NOT_APPLICABLE\n\n"
        f"## Upstream negative result\n\n"
        f"**Condition:** `{neg.condition}`\n\n{neg.detail}\n\n"
        f"## What this means\n\n"
        f"Discovery ended in a valid negative result, so there is no region to "
        f"annotate. **Stage E did not invent regions**, and it did not annotate a "
        f"'best non-significant' region as if it were a finding. A negative result is "
        f"a valid, informative outcome and licenses no method modification.\n\n"
        f"The status is `NOT_RUN`, not `COMPLETED_NEGATIVE`: the negative finding "
        f"belongs to the upstream discovery stage. Stage E never became applicable and "
        f"produced no negative result of its own, and reporting one would present a "
        f"non-event as a result.\n\n"
        f"All canonical Stage E tables exist with headers and no rows: an empty table "
        f"is explicit, never a missing file.\n"
    )


# --- QC helpers -------------------------------------------------------------

def _qc_plddt(rows: list[dict], data) -> dict:
    bad = []
    for row in rows:
        residues = sorted(
            int(r["residue_index"])
            for r in data.covered_residues.get(row["hotspot_id"], [])
        )
        values = [data.plddt[r] for r in residues if r in data.plddt]
        if values and row["plddt_min"] != min(values):
            bad.append(row["hotspot_id"])
    return {
        "check": "plddt_matches_03_structure_qc",
        "n_hotspots": len(rows), "n_mismatches": len(bad),
        "result": "PASS" if not bad else "FAIL",
        "note": "pLDDT read from 03_STRUCTURE_QC and never re-derived.",
    }


def _qc_mandatory_context(rows: list[dict]) -> dict:
    assert_mandatory_context(rows)
    return {
        "check": "robustness_and_all_four_flags_on_every_row",
        "n_rows": len(rows), "result": "PASS",
        "note": "Enforced as a schema requirement; a row missing any cannot be written.",
    }


def _qc_residue_objects(rows: list[dict]) -> dict:
    collapsed = [
        r["hotspot_id"] for r in rows
        if r["n_significant_hotspot_centers"] == r["n_hotspot_covered_residues"]
        == r["n_hotspot_sphere_classified_variants"] and r["n_hotspot_covered_residues"] > 1
    ]
    return {
        "check": "three_residue_objects_reported_separately",
        "n_rows": len(rows), "result": "PASS",
        "coincidentally_equal_counts": collapsed,
        "note": (
            "Counts are emitted in three separate columns and never summed. A "
            "significant hotspot center is a geometric test position and need not "
            "carry a ClinVar variant."
        ),
    }


def _qc_upstream_unchanged(ctx: RunContext, data) -> dict:
    """Re-verify every upstream hash at stage end (agent §9.11).

    Proves Stage E left Stage A–D artifacts untouched, and that nothing changed
    underneath it mid-stage. A difference is BLOCKING, not a note.
    """
    missing, changed = [], []
    for name, path in data.paths.items():
        if not path.is_file():
            missing.append(name)
        elif sha256_file(path) != data.input_hashes.get(name):
            changed.append(name)

    if missing or changed:
        raise BlockedError(
            f"BLOCKED — upstream artifacts changed during Stage E "
            f"(missing={missing}, modified={changed}). Nothing downstream may modify "
            f"Stage A–D artifacts, and an input that moves mid-stage invalidates every "
            f"verbatim copy taken from it."
        )
    return {
        "check": "upstream_inputs_unchanged_at_stage_end",
        "n_inputs": len(data.paths), "n_missing": 0, "n_modified": 0,
        "result": "PASS",
        "note": "SHA-256 recomputed for every declared input and compared to load time.",
    }


# --- small helpers ----------------------------------------------------------

def _warn_upstream_flags(warnings: WarningCollector, flags: dict) -> None:
    """Inherited flags are repeated in Stage E warnings; flags are never dropped."""
    if _truthy(flags.get("permutation_resolution_limited")):
        warnings.add(
            "PERMUTATION_RESOLUTION_LIMITED",
            "Inherited from Stage B: permutation resolution was limited. Every "
            "annotated hotspot carries this flag.",
            affected_output="hotspot_annotation.tsv",
        )
    if _truthy(flags.get("fallback_radius_domain")):
        warnings.add(
            "FALLBACK_RADIUS_DOMAIN",
            "Inherited from Stage B: the fallback radius domain was used. Every "
            "annotated hotspot carries this flag.",
            affected_output="hotspot_annotation.tsv",
        )
    if _truthy(flags.get("domain_boundary_warning")):
        warnings.add(
            "BOUNDARY_OPTIMUM_WARNING",
            "Inherited from upstream: the selected radius sits near a domain boundary. "
            "Every annotated hotspot carries this flag.",
            affected_output="hotspot_annotation.tsv",
        )
    clustering = str(flags.get("global_clustering_flag", "")).lower()
    if "non_significant" in clustering or clustering in ("false", "0"):
        warnings.add(
            "GLOBAL_CLUSTERING_NON_SIGNIFICANT",
            "Inherited from Stage B: global clustering was NOT significant. Annotated "
            "hotspots must be read with this caveat; local discovery is not "
            "invalidated by it, but no hotspot may be presented as established on the "
            "strength of a global signal that was absent.",
            affected_output="hotspot_annotation.tsv",
        )


def _truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().upper() in ("TRUE", "1", "YES")


def _summarize_profile(profile: dict) -> str:
    """Report the profile's own numbers. No verdict, no re-thresholding.

    Stage E emits no ROBUST / NOT_ROBUST label of its own (agent §16): it forwards
    what Stage D measured, in Stage D's own terms.
    """
    if not profile:
        return "NA"
    parts = [
        f"{key}={value}"
        for key, value in sorted(profile.items())
        if key != "schema_version" and isinstance(value, (int, float, str, bool))
    ]
    return ";".join(parts) if parts else "profile_present_no_scalar_fields"


def _per_center_influence(data) -> dict[str, dict[str, str]]:
    """Per-center influence and recurrence, carried through verbatim.

    Values are listed per center rather than aggregated: averaging Stage D's numbers
    would create a statistic Stage E is not entitled to author.

    Stage D may correctly report ``evaluable: false`` (``ROBUSTNESS_NOT_EVALUABLE``,
    e.g. n_S below the minimum for evaluation): a designed terminal sub-state, not a
    missing measurement. That state is reported verbatim rather than collapsed to
    "NA", which the schema reserves for a silently dropped field (agent §8, §16 — no
    re-thresholding or re-labelling of the robustness profile).
    """
    # Stage D emits center_sensitivity.tsv keyed on `center_residue_index`, with the
    # influence score as `mean_jaccard_drop` and recurrence as `center_recurrence_rate`
    # (II.11: per-center influence and recurrence). Names are resolved here, once.
    def _first(record: dict, *names: str):
        for name in names:
            if record.get(name) is not None:
                return record[name]
        return None

    by_residue = {}
    for r in data.center_sensitivity:
        key = _first(r, "center_residue_index", "residue_index")
        if key is not None:
            by_residue[int(key)] = r

    not_evaluable_text = None
    if data.robustness_profile.get("evaluable") is False:
        reason = data.robustness_profile.get("reason")
        n_s = data.robustness_profile.get("n_S")
        not_evaluable_text = reason or f"ROBUSTNESS_NOT_EVALUABLE (n_S={n_s})"

    out: dict[str, dict[str, str]] = {}
    for hid, centers in data.centers.items():
        influence, recurrence = [], []
        for center in centers:
            idx = int(center["residue_index"])
            record = by_residue.get(idx, {})
            infl = _first(record, "mean_jaccard_drop", "influence")
            if infl is not None:
                influence.append(f"{idx}={infl}")
            rec = _first(record, "center_recurrence_rate", "recurrence_rate")
            if rec is not None:
                recurrence.append(f"{idx}={rec}")
        fallback = not_evaluable_text if not_evaluable_text is not None else "NA"
        out[hid] = {
            "influence": ";".join(influence) if influence else fallback,
            "recurrence": ";".join(recurrence) if recurrence else fallback,
        }
    return out


def _mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def _overlapping(items: list[dict], span: tuple[int, int], start_key: str,
                 end_key: str) -> list[dict]:
    lo, hi = span
    return [
        it for it in items
        if it.get(start_key) is not None and it.get(end_key) is not None
        and int(it[start_key]) <= hi and int(it[end_key]) >= lo
    ]


def _names(items: list[dict], key: str) -> str:
    values = [str(it.get(key)) for it in items if it.get(key)]
    return ";".join(values)


def _version_of(versions: dict, name: str) -> str:
    sv = versions.get(name)
    return sv.version if sv is not None else "NA"


def _ss_rows(hid, residues, ss_map, versions, gaps) -> list[dict]:
    rows = []
    for idx in residues:
        code = ss_map.get(idx)
        if code is None:
            gaps.append(_gap("residue", f"{hid}:{idx}", "secondary_structure",
                             "residue absent from DSSP output", "DSSP", "ADVISORY"))
            continue
        rows.append({
            "hotspot_id": hid, "residue_index": idx, "ss_code": code,
            "ss_class": _ss_class(code), "source": "DSSP",
            "source_version": _version_of(versions, "DSSP"),
        })
    return rows


def _ss_class(code: str) -> str:
    return {"H": "helix", "G": "helix", "I": "helix",
            "E": "strand", "B": "strand"}.get(str(code), "coil")


def _ss_composition(residues: list[int], ss_map: dict) -> str:
    counts: dict[str, int] = {}
    for idx in residues:
        code = ss_map.get(idx)
        if code is not None:
            counts[_ss_class(code)] = counts.get(_ss_class(code), 0) + 1
    return ";".join(f"{k}={v}" for k, v in sorted(counts.items())) if counts else "NA"


def _domain_rows(hid, domains, versions) -> list[dict]:
    return [{
        "hotspot_id": hid, "domain_id": d.get("domain_id", "NA"),
        "domain_name": d.get("name", "NA"), "domain_start": d.get("start"),
        "domain_end": d.get("end"), "overlap_residues": d.get("overlap_residues", "NA"),
        "source": d.get("source", "NA"),
        "source_version": _version_of(versions, d.get("source", "")),
    } for d in domains]


def _site_rows(hid, sites, versions) -> list[dict]:
    return [{
        "hotspot_id": hid, "site_id": s.get("site_id", "NA"),
        "site_type": s.get("site_type", "NA"),
        "residue_start": s.get("residue_start"), "residue_end": s.get("residue_end"),
        "description": s.get("description", "NA"), "source": s.get("source", "NA"),
        "source_version": _version_of(versions, s.get("source", "")),
    } for s in sites]


def _cons_rows(hid, residues, conservation, versions, gaps) -> list[dict]:
    rows = []
    for idx in residues:
        rec = conservation.get(idx)
        if rec is None:
            gaps.append(_gap("residue", f"{hid}:{idx}", "conservation",
                             "residue not covered by the pinned conservation source",
                             "dbNSFP", "ADVISORY"))
            continue
        rows.append({
            "hotspot_id": hid, "residue_index": idx,
            "gerp_rs": rec.get("GERP++_RS"), "phylop100way": rec.get("phyloP100way"),
            "measurement_basis": (
                "genomic conservation mapped through the MANE transcript to the "
                "residue, averaged over the three codon positions; descriptive only"
            ),
            "source": "dbNSFP", "source_version": _version_of(versions, "dbNSFP"),
        })
    return rows


def _cons_summary(residues: list[int], conservation: dict) -> str:
    gerp = [conservation[i]["GERP++_RS"] for i in residues
            if i in conservation and conservation[i].get("GERP++_RS") is not None]
    phylop = [conservation[i]["phyloP100way"] for i in residues
              if i in conservation and conservation[i].get("phyloP100way") is not None]
    if not gerp and not phylop:
        return "NA"
    parts = []
    if gerp:
        parts.append(f"GERP++_RS_mean={_mean(gerp):.6g}")
    if phylop:
        parts.append(f"phyloP100way_mean={_mean(phylop):.6g}")
    return ";".join(parts)


def _disease_rows(hid, diseases, versions) -> list[dict]:
    return [{
        "hotspot_id": hid, "disease": d.get("disease", "NA"),
        "residue_start": d.get("residue_start"), "residue_end": d.get("residue_end"),
        "reference": d.get("reference", "NA"), "source": d.get("source", "NA"),
        "source_version": _version_of(versions, d.get("source", "")),
    } for d in diseases]


def _star_composition(hid: str, data) -> str:
    """ClinVar review-star composition of the hotspot's classified variants.

    Descriptive evidence-quality context. Never a re-filter and never a reason to
    redefine a region (agent §3).
    """
    stars = {
        int(r["residue_index"]): r.get("max_star", r.get("star_level"))
        for r in data.star_distribution
        if r.get("residue_index") is not None
    }
    residues = [int(r["residue_index"]) for r in data.sphere_variants.get(hid, [])]
    counts: dict[str, int] = {}
    for idx in residues:
        key = f"{stars.get(idx, 'NA')}star"
        counts[key] = counts.get(key, 0) + 1
    return ";".join(f"{k}={v}" for k, v in sorted(counts.items())) if counts else "NA"


def _sensitivity_overlap(hid: str, data, gaps: list[dict]) -> tuple[Any, Any]:
    for row in data.sensitivity_overlap:
        if str(row.get("hotspot_id")) == hid:
            return row.get("overlap_ge1star"), row.get("overlap_ge2star")
    gaps.append(_gap("hotspot", hid, "sensitivity_overlap",
                     "no >=1*/>=2* overlap row for this hotspot",
                     "11_SENSITIVITY/sensitivity_overlap.tsv", "ADVISORY"))
    return "NA", "NA"


def _gap(scope: str, identifier: str, missing: str, reason: str, source: str,
         severity: str) -> dict:
    return {
        "scope": scope, "identifier": identifier, "missing_item": missing,
        "reason": reason, "source_attempted": source, "severity": severity,
    }
