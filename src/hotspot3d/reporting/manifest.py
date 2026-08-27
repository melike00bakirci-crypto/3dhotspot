"""Run-level manifests, warnings and gene_summary — Lead-owned aggregation.

P1  REVIEW_PACK is derived, never recomputed.
P2  One writer per canonical file. Run-level aggregates are pure CONCATENATIONS
    or verbatim copies of agent-owned files — the Lead never authors a scientific
    number, so no value can diverge between inspection levels.
P3  Absence is explicit: the manifest lists EXPECTED files, so anything not
    produced appears as a NOT_CREATED row with a reason.
"""
from __future__ import annotations

from pathlib import Path

from ..utils.hashing import sha256_file
from ..utils.io import write_tsv
from ..utils.runctx import STAGE_OWNERS, RunContext
from ..utils.status import WARNING_COLUMNS, utc_now

MANIFEST_COLUMNS = [
    "relative_path", "stage", "agent_owner", "description", "file_type", "tier",
    "in_review_pack", "bytes", "sha256", "status", "created_utc",
]

# Expected canonical outputs per stage (Output Contract IX.3). A file listed here
# and absent from disk becomes a visible NOT_CREATED row rather than a silent gap.
EXPECTED: dict[str, list[tuple[str, str, str]]] = {
    "01_INPUT_RAW": [
        ("clinvar_raw.tsv", "PRIMARY", "immutable retrieved ClinVar payload"),
        ("alphafold_model.cif", "PRIMARY", "immutable retrieved structure"),
        ("uniprot_canonical.fasta", "PRIMARY", "immutable canonical sequence"),
        ("retrieval_log.json", "PRIMARY", "URLs, queries, release dates, HTTP status"),
    ],
    "02_CLINVAR": [
        ("variants_missense_all.tsv", "PRIMARY", "stratum 1 — all retrieved missense"),
        ("variants_residue_level.tsv", "PRIMARY", "stratum 2 — residue-level cohort"),
        ("variants_excluded_from_primary.tsv", "PRIMARY", "stratum 3 — exclusions with reasons"),
        ("review_star_distribution.tsv", "PRIMARY", "stratum 4 — star distribution"),
        ("cohort_summary.json", "PRIMARY", "cohort headline counts"),
        ("residue_class_conflicts.tsv", "PRIMARY", "F3 RESIDUE_CLASS_CONFLICT rows"),
    ],
    "03_STRUCTURE_QC": [
        ("structure_qc.json", "PRIMARY", "structural QC verdict"),
        ("residue_coordinates.tsv", "PRIMARY", "Ca coordinates + pLDDT"),
        ("positional_universe.tsv", "PRIMARY", "U_struct"),
        ("classified_cohort.tsv", "PRIMARY", "L"),
        ("plddt_profile.tsv", "SECONDARY", "per-residue pLDDT"),
        ("plddt_regions.tsv", "SECONDARY", "contiguous low-confidence regions"),
        ("mapping_report.json", "PRIMARY", "ClinVar <-> structure mapping"),
        ("handoff_01.json", "PRIMARY", "A -> B handoff"),
        ("stage_a_report.md", "PRIMARY", "Stage A narrative report"),
    ],
    "04_GLOBAL_CLUSTERING": [
        ("ripleys_k_plp.tsv", "PRIMARY", "Ripley K, P/LP cohort"),
        ("ripleys_k_blb.tsv", "PRIMARY", "Ripley K, B/LB cohort"),
        ("ripleys_k_summary.json", "PRIMARY", "global max|Z| test"),
        ("pair_correlation_plp.tsv", "PRIMARY", "g(r), P/LP"),
        ("pair_correlation_blb.tsv", "PRIMARY", "g(r), B/LB"),
        ("pair_correlation_peaks.json", "SECONDARY", "PCF peaks (informative only)"),
        ("candidate_radius_domain.json", "PRIMARY", "derived domain + fallback trigger"),
    ],
    "05_HOTSPOT_RADIUS": [
        ("r_hot_scan.tsv", "PRIMARY", "decision trace incl. rejected candidates"),
        ("loo_diagnostics.tsv", "PRIMARY", "sparse/zero-neighbour/tie reporting"),
        ("objective_matrix.json", "PRIMARY", "objectives + correlation diagnostic"),
        ("pareto_front.json", "PRIMARY", "non-dominated candidates"),
        ("pareto_dominated.json", "SECONDARY", "dominated candidates with dominators"),
        ("radius_decision.json", "PRIMARY", "normalization, utopia, distances, tie chain"),
        ("domain_boundary_diagnostic.json", "PRIMARY", "BW-1/2/3 evaluation"),
    ],
    "06_FINAL_HOTSPOTS": [
        ("all_residue_center_tests.tsv", "PRIMARY", "every tested candidate center"),
        ("significant_hotspot_centers.tsv", "PRIMARY", "S — BH survivors"),
        ("hotspot_classified_variants.tsv", "PRIMARY", "classified residues in spheres"),
        ("hotspot_covered_residues.tsv", "PRIMARY", "all covered residues"),
        ("hotspot_regions.tsv", "PRIMARY", "connected components of S"),
        ("bh_fdr_table.tsv", "PRIMARY", "BH ranks and boundary, BY comparison"),
        ("permutation_resolution_diagnostic.json", "PRIMARY", "resolution limitation check"),
        ("handoff_02.json", "PRIMARY", "B -> C handoff"),
        ("stage_b_report.md", "PRIMARY", "Stage B narrative report"),
    ],
    "07_FOOTPRINT_RADIUS": [
        ("r_fp_scan.tsv", "PRIMARY", "footprint decision trace"),
        ("footprint_domain.json", "PRIMARY", "MST merge scales and domain"),
        ("footprint_admissibility.tsv", "PRIMARY", "QC verdict per candidate"),
        ("merge_events.tsv", "SECONDARY", "component merge events"),
        ("footprint_objective_matrix.json", "PRIMARY", "objectives + correlation"),
        ("footprint_pareto.json", "PRIMARY", "non-dominated candidates"),
        ("footprint_pareto_dominated.json", "SECONDARY", "dominated with dominators"),
        ("footprint_decision.json", "PRIMARY", "selection trace"),
        ("footprint_domain_boundary_diagnostic.json", "PRIMARY", "BW-1/2/3 evaluation"),
    ],
    "08_FINAL_FOOTPRINT": [
        ("footprint_residues.tsv", "PRIMARY", "FP_original residues"),
        ("footprint_geometry.tsv", "PRIMARY", "final geometric descriptors"),
        ("footprint_components.tsv", "PRIMARY", "final component membership"),
        ("footprint_freeze.json", "PRIMARY", "FREEZE GATE hashes"),
        ("handoff_03.json", "PRIMARY", "C -> D handoff"),
        ("stage_c_report.md", "PRIMARY", "Stage C narrative report"),
    ],
    "09_ROBUSTNESS": [
        ("baseline_reproduction.json", "PRIMARY", "mandatory baseline check"),
        ("perturbation_universe.tsv", "PRIMARY", "S with variant-carrying flag"),
        ("perturbation_design.tsv", "PRIMARY", "per-level exhaustive/sampled design"),
        ("perturbation_results.tsv", "PRIMARY", "per-iteration metrics"),
        ("iteration_failures.tsv", "PRIMARY", "failed iterations with cause"),
        ("center_sensitivity.tsv", "PRIMARY", "per-center influence ranking"),
        ("robustness_summary.tsv", "PRIMARY", "aggregates with BCa CIs"),
        ("robustness_profile.json", "PRIMARY", "continuous robustness profile"),
        ("handoff_04.json", "PRIMARY", "D -> E handoff"),
        ("stage_d_report.md", "PRIMARY", "Stage D narrative report"),
    ],
    "10_ANNOTATION": [
        ("hotspot_annotation.tsv", "PRIMARY", "per-hotspot biological annotation"),
        ("functional_mechanism_variants.tsv", "PRIMARY", "GOF/LOF/DN evidence ledger"),
        ("mechanism_by_hotspot.tsv", "PRIMARY", "mechanism distribution"),
        ("annotation_gaps.tsv", "PRIMARY", "explicit annotation gaps"),
        ("literature_search_log.tsv", "PRIMARY", "reproducible search log"),
        ("handoff_05.json", "PRIMARY", "E -> Lead handoff"),
        ("stage_e_report.md", "PRIMARY", "Stage E narrative report"),
    ],
    "11_SENSITIVITY": [
        ("sensitivity_plddt70.json", "PRIMARY", "pLDDT >= 70 sensitivity"),
        ("sensitivity_review_status.json", "PRIMARY", ">=1* / >=2* sensitivity"),
        ("sensitivity_overlap.tsv", "PRIMARY", "overlap vs primary S"),
        ("SENSITIVITY_IS_NON_REDEFINING.txt", "PRIMARY", "boundary statement"),
    ],
    "12_REPRODUCIBILITY": [
        ("config.yaml", "PRIMARY",
         "effective merged configuration the run actually read (overlay included)"),
        ("config_base.yaml", "PRIMARY",
         "base config/pipeline.yaml verbatim, with its FROZEN annotations"),
        # config_overlay.yaml is deliberately unregistered: it exists only for an
        # overlaid run, and listing it here would emit a NOT_CREATED row on every
        # ordinary run — absence would read as a missing output rather than as the
        # normal case. When present it is picked up by the directory walk.
        ("run_metadata.json", "PRIMARY", "run identity and status"),
        ("software_versions.tsv", "PRIMARY", "package versions"),
        ("random_seeds.tsv", "PRIMARY", "every derived seed with its context"),
        ("checksums.tsv", "PRIMARY", "hash of every produced file"),
        ("INTERMEDIATE_RETENTION.md", "PRIMARY", "regeneration recipes for omitted data"),
    ],
}

REVIEW_PACK_MEMBERS = {
    "r_hot_scan.tsv", "r_fp_scan.tsv", "significant_hotspot_centers.tsv",
    "hotspot_regions.tsv", "footprint_residues.tsv", "robustness_profile.json",
    "hotspot_annotation.tsv",
}

# Stage directories that are always created so no mysteriously empty folder exists (P3).
ALL_STAGES = [s for s in STAGE_OWNERS if s != "00_RUN_SUMMARY"]


def _stage_of(rel: Path) -> str:
    parts = rel.parts
    if len(parts) >= 2 and parts[0] == "FULL_RESULTS":
        return parts[1]
    return "RUN_ROOT"


def build_manifest(ctx: RunContext, not_created_reasons: dict[str, str] | None = None
                   ) -> list[dict]:
    """One row per file — produced AND expected-but-absent."""
    reasons = not_created_reasons or {}
    rows: list[dict] = []
    seen: set[str] = set()

    for path in sorted(ctx.run_root.rglob("*")):
        if not path.is_file() or "__pycache__" in path.parts:
            continue
        rel = path.relative_to(ctx.run_root)
        if rel.parts and rel.parts[0] == "ARCHIVES":
            continue
        stage = _stage_of(rel)
        stat = path.stat()
        seen.add(path.name)
        rows.append({
            "relative_path": str(rel),
            "stage": stage,
            "agent_owner": STAGE_OWNERS.get(stage, "lead"),
            "description": _describe(stage, path.name),
            "file_type": path.suffix.lstrip(".") or "none",
            "tier": _tier(stage, path.name),
            "in_review_pack": path.name in REVIEW_PACK_MEMBERS or rel.parts[0] == "REVIEW_PACK",
            "bytes": stat.st_size,
            "sha256": sha256_file(path),
            "status": "CREATED" if stat.st_size > 0 else "EMPTY",
            "created_utc": utc_now(),
        })

    for stage, expected in EXPECTED.items():
        stage_dir = ctx.full_results / stage
        for name, tier, desc in expected:
            if (stage_dir / name).exists():
                continue
            reason = reasons.get(f"{stage}/{name}") or reasons.get(stage) or \
                "not produced — stage did not reach this output"
            rows.append({
                "relative_path": f"FULL_RESULTS/{stage}/{name}",
                "stage": stage, "agent_owner": STAGE_OWNERS.get(stage, "lead"),
                "description": f"{desc} [NOT CREATED: {reason}]",
                "file_type": Path(name).suffix.lstrip("."), "tier": tier,
                "in_review_pack": name in REVIEW_PACK_MEMBERS,
                "bytes": 0, "sha256": "NA", "status": "NOT_CREATED",
                "created_utc": utc_now(),
            })

    return sorted(rows, key=lambda r: r["relative_path"])


def _describe(stage: str, name: str) -> str:
    for fname, _tier, desc in EXPECTED.get(stage, []):
        if fname == name:
            return desc
    if name == "stage_status.json":
        return "stage status and outcome type"
    if name.startswith("warnings_"):
        return "per-stage warnings"
    if name == "stage_manifest.tsv":
        return "per-stage manifest"
    if name == "NOT_RUN.txt":
        return "human-readable explanation of a non-completed stage"
    return "pipeline output"


def _tier(stage: str, name: str) -> str:
    for fname, tier, _desc in EXPECTED.get(stage, []):
        if fname == name:
            return tier
    if "/" in name or name.endswith((".png", ".svg")):
        return "SECONDARY"
    return "INTERMEDIATE"


def write_run_manifest(ctx: RunContext, rows: list[dict]) -> Path:
    return write_tsv(ctx.run_root / "manifest.tsv", rows, MANIFEST_COLUMNS)


def aggregate_warnings(ctx: RunContext) -> list[dict]:
    """Concatenate every per-stage warnings file, severity-sorted (P2)."""
    from ..utils.io import read_tsv
    order = {"BLOCKING": 0, "MAJOR": 1, "ADVISORY": 2, "INFO": 3}
    rows: list[dict] = []
    for path in sorted(ctx.full_results.rglob("warnings_*.tsv")):
        for row in read_tsv(path):
            rows.append({k: row.get(k) for k in WARNING_COLUMNS})
    return sorted(rows, key=lambda r: (order.get(str(r.get("severity")), 9),
                                       str(r.get("stage")), str(r.get("warning_code"))))


def write_run_warnings(ctx: RunContext, rows: list[dict]) -> Path:
    return write_tsv(ctx.run_root / "warnings.tsv", rows, WARNING_COLUMNS)


GENE_SUMMARY_COLUMNS = [
    "gene", "uniprot_acc", "mane_transcript", "run_id", "run_status", "outcome_type",
    "clinvar_release", "alphafold_model_version", "n_missense_retrieved",
    "n_residues_plp", "n_residues_blb", "n_residue_class_conflict",
    "n_positional_universe", "global_clustering_plp_p", "global_clustering_blb_p",
    "r_hot", "r_hot_domain_source", "fallback_radius_domain", "loo_mcc",
    "perm_evidence_z", "fold_enrichment", "neighbor_stability", "distance_to_ideal_hot",
    "n_pareto_hot", "permutation_count", "permutation_resolution_limited",
    "b_recommended", "fdr_method", "fdr_q", "n_significant_centers",
    "n_hotspot_regions", "n_centers_without_variant", "r_fp", "footprint_volume_A3",
    "footprint_surface_A2", "footprint_n_residues", "footprint_coverage",
    "footprint_n_components", "distance_to_ideal_fp", "n_perturbations_possible",
    "n_perturbations_evaluated", "perturbation_mode", "jaccard_median",
    "jaccard_iqr_lo", "jaccard_iqr_hi", "preservation_freq_050",
    "preservation_freq_070", "n_warnings_major", "n_warnings_blocking",
    "review_pack_bytes", "full_results_bytes",
]


def write_gene_summary(ctx: RunContext, metrics: dict) -> Path:
    """Exactly one data row per gene/run — written even for negative runs (IX.8.5)."""
    row = {c: metrics.get(c) for c in GENE_SUMMARY_COLUMNS}
    row["gene"] = ctx.gene
    row["run_id"] = ctx.run_id
    return write_tsv(ctx.run_root / "gene_summary.tsv", [row], GENE_SUMMARY_COLUMNS)


def dir_size(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())
