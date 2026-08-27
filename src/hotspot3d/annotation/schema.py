"""Frozen output schemas and the structural guardrails — Output Contract Part IX.

Two things are enforced here that agent §8 states as scientific requirements, so
that violating them is a crash rather than a review finding:

  * **Robustness travels with the claim.** Every hotspot row must carry the
    robustness-profile summary, per-center influence, center recurrence and all four
    upstream flags. :func:`assert_mandatory_context` makes this a schema requirement,
    not a convention — a row missing a flag cannot be written.
  * **One writer per path.** :func:`assert_write_target` refuses any path outside
    ``10_ANNOTATION/`` (plus the single owned provenance file), so a coding slip
    cannot write into an upstream stage directory.
"""
from __future__ import annotations

import re
from pathlib import Path

from ..utils.errors import BlockedError, LeakageError

AGENT = "biological-annotation"
STAGE = "10_annotation"
STAGE_DIR = "10_ANNOTATION"

#: The four upstream flags. Flags are never dropped (agent §8) and never defaulted.
UPSTREAM_FLAGS = (
    "global_clustering_flag",
    "permutation_resolution_limited",
    "fallback_radius_domain",
    "domain_boundary_warning",
)

#: Robustness context that must accompany every hotspot-level statement (agent §8).
ROBUSTNESS_CONTEXT = (
    "robustness_profile_summary",
    "per_center_influence_summary",
    "center_recurrence_rate",
)

MANDATORY_ROW_CONTEXT = UPSTREAM_FLAGS + ROBUSTNESS_CONTEXT


# --- frozen column orders ---------------------------------------------------

HOTSPOT_ANNOTATION_COLUMNS = [
    "hotspot_id",
    "residue_start",
    "residue_end",
    "covered_interval",
    # --- copied verbatim from 06_FINAL_HOTSPOTS ---
    "n_plp",
    "n_blb",
    "fold_enrichment",
    "p_emp",
    "q_bh",
    # --- the three residue objects, reported separately, never collapsed ---
    "n_significant_hotspot_centers",
    "n_hotspot_sphere_classified_variants",
    "n_hotspot_covered_residues",
    # --- structural annotation ---
    "plddt_mean",
    "plddt_min",
    "secondary_structure_composition",
    "overlapping_domains",
    "functional_regions",
    "ligand_sites",
    "protein_interaction_sites",
    "conservation_summary",
    "disease_associations",
    # --- descriptive evidence quality (never a re-filter) ---
    "star_composition",
    "overlap_ge1star",
    "overlap_ge2star",
    # --- robustness + flags: mandatory on every row ---
    "center_recurrence_rate",
    "robustness_profile_summary",
    "per_center_influence_summary",
    "global_clustering_flag",
    "permutation_resolution_limited",
    "fallback_radius_domain",
    "domain_boundary_warning",
]

SECONDARY_STRUCTURE_COLUMNS = [
    "hotspot_id", "residue_index", "ss_code", "ss_class", "source", "source_version",
]

DOMAINS_COLUMNS = [
    "hotspot_id", "domain_id", "domain_name", "domain_start", "domain_end",
    "overlap_residues", "source", "source_version",
]

FUNCTIONAL_SITES_COLUMNS = [
    "hotspot_id", "site_id", "site_type", "residue_start", "residue_end",
    "description", "source", "source_version",
]

CONSERVATION_COLUMNS = [
    "hotspot_id", "residue_index", "gerp_rs", "phylop100way",
    "measurement_basis", "source", "source_version",
]

DISEASE_ASSOCIATIONS_COLUMNS = [
    "hotspot_id", "disease", "residue_start", "residue_end", "reference",
    "source", "source_version",
]

FUNCTIONAL_MECHANISM_COLUMNS = [
    "residue_index", "aa_change", "mechanism", "reference", "experimental_system",
    "assay_type", "principal_finding", "evidence_strength", "curator_note",
]

MECHANISM_BY_HOTSPOT_COLUMNS = [
    "hotspot_id", "mechanism", "n_variants", "n_at_or_above_moderate",
    "residues", "robustness_profile_summary", "global_clustering_flag",
]

ANNOTATION_GAPS_COLUMNS = [
    "scope", "identifier", "missing_item", "reason", "source_attempted", "severity",
]

LITERATURE_LOG_COLUMNS = [
    "query", "database", "search_date", "hit_count", "record_reference",
    "residue_index", "decision", "decision_reason",
]

STAGE_MANIFEST_COLUMNS = [
    "file", "kind", "status", "reason",
]


# --- structural guardrails --------------------------------------------------

def assert_mandatory_context(rows: list[dict], *, where: str = "hotspot_annotation.tsv") -> None:
    """Every hotspot row carries robustness context and all four upstream flags.

    ``NA`` is a legal value for a *measurement* that is genuinely unavailable, but not
    for these fields: a dropped flag is indistinguishable from an absent one to a
    reader, which is exactly the failure agent §8 forbids.
    """
    for row in rows:
        hid = row.get("hotspot_id", "<unknown>")
        missing = [k for k in MANDATORY_ROW_CONTEXT if k not in row]
        if missing:
            raise BlockedError(
                f"SCHEMA violation in {where}: hotspot {hid!r} is missing mandatory "
                f"context {missing}. Robustness and upstream flags travel with every "
                f"claim (agent §8); a row cannot be written without them."
            )
        empty = [
            k for k in MANDATORY_ROW_CONTEXT
            if row[k] is None or (isinstance(row[k], str) and row[k].strip() in ("", "NA"))
        ]
        if empty:
            raise BlockedError(
                f"SCHEMA violation in {where}: hotspot {hid!r} has empty/NA mandatory "
                f"context {empty}. Flags are never dropped and never defaulted."
            )


_STAGE_DIR_RE = re.compile(r"(?:^|/)(0[1-9]|1[01])_[A-Z_]+(?:/|$)")


def assert_write_target(run_root: str | Path, path: str | Path) -> Path:
    """Refuse to write anywhere Stage E does not own.

    Permitted: ``FULL_RESULTS/10_ANNOTATION/**`` and the single owned provenance file
    ``FULL_RESULTS/12_REPRODUCIBILITY/provenance/provenance_10_annotation.json``.
    Anything under ``01_*``–``09_*`` or ``11_*`` is a leakage event (agent §6.6).
    """
    run_root = Path(run_root).resolve()
    target = Path(path).resolve()

    try:
        rel = target.relative_to(run_root)
    except ValueError as exc:
        raise LeakageError(
            f"OWNERSHIP violation: {target} is outside the run root {run_root}. "
            f"Stage E writes only inside its own run."
        ) from exc

    rel_str = str(rel).replace("\\", "/")

    if rel_str.startswith(f"FULL_RESULTS/{STAGE_DIR}/") or rel_str == f"FULL_RESULTS/{STAGE_DIR}":
        return target

    if rel_str == (
        "FULL_RESULTS/12_REPRODUCIBILITY/provenance/provenance_10_annotation.json"
    ):
        return target

    hit = _STAGE_DIR_RE.search(rel_str)
    if hit:
        raise LeakageError(
            f"OWNERSHIP violation: Stage E attempted to write {rel_str}, which belongs "
            f"to stage {hit.group(0).strip('/')!r}. Nothing downstream may modify "
            f"Stage A–D artifacts (agent §11)."
        )

    raise LeakageError(
        f"OWNERSHIP violation: Stage E may write only FULL_RESULTS/{STAGE_DIR}/** and "
        f"its own provenance file; refused {rel_str}."
    )


#: handoff_05 payload keys, frozen. Anything else is a contract change.
HANDOFF_05_KEYS = (
    "annotation_sources",
    "n_hotspots_annotated",
    "mechanism_counts",
    "posthoc_gate_passed",
    "posthoc_gate_numbers",
    "robustness_profile",
    "global_clustering_flag",
    "interpretation_caveats",
    # Stage E additions required by agent §11 (upstream hashes + terminality proof).
    "upstream_handoff_hashes",
    "annotation_source_versions",
    "residue_objects",
    "terminal_stage",
    "may_influence_discovery",
)

#: Keys whose contents are a byte-identical echo of upstream and are therefore
#: exempt from the radius-name scan below.
HANDOFF_05_PASSTHROUGH = ("robustness_profile", "upstream_handoff_hashes")

#: Names that would let a downstream reader set an upstream parameter. Stage E is
#: terminal: none of these may appear in anything it authors.
FORBIDDEN_UPSTREAM_KEYS = (
    "r_hot", "r_fp", "hotspot_radius", "footprint_radius", "radius",
    "preservation_threshold", "robustness_threshold", "q", "fdr",
    "recommended_radius", "suggested_radius", "revised_radius",
)


def assert_no_upstream_influence(payload: dict) -> None:
    """Prove Stage E's handoff cannot carry an upstream directive (agent §16).

    Scans everything Stage E *authors*. Upstream passthrough blocks are exempt because
    they are echoes of another stage's own output, verified elsewhere to be identical.
    """
    unknown = sorted(set(payload) - set(HANDOFF_05_KEYS))
    if unknown:
        raise BlockedError(
            f"BLOCKED — handoff_05 carries undeclared key(s) {unknown}. The payload "
            f"schema is frozen; a new key is a contract change, not an implementation "
            f"detail."
        )

    offending: list[str] = []

    def walk(node, trail: str) -> None:
        if isinstance(node, dict):
            for key, value in node.items():
                if str(key).lower() in FORBIDDEN_UPSTREAM_KEYS:
                    offending.append(f"{trail}.{key}" if trail else str(key))
                walk(value, f"{trail}.{key}" if trail else str(key))
        elif isinstance(node, (list, tuple)):
            for i, value in enumerate(node):
                walk(value, f"{trail}[{i}]")

    for key, value in payload.items():
        if key in HANDOFF_05_PASSTHROUGH:
            continue
        if str(key).lower() in FORBIDDEN_UPSTREAM_KEYS:
            offending.append(str(key))
        walk(value, str(key))

    if offending:
        raise LeakageError(
            f"BACK-PROPAGATION violation: handoff_05 contains upstream-parameter "
            f"key(s) {sorted(set(offending))}. Stage E is terminal — no output may "
            f"reach hotspot discovery, r_hot, r_fp, footprint geometry or any "
            f"robustness threshold (agent §2, §16)."
        )
