"""``stage_e_report.md`` — the narrative, with the barrier kept visible.

Every section keeps mechanism-blind discovery (Stages B–D) and this overlay
separate, because the entire scientific value of ¶55 rests on the mechanism labels
having been unavailable to Stages B, C and D. The report therefore never says a
hotspot is "supported by" or "confirmed by" biology; it says what was found without
biology, and then what biology was laid over it afterwards.
"""
from __future__ import annotations

from typing import Any

from .rubric import RUBRIC_DEFINITIONS, Evidence

#: Caveats that travel with every Stage E result, in handoff_05 and in the report.
STANDING_CAVEATS = [
    "Hotspots and the footprint were discovered and validated WITHOUT any biological "
    "input. Every annotation in this stage is a downstream overlay and did not "
    "influence discovery, r_hot, r_fp, footprint geometry or any robustness threshold.",
    "A significant hotspot center is a GEOMETRIC TEST POSITION and need not carry a "
    "ClinVar variant. SIGNIFICANT_HOTSPOT_CENTERS, HOTSPOT_SPHERE_CLASSIFIED_VARIANTS "
    "and HOTSPOT_COVERED_RESIDUES are three distinct objects and are never collapsed.",
    "The upstream robustness result is GEOMETRIC FOOTPRINT ROBUSTNESS ONLY. A "
    "preserved footprint is not confirmation that the underlying hotspots are "
    "validated, and center-set perturbation is not validation of hotspot discovery.",
    "Every statistic in this stage is copied verbatim from 06_FINAL_HOTSPOTS/ and was "
    "verified byte-for-byte. Stage E recomputes no upstream statistic.",
    "Mechanism labels were assigned from literature evidence alone, under the frozen "
    "rubric, with no access to hotspot membership. Hotspot membership is never "
    "evidence for a functional claim, and a functional claim is never evidence that a "
    "hotspot is real.",
    "ClinVar review-star composition and the >=1*/>=2* sensitivity overlap are "
    "DESCRIPTIVE evidence-quality context only — never a re-filter and never a reason "
    "to redefine a region.",
    "Conservation is genomic conservation mapped through the MANE transcript to a "
    "residue, which is what it is. It is purely descriptive.",
]


def build_report(
    *,
    gene: str,
    run_id: str,
    annotation_rows: list[dict],
    curated: list[dict],
    counts: dict[str, int],
    gate_numbers: dict[str, Any],
    posthoc: dict[str, Any] | None,
    qc: list[dict],
    source_versions: dict,
    gaps: list[dict],
    warnings: list,
    synthetic: bool,
) -> str:
    """Assemble the Stage E report."""
    lines: list[str] = []
    add = lines.append

    add(f"# Stage E — Biological Annotation: {gene}")
    add("")
    add(f"- **Run ID:** `{run_id}`")
    add(f"- **Stage:** `10_ANNOTATION` (terminal; nothing flows backwards)")
    add(f"- **Hotspots annotated:** {len(annotation_rows)}")
    add(f"- **Mechanism records curated:** {len(curated)}")
    if synthetic:
        add("- **INPUT MODE: SYNTHETIC/MOCK.** Annotation content is not real "
            "biological data and must not be interpreted as such.")
    add("")

    add("## 1. What this stage did and did not do")
    add("")
    add("Stages B, C and D discovered the hotspots and evaluated footprint robustness "
        "**without any biological information whatsoever** — no domains, no "
        "conservation, no literature, no mechanism labels. This stage runs strictly "
        "afterwards and only describes what those stages found.")
    add("")
    add("The separation is not a convention here; it is enforced structurally. "
        "Mechanism labels are assigned by a rubric function that receives literature "
        "records and nothing else, and the handoff is scanned for upstream-parameter "
        "keys before it is written.")
    add("")
    for caveat in STANDING_CAVEATS:
        add(f"- {caveat}")
    add("")

    add("## 2. Mechanism-blind discovery (Stages B–D) — reported, not re-derived")
    add("")
    add("| Hotspot | Interval | n_P/LP | n_B/LB | Fold enr. | p_emp | q_BH |")
    add("|---|---|---|---|---|---|---|")
    for row in annotation_rows:
        add(
            f"| {row['hotspot_id']} | {row['covered_interval']} | {row['n_plp']} | "
            f"{row['n_blb']} | {row['fold_enrichment']} | {row['p_emp']} | "
            f"{row['q_bh']} |"
        )
    add("")
    add("Every value above was copied verbatim from `06_FINAL_HOTSPOTS/` as raw text "
        "and verified byte-for-byte against the source file. None was recomputed.")
    add("")

    add("### 2.1 The three residue objects, kept separate")
    add("")
    add("| Hotspot | Significant centers | Sphere-classified variants | Covered residues |")
    add("|---|---|---|---|")
    for row in annotation_rows:
        add(
            f"| {row['hotspot_id']} | {row['n_significant_hotspot_centers']} | "
            f"{row['n_hotspot_sphere_classified_variants']} | "
            f"{row['n_hotspot_covered_residues']} |"
        )
    add("")
    add("**A significant hotspot center is a geometric test position and need not "
        "carry a ClinVar variant.** These three counts answer three different "
        "questions and are never added together or substituted for one another.")
    add("")

    add("### 2.2 Robustness and upstream flags")
    add("")
    add("Every annotated row carries the robustness profile summary, the per-center "
        "influence summary, the center recurrence rate and all four upstream flags. "
        "This is a schema requirement — a row missing any of them cannot be written.")
    add("")
    add("| Hotspot | Robustness profile | Per-center influence | Center recurrence | "
        "Global clustering | Perm. res. limited | Fallback radius | Domain boundary |")
    add("|---|---|---|---|---|---|---|---|")
    for row in annotation_rows:
        add(
            f"| {row['hotspot_id']} | {row['robustness_profile_summary']} | "
            f"{row['per_center_influence_summary']} | {row['center_recurrence_rate']} | "
            f"{row['global_clustering_flag']} | {row['permutation_resolution_limited']} | "
            f"{row['fallback_radius_domain']} | {row['domain_boundary_warning']} |"
        )
    add("")
    add("The robustness profile describes **geometric footprint preservation under "
        "center-set perturbation**. It is not a verdict on whether the hotspots are "
        "real, and Stage E emits no ROBUST / NOT_ROBUST verdict of its own.")
    add("")

    add("## 3. Biological overlay")
    add("")
    add("### 3.1 Annotation sources (all version-pinned)")
    add("")
    add("| Source | Version | URL | Accession |")
    add("|---|---|---|---|")
    for name in sorted(source_versions):
        sv = source_versions[name]
        add(f"| {sv.name} | `{sv.version}` | {sv.url} | {sv.accession} |")
    add("")
    add("An unversioned source may not be used; that condition is BLOCKING.")
    add("")

    add("### 3.2 Functional mechanism curation (¶55)")
    add("")
    add("Applied rubric, reproduced as executed:")
    add("")
    add("| Grade | Definition |")
    add("|---|---|")
    for grade in (Evidence.STRONG, Evidence.MODERATE, Evidence.WEAK, Evidence.NONE):
        add(f"| `{grade.value}` | {RUBRIC_DEFINITIONS[grade]} |")
    add("")
    add("A mechanism label requires **>= MODERATE**. Conflicting STRONG sources resolve "
        "to **Mixed** with both references retained. WEAK-only or conflicting-WEAK "
        "resolve to **Unclear**. No evidence gives "
        "**Not_Experimentally_Characterized** — a real category, never filled by "
        "inference.")
    add("")
    add("| Mechanism | Variants |")
    add("|---|---|")
    for name, n in counts.items():
        add(f"| {name} | {n} |")
    add("")

    add("### 3.3 Post hoc mechanism-spatial gate")
    add("")
    add(f"- Gate: at least **{gate_numbers['threshold_min_variants_per_category']}** "
        f"variants at **>= {gate_numbers['threshold_min_evidence']}** evidence in each "
        f"of at least **{gate_numbers['threshold_min_categories']}** categories.")
    add(f"- Observed counts: `{gate_numbers['counts_by_category_at_or_above_min_evidence']}`")
    add(f"- Qualifying categories: {gate_numbers['n_qualifying_categories']} "
        f"({', '.join(gate_numbers['qualifying_categories']) or 'none'})")
    add(f"- Eligible variants: {gate_numbers['n_eligible_variants']}")
    add(f"- **Gate passed: {gate_numbers['passed']}**")
    add("")
    add(gate_numbers["reason"])
    add("")
    if posthoc and posthoc.get("status") != "NOT_RUN":
        add("The analysis ran and is labelled `analysis_type: post_hoc`, corrected "
            "separately from Stage B, and is **hypothesis-generating only**. It feeds "
            "nothing upstream.")
        add("")
        add("| Mechanism | n | Mean intra-category distance (A) | p_emp | q_BH | Significant |")
        add("|---|---|---|---|---|---|")
        for cat in posthoc.get("categories", []):
            add(
                f"| {cat['mechanism']} | {cat['n_variants']} | "
                f"{cat['observed_mean_intra_distance_A']} | {cat['p_emp']} | "
                f"{cat['q_bh']} | {cat['significant']} |"
            )
        add("")
    else:
        add("The analysis was **NOT RUN**. Below the gate the test has no power. This "
            "is recorded as a gate failure with its numbers, not as an omission, and a "
            "near-miss is not a pass.")
        add("")

    add("## 4. Annotation gaps")
    add("")
    if gaps:
        add(f"{len(gaps)} gap(s) recorded in `annotation_gaps.tsv`. A missing "
            "annotation is recorded as a gap and is never filled by inference.")
        add("")
        add("| Scope | Identifier | Missing | Reason |")
        add("|---|---|---|---|")
        for gap in gaps[:50]:
            add(f"| {gap['scope']} | {gap['identifier']} | {gap['missing_item']} | "
                f"{gap['reason']} |")
    else:
        add("**No gaps were recorded, which is itself suspicious and requires "
            "justification.** Complete annotation coverage of every hotspot from every "
            "source is unusual for real data; it is expected only when running against "
            "a synthetic source constructed to be complete. Justification for this run: "
            + ("synthetic/mock annotation source with deliberately complete coverage."
               if synthetic else
               "NONE RECORDED — the Lead must confirm this before publication."))
    add("")

    add("## 5. QC results")
    add("")
    add("| Check | Result | Detail |")
    add("|---|---|---|")
    for item in qc:
        detail = {k: v for k, v in item.items() if k not in ("check", "result")}
        add(f"| `{item['check']}` | **{item['result']}** | `{detail}` |")
    add("")

    add("## 6. Warnings")
    add("")
    if warnings:
        add("| Severity | Code | Message |")
        add("|---|---|---|")
        for w in warnings:
            add(f"| {w.severity.value} | `{w.warning_code}` | {w.message} |")
    else:
        add("No warnings raised.")
    add("")

    add("## 7. Non-actionable observations for the Lead")
    add("")
    add("If the literature suggests a region this analysis did not identify, that is a "
        "**reported observation only**, explicitly non-actionable for the current run. "
        "It is never a modification request and never a pipeline change. No output of "
        "this stage may reach hotspot discovery, `r_hot`, `r_fp`, footprint geometry "
        "or any robustness threshold — in this run or through any side channel.")
    add("")

    return "\n".join(lines) + "\n"
