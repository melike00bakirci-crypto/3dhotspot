"""``stage_d_report.md`` — the narrative record of the robustness measurement.

Instability is stated plainly. There is no summary label, no verdict, and no
sentence that converts a low preservation frequency into a reassuring adjective.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..utils.io import write_text


def _fmt(value, spec: str = ".4g") -> str:
    if value is None:
        return "NA"
    if isinstance(value, float) and not np.isfinite(value):
        return "NA"
    return format(value, spec) if isinstance(value, (int, float)) else str(value)


def write_stage_d_report(path: Path, profile: dict, design, params, results,
                         sensitivity: list[dict], phase_c, warnings: list,
                         baseline: dict, run_id: str, gene: str) -> Path:
    lines: list[str] = []
    add = lines.append
    iterations = profile["iterations"]

    add(f"# Stage D — geometric footprint robustness ({gene})")
    add("")
    add(f"- run_id: `{run_id}`")
    add(f"- footprint code_version: `{phase_c.code_version}`")
    add(f"- validation target: **FP_original at r_fp = "
        f"{_fmt(phase_c.r_fp)} A** (frozen before this phase began)")
    add("")
    add("> This result is **geometric footprint robustness only**. It is not "
        "validation of hotspot discovery and must never be presented as hotspot "
        "leave-one-out validation. No Ripley's K, pair correlation, radius scan, "
        "`r_hot` selection, permutation test or FDR was recomputed at any point, "
        "including inside an iteration.")
    add("")

    add("## QC RESULTS")
    add("")
    add(f"- baseline reproduction: "
        f"**{'PASSED' if baseline['byte_identical'] else 'FAILED'}** — the frozen "
        f"footprint API reproduced `FP_original` byte-identically from the original "
        f"inputs before iteration 1")
    add(f"- planned / completed / failed: **{iterations['planned']} / "
        f"{iterations['completed']} / {iterations['failed']}**")
    add(f"- `planned = completed + failed` asserted; no iteration was removed from "
        f"any denominator")
    add("- freeze hashes for `07_FOOTPRINT_RADIUS/` and `08_FINAL_FOOTPRINT/` "
        "re-verified at stage entry and at stage end")
    add("- classified cohort fingerprint re-verified at stage end: unchanged")
    flags = profile.get("evidence_flags", {})
    if flags.get("minimal_design"):
        add(f"- **minimal design**: {flags['minimal_design_detail']}")
    if flags.get("n_iterations_failed"):
        marker = "**" if flags.get("iterations_failed_en_masse") else ""
        add(f"- {marker}iteration failures{marker}: {flags['failure_detail']}")
    add("")

    add("## SCIENTIFIC DECISIONS")
    add("")
    add(f"- `n_S` = **{design.n_S}** significant hotspot centers")
    add(f"- `K_MAX` = **{design.K_MAX}** (= min(n_S - 1, floor(n_S / 2)))")
    add(f"- total subset space `sum_k T(k)` = **{design.total_subset_space:,}**; "
        f"evaluated = **{design.total_evaluated:,}** (N_CAP = {params.n_cap:,})")
    add("")
    add("| k | T(k) | evaluated | mode | sampling fraction | seed context |")
    add("|---|---|---|---|---|---|")
    for level in design.levels:
        add(f"| {level.k} | {level.total_subsets:,} | {level.n_evaluated:,} | "
            f"{level.mode} | {_fmt(level.sampling_fraction)} | "
            f"`{level.seed_context}` |")
    add("")
    add("The exhaustive-versus-sampled split was decided by the combinatorial "
        "arithmetic against `N_CAP` alone, recorded in full before iteration 1, and "
        "never re-drawn, re-seeded or truncated in response to results.")
    add("")
    add("**Perturbation acted on geometric center positions.** `S_iter = S \\ T`. "
        "The clinical dataset was bit-identical across every iteration: no ClinVar "
        "record was removed and no class label was altered.")
    add("")

    inv = profile["robustness_r_fp_invariant"]
    add("### Footprint reconstruction radius")
    add("")
    add(f"- DECISION-STAGE-D-FIXED-RFP-0001: every iteration reconstructs at the "
        f"frozen `r_fp` = {_fmt(inv['r_fp_original'])} A selected once by Stage C. "
        f"Radius selection (domain re-derivation, multi-scale sweep, QC, Pareto, "
        f"distance-to-ideal) is never re-run inside Phase D.")
    add(f"- iterations completed: {inv['n_iterations_completed']}; equal to the "
        f"original `r_fp`: {inv['n_equal_to_original']} "
        f"(invariant holds: {inv['invariant_holds']})")
    add("")

    add("## ROBUSTNESS PROFILE (the primary result)")
    add("")
    j = profile["jaccard_residues"]
    d = profile["dice_residues"]
    add("| metric | n | median | IQR | mean | SD |")
    add("|---|---|---|---|---|---|")
    for name, block in (("covered-residue Jaccard", j),
                        ("covered-residue Dice", d)):
        iqr = (None if block["q3"] is None or block["q1"] is None
               else block["q3"] - block["q1"])
        add(f"| {name} | {block['n']} | {_fmt(block['median'])} | {_fmt(iqr)} | "
            f"{_fmt(block['mean'])} | {_fmt(block['sd'])} |")
    add("")
    pres = profile["preservation"]
    add(f"- `P(J >= {pres['threshold_primary']:g})` = "
        f"**{_fmt(pres['P_jaccard_ge_primary'])}** "
        f"({pres['n_iterations_ge_primary']}/{pres['denominator']})")
    add(f"- `P(J >= {pres['threshold_secondary']:g})` = "
        f"**{_fmt(pres['P_jaccard_ge_secondary'])}** "
        f"({pres['n_iterations_ge_secondary']}/{pres['denominator']})")
    add("- both thresholds are pre-registered in config and exist **solely** to "
        "compute the preservation proportion. **No ROBUST / NOT_ROBUST verdict is "
        "emitted.**")
    add("")
    add(f"- component preservation rate: "
        f"{_fmt(profile['component_preservation_rate']['rate'])} "
        f"({profile['component_preservation_rate']['n_true']}/"
        f"{profile['component_preservation_rate']['denominator']})")
    add(f"- geometric continuity rate: "
        f"{_fmt(profile['geometric_continuity_rate']['rate'])}")
    add(f"- centroid shift: median "
        f"{_fmt(profile['centroid_shift_A']['median'])} A; "
        f"95th-percentile Hausdorff: median "
        f"{_fmt(profile['hausdorff95_A']['median'])} A")
    add("")
    add("### Classification performance, both families")
    add("")
    add("| family | n | median MCC | mean MCC | SD |")
    add("|---|---|---|---|---|")
    for name, block in (("against FP_iter (footprint)",
                         profile["mcc_footprint_family"]),
                        ("within r_hot of S_iter (`_hs`)",
                         profile["mcc_hotspot_coverage_family"])):
        add(f"| {name} | {block['n']} | {_fmt(block['median'])} | "
            f"{_fmt(block['mean'])} | {_fmt(block['sd'])} |")
    add("")
    add("Both families are reported because F10 permits `r_fp < r_hot`; in that case "
        "a footprint-only confusion matrix deflates MCC by construction. The `_hs` "
        "family is a pure geometric coverage query over the unchanged cohort `L` — "
        "it is not hotspot re-testing.")
    add("")

    add("### Center-loss sensitivity and recurrence")
    add("")
    add("| rank | center | n removing | mean Jaccard drop | SD | recurrence |")
    add("|---|---|---|---|---|---|")
    for row in sensitivity[:15]:
        add(f"| {row['rank']} | {row['center_residue_index']} | "
            f"{row['n_iterations_removing']} | "
            f"{_fmt(row['mean_jaccard_drop'])} | {_fmt(row['sd_jaccard_drop'])} | "
            f"{_fmt(row['center_recurrence_rate'])} "
            f"({row['n_recurrent_iterations']}/{row['n_iterations_removing']}) |")
    add("")
    add("**Influence** is the mean drop in covered-residue Jaccard when a center is "
        "removed — how much the whole footprint degrades without it. **Recurrence** "
        "is the fraction of those same iterations in which the reconstructed "
        "footprint still covers that center's position anyway — whether the "
        "surviving centers reach it. The two come apart: a footprint can degrade "
        "badly overall while still covering a removed center, or barely degrade "
        "while losing it.")
    add("")
    add("Failed iterations stay in the recurrence denominator and count as "
        "non-recurrence. A center that no evaluated iteration removed has no "
        "recurrence evidence and reports NA, never 0.")
    add("")

    add("## WARNINGS")
    add("")
    if not warnings:
        add("- none")
    for w in warnings:
        add(f"- **{w.severity.value}** `{w.warning_code}` — {w.message}")
    add("")
    add("## Interpretation")
    add("")
    add(_interpretation(profile, params))
    add("")
    return write_text(path, "\n".join(lines))


def _interpretation(profile: dict, params) -> str:
    p = profile["preservation"]["P_jaccard_ge_primary"]
    median = profile["jaccard_residues"]["median"]
    if p is None or not np.isfinite(p):
        return ("The preservation proportion could not be computed. See the "
                "iteration table for the reason.")
    body = (
        f"Across {profile['iterations']['planned']} evaluated iterations the median "
        f"covered-residue Jaccard against `FP_original` was {median:.3f}, and the "
        f"footprint met the pre-registered preservation threshold of "
        f"{params.preservation_primary:g} in {p:.1%} of them "
        f"(failed iterations included, contributing J = 0)."
    )
    if p < params.preservation_primary:
        body += (
            "\n\nThis is **limited geometric footprint robustness**: the spatial "
            "footprint changes substantially when significant centers are removed, "
            "so it rests appreciably on the specific set of centers present. That is "
            "a finding, reported as measured. Nothing was re-tuned, no threshold was "
            "moved and no iteration was dropped in response to it. Every downstream "
            "statement about this footprint must carry this profile."
        )
    else:
        body += (
            "\n\nThe footprint is largely preserved under center removal at the "
            "pre-registered threshold. This is a statement about geometric stability "
            "under perturbation of the center set, and about nothing else: it does "
            "not corroborate hotspot significance, which was decided upstream and "
            "was not re-tested here."
        )
    return body


def write_not_evaluable_report(path: Path, reason: str, phase_c, baseline: dict,
                               run_id: str, gene: str) -> Path:
    return write_text(path, "\n".join([
        f"# Stage D — ROBUSTNESS NOT EVALUABLE ({gene})",
        "",
        f"- run_id: `{run_id}`",
        f"- footprint code_version: `{phase_c.code_version}`",
        f"- r_fp: {phase_c.r_fp}",
        "",
        "## Reason",
        "",
        reason,
        "",
        "## How to read this",
        "",
        "This is an honest statement of non-evaluability, not a failure and not a "
        "result. No surrogate perturbation was invented, no weaker design was "
        "substituted, and no robustness number was fabricated to fill the gap.",
        "",
        f"- baseline reproduction: "
        f"{'PASSED' if baseline.get('byte_identical') else 'NOT PERFORMED'}",
        "",
        "The footprint itself remains valid and is reported downstream with an "
        "explicit not-evaluable robustness statement attached.",
        "",
    ]))
