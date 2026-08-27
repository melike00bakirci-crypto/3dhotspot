"""Terminal completion summary — Output Contract IX.10 (Lead-owned).

The identical text is saved to 00_RUN_SUMMARY/terminal_summary.txt.
For negative or partial runs, metric lines that do not exist print
"— not computed (see NOT_RUN.txt)" rather than being omitted.
"""
from __future__ import annotations

from ..utils.status import BANNERS, EXIT_CODES, RunStatus

_NOT_COMPUTED = "— not computed (see NOT_RUN.txt)"


def _v(value, suffix: str = "") -> str:
    if value is None or value == "NA":
        return _NOT_COMPUTED
    # Round floats for display only. The unrounded value is always available in the
    # canonical artifact this line was copied from; a terminal banner reading
    # "564.8922092825686 A^2" is noise, not precision.
    if isinstance(value, float):
        return f"{value:.6g}{suffix}"
    return f"{value}{suffix}"


def render_terminal_summary(*, gene: str, run_id: str, run_status: RunStatus,
                            outcome_type: str, identity: dict, metrics: dict,
                            warnings_rows: list[dict], paths: dict) -> str:
    banner = BANNERS[run_status].format(gene=gene)
    rule = "=" * max(60, len(banner))

    def m(key, suffix=""):
        entry = metrics.get(key)
        return _v(entry.get("value") if isinstance(entry, dict) else entry, suffix)

    lines = [
        rule, banner, rule, "",
        f"  RUN_ID              {run_id}",
        f"  outcome_type        {outcome_type}",
        f"  gene / UniProt      {identity.get('Gene','NA')} / {identity.get('UniProt','NA')}",
        f"  ClinVar release     {identity.get('ClinVar release','NA')}",
        f"  Structure           {identity.get('Structure','NA')}",
        "",
        "  COHORT",
        f"    missense retrieved  {m('n_missense_retrieved')}",
        f"    P/LP residues       {m('n_residues_plp')}",
        f"    B/LB residues       {m('n_residues_blb')}",
        f"    class conflicts     {m('n_residue_class_conflict')}",
        f"    positional universe {m('n_positional_universe')}",
        "",
        "  HOTSPOT DISCOVERY",
        f"    r_hot               {m('r_hot', ' A')}",
        f"    domain source       {m('r_hot_domain_source')}"
        f"   (fallback={m('fallback_radius_domain')})",
        f"    permutations B      {m('permutation_count')}",
        f"    resolution limited  {m('permutation_resolution_limited')}"
        f"   (B_rec={m('b_recommended')})",
        f"    FDR                 {m('fdr_method')} at q={m('fdr_q')}",
        f"    significant centers {m('n_significant_centers')}"
        f"   in {m('n_hotspot_regions')} region(s)",
        f"    centers w/o variant {m('n_centers_without_variant')}",
        "",
        "  FOOTPRINT",
        f"    r_fp                {m('r_fp', ' A')}",
        f"    volume / surface    {m('footprint_volume_A3', ' A^3')}"
        f" / {m('footprint_surface_A2', ' A^2')}",
        f"    residues / coverage {m('footprint_n_residues')} / {m('footprint_coverage')}",
        f"    components          {m('footprint_n_components')}",
        "",
        "  ROBUSTNESS (geometric footprint robustness only)",
        f"    subsets evaluated   {m('n_perturbations_evaluated')}"
        f" of {m('n_perturbations_possible')} possible",
        f"    median Jaccard      {m('jaccard_median')}"
        f"   IQR [{m('jaccard_iqr_lo')}, {m('jaccard_iqr_hi')}]",
        f"    P(J>=0.50)          {m('preservation_freq_050')}"
        f"   P(J>=0.70) {m('preservation_freq_070')}",
        "    (no ROBUST/NOT_ROBUST verdict is emitted by design)",
        "",
    ]

    flagged = [w for w in warnings_rows
               if str(w.get("severity")) in ("MAJOR", "BLOCKING")]
    lines.append("  WARNINGS (MAJOR / BLOCKING)")
    if flagged:
        seen = set()
        for w in flagged:
            code = str(w.get("warning_code"))
            if code in seen:
                continue
            seen.add(code)
            lines.append(f"    [{w.get('severity')}] {code}")
            lines.append(f"        {w.get('message')}")
    else:
        lines.append("    none")
    lines += [
        "",
        "  OUTPUTS",
        f"    review pack   {paths.get('index_html', 'NA')}",
        f"    full results  {paths.get('full_results', 'NA')}",
        f"    archives      {paths.get('review_pack_zip', 'NA')}",
        f"                  {paths.get('full_results_tar', 'NA')}",
        "", rule,
        f"exit code {EXIT_CODES[run_status]}", "",
    ]
    return "\n".join(lines)
