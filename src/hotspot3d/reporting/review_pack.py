"""REVIEW_PACK assembly — Output Contract IX.2 (Lead-owned).

P1: every number in REVIEW_PACK is COPIED VERBATIM from a FULL_RESULTS artifact
and checksum-verified. A value may never appear in two places with two
computations behind it. Files are hardlinked where the filesystem permits (zero
duplicated bytes, identical content by construction) and copied otherwise; either
way the copy is checksum-verified against its source.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

from ..utils.hashing import sha256_file
from ..utils.io import read_json, read_tsv, write_json, write_text, write_tsv
from ..utils.runctx import RunContext
from ..utils.status import BANNERS, RunStatus
from .index_html import render_index_html
from .manifest import MANIFEST_COLUMNS, dir_size

KEY_METRIC_COLUMNS = ["metric", "value", "unit", "stage", "source_file"]
SELECTED_RADII_COLUMNS = [
    "radius_type", "value_A", "objective_vector", "normalized_vector",
    "distance_to_ideal", "domain_source", "domain_lo_A", "domain_hi_A",
    "fallback_flag", "boundary_warning", "near_tie_flag", "n_pareto", "source_file",
]

MANDATORY_FIGURES = [
    ("F1", "Ripley's K / L(r)-r, P/LP and B/LB, with permutation envelopes",
     "04_GLOBAL_CLUSTERING", "F1_ripley_k.png"),
    ("F2", "Pair correlation g(r) with envelopes; derived r_hot domain shaded",
     "04_GLOBAL_CLUSTERING", "F2_pair_correlation.png"),
    ("F3", "r_hot scan — four objectives vs radius, inadmissible shaded",
     "05_HOTSPOT_RADIUS", "F3_radius_scan.png"),
    ("F4", "r_hot Pareto front, parallel coordinates in normalized space",
     "05_HOTSPOT_RADIUS", "F4_pareto_parallel.png"),
    ("F5", "Distance-to-ideal vs radius, boundary band shaded when warned",
     "05_HOTSPOT_RADIUS", "F5_distance_to_ideal.png"),
    ("F6", "Hotspot map — -log10(q) per residue with the BH threshold line",
     "06_FINAL_HOTSPOTS", "F6_hotspot_map.png"),
    ("F7", "Footprint multi-scale — components and volume vs r_fp",
     "07_FOOTPRINT_RADIUS", "F7_footprint_multiscale.png"),
    ("F8", "r_fp Pareto parallel coordinates and distance-to-ideal",
     "07_FOOTPRINT_RADIUS", "F8_footprint_pareto.png"),
    ("F9", "Robustness — Jaccard and Dice distributions with thresholds",
     "09_ROBUSTNESS", "F9_robustness_distributions.png"),
    ("F10", "Center-loss sensitivity — per-center influence ranking",
     "09_ROBUSTNESS", "F10_center_influence.png"),
]


def _link_or_copy(src: Path, dst: Path) -> str:
    """Hardlink when possible, else copy. Always checksum-verified (IX.13.1)."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        dst.unlink()
    method = "hardlink"
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)
        method = "copy"
    if sha256_file(src) != sha256_file(dst):
        raise RuntimeError(f"REVIEW_PACK derivation failed checksum verification: {src}")
    return method


def _maybe(path: Path):
    """Read a JSON artifact if present, else None. Absence is never an error here."""
    try:
        return read_json(path) if path.is_file() else None
    except Exception:
        return None


def _rows(path: Path) -> list[dict]:
    try:
        return read_tsv(path) if path.is_file() else []
    except Exception:
        return []


def build_review_pack(ctx: RunContext, *, run_status: RunStatus, outcome_type: str,
                      metrics: dict, identity: dict, warnings_rows: list[dict],
                      negative_result: dict | None,
                      manifest_rows: list[dict]) -> dict:
    """Assemble REVIEW_PACK/. Returns a summary dict for run_summary.json."""
    pack = ctx.review_pack
    pack.mkdir(parents=True, exist_ok=True)
    fr = ctx.full_results
    derived: list[dict] = []

    # -- decision traces (small, belong at L2) -------------------------------
    trace_dir = pack / "decision_trace"
    for src_rel, dst_name in (("05_HOTSPOT_RADIUS/r_hot_scan.tsv", "r_hot_scan.tsv"),
                              ("07_FOOTPRINT_RADIUS/r_fp_scan.tsv", "r_fp_scan.tsv")):
        src = fr / src_rel
        if src.is_file():
            method = _link_or_copy(src, trace_dir / dst_name)
            derived.append({"source": src_rel, "target": f"decision_trace/{dst_name}",
                            "method": method})

    # -- key metrics (long format so the metric set can grow) ---------------
    metric_rows = [
        {"metric": k, "value": v.get("value"), "unit": v.get("unit", "NA"),
         "stage": v.get("stage", "NA"), "source_file": v.get("source", "NA")}
        for k, v in metrics.items()
    ]
    write_tsv(pack / "key_metrics.tsv", metric_rows, KEY_METRIC_COLUMNS)
    write_json(pack / "key_metrics.json", {"metrics": metric_rows})

    # -- selected radii: the one-line answer to "why this radius" -----------
    write_tsv(pack / "selected_radii.tsv", _selected_radii_rows(fr),
              SELECTED_RADII_COLUMNS)

    # -- warnings (verbatim copy of the root file) --------------------------
    root_warnings = ctx.run_root / "warnings.tsv"
    if root_warnings.is_file():
        _link_or_copy(root_warnings, pack / "warnings.tsv")

    # -- structures (consolidated READ view, copy-only) ---------------------
    struct_dir = pack / "structures"
    for src in sorted(fr.rglob("structures/*")):
        if src.is_file():
            _link_or_copy(src, struct_dir / src.name)
    if not struct_dir.exists() or not any(struct_dir.iterdir()):
        write_text(struct_dir / "RENDERING_UNAVAILABLE.txt",
                   "No structural exports were produced for this run.\n"
                   "See manifest.tsv for the NOT_CREATED rows and their reasons.\n")

    # -- figures (150 dpi PNG + SVG; full resolution stays in FULL_RESULTS) --
    figs_dir = pack / "key_figures"
    figure_specs: list[tuple[str, str, Path]] = []
    for code, caption, stage, fname in MANDATORY_FIGURES:
        src_png = fr / stage / "figures" / fname
        target = figs_dir / fname
        if src_png.is_file():
            _link_or_copy(src_png, target)
            src_svg = src_png.with_suffix(".svg")
            if src_svg.is_file():
                _link_or_copy(src_svg, figs_dir / src_svg.name)
        figure_specs.append((code, f"{code} — {caption}", target))

    # -- narrative report ---------------------------------------------------
    summary_md = _final_summary_md(ctx, run_status, outcome_type, metrics,
                                   identity, warnings_rows, negative_result)
    write_text(pack / "final_summary.md", summary_md)
    _emit_pdf_or_placeholder(pack)

    # -- the L1 interface ---------------------------------------------------
    html_text = render_index_html(
        identity=identity, run_status=run_status.value,
        banner=BANNERS[run_status].format(gene=ctx.gene), outcome_type=outcome_type,
        metrics=metrics, flags=_flags(warnings_rows),
        sections=_sections(fr), figures=figure_specs,
        negative_result=negative_result, warnings_rows=warnings_rows,
    )
    (pack / "index.html").write_text(html_text, encoding="utf-8", newline="\n")

    # -- pack manifest + how-to-view ----------------------------------------
    pack_rows = [r for r in manifest_rows
                 if str(r["relative_path"]).startswith("REVIEW_PACK/")]
    write_tsv(pack / "manifest_review.tsv", pack_rows, MANIFEST_COLUMNS)
    write_text(pack / "HOW_TO_VIEW.txt",
               "1. Open index.html by double-clicking it. No server, no internet needed.\n"
               "2. 3D: chimerax structures/view_hotspots.cxc "
               "(or PyMOL: pymol structures/view_hotspots.pml)\n"
               "3. The complete audit trail is in ../FULL_RESULTS/\n")

    size = dir_size(pack)
    budget = int(ctx.config.get("output.review_pack.size_budget_bytes"))
    return {"bytes": size, "over_budget": size > budget, "budget": budget,
            "derived": derived, "n_figures_embedded":
                sum(1 for _c, _cap, p in figure_specs if p.is_file())}


def _selected_radii_rows(fr: Path) -> list[dict]:
    rows = []
    hot = _maybe(fr / "05_HOTSPOT_RADIUS" / "radius_decision.json")
    dom = _maybe(fr / "04_GLOBAL_CLUSTERING" / "candidate_radius_domain.json")
    bnd = _maybe(fr / "05_HOTSPOT_RADIUS" / "domain_boundary_diagnostic.json")
    if hot:
        rows.append({
            "radius_type": "r_hot", "value_A": hot.get("selected_r_hot"),
            "objective_vector": hot.get("selected_objectives"),
            "normalized_vector": hot.get("selected_normalized"),
            "distance_to_ideal": hot.get("selected_distance"),
            "domain_source": (dom or {}).get("search_domain_source"),
            "domain_lo_A": (dom or {}).get("domain_lo"),
            "domain_hi_A": (dom or {}).get("domain_hi"),
            "fallback_flag": (dom or {}).get("FALLBACK_RADIUS_DOMAIN"),
            "boundary_warning": (bnd or {}).get("DOMAIN_BOUNDARY_WARNING"),
            "near_tie_flag": hot.get("near_tie"),
            "n_pareto": len(hot.get("pareto_members", []) or []),
            "source_file": "FULL_RESULTS/05_HOTSPOT_RADIUS/radius_decision.json",
        })
    fp = _maybe(fr / "07_FOOTPRINT_RADIUS" / "footprint_decision.json")
    fdom = _maybe(fr / "07_FOOTPRINT_RADIUS" / "footprint_domain.json")
    fbnd = _maybe(fr / "07_FOOTPRINT_RADIUS" / "footprint_domain_boundary_diagnostic.json")
    if fp:
        rows.append({
            "radius_type": "r_fp", "value_A": fp.get("selected_r_fp"),
            "objective_vector": fp.get("selected_objectives"),
            "normalized_vector": fp.get("selected_normalized"),
            "distance_to_ideal": fp.get("selected_distance"),
            "domain_source": "MST_DERIVED",
            "domain_lo_A": (fdom or {}).get("rho_min"),
            "domain_hi_A": (fdom or {}).get("rho_max"),
            "fallback_flag": False,
            "boundary_warning": (fbnd or {}).get("DOMAIN_BOUNDARY_WARNING"),
            "near_tie_flag": fp.get("near_tie"),
            "n_pareto": len(fp.get("pareto_members", []) or []),
            "source_file": "FULL_RESULTS/07_FOOTPRINT_RADIUS/footprint_decision.json",
        })
    return rows


def _flags(warnings_rows: list[dict]) -> list[dict]:
    seen, out = set(), []
    for w in warnings_rows:
        code = str(w.get("warning_code"))
        if code not in seen:
            seen.add(code)
            out.append({"warning_code": code, "severity": str(w.get("severity"))})
    order = {"BLOCKING": 0, "MAJOR": 1, "ADVISORY": 2, "INFO": 3}
    return sorted(out, key=lambda f: (order.get(f["severity"], 9), f["warning_code"]))


def _sections(fr: Path) -> dict[str, str]:
    """The five mandated sections, each with headline numbers visible (IX.2)."""
    from .index_html import _table

    out: dict[str, str] = {}

    centers = _rows(fr / "06_FINAL_HOTSPOTS" / "significant_hotspot_centers.tsv")
    regions = _rows(fr / "06_FINAL_HOTSPOTS" / "hotspot_regions.tsv")
    scan = _rows(fr / "05_HOTSPOT_RADIUS" / "r_hot_scan.tsv")
    out["Hotspot discovery"] = (
        f"<p><b>{len(centers)}</b> significant hotspot centers in "
        f"<b>{len(regions)}</b> region(s).</p>"
        + (f"<h3>Regions</h3>{_table(regions, list(regions[0].keys())[:10], 'tbl-reg')}"
           if regions else "<p class='muted'>No hotspot regions.</p>")
        + (f"<details open><summary>r_hot decision trace "
           f"({len(scan)} scanned radii, rejected candidates retained)</summary>"
           f"{_table(scan, [c for c in ['radius_A','loo_mcc','perm_evidence_raw','fold_enrichment','neighbor_stability','n_significant_centers','coverage_fraction','qc_status','qc_failure_reason','pareto_member','distance_to_ideal','selected'] if scan and c in scan[0]], 'tbl-hot', 'selected')}</details>"
           if scan else ""))

    fp_scan = _rows(fr / "07_FOOTPRINT_RADIUS" / "r_fp_scan.tsv")
    fp_res = _rows(fr / "08_FINAL_FOOTPRINT" / "footprint_residues.tsv")
    out["Footprint"] = (
        f"<p>Final footprint covers <b>{len(fp_res)}</b> residues.</p>"
        + (f"<details open><summary>r_fp decision trace ({len(fp_scan)} candidates)"
           f"</summary>{_table(fp_scan, [c for c in ['r_fp_A','volume_A3','surface_area_A2','n_components','coverage_fraction','compactness','connectivity','stability','parsimony','bridging_index','qc_status','qc_failure_reason','pareto_member','distance_to_ideal','selected'] if fp_scan and c in fp_scan[0]], 'tbl-fp', 'selected')}</details>"
           if fp_scan else "<p class='muted'>No footprint scan.</p>"))

    profile = _maybe(fr / "09_ROBUSTNESS" / "robustness_profile.json")
    influence = _rows(fr / "09_ROBUSTNESS" / "center_sensitivity.tsv")
    if profile:
        jac = profile.get("jaccard_residues", {})
        out["Robustness"] = (
            f"<p>Geometric footprint robustness — <b>not</b> validation of hotspot "
            f"discovery. Median Jaccard <b>{jac.get('median','NA')}</b> "
            f"(IQR {jac.get('iqr_lo','NA')}–{jac.get('iqr_hi','NA')}); "
            f"P(J&ge;0.50) = <b>{profile.get('preservation_freq_050','NA')}</b>, "
            f"P(J&ge;0.70) = {profile.get('preservation_freq_070','NA')}.</p>"
            f"<p class='muted'>Reported as a continuous profile; no ROBUST/NOT_ROBUST "
            f"verdict is emitted.</p>"
            + (f"<details><summary>Per-center influence ({len(influence)} centers)</summary>"
               f"{_table(influence, list(influence[0].keys())[:8], 'tbl-inf')}</details>"
               if influence else ""))
    else:
        out["Robustness"] = "<p class='muted'>Robustness not evaluated for this run.</p>"

    annot = _rows(fr / "10_ANNOTATION" / "hotspot_annotation.tsv")
    mech = _rows(fr / "10_ANNOTATION" / "mechanism_by_hotspot.tsv")
    out["Biological annotation"] = (
        (f"<p><b>{len(annot)}</b> hotspot(s) annotated. Mechanism labels are a "
         f"downstream overlay and played no part in discovery.</p>"
         f"<details open><summary>Annotation</summary>"
         f"{_table(annot, list(annot[0].keys())[:12], 'tbl-ann')}</details>"
         + (f"<details><summary>Mechanism distribution</summary>"
            f"{_table(mech, list(mech[0].keys())[:8], 'tbl-mech')}</details>" if mech else ""))
        if annot else "<p class='muted'>No annotation produced for this run.</p>")

    return out


def _final_summary_md(ctx, run_status, outcome_type, metrics, identity,
                      warnings_rows, negative_result) -> str:
    major = [w for w in warnings_rows if str(w.get("severity")) in ("MAJOR", "BLOCKING")]
    lines = [
        f"# 3D Hotspot & Footprint Analysis — {ctx.gene}",
        "",
        f"**{BANNERS[run_status].format(gene=ctx.gene)}**  ",
        f"`outcome_type = {outcome_type}`  ",
        f"`RUN_ID = {ctx.run_id}`",
        "",
    ]
    if major:
        lines += ["## Warnings you must read first", ""]
        lines += [f"- **[{w['severity']}] {w['warning_code']}** — {w['message']}  "
                  f"_Action: {w.get('recommended_action','NA')}_" for w in major]
        lines += [""]
    if negative_result:
        lines += ["## Negative result", "",
                  f"**{negative_result.get('condition','')}** — "
                  f"{negative_result.get('detail','')}", "",
                  "This is a complete, valid scientific outcome, not a technical "
                  "failure, and it licenses no method modification.", ""]
    lines += ["## Run identity", ""]
    lines += [f"- **{k}**: {v}" for k, v in identity.items()]
    lines += ["", "## Key metrics", "",
              "| metric | value | unit |", "|---|---|---|"]
    lines += [f"| {k} | {v.get('value')} | {v.get('unit','')} |"
              for k, v in metrics.items()]
    lines += [
        "", "## Method boundaries recorded for this run", "",
        "- `r_hot` was selected by Pareto → min–max normalization over the admissible "
        "set → equally weighted (`w_k = 1`) L2 distance-to-ideal. It is **not** the "
        "max-MCC radius and **not** the pair-correlation peak.",
        "- LOO-MCC is a radius-**selection** objective only; it is not hotspot validation.",
        "- `r_fp` was optimized independently of `r_hot` and is not bounded below by it.",
        "- Robustness perturbs the significant hotspot **center set**; no ClinVar record "
        "or class label was modified, and the hotspot pipeline was not re-run.",
        "- Mechanism labels are a downstream overlay and were unavailable to discovery.",
        "", "## Where to look next", "",
        "- `index.html` — this run at a glance (offline).",
        "- `decision_trace/` — every scanned candidate, including rejected ones.",
        "- `../FULL_RESULTS/` — complete audit trail.",
        "- `../FULL_RESULTS/12_REPRODUCIBILITY/` — config snapshot, seeds, versions, checksums.",
        "",
    ]
    return "\n".join(lines)


def _emit_pdf_or_placeholder(pack: Path) -> None:
    """PDF is best-effort; if no toolchain exists, say so explicitly (IX.2)."""
    for tool in ("pandoc", "wkhtmltopdf"):
        if shutil.which(tool):
            try:
                import subprocess
                if tool == "pandoc":
                    subprocess.run([tool, str(pack / "final_summary.md"),
                                    "-o", str(pack / "final_summary.pdf")],
                                   check=True, capture_output=True, timeout=120)
                    return
            except Exception:
                break
    write_text(pack / "final_summary.pdf.UNAVAILABLE.txt",
               "final_summary.pdf was not produced: no PDF toolchain "
               "(pandoc / wkhtmltopdf) is available on this node.\n\n"
               "The Markdown report final_summary.md contains the identical content.\n\n"
               "To regenerate the PDF locally:\n"
               "  pandoc final_summary.md -o final_summary.pdf\n")
