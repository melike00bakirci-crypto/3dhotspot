"""Mandatory Stage B figures F1-F6 (Output Contract IX; methodology paragraph 57).

PNG at 300 dpi for FULL_RESULTS plus SVG vector companions. The Lead derives the
150 dpi REVIEW_PACK copies; this module never writes there.

Figures are written **deterministically**: a fixed hash salt and suppressed
timestamp metadata mean two runs produce byte-identical image files, so figures can
participate in the re-run determinism assertion instead of being excluded from it.
The settings live in ``hotspot3d.utils.figures`` because rcParams are process-wide
and any module may reset them — setting them here alone was not enough.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

from ..utils.figures import SVG_METADATA, apply_deterministic_rcparams

matplotlib.use("Agg")
apply_deterministic_rcparams()

import matplotlib.pyplot as plt                                        # noqa: E402
import numpy as np                                                     # noqa: E402

DPI = 300
_INADMISSIBLE = dict(color="0.85", alpha=0.55, zorder=0)


def _save(fig, stage_dir: Path, name: str) -> list[Path]:
    figures = Path(stage_dir) / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    png = figures / f"{name}.png"
    svg = figures / f"{name}.svg"
    fig.savefig(png, dpi=DPI, format="png", metadata={"Software": "hotspot3d"})
    fig.savefig(svg, format="svg", metadata=SVG_METADATA)
    plt.close(fig)
    return [png, svg]


def _shade_inadmissible(ax, radii, admissible, step):
    first = True
    for r, ok in zip(radii, admissible):
        if ok:
            continue
        ax.axvspan(r - step / 2, r + step / 2,
                   label="inadmissible (QC)" if first else None, **_INADMISSIBLE)
        first = False


# --- F1 ----------------------------------------------------------------------

def figure_f1_ripley(stage_dir: Path, stats: dict) -> list[Path]:
    """K(r) and L(r) - r with 2.5/97.5 positional-null envelopes, both cohorts."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for row, cohort in enumerate(("PLP", "BLB")):
        s = stats[cohort]
        r = s.grid
        ax = axes[row, 0]
        ax.fill_between(r, s.k_env_lo, s.k_env_hi, color="0.80",
                        label="null envelope 2.5-97.5%")
        ax.plot(r, s.k_null_mean, color="0.45", lw=1.0, ls="--", label="null mean")
        ax.plot(r, s.k_obs, color="firebrick" if cohort == "PLP" else "steelblue",
                lw=1.8, label=f"observed {cohort}")
        ax.set_xlabel("r (A)")
        ax.set_ylabel("K(r)")
        ax.set_title(f"{cohort}: Ripley K  (n = {s.n_points}, "
                     f"global max|Z| = {s.t_obs:.3g}, p = {s.p_global:.4g})")
        ax.legend(fontsize=7, loc="upper left")

        ax = axes[row, 1]
        ax.fill_between(r, s.l_env_lo - r, s.l_env_hi - r, color="0.80")
        ax.axhline(0.0, color="0.35", lw=1.0, ls=":")
        ax.plot(r, s.l_obs - r, color="firebrick" if cohort == "PLP" else "steelblue",
                lw=1.8)
        ax.set_xlabel("r (A)")
        ax.set_ylabel("L(r) - r")
        ax.set_title(f"{cohort}: L(r) - r  (>0 = clustering)")
    fig.suptitle("F1 — global spatial statistics under the structure-aware positional null\n"
                 "Global assessment only: Ripley's K never selects r_hot (F5/F6).",
                 fontsize=10)
    return _save(fig, stage_dir, "F1_ripleys_k")


# --- F2 ----------------------------------------------------------------------

def figure_f2_pcf(stage_dir: Path, stats: dict, peaks: dict, domain) -> list[Path]:
    """g(r) with envelopes, peaks marked, the derived domain shaded."""
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
    for ax, cohort in zip(axes, ("PLP", "BLB")):
        s = stats[cohort]
        r = s.grid
        ax.fill_between(r, s.g_env_lo, s.g_env_hi, color="0.80",
                        label="null envelope 2.5-97.5%")
        ax.plot(r, s.g_null_p95, color="0.45", lw=1.0, ls="--", label="null 95th pct")
        ax.axhline(1.0, color="0.30", lw=1.0, ls=":", label="g = 1 (CSR)")
        ax.plot(r, s.g_obs, color="firebrick" if cohort == "PLP" else "steelblue",
                lw=1.8, label=f"observed {cohort}")
        if cohort == "PLP":
            for p in peaks.get("peaks", []):
                ax.plot(p["radius_A"], p["g"], marker="v", color="darkgreen", ms=6,
                        ls="none")
            principal = peaks.get("principal_peak")
            if principal:
                ax.annotate(f"principal PCF peak {principal['radius_A']:g} A\n"
                            f"(NEVER used as r_hot)",
                            xy=(principal["radius_A"], principal["g"]),
                            xytext=(6, 14), textcoords="offset points", fontsize=7,
                            color="darkgreen")
            if len(domain.grid):
                ax.axvspan(float(domain.grid[0]), float(domain.grid[-1]),
                           color="gold", alpha=0.25, zorder=0,
                           label=f"candidate r_hot domain ({domain.source})")
            if domain.fallback:
                ax.text(0.02, 0.92,
                        f"FALLBACK_RADIUS_DOMAIN = TRUE  (trigger {domain.trigger_id})",
                        transform=ax.transAxes, fontsize=8, color="darkred")
        ax.set_ylabel("g(r)")
        ax.set_title(f"{cohort}: pair correlation")
        ax.legend(fontsize=7, loc="upper right")
    axes[-1].set_xlabel("r (A)")
    fig.suptitle("F2 — pair correlation and the derived candidate radius domain\n"
                 "The PCF peak is never r_hot and never defines the interval alone (F6/F14).",
                 fontsize=10)
    return _save(fig, stage_dir, "F2_pair_correlation")


# --- F3 ----------------------------------------------------------------------

def figure_f3_scan(stage_dir: Path, rows: list[dict], selected: float | None,
                   step: float) -> list[Path]:
    """The four Pareto objectives against radius; inadmissible radii shaded."""
    radii = [r["radius_A"] for r in rows]
    admissible = [r["qc_status"] == "PASS" for r in rows]
    panels = [
        ("loo_mcc", "LOO-MCC\n(selection objective only)"),
        ("perm_evidence_raw", "permutation evidence Zg"),
        ("fold_enrichment", "fold enrichment"),
        ("neighbor_stability", "neighbouring-radius stability"),
    ]
    fig, axes = plt.subplots(len(panels), 1, figsize=(10, 11), sharex=True,
                             constrained_layout=True)
    for ax, (key, label) in zip(axes, panels):
        _shade_inadmissible(ax, radii, admissible, step)
        values = [r[key] if r[key] is not None else np.nan for r in rows]
        ax.plot(radii, values, marker="o", ms=4, lw=1.4, color="navy")
        if selected is not None:
            ax.axvline(selected, color="firebrick", lw=1.6,
                       label=f"selected r_hot = {selected:g} A")
        ax.set_ylabel(label, fontsize=8)
        ax.legend(fontsize=7, loc="best")
    axes[-1].set_xlabel("candidate radius (A)")
    fig.suptitle("F3 — radius scan: four objectives, all maximizing, none decisive alone\n"
                 "Max-MCC has no privileged status (F7/F8); QC rules are constraints, "
                 "not objectives.", fontsize=10)
    return _save(fig, stage_dir, "F3_radius_scan")


# --- F4 ----------------------------------------------------------------------

def figure_f4_pareto(stage_dir: Path, normalized: dict, pareto: list[float],
                     names: list[str], selected: float | None) -> list[Path]:
    """Pareto front as PARALLEL COORDINATES in normalized objective space."""
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    xs = np.arange(len(names))
    if not pareto:
        ax.text(0.5, 0.5, "no admissible candidate — Pareto front is empty",
                ha="center", va="center", transform=ax.transAxes)
    for key in sorted(pareto):
        values = [normalized[key][n] for n in names]
        is_sel = selected is not None and key == selected
        ax.plot(xs, values, marker="o", ms=5,
                lw=2.6 if is_sel else 1.1,
                color="firebrick" if is_sel else "0.55",
                zorder=3 if is_sel else 2,
                label=f"r = {key:g} A" + (" (SELECTED)" if is_sel else ""))
    ax.plot(xs, np.ones(len(names)), ls="--", lw=1.4, color="darkgreen",
            label="utopia z* = (1,1,1,1)")
    ax.set_xticks(xs)
    ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel("min-max normalized value (over the admissible set)")
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=7, ncol=2)
    fig.suptitle("F4 — Pareto front in normalized objective space (equal weights, w_k = 1)",
                 fontsize=10)
    return _save(fig, stage_dir, "F4_pareto_parallel_coordinates")


# --- F5 ----------------------------------------------------------------------

def figure_f5_distance(stage_dir: Path, distances: dict, pareto: list[float],
                       selected: float | None, boundary: dict, step: float) -> list[Path]:
    """L2 distance-to-ideal against radius, with the boundary band shaded if it fired."""
    fig, ax = plt.subplots(figsize=(10, 5.5), constrained_layout=True)
    if boundary.get("DOMAIN_BOUNDARY_WARNING") and boundary.get("band_members"):
        first = True
        for r in boundary["band_members"]:
            ax.axvspan(r - step / 2, r + step / 2, color="moccasin", alpha=0.7,
                       zorder=0, label="boundary band" if first else None)
            first = False
    keys = sorted(distances)
    ax.plot(keys, [distances[k] for k in keys], marker="o", ms=4, lw=1.3,
            color="navy", label="all admissible candidates")
    if pareto:
        ax.plot(sorted(pareto), [distances[k] for k in sorted(pareto)], ls="none",
                marker="s", ms=9, mfc="none", mec="darkgreen", label="Pareto members")
    if selected is not None:
        ax.axvline(selected, color="firebrick", lw=1.6,
                   label=f"selected r_hot = {selected:g} A")
    ax.set_xlabel("candidate radius (A)")
    ax.set_ylabel("L2 distance to utopia (lower is better)")
    ax.legend(fontsize=7)
    title = "F5 — distance-to-ideal over the admissible set"
    if boundary.get("DOMAIN_BOUNDARY_WARNING"):
        title += "\nDOMAIN_BOUNDARY_WARNING fired — reported, never self-corrected"
    fig.suptitle(title, fontsize=10)
    return _save(fig, stage_dir, "F5_distance_to_ideal")


# --- F6 ----------------------------------------------------------------------

def figure_f6_hotspot_map(stage_dir: Path, rows: list[dict], q: float,
                          boundary_p: float | None, radius: float) -> list[Path]:
    """-log10(q_bh) per residue along the sequence, with the BH decision line."""
    fig, ax = plt.subplots(figsize=(11, 5), constrained_layout=True)
    tested = [r for r in rows if r["in_test_family"]]
    if tested:
        idx = np.array([r["center_residue_index"] for r in tested])
        qv = np.array([r["q_bh"] for r in tested], dtype=np.float64)
        sig = np.array([r["significant"] for r in tested], dtype=bool)
        y = -np.log10(np.clip(qv, 1e-300, None))
        ax.vlines(idx[~sig], 0, y[~sig], color="0.72", lw=0.9)
        ax.vlines(idx[sig], 0, y[sig], color="firebrick", lw=1.6)
        ax.plot(idx[sig], y[sig], ls="none", marker="o", ms=3.5, color="firebrick",
                label=f"significant (BH q <= {q:g}): {int(sig.sum())}")
    ax.axhline(-np.log10(q), color="darkgreen", ls="--", lw=1.3,
               label=f"BH threshold q = {q:g}")
    excluded = [r for r in rows if not r["in_test_family"]]
    if excluded:
        ax.plot([r["center_residue_index"] for r in excluded],
                np.zeros(len(excluded)), ls="none", marker="x", ms=3, color="0.55",
                label=f"excluded from the family: {len(excluded)}")
    ax.set_xlabel("residue index")
    ax.set_ylabel("-log10(q_BH)")
    ax.legend(fontsize=7, loc="upper right")
    sub = f"r_hot = {radius:g} A"
    if boundary_p is not None:
        sub += f"; BH boundary p = {boundary_p:.6g}"
    fig.suptitle(f"F6 — hotspot map: every residue of U_struct is a candidate center\n{sub}",
                 fontsize=10)
    return _save(fig, stage_dir, "F6_hotspot_map")
