"""Mandatory Phase C figures F7 and F8 (PNG 300 dpi + SVG).

Figures are *reports* of the decision, never part of it: nothing here feeds back
into selection, and no footprint is ever chosen by eye (§8, "no visual selection
at any point").
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

from ..utils.figures import SVG_METADATA, apply_deterministic_rcparams

matplotlib.use("Agg")
apply_deterministic_rcparams()

import matplotlib.pyplot as plt        # noqa: E402
import numpy as np                     # noqa: E402

from .api import FootprintSolution     # noqa: E402
from .params import FootprintParams    # noqa: E402

DPI = 300


def _save(fig, out_dir: Path, stem: str) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("png", "svg"):
        path = out_dir / f"{stem}.{ext}"
        fig.savefig(path, dpi=DPI, bbox_inches="tight",
                    metadata=SVG_METADATA if ext == "svg" else None)
        paths.append(path)
    plt.close(fig)
    return paths


def figure_f7(out_dir: Path, solution: FootprintSolution,
              params: FootprintParams) -> list[Path]:
    """F7 — component count and volume vs r_fp, merges marked, coverage>0.50 shaded."""
    rows = solution.sweep.rows
    rho = np.array([r.rho for r in rows])
    vol = np.array([r.volume for r in rows])
    ncomp = np.array([r.n_components for r in rows])
    cov = np.array([r.coverage for r in rows])

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(rho, vol, color="#1f4e79", lw=2, marker="o", ms=3, label="volume $V$")
    ax.set_xlabel(r"footprint radius $\rho$ ($\AA$)")
    ax.set_ylabel(r"volume ($\AA^3$)", color="#1f4e79")
    ax.tick_params(axis="y", labelcolor="#1f4e79")

    ax2 = ax.twinx()
    ax2.step(rho, ncomp, where="post", color="#b03a2e", lw=1.8,
             label="components $n_{comp}$")
    ax2.set_ylabel(r"connected components", color="#b03a2e")
    ax2.tick_params(axis="y", labelcolor="#b03a2e")

    excessive = cov > params.qc_f1_max_coverage
    if excessive.any():
        ax.fill_between(rho, 0, 1, where=excessive, transform=ax.get_xaxis_transform(),
                        color="#d5d8dc", alpha=0.7, step="mid",
                        label=f"coverage > {params.qc_f1_max_coverage:g} (QC-F1 fail)")

    for event in solution.sweep.merge_events:
        colour = {"natural": "#117864", "artificial": "#d68910"}.get(
            event.connection_type, "#7d3c98")
        ax.axvline(event.rho, color=colour, ls="--", lw=1.0, alpha=0.85)

    if solution.r_fp is not None:
        ax.axvline(solution.r_fp, color="black", lw=2.2)
        ax.annotate(f"selected $r_{{fp}}$ = {solution.r_fp:g} $\\AA$",
                    xy=(solution.r_fp, ax.get_ylim()[1]),
                    xytext=(4, -12), textcoords="offset points", fontsize=9)

    handles, labels = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(handles + h2, labels + l2, loc="upper left", fontsize=8, frameon=False)
    ax.set_title("F7 — multi-scale footprint evolution\n"
                 "dashed: merge events (teal natural, amber artificial); "
                 "shaded: QC-F1 inadmissible", fontsize=10)
    return _save(fig, out_dir, "F7_footprint_multiscale_evolution")


def figure_f8(out_dir: Path, solution: FootprintSolution,
              params: FootprintParams) -> list[Path]:
    """F8 — Pareto parallel coordinates with distance-to-ideal, selection marked."""
    result = solution.selection
    names = list(params.objective_names)
    members = list(result.pareto_members)

    fig, (ax, axd) = plt.subplots(
        1, 2, figsize=(11, 5), gridspec_kw={"width_ratios": [2.0, 1.0]})

    if members:
        dmax = max(result.distances[k] for k in members) or 1.0
        cmap = plt.get_cmap("viridis")
        for key in members:
            values = [result.normalized[key][n] for n in names]
            selected = (solution.r_fp is not None and abs(key - solution.r_fp) < 1e-9)
            ax.plot(range(len(names)), values,
                    color="black" if selected else cmap(result.distances[key] / dmax),
                    lw=2.6 if selected else 1.1,
                    alpha=1.0 if selected else 0.75,
                    marker="o", ms=5 if selected else 3, zorder=5 if selected else 2,
                    label=f"selected {key:g} $\\AA$" if selected else None)
        ax.plot(range(len(names)), [1.0] * len(names), color="#c0392b", ls=":",
                lw=1.6, label="utopia (1,1,1,1)")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=15)
    ax.set_ylabel("min-max normalized objective (admissible set)")
    ax.set_ylim(-0.05, 1.08)
    if ax.get_legend_handles_labels()[1]:
        ax.legend(fontsize=8, frameon=False, loc="lower left")
    ax.set_title("F8a — Pareto front, parallel coordinates", fontsize=10)

    keys = sorted(result.distances)
    if keys:
        axd.plot(keys, [result.distances[k] for k in keys], color="#1f4e79",
                 lw=1.4, marker="o", ms=3, label="all admissible")
        axd.scatter(members, [result.distances[k] for k in members],
                    color="#e67e22", zorder=4, s=28, label="Pareto member")
        if solution.r_fp is not None:
            axd.scatter([solution.r_fp], [result.distances[solution.r_fp]],
                        color="black", zorder=6, s=70, marker="*",
                        label="selected $r_{fp}$")
    axd.set_xlabel(r"$\rho$ ($\AA$)")
    axd.set_ylabel("L2 distance to utopia")
    if axd.get_legend_handles_labels()[1]:
        axd.legend(fontsize=8, frameon=False)
    axd.set_title("F8b — distance-to-ideal", fontsize=10)

    fig.tight_layout()
    return _save(fig, out_dir, "F8_footprint_pareto_selection")
