"""Mandatory Phase D figures F9 and F10 (PNG 300 dpi + SVG).

Both show the distribution, not a headline. The preservation thresholds are drawn
as reference lines because that is all they are: pre-registered marks used to
compute the ¶51 proportion. Neither figure carries a verdict.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

from ..utils.figures import SVG_METADATA, apply_deterministic_rcparams

matplotlib.use("Agg")
apply_deterministic_rcparams()

import matplotlib.pyplot as plt        # noqa: E402
import numpy as np                     # noqa: E402

DPI = 300


def _save(fig, out_dir: Path, stem: str) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("png", "svg"):
        path = out_dir / f"{stem}.{ext}"
        # The SVG date was previously left at matplotlib's default, so F9 and F10
        # embedded a wall-clock <dc:date> and could never be byte-reproducible.
        # No Phase D assertion compared them, so it went unnoticed.
        fig.savefig(path, dpi=DPI, bbox_inches="tight",
                    metadata=SVG_METADATA if ext == "svg" else None)
        paths.append(path)
    plt.close(fig)
    return paths


def figure_f9(out_dir: Path, jaccards: list, dices: list,
              thresholds: tuple[float, float]) -> list[Path]:
    """F9 — Jaccard and Dice distributions with the preservation thresholds marked."""
    j = np.asarray([0.0 if v is None else float(v) for v in jaccards],
                   dtype=np.float64)
    d = np.asarray([0.0 if v is None else float(v) for v in dices], dtype=np.float64)

    fig, (ax_j, ax_d) = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    bins = np.linspace(0.0, 1.0, 21)
    for ax, values, name, colour in ((ax_j, j, "Jaccard", "#1f4e79"),
                                     (ax_d, d, "Dice", "#117864")):
        ax.hist(values, bins=bins, color=colour, alpha=0.85, edgecolor="white")
        for threshold, style in zip(thresholds, ("-", "--")):
            ax.axvline(threshold, color="#c0392b", ls=style, lw=1.6,
                       label=f"preservation threshold {threshold:g}")
        if values.size:
            ax.axvline(float(np.median(values)), color="black", lw=1.8,
                       label=f"median = {np.median(values):.3f}")
        ax.set_xlabel(f"covered-residue {name} vs FP_original")
        ax.set_xlim(0, 1)
        ax.legend(fontsize=8, frameon=False)
    ax_j.set_ylabel("iterations")
    fig.suptitle("F9 — geometric footprint robustness: similarity distributions\n"
                 "failed iterations are included and contribute J = 0", fontsize=10)
    fig.tight_layout()
    return _save(fig, out_dir, "F9_robustness_similarity_distributions")


def figure_f10(out_dir: Path, sensitivity_rows: list[dict]) -> list[Path]:
    """F10 — per-center influence ranking with contributing iteration counts."""
    rows = [r for r in sensitivity_rows if r["mean_jaccard_drop"] is not None]
    fig, ax = plt.subplots(figsize=(9, max(3.0, 0.32 * max(len(rows), 1) + 1.6)))
    if rows:
        labels = [f"{r['center_residue_index']} (n={r['n_iterations_removing']})"
                  for r in rows]
        values = [r["mean_jaccard_drop"] for r in rows]
        errors = [r["sd_jaccard_drop"] or 0.0 for r in rows]
        y = np.arange(len(rows))
        ax.barh(y, values, xerr=errors, color="#7d3c98", alpha=0.85,
                error_kw={"lw": 1.0, "ecolor": "#4a235a"})
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()
    ax.set_xlabel("mean drop in covered-residue Jaccard when this center is removed")
    ax.set_ylabel("significant hotspot center (contributing iterations)")
    ax.set_title("F10 — center-loss sensitivity ranking", fontsize=10)
    fig.tight_layout()
    return _save(fig, out_dir, "F10_center_influence_ranking")
