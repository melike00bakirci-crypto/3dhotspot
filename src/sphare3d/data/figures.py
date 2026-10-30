"""Stage A figures — 300 dpi PNG in FULL_RESULTS (Output Contract IX).

The Lead derives the 150 dpi and SVG REVIEW_PACK copies; Stage A emits the
high-resolution originals only. Rendering is deterministic (fixed colours, fixed
ordering, Agg backend) and never blocks the stage: if matplotlib is unavailable
the figure is recorded as NOT_CREATED with a reason, which is what
``RENDERING_UNAVAILABLE`` exists for.
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence

DPI = 300
FIGSIZE_WIDE = (10.0, 4.0)
FIGSIZE_SQUARE = (7.0, 5.0)

CLASS_COLOURS = {
    "PLP": "#b2182b", "BLB": "#2166ac", "CONFLICT": "#7f7f7f",
    "VUS": "#bdbdbd", "CONFLICTING": "#969696", "OTHER": "#d9d9d9",
    "RECORD_SIGNIFICANCE_CONFLICT": "#525252",
}
BAND_COLOURS = ["#d73027", "#fc8d59", "#91bfdb", "#4575b4", "#cccccc"]


class FigureUnavailable(RuntimeError):
    """matplotlib is not importable; the figure is recorded, not fabricated."""


def _pyplot():
    try:
        import matplotlib
    except ImportError as exc:                                  # pragma: no cover
        raise FigureUnavailable(f"matplotlib unavailable: {exc}") from exc
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    from ..utils.figures import apply_deterministic_rcparams

    # rcdefaults() gives Stage A a clean baseline so a stray rcParam set elsewhere
    # in the process cannot alter its figures. It resets EVERY parameter though,
    # including the two that make vector output reproducible — and because Stage A
    # draws before every other stage, dropping them here left every later SVG in
    # the run non-reproducible. The reset is kept; the guarantee is put back.
    plt.rcdefaults()
    apply_deterministic_rcparams()
    return plt


def _save(fig, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI, format="png")
    import matplotlib.pyplot as plt

    plt.close(fig)
    return path


# --- 02_CLINVAR -------------------------------------------------------------

def figure_significance_composition(rows: Sequence[dict], path: Path, gene: str) -> Path:
    """Record counts by significance class — what was retrieved vs what qualifies."""
    plt = _pyplot()
    tally: dict[str, int] = {}
    for row in rows:
        key = str(row.get("significance_class", "NA"))
        tally[key] = tally.get(key, 0) + 1
    labels = sorted(tally, key=lambda k: (-tally[k], k))
    values = [tally[k] for k in labels]

    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    ax.bar(labels, values, color=[CLASS_COLOURS.get(k, "#999999") for k in labels])
    for i, v in enumerate(values):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("ClinVar records")
    ax.set_title(f"{gene} — retrieved missense records by significance class\n"
                 f"(all review star levels retained; F12)", fontsize=10)
    ax.tick_params(axis="x", rotation=30)
    return _save(fig, path)


def figure_star_distribution(rows: Sequence[dict], path: Path, gene: str) -> Path:
    """Review-star composition. Descriptive only — no star ever excluded a record."""
    plt = _pyplot()
    overall = [r for r in rows if r["class"] == "ALL" and r["inclusion_stratum"] == "ALL"]
    tally: dict[str, int] = {}
    for row in overall:
        key = str(row["star_level"])
        tally[key] = tally.get(key, 0) + int(row["n_records"])
    labels = sorted(tally, key=lambda k: (k == "NA", k))
    values = [tally[k] for k in labels]

    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    ax.bar(labels, values, color="#4d4d4d")
    for i, v in enumerate(values):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=8)
    ax.set_xlabel("ClinVar review stars")
    ax.set_ylabel("records")
    ax.set_title(f"{gene} — review-star distribution (recorded, never used to filter)",
                 fontsize=10)
    return _save(fig, path)


def figure_cohort_along_sequence(residue_rows: Sequence[dict], path: Path, gene: str,
                                 sequence_length: int) -> Path:
    """Where the classified residues sit along the sequence."""
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=FIGSIZE_WIDE)
    lanes = {"PLP": 2.0, "BLB": 1.0, "CONFLICT": 0.0}
    for klass, y in lanes.items():
        xs = [r["residue_index"] for r in residue_rows if r["class"] == klass]
        ax.plot(xs, [y] * len(xs), "|", color=CLASS_COLOURS[klass], markersize=12,
                markeredgewidth=1.2, label=f"{klass} (n={len(xs)})")
    ax.set_yticks(list(lanes.values()))
    ax.set_yticklabels(list(lanes))
    ax.set_ylim(-0.6, 2.6)
    ax.set_xlim(0, max(sequence_length, 1))
    ax.set_xlabel("residue index (UniProt canonical numbering)")
    ax.set_title(f"{gene} — residue-level cohort along the sequence", fontsize=10)
    ax.legend(loc="upper right", fontsize=8, ncol=3)
    return _save(fig, path)


# --- 03_STRUCTURE_QC --------------------------------------------------------

def figure_plddt_profile(profile: Sequence[dict], path: Path, gene: str,
                         bands: Sequence[float], cohort: Sequence[int] = ()) -> Path:
    """Per-residue pLDDT with the band cut points drawn, and the cohort marked."""
    plt = _pyplot()
    ordered = sorted(profile, key=lambda r: r["residue_index"])
    xs = [r["residue_index"] for r in ordered]
    ys = [r["plddt"] for r in ordered]

    fig, ax = plt.subplots(figsize=FIGSIZE_WIDE)
    ax.plot(xs, ys, color="#333333", linewidth=0.9)
    for edge, colour in zip(sorted(bands), BAND_COLOURS):
        ax.axhline(edge, color=colour, linestyle="--", linewidth=0.8,
                   label=f"band edge {edge:g}")
    marked = sorted(set(cohort))
    if marked:
        lookup = {r["residue_index"]: r["plddt"] for r in ordered}
        ax.plot(marked, [lookup.get(i) for i in marked], ".", color="#b2182b",
                markersize=3, label=f"classified residues (n={len(marked)})")
    ax.set_ylim(0, 100)
    ax.set_xlabel("residue index")
    ax.set_ylabel("pLDDT")
    ax.set_title(f"{gene} — AlphaFold pLDDT profile (recorded, never a filter; F2)",
                 fontsize=10)
    ax.legend(loc="lower right", fontsize=7, ncol=2)
    return _save(fig, path)


def figure_plddt_bands(summary: dict, path: Path, gene: str) -> Path:
    """Band composition of the model."""
    plt = _pyplot()
    counts = summary.get("band_counts", {})
    labels = list(counts)
    values = [counts[k] for k in labels]

    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    ax.bar(labels, values, color=BAND_COLOURS[:len(labels)])
    for i, v in enumerate(values):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=8)
    ax.set_xlabel("pLDDT band")
    ax.set_ylabel("residues")
    ax.set_title(f"{gene} — pLDDT band distribution (no residue excluded)", fontsize=10)
    return _save(fig, path)


def figure_universe_projection(universe: Sequence[dict], cohort: Sequence[dict],
                               path: Path, gene: str) -> Path:
    """A plain XY projection of U_struct with L overlaid — a sanity view, not analysis."""
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=FIGSIZE_SQUARE)
    ax.scatter([r["x_ca"] for r in universe], [r["y_ca"] for r in universe],
               s=6, c="#cccccc", label=f"U_struct (M={len(universe)})")
    for klass in ("PLP", "BLB"):
        pts = [r for r in cohort if r["class"] == klass]
        ax.scatter([r["x_ca"] for r in pts], [r["y_ca"] for r in pts], s=18,
                   c=CLASS_COLOURS[klass], label=f"{klass} (n={len(pts)})")
    ax.set_xlabel("x (A)")
    ax.set_ylabel("y (A)")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_title(f"{gene} — CA positional universe with the classified cohort",
                 fontsize=10)
    ax.legend(loc="best", fontsize=8)
    return _save(fig, path)
