"""pLDDT extraction, banding and low-confidence regions — F2: RECORD, NEVER FILTER.

``plddt.primary_filtering_enabled`` is FROZEN at FALSE and ``plddt.pae_threshold``
is FROZEN at null. Nothing in this module removes a residue: it produces the
per-residue profile, the four bands (``<50``, ``50-70``, ``70-90``, ``>90``) and
the contiguous low-confidence runs, all of which are *reported* so that a reader
can judge structural confidence for themselves.

The pLDDT >= 70 sensitivity analysis of METHOD_SPEC II.7 belongs to Stage B and
is post-primary and non-redefining. Stage A supplies the numbers; it draws no
line with them.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .sources import StructureModel, StructureResidue

BAND_BELOW = "<{low:g}"
BAND_ABOVE = ">{high:g}"
BAND_MIDDLE = "{low:g}-{high:g}"
BAND_UNKNOWN = "NA"

PLDDT_PROFILE_COLUMNS: tuple[str, ...] = (
    "residue_index", "aa", "chain_id", "plddt", "plddt_band", "ca_usable",
)
PLDDT_REGION_COLUMNS: tuple[str, ...] = (
    "region_id", "threshold", "start_residue", "end_residue", "length",
    "mean_plddt", "min_plddt", "max_plddt", "n_cohort_residues", "excluded",
)


def band_labels(bands: Sequence[float]) -> list[str]:
    """Frozen band labels for the configured cut points (``plddt.bands``)."""
    edges = sorted(float(b) for b in bands)
    labels = [BAND_BELOW.format(low=edges[0])]
    labels += [BAND_MIDDLE.format(low=lo, high=hi) for lo, hi in zip(edges, edges[1:])]
    labels.append(BAND_ABOVE.format(high=edges[-1]))
    return labels


def band_for(value: float | None, bands: Sequence[float]) -> str:
    """Assign one pLDDT to its band. Half-open ``[lo, hi)``; ``None`` -> ``NA``."""
    if value is None:
        return BAND_UNKNOWN
    edges = sorted(float(b) for b in bands)
    labels = band_labels(edges)
    for label, edge in zip(labels, edges):
        if float(value) < edge:
            return label
    return labels[-1]


@dataclass(frozen=True)
class PlddtSummary:
    n_residues: int
    n_with_plddt: int
    mean: float | None
    median: float | None
    minimum: float | None
    maximum: float | None
    band_counts: dict[str, int]
    band_fractions: dict[str, float]

    def as_dict(self) -> dict:
        return {
            "n_residues": self.n_residues, "n_with_plddt": self.n_with_plddt,
            "mean": self.mean, "median": self.median,
            "min": self.minimum, "max": self.maximum,
            "band_counts": self.band_counts, "band_fractions": self.band_fractions,
            "used_as_filter": False,
            "frozen_decision": "F2 — pLDDT is recorded, never a primary filter; "
                               "no PAE threshold exists.",
        }


def plddt_profile(model: StructureModel, bands: Sequence[float]) -> list[dict]:
    """Per-residue pLDDT profile for every modelled residue (no exclusions)."""
    rows: list[dict] = []
    for residue in sorted(model.residues, key=lambda r: (r.chain_id, r.residue_index)):
        value = residue.plddt
        rows.append({
            "residue_index": residue.residue_index,
            "aa": residue.aa,
            "chain_id": residue.chain_id,
            "plddt": value,
            "plddt_band": band_for(value, bands),
            "ca_usable": residue.ca_usable,
        })
    return rows


def summarize(profile: Sequence[dict], bands: Sequence[float]) -> PlddtSummary:
    values = sorted(r["plddt"] for r in profile if r["plddt"] is not None)
    counts = {label: 0 for label in band_labels(bands)}
    counts[BAND_UNKNOWN] = 0
    for row in profile:
        counts[row["plddt_band"]] = counts.get(row["plddt_band"], 0) + 1
    total = float(len(profile)) or 1.0
    n = len(values)
    median = None
    if n:
        median = (values[n // 2] if n % 2 else 0.5 * (values[n // 2 - 1] + values[n // 2]))
    return PlddtSummary(
        n_residues=len(profile), n_with_plddt=n,
        mean=(sum(values) / n if n else None), median=median,
        minimum=(values[0] if n else None), maximum=(values[-1] if n else None),
        band_counts={k: v for k, v in counts.items() if v or k != BAND_UNKNOWN},
        band_fractions={k: round(v / total, 6) for k, v in counts.items()
                        if v or k != BAND_UNKNOWN},
    )


def low_confidence_regions(profile: Sequence[dict], thresholds: Sequence[float],
                           cohort_residues: Sequence[int] = ()) -> list[dict]:
    """Maximal runs of consecutive residues below each threshold.

    Reported for the reader's judgement. ``excluded`` is FALSE on every row and
    is written out explicitly, so the file itself states that identifying a
    low-confidence region never removed anything (F2).
    """
    cohort = set(cohort_residues)
    ordered = sorted(profile, key=lambda r: r["residue_index"])
    regions: list[dict] = []
    region_id = 0

    for threshold in sorted(float(t) for t in thresholds):
        run: list[dict] = []

        def flush(current: list[dict]) -> None:
            nonlocal region_id
            if not current:
                return
            region_id += 1
            values = [r["plddt"] for r in current]
            indices = [r["residue_index"] for r in current]
            regions.append({
                "region_id": f"LC{region_id:04d}",
                "threshold": threshold,
                "start_residue": indices[0], "end_residue": indices[-1],
                "length": len(current),
                "mean_plddt": sum(values) / len(values),
                "min_plddt": min(values), "max_plddt": max(values),
                "n_cohort_residues": sum(1 for i in indices if i in cohort),
                "excluded": False,
            })

        previous_index: int | None = None
        for row in ordered:
            value = row["plddt"]
            below = value is not None and value < threshold
            contiguous = previous_index is not None and row["residue_index"] == previous_index + 1
            if below and contiguous:
                run.append(row)
            elif below:
                flush(run)
                run = [row]
            else:
                flush(run)
                run = []
            previous_index = row["residue_index"]
        flush(run)

    return regions
