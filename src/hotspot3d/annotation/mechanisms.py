"""Functional-mechanism curation — ¶55, METHOD_SPEC II.13.

The order of operations in this module is the scientific point of the whole stage.
Curation happens **first and blind**: :func:`curate` receives literature records and
the frozen rubric, and nothing else. It has no access to hotspot intervals, so a
mechanism *cannot* be inferred from hotspot membership. Only afterwards does
:func:`mechanism_by_hotspot` overlay the finished labels onto the independently
discovered regions.

:func:`assert_no_membership_inference` turns that argument into a test: it re-derives
every label from the evidence alone and requires an exact match, which holds only if
the label is a pure function of the evidence.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from ..utils.errors import BlockedError
from .rubric import (
    ASSERTABLE,
    Evidence,
    Mechanism,
    MechanismRecord,
    RANK,
    Resolution,
    Rubric,
    resolve_mechanism,
)

VariantKey = tuple[int, str]


@dataclass
class CuratedVariant:
    """One variant, its evidence, and the rubric's verdict."""

    residue_index: int
    aa_change: str
    resolution: Resolution
    records: list[MechanismRecord] = field(default_factory=list)

    @property
    def key(self) -> VariantKey:
        return (self.residue_index, self.aa_change)

    @property
    def mechanism(self) -> Mechanism:
        return self.resolution.mechanism

    @property
    def evidence_strength(self) -> Evidence:
        return self.resolution.evidence_strength

    def meets(self, minimum: Evidence) -> bool:
        return RANK[self.evidence_strength] >= RANK[minimum]

    def as_row(self) -> dict:
        res = self.resolution
        return {
            "residue_index": self.residue_index,
            "aa_change": self.aa_change,
            "mechanism": res.mechanism.value,
            "reference": ";".join(res.references) if res.references else "NA",
            "experimental_system": _join(self.records, "experimental_system"),
            "assay_type": _join(self.records, "assay_type"),
            "principal_finding": _join(self.records, "principal_finding"),
            "evidence_strength": res.evidence_strength.value,
            "curator_note": res.curator_note,
        }


def _join(records: list[MechanismRecord], attr: str) -> str:
    """Retain every source's value; contradicting sources are never dropped."""
    values: list[str] = []
    for rec in records:
        if not rec.included:
            continue
        value = str(getattr(rec, attr) or "").strip()
        if value and value != "NA" and value not in values:
            values.append(value)
    return "; ".join(values) if values else "NA"


# --- curation (blind to hotspots) -------------------------------------------

def curate(records: list[MechanismRecord], rubric: Rubric) -> list[CuratedVariant]:
    """Group evidence per variant and apply the frozen rubric.

    Deterministic: output is sorted by residue then amino-acid change, and the rubric
    itself is order-independent, so the same evidence always yields the same labels.
    """
    grouped: dict[VariantKey, list[MechanismRecord]] = {}
    for rec in records:
        grouped.setdefault((int(rec.residue_index), str(rec.aa_change)), []).append(rec)

    curated = [
        CuratedVariant(
            residue_index=key[0],
            aa_change=key[1],
            resolution=resolve_mechanism(recs, rubric),
            records=recs,
        )
        for key, recs in grouped.items()
    ]
    curated.sort(key=lambda cv: (cv.residue_index, cv.aa_change))
    return curated


def assert_no_membership_inference(
    curated: list[CuratedVariant], rubric: Rubric
) -> dict:
    """QC §9.6 — prove no label depended on anything but the evidence.

    Re-resolves each variant from its records in isolation. Any divergence means some
    other input reached the labelling decision, which is the circularity failure the
    architecture exists to prevent.
    """
    divergent: list[dict] = []
    for cv in curated:
        replay = resolve_mechanism(cv.records, rubric)
        if (replay.mechanism, replay.evidence_strength) != (
            cv.mechanism, cv.evidence_strength
        ):
            divergent.append(
                {
                    "residue_index": cv.residue_index,
                    "aa_change": cv.aa_change,
                    "assigned": cv.mechanism.value,
                    "evidence_only": replay.mechanism.value,
                }
            )
    if divergent:
        raise BlockedError(
            f"CIRCULARITY violation: {len(divergent)} mechanism label(s) are not a "
            f"pure function of their evidence: {divergent[:3]}. A mechanism may never "
            f"be inferred from hotspot membership (agent §8)."
        )
    return {
        "check": "no_mechanism_inferred_from_hotspot_membership",
        "n_variants_replayed": len(curated),
        "n_divergent": 0,
        "result": "PASS",
        "method": (
            "Each label re-derived from its literature records alone, with no hotspot "
            "context available to the rubric; results compared for exact equality."
        ),
    }


# --- counts and the overlay -------------------------------------------------

def mechanism_counts(curated: list[CuratedVariant]) -> dict[str, int]:
    """Counts over every frozen category, including the zeros (P3: absence explicit)."""
    counts = {m.value: 0 for m in Mechanism}
    for cv in curated:
        counts[cv.mechanism.value] += 1
    return counts


def counts_at_or_above(
    curated: list[CuratedVariant], minimum: Evidence, *, assertable_only: bool = True
) -> dict[str, int]:
    """Per-category counts restricted to variants meeting an evidence bar.

    ``assertable_only`` excludes ``Unclear`` and ``Not_Experimentally_Characterized``:
    they are rubric outcomes, not mechanisms, and a spatial test over them would be
    testing the absence of evidence rather than a biological hypothesis.
    """
    counts: dict[str, int] = {
        m.value: 0 for m in (ASSERTABLE if assertable_only else set(Mechanism))
    }
    for cv in curated:
        if assertable_only and cv.mechanism not in ASSERTABLE:
            continue
        if cv.meets(minimum):
            counts[cv.mechanism.value] += 1
    return dict(sorted(counts.items()))


def mechanism_by_hotspot(
    curated: list[CuratedVariant],
    covered_residues: dict[str, list[dict]],
    *,
    minimum: Evidence,
    row_context: dict[str, dict],
) -> list[dict]:
    """Overlay finished labels onto the independently discovered regions.

    This runs strictly after :func:`curate`. ``row_context`` supplies the robustness
    summary and clustering flag that must travel with every hotspot-level statement.
    """
    by_hotspot: dict[str, list[int]] = {
        hid: sorted({int(r["residue_index"]) for r in rows})
        for hid, rows in covered_residues.items()
    }

    rows: list[dict] = []
    for hid in sorted(by_hotspot):
        residues = set(by_hotspot[hid])
        inside = [cv for cv in curated if cv.residue_index in residues]
        context = row_context.get(hid, {})
        for mech in Mechanism:
            members = [cv for cv in inside if cv.mechanism is mech]
            if not members:
                continue
            rows.append(
                {
                    "hotspot_id": hid,
                    "mechanism": mech.value,
                    "n_variants": len(members),
                    "n_at_or_above_moderate": sum(
                        1 for cv in members if cv.meets(minimum)
                    ),
                    "residues": ";".join(str(cv.residue_index) for cv in members),
                    "robustness_profile_summary": context.get(
                        "robustness_profile_summary", "NA"
                    ),
                    "global_clustering_flag": context.get("global_clustering_flag", "NA"),
                }
            )
    return rows


# --- literature search log --------------------------------------------------

def build_search_log(records: list[MechanismRecord]) -> list[dict]:
    """Reproducible log: every query, database, date, hit count and decision.

    Excluded records appear alongside included ones — an exclusion that leaves no
    trace is indistinguishable from a record that was never retrieved (agent §9.7).
    """
    rows = [
        {
            "query": rec.query,
            "database": rec.database,
            "search_date": rec.search_date,
            "hit_count": rec.hit_count,
            "record_reference": rec.reference,
            "residue_index": rec.residue_index,
            "decision": "INCLUDED" if rec.included else "EXCLUDED",
            "decision_reason": (
                "Meets curation criteria; graded "
                f"{rec.evidence_strength.value} under the frozen rubric."
                if rec.included
                else rec.exclusion_reason
            ),
        }
        for rec in records
    ]
    rows.sort(
        key=lambda r: (
            str(r["database"]), str(r["query"]), int(r["residue_index"]),
            str(r["record_reference"]),
        )
    )
    return rows
