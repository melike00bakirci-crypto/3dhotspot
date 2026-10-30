"""Structural quality control — METHOD_SPEC II.13, agent §3.8.

Checks model integrity, chain composition, missing coordinates, non-standard
residues, residue count against the UniProt canonical length, and the numbering
offset established by an **explicit full-sequence alignment**.

The offset is *reported, never silently applied*. An AlphaFold model is numbered
against the UniProt canonical sequence, so the expected offset is 0; a non-zero
offset means the two sequences disagree about what residue *n* is, and silently
shifting one to match the other would fabricate the correspondence that this
stage exists to verify. It is escalated instead.
"""
from __future__ import annotations

import difflib
from dataclasses import dataclass, field
from typing import Sequence

from ..utils.errors import BlockedError, EscalationRequired
from .sources import StructureModel

EXPECTED_NUMBERING_OFFSET = 0


@dataclass(frozen=True)
class Alignment:
    """Result of the explicit model-vs-canonical sequence alignment."""

    offset: int
    n_overlap: int
    n_matched: int
    identity: float
    longest_block: int
    model_length: int
    canonical_length: int
    applied: bool = False

    def as_dict(self) -> dict:
        return {
            "method": "difflib.SequenceMatcher longest matching block",
            "numbering_offset": self.offset,
            "expected_offset": EXPECTED_NUMBERING_OFFSET,
            "n_overlap": self.n_overlap,
            "n_matched": self.n_matched,
            "identity": round(self.identity, 6),
            "longest_matching_block": self.longest_block,
            "model_length": self.model_length,
            "canonical_length": self.canonical_length,
            "offset_applied": self.applied,
            "policy": "reported, never silently applied",
        }


def align_numbering(model_sequence: str, canonical_sequence: str) -> Alignment:
    """Full-sequence alignment giving the model -> canonical numbering offset.

    ``offset`` is defined so that ``canonical_position = model_position + offset``.
    Deterministic: the longest matching block of a ``SequenceMatcher`` with
    autojunk disabled, which for two near-identical sequences is the whole span.
    """
    matcher = difflib.SequenceMatcher(None, model_sequence, canonical_sequence,
                                      autojunk=False)
    blocks = [b for b in matcher.get_matching_blocks() if b.size > 0]
    if not blocks:
        return Alignment(offset=0, n_overlap=0, n_matched=0, identity=0.0,
                         longest_block=0, model_length=len(model_sequence),
                         canonical_length=len(canonical_sequence))

    best = max(blocks, key=lambda b: (b.size, -b.a, -b.b))
    offset = best.b - best.a
    n_matched = 0
    n_overlap = 0
    for i, aa in enumerate(model_sequence):
        j = i + offset
        if 0 <= j < len(canonical_sequence):
            n_overlap += 1
            n_matched += int(aa == canonical_sequence[j])
    identity = (n_matched / n_overlap) if n_overlap else 0.0
    return Alignment(offset=offset, n_overlap=n_overlap, n_matched=n_matched,
                     identity=identity, longest_block=best.size,
                     model_length=len(model_sequence),
                     canonical_length=len(canonical_sequence))


@dataclass
class QCCheck:
    """One recorded QC check. Recorded whether it passes or fails (agent §9)."""

    check_id: str
    description: str
    passed: bool
    detail: str = "NA"
    severity: str = "MAJOR"

    def as_dict(self) -> dict:
        return {"check_id": self.check_id, "description": self.description,
                "passed": self.passed, "detail": self.detail,
                "severity_on_failure": self.severity}


@dataclass
class StructureQCResult:
    checks: list[QCCheck] = field(default_factory=list)
    alignment: Alignment | None = None
    payload: dict = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        return all(c.passed for c in self.checks)

    @property
    def failures(self) -> list[QCCheck]:
        return [c for c in self.checks if not c.passed]

    def as_dict(self) -> dict:
        return {
            **self.payload,
            "alignment": self.alignment.as_dict() if self.alignment else None,
            "checks": [c.as_dict() for c in self.checks],
            "n_checks": len(self.checks),
            "n_failed": len(self.failures),
            "verdict": "PASS" if self.passed else "FAIL",
        }


def structure_qc(model: StructureModel, *, canonical_sequence: str,
                 expected_chains: Sequence[str] = ("A",)) -> StructureQCResult:
    """Run and record every structural QC check. Never repairs, never trims."""
    residues = sorted(model.residues, key=lambda r: (r.chain_id, r.residue_index))
    indices = [r.residue_index for r in residues]
    chains = model.chain_ids
    non_standard = [{"residue_index": r.residue_index, "comp_id": r.aa3}
                    for r in residues if not r.is_standard]
    missing_ca = [r.residue_index for r in residues if r.ca is None]
    zero_occupancy = [r.residue_index for r in residues
                      if r.ca is not None and r.ca.occupancy <= 0.0]
    non_finite = [r.residue_index for r in residues
                  if r.ca is not None and not r.ca.is_finite]
    missing_plddt = [r.residue_index for r in residues if r.plddt is None]
    duplicates = sorted({i for i in indices if indices.count(i) > 1})
    gaps = [(a, b) for a, b in zip(indices, indices[1:]) if b != a + 1]

    alignment = align_numbering(model.sequence, canonical_sequence)

    checks = [
        QCCheck("SQC-1", "model contains at least one residue",
                bool(residues), f"n_residues={len(residues)}", "BLOCKING"),
        QCCheck("SQC-2", "chain composition matches expectation",
                chains == list(expected_chains),
                f"chains={chains};expected={list(expected_chains)}", "BLOCKING"),
        QCCheck("SQC-3", "residue numbering is unique",
                not duplicates, f"duplicate_indices={duplicates[:20]}", "BLOCKING"),
        QCCheck("SQC-4", "residue numbering is contiguous",
                not gaps, f"n_gaps={len(gaps)};first={gaps[:5]}", "MAJOR"),
        QCCheck("SQC-5", "modelled residue count equals UniProt canonical length",
                len(residues) == len(canonical_sequence),
                f"n_residues={len(residues)};uniprot_length={len(canonical_sequence)}",
                "BLOCKING"),
        QCCheck("SQC-6", "every residue is a standard amino acid",
                not non_standard, f"n_non_standard={len(non_standard)}", "MAJOR"),
        QCCheck("SQC-7", "every residue has a CA atom",
                not missing_ca, f"n_missing_ca={len(missing_ca)}", "MAJOR"),
        QCCheck("SQC-8", "every CA has non-zero occupancy and finite coordinates",
                not zero_occupancy and not non_finite,
                f"n_zero_occupancy={len(zero_occupancy)};n_non_finite={len(non_finite)}",
                "MAJOR"),
        QCCheck("SQC-9", "pLDDT present for every modelled residue",
                not missing_plddt, f"n_missing_plddt={len(missing_plddt)}", "BLOCKING"),
        QCCheck("SQC-10", "sequence alignment identity is total over the overlap",
                alignment.identity >= 1.0,
                f"identity={alignment.identity:.6f};n_overlap={alignment.n_overlap}",
                "BLOCKING"),
        QCCheck("SQC-11", "numbering offset is zero (AlphaFold numbers on UniProt)",
                alignment.offset == EXPECTED_NUMBERING_OFFSET,
                f"offset={alignment.offset}", "BLOCKING"),
        QCCheck("SQC-12", "AFDB entry is a single fragment",
                model.n_fragments == 1, f"n_fragments={model.n_fragments}", "BLOCKING"),
    ]

    payload = {
        "accession": model.accession,
        "model_id": model.model_id(),
        "model_version": model.model_version,
        "database_version": model.db_version,
        "n_fragments": model.n_fragments,
        "n_residues_modelled": len(residues),
        "n_atoms": model.n_atoms,
        "uniprot_canonical_length": len(canonical_sequence),
        "chains": chains,
        "residue_index_min": min(indices) if indices else None,
        "residue_index_max": max(indices) if indices else None,
        "n_numbering_gaps": len(gaps),
        "numbering_gaps": [{"after": a, "next": b} for a, b in gaps[:50]],
        "duplicate_residue_indices": duplicates,
        "non_standard_residues": non_standard,
        "residues_missing_ca": missing_ca,
        "residues_zero_occupancy_ca": zero_occupancy,
        "residues_non_finite_ca": non_finite,
        "residues_missing_plddt": missing_plddt,
        "structure_edited": False,
        "sequence_trimmed": False,
        "numbering_repaired": False,
    }
    return StructureQCResult(checks=checks, alignment=alignment, payload=payload)


def assert_numbering_usable(result: StructureQCResult) -> None:
    """Escalate a non-zero offset rather than applying it (agent §8, §15)."""
    alignment = result.alignment
    if alignment is None:
        raise BlockedError("STRUCTURE_MAPPING_FAILURE — no alignment was computed.")
    if alignment.offset != EXPECTED_NUMBERING_OFFSET:
        raise EscalationRequired(
            ambiguity=(
                f"The AlphaFold model and the UniProt canonical sequence disagree about "
                f"residue numbering by {alignment.offset} positions "
                f"(alignment identity {alignment.identity:.4f} over "
                f"{alignment.n_overlap} residues)."
            ),
            options=[
                "Escalate and stop",
                f"Apply the {alignment.offset}-residue offset as a recorded transformation",
                "Re-resolve the accession/isoform and retry",
            ],
            consequences=[
                "Stage A stops; no cohort is produced and no coordinate is claimed.",
                "Every ClinVar position would be shifted; if the offset is wrong the "
                "entire spatial analysis is silently wrong, and silent numbering repair "
                "is prohibited.",
                "The likeliest root cause is an isoform mismatch, which is a Lead decision.",
            ],
            recommendation=(
                "Escalate and stop. Verify the canonical isoform before any offset is "
                "applied; an offset is never applied by Stage A on its own authority."
            ),
        )
