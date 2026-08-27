"""ClinVar-to-structure mapping, ``U_struct`` and ``L`` — METHOD_SPEC II.0, F1.

``U_struct`` is every residue of the model with a **usable CA**: the atom is
present, its occupancy is non-zero and its coordinates are finite. ``L`` is the
classified cohort, and ``L`` is a subset of ``U_struct`` by construction — an
eligible residue with no usable CA is excluded with the reason ``no_ca_coordinate``
and appears in the unmapped list, never silently dropped.

F1 is frozen: the residue representation is CA and the distance is CA-CA
Euclidean. CB and side-chain-centroid columns are emitted for future flexibility
and carry ``unused_under_F1 = TRUE`` on every row.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ..utils.geometry import protein_diameter
from .plddt import band_for
from .sources import StructureModel, StructureResidue

RESIDUE_COORDINATE_COLUMNS: tuple[str, ...] = (
    "residue_index", "aa", "x_ca", "y_ca", "z_ca", "plddt", "plddt_band", "ca_usable",
    "x_cb", "y_cb", "z_cb", "x_sc", "y_sc", "z_sc", "unused_under_F1",
)
#: ``U_struct``. Declared to Stage B in handoff_01 as ``stage_b_universe_columns``.
POSITIONAL_UNIVERSE_COLUMNS: tuple[str, ...] = (
    "residue_index", "aa", "x_ca", "y_ca", "z_ca", "plddt", "plddt_band", "ca_usable",
)
#: ``L``. Exactly ``STAGE_B_PERMITTED_COLUMNS`` — no free-text field is present.
CLASSIFIED_COHORT_COLUMNS: tuple[str, ...] = (
    "residue_index", "class", "x_ca", "y_ca", "z_ca", "plddt", "ca_usable",
)
UNMAPPED_COLUMNS: tuple[str, ...] = (
    "residue_index", "aa_ref", "class", "reason", "detail",
)


def _coordinate_row(residue: StructureResidue, bands: Sequence[float]) -> dict:
    ca = residue.ca
    cb = residue.cb
    sc = residue.sidechain_centroid()
    return {
        "residue_index": residue.residue_index,
        "aa": residue.aa,
        "x_ca": ca.x if ca else None,
        "y_ca": ca.y if ca else None,
        "z_ca": ca.z if ca else None,
        "plddt": residue.plddt,
        "plddt_band": band_for(residue.plddt, bands),
        "ca_usable": residue.ca_usable,
        "x_cb": cb.x if cb else None,
        "y_cb": cb.y if cb else None,
        "z_cb": cb.z if cb else None,
        "x_sc": sc[0] if sc else None,
        "y_sc": sc[1] if sc else None,
        "z_sc": sc[2] if sc else None,
        "unused_under_F1": True,
    }


def residue_coordinates(model: StructureModel, bands: Sequence[float]) -> list[dict]:
    """Every modelled residue's coordinates, usable or not (nothing is hidden)."""
    return [_coordinate_row(r, bands)
            for r in sorted(model.residues, key=lambda r: (r.chain_id, r.residue_index))]


def positional_universe(coordinate_rows: Sequence[dict]) -> list[dict]:
    """``U_struct`` — every residue with a usable CA (METHOD_SPEC II.0, F9)."""
    return [{c: row[c] for c in POSITIONAL_UNIVERSE_COLUMNS}
            for row in coordinate_rows if row["ca_usable"]]


def unusable_residue_indices(coordinate_rows: Sequence[dict]) -> list[int]:
    return [row["residue_index"] for row in coordinate_rows if not row["ca_usable"]]


@dataclass
class MappingResult:
    universe: list[dict]
    cohort: list[dict]
    unmapped: list[dict]
    report: dict

    @property
    def M(self) -> int:
        return len(self.universe)

    @property
    def N(self) -> int:
        return len(self.cohort)


def map_cohort(residue_rows, coordinate_rows: Sequence[dict]) -> MappingResult:
    """Place the classified cohort onto the structure and build ``U_struct``/``L``.

    ``residue_rows`` are :class:`hotspot3d.data.cohort.ResidueRow` objects whose
    structural exclusions have already been applied, so an ineligible residue
    here is one the strata already account for.
    """
    universe = positional_universe(coordinate_rows)
    by_index = {row["residue_index"]: row for row in coordinate_rows}
    universe_index = {row["residue_index"] for row in universe}

    cohort: list[dict] = []
    unmapped: list[dict] = []
    for residue in sorted(residue_rows, key=lambda r: r.residue_index):
        coord = by_index.get(residue.residue_index)
        if coord is None:
            unmapped.append({
                "residue_index": residue.residue_index, "aa_ref": residue.aa_ref,
                "class": residue.residue_class, "reason": "not_in_model",
                "detail": "residue index absent from the AlphaFold model",
            })
            continue
        if residue.residue_index not in universe_index:
            unmapped.append({
                "residue_index": residue.residue_index, "aa_ref": residue.aa_ref,
                "class": residue.residue_class, "reason": "no_ca_coordinate",
                "detail": "CA absent, zero occupancy or non-finite coordinates",
            })
            continue
        if not residue.eligible_primary:
            continue
        cohort.append({
            "residue_index": residue.residue_index,
            "class": residue.residue_class,
            "x_ca": coord["x_ca"], "y_ca": coord["y_ca"], "z_ca": coord["z_ca"],
            "plddt": coord["plddt"], "ca_usable": coord["ca_usable"],
        })

    aa_agreement = [
        {"residue_index": r.residue_index, "aa_ref_uniprot": r.aa_ref,
         "aa_model": by_index[r.residue_index]["aa"]}
        for r in residue_rows if r.residue_index in by_index
        and by_index[r.residue_index]["aa"] != r.aa_ref
    ]

    coords = [(row["x_ca"], row["y_ca"], row["z_ca"]) for row in universe]
    report = {
        "representation": "CA",
        "distance_metric": "euclidean",
        "frozen_decision": "F1 — CB and side-chain centroid columns are emitted "
                           "for future flexibility and are unused under F1.",
        "M_positional_universe": len(universe),
        "n_modelled_residues": len(coordinate_rows),
        "n_residues_without_usable_ca": len(coordinate_rows) - len(universe),
        "n_cohort_residues_considered": len(residue_rows),
        "N_classified_cohort": len(cohort),
        "n_unmapped": len(unmapped),
        "unmapped_residues": unmapped,
        "position_correspondence_verified": not aa_agreement,
        "n_position_correspondence_mismatches": len(aa_agreement),
        "position_correspondence_mismatches": aa_agreement[:50],
        "cohort_subset_of_universe": all(
            row["residue_index"] in universe_index for row in cohort),
        "D_max": protein_diameter(coords) if len(coords) > 1 else 0.0,
    }
    return MappingResult(universe=universe, cohort=cohort, unmapped=unmapped,
                         report=report)
