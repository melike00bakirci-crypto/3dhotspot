"""Stage A — structure acquisition, CA extraction, pLDDT and structural QC.

Owned by agent ``data-structure``. F1 is frozen: the residue representation is CA
and the distance is CA-CA Euclidean; CB and side-chain-centroid columns are
emitted for future flexibility and marked ``unused_under_F1``.
"""
from __future__ import annotations

from .mapping import map_cohort, positional_universe, residue_coordinates
from .plddt import band_for, low_confidence_regions, plddt_profile
from .qc import align_numbering, structure_qc
from .sources import (
    Atom,
    StructureModel,
    StructureResidue,
    StructureSource,
    assert_model_admissible,
)

__all__ = [
    "Atom", "StructureModel", "StructureResidue", "StructureSource",
    "assert_model_admissible", "structure_qc", "align_numbering",
    "plddt_profile", "band_for", "low_confidence_regions",
    "residue_coordinates", "positional_universe", "map_cohort",
]
