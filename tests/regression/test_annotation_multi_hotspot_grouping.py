"""Regression pin for a real defect found running Stage E against real KCNA2
data for the first time (LiveAnnotationSource work, 2026-08-18).

``hotspot_covered_residues.tsv`` and ``hotspot_classified_variants.tsv`` (Stage
B's real output) name their grouping column ``hotspot_ids`` (plural,
``;``-separated) -- a residue or variant can legitimately fall inside more
than one overlapping hotspot's sphere. ``significant_hotspot_centers.tsv``
names it ``hotspot_id`` (singular) -- a center belongs to exactly one hotspot
by construction. ``annotation/inputs.py::_group_by_hotspot`` only ever looked
up ``hotspot_id`` (singular), so on real data every row silently grouped under
the string ``"None"`` instead of a real hotspot id -- ``data.covered_residues``
and ``data.sphere_variants`` were empty for every real hotspot, so
secondary_structure.tsv, conservation.tsv and plddt_min/plddt_mean all came out
empty/None, while the stage still reported COMPLETED. No existing synthetic
test caught this: ``tests/fixtures/synthetic_annotation.py`` used the same
(wrong) singular column name Stage E's code expected, so the mismatch was
invisible until Stage E ran against real Stage B output for the first time.

This never ran against N_CAP, K_MAX, B, any threshold or any cohort
definition -- it is pure cross-stage schema wiring.
"""
from __future__ import annotations

import pytest

from hotspot3d.annotation.inputs import _group_by_hotspot

pytestmark = pytest.mark.unit


def test_plural_hotspot_ids_column_is_grouped_correctly():
    """The real Stage B schema for covered_residues/sphere_variants."""
    rows = [
        {"hotspot_ids": "H1", "covered_residue_index": "10"},
        {"hotspot_ids": "H2", "covered_residue_index": "20"},
    ]
    grouped = _group_by_hotspot(rows, "covered_residue_index")
    assert set(grouped) == {"H1", "H2"}
    assert grouped["H1"][0]["residue_index"] == "10"
    assert grouped["H2"][0]["residue_index"] == "20"
    assert "None" not in grouped, (
        "a real hotspot_ids value must never fall through to the None bucket")


def test_a_residue_covered_by_two_hotspots_appears_under_both():
    """The exact scenario the plural column name exists for: real geometric
    overlap between two adjacent hotspot spheres."""
    rows = [{"hotspot_ids": "H1;H2", "covered_residue_index": "42"}]
    grouped = _group_by_hotspot(rows, "covered_residue_index")
    assert set(grouped) == {"H1", "H2"}
    assert grouped["H1"][0]["covered_residue_index"] == "42"
    assert grouped["H2"][0]["covered_residue_index"] == "42"


def test_singular_hotspot_id_column_still_works():
    """The real Stage B schema for significant_hotspot_centers.tsv — a
    center is singular by construction; backward compatibility is required."""
    rows = [{"hotspot_id": "H1", "center_residue_index": "258"}]
    grouped = _group_by_hotspot(rows, "center_residue_index")
    assert set(grouped) == {"H1"}
    assert grouped["H1"][0]["residue_index"] == "258"


def test_a_row_with_neither_column_present_does_not_crash():
    """Defensive, pre-existing behavior preserved as-is (out of scope to
    change here): with neither hotspot_ids nor hotspot_id present, the row
    groups under the literal string "None" rather than raising. Real Stage B
    output always has one of the two columns; this only pins that the
    plural-column fix did not introduce a new crash on malformed input."""
    rows = [{"something_else": "x"}]
    grouped = _group_by_hotspot(rows, None)
    assert grouped == {"None": [{"something_else": "x"}]}
