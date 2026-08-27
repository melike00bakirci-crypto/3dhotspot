"""Verbatim copying: values pass through as raw text, and mismatches escalate."""
from __future__ import annotations

import pytest

from hotspot3d.annotation.verbatim import (
    VerbatimAudit,
    load_verbatim_table,
    verbatim_copy,
)
from hotspot3d.utils.errors import BlockedError, EscalationRequired
from hotspot3d.utils.io import write_tsv

pytestmark = pytest.mark.unit

COLUMNS = ["hotspot_id", "n_plp", "fold_enrichment", "p_emp", "q_bh", "Z"]


def make_table(tmp_path, rows=None):
    rows = rows if rows is not None else [
        {"hotspot_id": "HS1", "n_plp": 12, "fold_enrichment": 3.25,
         "p_emp": 0.0001, "q_bh": 0.002, "Z": 4.12},
    ]
    path = write_tsv(tmp_path / "hotspot_regions.tsv", rows, COLUMNS)
    return load_verbatim_table(path, "hotspot_id")


def test_values_are_returned_as_raw_text(tmp_path):
    """Never parsed to float — so there is no representation that could be recomputed."""
    table = make_table(tmp_path)
    value = verbatim_copy(table, "HS1", "p_emp")
    assert isinstance(value, str)
    assert value == "0.0001"


def test_raw_text_survives_a_round_trip_unchanged(tmp_path):
    table = make_table(tmp_path)
    audit = VerbatimAudit(source=table)
    copied = audit.copy_all("HS1", ("n_plp", "fold_enrichment", "p_emp", "q_bh", "Z"))

    out = write_tsv(
        tmp_path / "hotspot_annotation.tsv",
        [{"hotspot_id": "HS1", **copied}], COLUMNS,
    )
    record = audit.verify_written(out)
    assert record["result"] == "PASS"
    assert record["n_mismatches"] == 0
    assert record["n_cells_compared"] == 5


def test_missing_hotspot_escalates(tmp_path):
    table = make_table(tmp_path)
    with pytest.raises(EscalationRequired) as exc:
        verbatim_copy(table, "HS_ABSENT", "p_emp")
    assert "no row" in exc.value.ambiguity


def test_missing_column_escalates_and_refuses_recomputation(tmp_path):
    table = make_table(tmp_path)
    with pytest.raises(EscalationRequired) as exc:
        verbatim_copy(table, "HS1", "not_a_column")
    assert "PROHIBITED" in " ".join(exc.value.consequences)


def test_written_mismatch_raises_escalation_required(tmp_path):
    """The core QC guarantee: a changed statistic cannot reach the output silently."""
    table = make_table(tmp_path)
    audit = VerbatimAudit(source=table)
    copied = audit.copy_all("HS1", ("n_plp", "p_emp"))

    tampered = dict(copied)
    tampered["p_emp"] = "0.05"          # a "tidied" p-value
    out = write_tsv(
        tmp_path / "hotspot_annotation.tsv",
        [{"hotspot_id": "HS1", **tampered}], ["hotspot_id", "n_plp", "p_emp"],
    )

    with pytest.raises(EscalationRequired) as exc:
        audit.verify_written(out)
    assert "byte-for-byte" in exc.value.ambiguity or "do not match" in exc.value.ambiguity
    assert "never self-resolved" in exc.value.recommendation


def test_absent_row_in_output_is_also_a_mismatch(tmp_path):
    table = make_table(tmp_path)
    audit = VerbatimAudit(source=table)
    audit.copy_all("HS1", ("n_plp",))
    out = write_tsv(tmp_path / "hotspot_annotation.tsv", [], ["hotspot_id", "n_plp"])
    with pytest.raises(EscalationRequired):
        audit.verify_written(out)


def test_duplicate_key_is_blocking(tmp_path):
    rows = [
        {"hotspot_id": "HS1", "n_plp": 12, "fold_enrichment": 3.0, "p_emp": 0.01,
         "q_bh": 0.02, "Z": 3.0},
        {"hotspot_id": "HS1", "n_plp": 13, "fold_enrichment": 3.1, "p_emp": 0.02,
         "q_bh": 0.03, "Z": 3.1},
    ]
    path = write_tsv(tmp_path / "dup.tsv", rows, COLUMNS)
    with pytest.raises(BlockedError, match="duplicate"):
        load_verbatim_table(path, "hotspot_id")


def test_missing_file_is_blocking(tmp_path):
    with pytest.raises(BlockedError, match="not found"):
        load_verbatim_table(tmp_path / "nope.tsv", "hotspot_id")
