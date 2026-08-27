"""Workflow v2 §2 — oligomeric assembly determination and its three branches.

Unit-level coverage for the classifier and the admissibility gate, plus one
full-``run_stage_a`` test per branch named in the Lead's task:
  (a) confident oligomer, monomer exception authorized -> declared, ADVISORY,
      ``intra_subunit_only`` flag set true, handoff carries it;
  (b) confident oligomer, monomer exception NOT authorized -> BLOCKED with
      ``OLIGOMERIC_ASSEMBLY_UNAVAILABLE``;
  (c) no subunit signal at all -> recorded as ``unknown`` (never silently
      "monomer"), no block, no declaration.

No network: everything runs through the synthetic fixture and a local config
overlay (never ``config/pipeline.yaml`` itself, which stays untouched on disk).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from hotspot3d.data.sources import (
    OLIGOMERIC_STATE_MONOMER,
    OLIGOMERIC_STATE_OLIGOMERIC,
    OLIGOMERIC_STATE_UNKNOWN,
    classify_oligomeric_state,
    extract_subunit_comment,
)
from hotspot3d.data.stage import run_stage_a
from hotspot3d.structure.sources import assert_oligomeric_assembly_permitted
from hotspot3d.utils.config import load_config
from hotspot3d.utils.errors import BlockedError
from hotspot3d.utils.runctx import RunContext

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
PIPELINE_CONFIG = REPO_ROOT / "config" / "pipeline.yaml"


# --- classify_oligomeric_state (pure) ---------------------------------------

@pytest.mark.parametrize("text,expected_state,expects_phrase", [
    ("Homodimer.", OLIGOMERIC_STATE_OLIGOMERIC, True),
    ("Homotetramer.", OLIGOMERIC_STATE_OLIGOMERIC, True),
    ("Homooligomer.", OLIGOMERIC_STATE_OLIGOMERIC, True),
    ("Heterodimer of ATP1A1 and ATP1B1.", OLIGOMERIC_STATE_OLIGOMERIC, True),
    ("Heterotetramer of 2 alpha and 2 beta chains.", OLIGOMERIC_STATE_OLIGOMERIC, True),
    ("Monomer.", OLIGOMERIC_STATE_MONOMER, False),
    ("Monomer in solution.", OLIGOMERIC_STATE_MONOMER, False),
    ("Interacts with CALM1 in a calcium-dependent manner.", OLIGOMERIC_STATE_MONOMER, False),
    ("Binds DNA as part of a larger transcription complex.", OLIGOMERIC_STATE_MONOMER, False),
])
def test_classify_oligomeric_state_confident_matches(text, expected_state, expects_phrase):
    result = classify_oligomeric_state(text)
    assert result["oligomeric_state"] == expected_state
    assert (result["matched_phrase"] is not None) is expects_phrase
    assert result["subunit_text"] == text


@pytest.mark.parametrize("text", [None, "", "   "])
def test_classify_oligomeric_state_absent_comment_is_unknown_not_monomer(text):
    result = classify_oligomeric_state(text)
    assert result["oligomeric_state"] == OLIGOMERIC_STATE_UNKNOWN
    assert result["subunit_text"] is None
    assert result["matched_phrase"] is None


def test_classify_oligomeric_state_never_over_declares_on_bare_interaction_mentions():
    # A monomeric enzyme that merely "interacts with" a partner must not be
    # read as an oligomer; only homo-/hetero-<n>mer language is confident enough.
    result = classify_oligomeric_state(
        "Interacts with the homologous recombination machinery."
    )
    assert result["oligomeric_state"] == OLIGOMERIC_STATE_MONOMER


def test_extract_subunit_comment_reads_only_subunit_comment_type():
    entry = {
        "comments": [
            {"commentType": "FUNCTION", "texts": [{"value": "Catalyzes X."}]},
            {"commentType": "SUBUNIT", "texts": [{"value": "Homotetramer."}]},
        ]
    }
    assert extract_subunit_comment(entry) == "Homotetramer."


def test_extract_subunit_comment_absent_returns_none():
    assert extract_subunit_comment({}) is None
    assert extract_subunit_comment({"comments": []}) is None
    assert extract_subunit_comment(
        {"comments": [{"commentType": "FUNCTION", "texts": [{"value": "x"}]}]}
    ) is None


# --- assert_oligomeric_assembly_permitted (the gate) ------------------------

def test_gate_blocks_confident_oligomer_when_monomer_not_authorized():
    with pytest.raises(BlockedError, match="OLIGOMERIC_ASSEMBLY_UNAVAILABLE"):
        assert_oligomeric_assembly_permitted(OLIGOMERIC_STATE_OLIGOMERIC, allow_monomer=False)


def test_gate_permits_confident_oligomer_when_monomer_authorized():
    assert_oligomeric_assembly_permitted(OLIGOMERIC_STATE_OLIGOMERIC, allow_monomer=True)


@pytest.mark.parametrize("state", [OLIGOMERIC_STATE_MONOMER, OLIGOMERIC_STATE_UNKNOWN])
def test_gate_never_blocks_monomer_or_unknown_regardless_of_config(state):
    assert_oligomeric_assembly_permitted(state, allow_monomer=False)
    assert_oligomeric_assembly_permitted(state, allow_monomer=True)


# --- full-stage branches (synthetic payloads, no network) -------------------

def _ctx_with_overlay(tmp_path, run_id: str, *, allow_monomer: bool) -> RunContext:
    overlay = tmp_path / "overlay.yaml"
    overlay.write_text(
        yaml.safe_dump({"structure": {"allow_monomer_for_obligate_oligomer": allow_monomer}}),
        encoding="utf-8",
    )
    cfg = load_config(PIPELINE_CONFIG, gene_overlay=overlay)
    return RunContext.create(gene="SYNGENE", config=cfg, results_root=tmp_path / "results",
                             synthetic=True, run_id=run_id)


def test_branch_confident_oligomer_monomer_allowed_is_declared_and_advisory(
        clustered, tmp_path):
    ctx = _ctx_with_overlay(tmp_path, "TSTOLIG001", allow_monomer=True)
    case = clustered.with_oligomeric_state(
        "oligomeric", evidence="Homotetramer.", matched_phrase="Homotetramer")
    handoff = run_stage_a(ctx, gene=case.gene, source=case.source)

    assembly = handoff.payload["oligomeric_assembly"]
    assert assembly["oligomeric_state"] == "oligomeric"
    assert assembly["intra_subunit_only"] is True
    assert assembly["allow_monomer_for_obligate_oligomer"] is True

    warnings_path = (ctx.full_results / "03_STRUCTURE_QC" / "warnings_03_structure_qc.tsv")
    text = warnings_path.read_text(encoding="utf-8")
    assert "OLIGOMERIC_ASSEMBLY_MONOMER_DECLARED" in text
    assert "ADVISORY" in text
    assert "OLIGOMERIC_ASSEMBLY_UNAVAILABLE" not in text

    qc = json.loads((ctx.full_results / "03_STRUCTURE_QC" / "structure_qc.json")
                     .read_text(encoding="utf-8"))
    assert qc["oligomeric_assembly"]["intra_subunit_only"] is True
    assert qc["oligomeric_assembly"]["oligomeric_state"] == "oligomeric"


def test_branch_confident_oligomer_monomer_disallowed_is_blocked(clustered, tmp_path):
    ctx = _ctx_with_overlay(tmp_path, "TSTOLIG002", allow_monomer=False)
    case = clustered.with_oligomeric_state(
        "oligomeric", evidence="Heterodimer of A and B.", matched_phrase="Heterodimer")
    with pytest.raises(BlockedError, match="OLIGOMERIC_ASSEMBLY_UNAVAILABLE"):
        run_stage_a(ctx, gene=case.gene, source=case.source)

    warnings_path = (ctx.full_results / "03_STRUCTURE_QC" / "warnings_03_structure_qc.tsv")
    text = warnings_path.read_text(encoding="utf-8")
    assert "OLIGOMERIC_ASSEMBLY_UNAVAILABLE" in text
    assert "BLOCKING" in text
    assert "OLIGOMERIC_ASSEMBLY_MONOMER_DECLARED" not in text

    status = json.loads((ctx.full_results / "03_STRUCTURE_QC" / "stage_status.json")
                        .read_text(encoding="utf-8"))
    assert status["status"] == "BLOCKED"


def test_branch_no_signal_is_recorded_unknown_not_blocked_not_declared(clustered, tmp_path):
    ctx = _ctx_with_overlay(tmp_path, "TSTOLIG003", allow_monomer=False)
    case = clustered.with_oligomeric_state("unknown")   # no SUBUNIT comment at all
    handoff = run_stage_a(ctx, gene=case.gene, source=case.source)

    assembly = handoff.payload["oligomeric_assembly"]
    assert assembly["oligomeric_state"] == "unknown"
    assert assembly["intra_subunit_only"] is False

    warnings_path = (ctx.full_results / "03_STRUCTURE_QC" / "warnings_03_structure_qc.tsv")
    text = warnings_path.read_text(encoding="utf-8")
    assert "OLIGOMERIC_ASSEMBLY_UNAVAILABLE" not in text
    assert "OLIGOMERIC_ASSEMBLY_MONOMER_DECLARED" not in text


def test_branch_confirmed_monomer_is_recorded_distinctly_from_unknown(clustered, tmp_path):
    ctx = _ctx_with_overlay(tmp_path, "TSTOLIG004", allow_monomer=False)
    case = clustered.with_oligomeric_state("monomer", evidence="Monomer.")
    handoff = run_stage_a(ctx, gene=case.gene, source=case.source)

    assembly = handoff.payload["oligomeric_assembly"]
    assert assembly["oligomeric_state"] == "monomer"
    assert assembly["intra_subunit_only"] is False
