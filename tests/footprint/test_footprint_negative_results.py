"""Negative results are valid, complete, terminating outcomes — never failures (P4)."""
from __future__ import annotations

import json

import pytest

from synthetic_footprint import lattice, residues_beyond, synthetic_run

from hotspot3d.footprint.stage import run_stage_cd
from hotspot3d.utils.io import read_tsv

pytestmark = pytest.mark.unit


NEGATIVE = "NO_EVIDENCE_SET"


def _no_evidence_run(tmp_path):
    """The containment-era negative: an evidence set with no members.

    [DECISION-RFP-CONTAINMENT-0001] The retired trigger was an empty search domain
    (`rho_min > rho_max`). Containment has no search and, by L2, `r_fp <= r_hot`
    always, so the domain is never empty — the one thing that can leave `r_contain`
    undefined is no P/LP residue lying within `r_hot` of any center. L1 says a
    significant center cannot be P/LP-empty, so this is the lemma's own guard, and a
    test must state the cohort to reach it.
    """
    ids, coords = lattice(6, 2.0)
    centers = [1, 2]
    far = residues_beyond(ids, coords, centers, radius_A=8.0, k=20)
    return synthetic_run(tmp_path, n_side=6, spacing=2.0, center_ids=centers,
                         plp_ids=far[:10], blb_ids=far[10:])


def test_an_empty_evidence_set_terminates_the_chain_honestly(tmp_path):
    """No footprint exists and nothing is relaxed to invent one."""
    ctx, upstream = _no_evidence_run(tmp_path)
    handoff_03, handoff_04 = run_stage_cd(ctx, upstream=upstream)

    assert handoff_03.negative_result["condition"] == NEGATIVE

    # [DECISION-REPORTING-CONTEXT-0001] The flag the orchestrator's declared-negative
    # hook keys on. Stage C does not raise to the orchestrator — it records the
    # negative, writes NOT_RUN for 08/09 and returns — so without this flag the hook
    # skips the run, Stage E is invoked, blocks on the absent footprint_residues.tsv,
    # and a legitimate scientific negative is reported BLOCKED/TECHNICAL_FAILURE with
    # negative_result: null. Found live in the v3 batch: CACNA1S and CASQ2 both ended
    # 07/08 = COMPLETED_NEGATIVE and came back BLOCKED. That the hook then honours the
    # flag is covered by
    # test_end_to_end.py::test_unanalysable_cohort_is_a_negative_not_a_technical_failure.
    assert handoff_03.negative_result["settled_negative"] is True, (
        "a Stage C negative is SETTLED; without the flag the orchestrator reports it "
        "as a technical failure")
    assert handoff_03.payload["footprint_radius"] is None
    assert handoff_03.payload["hotspot_radius"] == 8.0      # passthrough survives
    assert handoff_03.qc_status != "FAIL", "a negative is not a QC failure"

    for stage in ("07_FOOTPRINT_RADIUS", "08_FINAL_FOOTPRINT"):
        status = json.loads(
            (ctx.full_results / stage / "stage_status.json").read_text())
        assert status["status"] == "COMPLETED_NEGATIVE"
        assert status["outcome_type"] == "SCIENTIFIC_NEGATIVE"
        assert (ctx.full_results / stage / "NOT_RUN.txt").is_file()

    # 09 exists and names the upstream cause: no mysteriously empty folder
    status_09 = json.loads(
        (ctx.full_results / "09_ROBUSTNESS" / "stage_status.json").read_text())
    assert status_09["status"] == "NOT_RUN"
    assert status_09["outcome_type"] == "NOT_APPLICABLE"
    assert NEGATIVE in status_09["reason"]

    profile = handoff_04.payload["robustness_profile"]
    assert profile["evaluable"] is False
    assert profile["upstream_condition"] == NEGATIVE

    report = (ctx.full_results / "08_FINAL_FOOTPRINT" / "stage_c_report.md").read_text()
    assert "TRUE SCIENTIFIC NEGATIVE" in report
    assert "No threshold was relaxed" in report


def test_negative_stage_still_writes_its_manifests_and_warnings(tmp_path):
    ctx, upstream = _no_evidence_run(tmp_path)
    run_stage_cd(ctx, upstream=upstream)

    for stage in ("07_FOOTPRINT_RADIUS", "08_FINAL_FOOTPRINT"):
        rows = read_tsv(ctx.full_results / stage / "stage_manifest.tsv")
        assert rows
        assert any(r["status"] == "NOT_CREATED" for r in rows), (
            "absence is explicit: files that could not be produced are named")
        assert (ctx.full_results / stage / f"warnings_{stage}.tsv").is_file()
