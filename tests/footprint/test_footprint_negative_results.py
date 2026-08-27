"""Negative results are valid, complete, terminating outcomes — never failures (P4)."""
from __future__ import annotations

import json

import pytest

from synthetic_footprint import synthetic_run

from hotspot3d.footprint.stage import run_stage_cd
from hotspot3d.utils.io import read_tsv

pytestmark = pytest.mark.unit


def test_empty_domain_terminates_the_chain_honestly(tmp_path):
    """rho_min > rho_max: no footprint exists and nothing is relaxed to invent one."""
    # two adjacent lattice points 2.0 A apart -> rho_all = 1.0 -> rho_max = 1.25
    ctx, upstream = synthetic_run(tmp_path, n_side=6, spacing=2.0,
                                  center_ids=[1, 2], n_plp=20, n_blb=20)
    handoff_03, handoff_04 = run_stage_cd(ctx, upstream=upstream)

    assert handoff_03.negative_result["condition"] == "NO_ADMISSIBLE_FOOTPRINT_DOMAIN"
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
    assert "NO_ADMISSIBLE_FOOTPRINT_DOMAIN" in status_09["reason"]

    profile = handoff_04.payload["robustness_profile"]
    assert profile["evaluable"] is False
    assert profile["upstream_condition"] == "NO_ADMISSIBLE_FOOTPRINT_DOMAIN"

    report = (ctx.full_results / "08_FINAL_FOOTPRINT" / "stage_c_report.md").read_text()
    assert "TRUE SCIENTIFIC NEGATIVE" in report
    assert "No threshold was relaxed" in report


def test_negative_stage_still_writes_its_manifests_and_warnings(tmp_path):
    ctx, upstream = synthetic_run(tmp_path, n_side=6, spacing=2.0,
                                  center_ids=[1, 2], n_plp=20, n_blb=20)
    run_stage_cd(ctx, upstream=upstream)

    for stage in ("07_FOOTPRINT_RADIUS", "08_FINAL_FOOTPRINT"):
        rows = read_tsv(ctx.full_results / stage / "stage_manifest.tsv")
        assert rows
        assert any(r["status"] == "NOT_CREATED" for r in rows), (
            "absence is explicit: files that could not be produced are named")
        assert (ctx.full_results / stage / f"warnings_{stage}.tsv").is_file()
