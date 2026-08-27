"""Workflow v2 §5.5 — scan/detection consistency.

Three requirements, and the third is the one that was violated:

  1. the radius scan and the final detection draw from **independent seed contexts**;
  2. the final detection is **authoritative** — a provisional scan center is never
     promoted into the published set;
  3. when the §5.4 certificate has FAILED, the disagreement must **not** be described
     as the evidence sitting near the significance boundary. It is evidence that the
     test is operating at its floor, which is a different claim with a different
     remedy.
"""
from __future__ import annotations

import json

import pytest

from hotspot3d.hotspot.stage import (
    SCAN_FINAL_NOTE_CERTIFIED,
    SCAN_FINAL_NOTE_UNCERTIFIED,
    run_stage_b,
)
from hotspot3d.utils.io import read_tsv
from hotspot3d.utils.seeds import ctx_final_detection, ctx_scan, derive_seed
from synthetic_hotspot import full_synthetic_run


@pytest.fixture(scope="module")
def positive_run(tmp_path_factory):
    """One complete Stage B run at the FROZEN default B = 10000, with a real hotspot."""
    ctx, handoff, geo = full_synthetic_run(tmp_path_factory.mktemp("consistency"),
                                           B=10_000)
    return ctx, run_stage_b(ctx, upstream=handoff), geo


@pytest.mark.unit
def test_the_two_readings_of_a_disagreement_are_kept_apart():
    """The certified note may invoke the BH boundary; the uncertified one may not."""
    assert "near the BH boundary" in SCAN_FINAL_NOTE_CERTIFIED
    assert "PASSED" in SCAN_FINAL_NOTE_CERTIFIED

    assert "FAILED" in SCAN_FINAL_NOTE_UNCERTIFIED
    assert "must NOT be described as the evidence sitting near the significance " \
           "boundary" in SCAN_FINAL_NOTE_UNCERTIFIED
    assert "operating at its floor" in SCAN_FINAL_NOTE_UNCERTIFIED

    for note in (SCAN_FINAL_NOTE_CERTIFIED, SCAN_FINAL_NOTE_UNCERTIFIED):
        assert "FINAL detection is authoritative" in note


@pytest.mark.unit
def test_scan_and_final_draw_from_independent_seed_contexts():
    """Different contexts, and therefore genuinely different draws."""
    contexts = {ctx_scan(7.5), ctx_final_detection(),
                f"{ctx_scan(7.5)}|secondary_label_null",
                f"{ctx_final_detection()}|secondary_label_null"}
    assert len(contexts) == 4
    seeds = {derive_seed(20250101, "run", c) for c in contexts}
    assert len(seeds) == 4


@pytest.mark.integration
@pytest.mark.slow
def test_the_final_detection_is_authoritative(positive_run):
    """The published set is the final one, whatever the scan provisionally found."""
    ctx, handoff, _ = positive_run
    r_hot = handoff.payload["hotspot_radius"]

    published = {r["center_residue_index"] for r in read_tsv(
        ctx.full_results / "06_FINAL_HOTSPOTS" / "significant_hotspot_centers.tsv")}
    final_rows = read_tsv(ctx.full_results / "06_FINAL_HOTSPOTS" /
                          "all_residue_center_tests.tsv")
    assert published == {r["center_residue_index"] for r in final_rows
                         if r["significant"]}

    provisional = {r["center_residue_index"] for r in read_tsv(
        ctx.full_results / "05_HOTSPOT_RADIUS" / "radius_scan_centers" /
        f"r_{r_hot:.1f}.tsv")}
    assert handoff.payload["n_provisional_centers_at_r_hot"] == len(provisional)
    assert handoff.payload["scan_final_disagreement"] is (published != provisional)


@pytest.mark.integration
@pytest.mark.slow
def test_a_certified_run_uses_the_certified_wording(positive_run):
    ctx, handoff, _ = positive_run
    assert handoff.payload["power_certificate_passed"] is True
    assert handoff.payload["scan_final_disagreement_note"] == SCAN_FINAL_NOTE_CERTIFIED

    # And the certificate that licenses that wording is on disk beside the result.
    certificate = json.loads((ctx.full_results / "06_FINAL_HOTSPOTS" /
                              "power_certificate.json").read_text())
    assert certificate["PASSES"] is True
    assert certificate["posthoc"]["radius_A"] == handoff.payload["hotspot_radius"]
