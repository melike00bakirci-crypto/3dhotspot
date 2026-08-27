"""Workflow v2 §8 [NEW] — "UNDERPOWERED iteration" investigated, found INAPPLICABLE.

v2 §8's final [NEW] sentence: "Re-running the analysis inside a validation
iteration must re-run §5.4. If removing a subset of centers causes the power
certificate to fail, that iteration is recorded as UNDERPOWERED, not as a
footprint that failed to reproduce."

§5.4's power certificate is a property of the per-center significance TEST
(``p_floor`` vs ``c_1 = q/m``, needing ``N_P``/``N_B``/``m``/the correction level),
computed only when the primary hotspot detection is (re-)run. This repo's FROZEN
config fixes Phase D to pure footprint-geometry reconstruction on an already-frozen
center subset:

  * ``robustness.rerun_hotspot_pipeline: false``   # FROZEN (F11/A7)
  * ``robustness.perturbation_universe: SIGNIFICANT_HOTSPOT_CENTERS``  # FROZEN (A21)

No p-value is ever computed inside a Phase D iteration, so there is no power
certificate to pass or fail, and nothing an ``UNDERPOWERED`` verdict could report.
The antecedent action the [NEW] sentence presupposes IS the withdrawn full-pipeline
re-execution design A7/F11 prohibits under any name. The registered
``ROBUSTNESS_ITERATION_UNDERPOWERED`` warning code (ADVISORY) therefore has no
operational referent under the frozen design and is never emitted here — pinned
below rather than left to be mistaken for an oversight. See the finding recorded
in ``hotspot3d.robustness.perturb``'s module docstring; the Lead makes the final
ruling on this investigation.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "footprint"))

from synthetic_footprint import synthetic_run, two_cluster_centers      # noqa: E402

from hotspot3d.footprint.stage import run_phase_c                       # noqa: E402
from hotspot3d.robustness.perturb import (STATUS_GEOMETRY_FAULT,        # noqa: E402
                                          STATUS_OK)
from hotspot3d.robustness.stage import (V2_SECTION8_UNDERPOWERED_APPLICABLE,  # noqa: E402
                                        assert_no_hotspot_statistics_import,
                                        run_phase_d)
from hotspot3d.utils.config import load_config                          # noqa: E402
from hotspot3d.utils.io import read_tsv                                 # noqa: E402
from hotspot3d.utils.status import WARNING_CODES, Severity              # noqa: E402

pytestmark = pytest.mark.unit

CONFIG_PATH = Path(__file__).resolve().parents[2] / "config" / "pipeline.yaml"


def _stage(ctx):
    return ctx.full_results / "09_ROBUSTNESS"


# --- the frozen config keys that make the clause inapplicable -----------------

def test_the_frozen_keys_cited_by_the_investigation_are_exactly_as_recorded():
    cfg = load_config(CONFIG_PATH)
    assert cfg.get("robustness.rerun_hotspot_pipeline") is False
    assert cfg.get("robustness.perturbation_universe") == "SIGNIFICANT_HOTSPOT_CENTERS"


# --- the mechanical barrier: §5.4's implementation is unreachable from Phase D -

def test_hotspot_power_certificate_module_is_unreachable_from_phase_d():
    """hotspot/power.py implements §5.4; hotspot3d.hotspot as a whole is forbidden."""
    from hotspot3d.robustness.stage import FORBIDDEN_IMPORTS

    assert "hotspot3d.hotspot" in FORBIDDEN_IMPORTS
    assert_no_hotspot_statistics_import()      # raises BlockedError on any violation


# --- the code registers the warning but never emits it -------------------------

def test_the_warning_code_is_registered_but_reserved_not_implemented():
    """Lead-registered for exactly this clause; ADVISORY severity, never emitted."""
    assert WARNING_CODES["ROBUSTNESS_ITERATION_UNDERPOWERED"] == Severity.ADVISORY


def test_no_footprint_or_robustness_source_emits_the_underpowered_code():
    """Source-level pin: the code exists in the vocabulary but is not wired to fire.

    Mirrors the existing barrier-test pattern (test_footprint_barriers.py) rather
    than trusting a single synthetic run to exercise every code path. Only the
    ``warn.add("ROBUSTNESS_ITERATION_UNDERPOWERED", ...)`` emission call is
    checked for — the code is named in prose (module docstrings, this file) to
    document the investigation, which is not an emission.
    """
    src_roots = [Path(__file__).resolve().parents[2] / "src" / "hotspot3d" / "footprint",
                 Path(__file__).resolve().parents[2] / "src" / "hotspot3d" / "robustness"]
    for root in src_roots:
        for path in sorted(root.rglob("*.py")):
            text = path.read_text(encoding="utf-8")
            assert 'warn.add("ROBUSTNESS_ITERATION_UNDERPOWERED"' not in text, path.name


# --- the existing failure vocabulary stays purely geometric --------------------

def test_iteration_status_vocabulary_has_no_statistical_member(phase_d_run):
    """Every iteration status is a GEOMETRY outcome, never a power/significance one.

    STATUS_OK / STATUS_GEOMETRY_FAULT are the complete failure vocabulary — since
    DECISION-STAGE-D-FIXED-RFP-0001 removed per-iteration radius re-selection, the
    only way an iteration can fail is a computational fault in the fixed-radius
    reconstruction. Categorically different from 'the design could not have
    rejected any center' — reusing this status to mean UNDERPOWERED would
    misrepresent what actually failed.
    """
    ctx, _, _, _ = phase_d_run
    rows = read_tsv(_stage(ctx) / "perturbation_results.tsv")
    assert rows

    allowed = {STATUS_OK, STATUS_GEOMETRY_FAULT}
    assert allowed == {"OK", "FAILED_GEOMETRY_FAULT"}
    statuses = {row["status"] for row in rows}
    assert statuses <= allowed
    assert "UNDERPOWERED" not in statuses


def test_no_warning_of_that_code_is_emitted_on_a_real_phase_d_run(phase_d_run):
    ctx, _, _, _ = phase_d_run
    warnings = read_tsv(_stage(ctx) / "warnings_09_ROBUSTNESS.tsv")
    codes = {w["warning_code"] for w in warnings}
    assert "ROBUSTNESS_ITERATION_UNDERPOWERED" not in codes


# --- the finding is surfaced machine-readably, unconditionally -----------------

def test_the_finding_is_a_first_class_field_in_robustness_profile_json(phase_d_run):
    """Lead ruling: always present, never conditional, with a pointer to the
    frozen config key so a reader never has to spelunk in a docstring."""
    ctx, _, _, handoff_04 = phase_d_run
    on_disk = json.loads((_stage(ctx) / "robustness_profile.json").read_text())

    for profile in (on_disk, handoff_04.payload["robustness_profile"]):
        assert profile["v2_section8_underpowered_iteration_applicable"] is \
            V2_SECTION8_UNDERPOWERED_APPLICABLE is False
        note = profile["v2_section8_underpowered_iteration_note"]
        assert "rerun_hotspot_pipeline" in note
        assert "UNDERPOWERED" in note


def test_the_finding_is_present_even_on_the_not_evaluable_path(tmp_path):
    """Unconditional means unconditional: also present when n_S = 1."""
    ctx, upstream = synthetic_run(tmp_path, center_ids=two_cluster_centers()[:1])
    handoff_03, phase_c = run_phase_c(ctx, upstream=upstream)
    handoff_04 = run_phase_d(ctx, phase_c=phase_c, handoff_03=handoff_03)

    profile = handoff_04.payload["robustness_profile"]
    assert profile["evaluable"] is False                # sanity: the right path ran
    assert profile["v2_section8_underpowered_iteration_applicable"] is False
    assert "rerun_hotspot_pipeline" in profile["v2_section8_underpowered_iteration_note"]
