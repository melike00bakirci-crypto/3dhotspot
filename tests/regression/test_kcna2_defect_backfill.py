"""Permanent regression coverage for implementation-defect classes this pipeline
has actually hit, per the QA architecture's Section 9 (regression promotion).

This module has two jobs:

1. Record, for each defect class named in the QA architecture spec, WHERE its
   regression coverage already lives — so a future reader (or QA itself) does
   not re-derive "is this covered?" from scratch, and so nobody duplicates
   protection that already exists (explicitly prohibited).
2. Add the ONE genuinely uncovered class: an invalid test fixture that could
   silently mix a FROZEN scientific parameter into a real pipeline run.

Fast, synthetic, no network, no real gene — every case here reproduces from a
tiny fixture, not the originating real-data run.
"""
from __future__ import annotations

import pytest

pytestmark = pytest.mark.integration

# --------------------------------------------------------------------------- #
# Coverage audit — 7 of 8 listed classes already have regression protection.  #
# Verified present as of this writing; each assertion below is a TRIPWIRE:    #
# if the referenced test is ever removed or renamed without a replacement,    #
# this module fails loudly instead of silently losing coverage.               #
# --------------------------------------------------------------------------- #


def test_cross_stage_schema_handoff_mismatch_is_covered():
    """Already covered: HANDOFF_REQUIRED_KEYS/assert_handoff_shape exercised in
    tests/data_structure/test_stage_a.py, tests/hotspot_stats/test_stage.py,
    tests/hotspot_stats/test_stage_a_to_b.py, tests/annotation/test_stage_e.py,
    tests/integration/test_lead_audit.py — plus the new deterministic validator
    in test_qa_architecture.py::test_case1_*. Not duplicated here."""
    import importlib

    for module in ("tests.hotspot_stats.test_stage_a_to_b",
                   "tests.integration.test_lead_audit"):
        importlib.import_module(module)  # fails loudly if the module is gone


def test_incorrect_terminal_state_propagation_is_covered():
    """Already covered: tests/integration/test_underpowered_terminal_state.py
    proves UNDERPOWERED reaches the run level, is never BLOCKED/TECHNICAL_FAILURE
    or COMPLETED_NEGATIVE, and downstream stages come back NOT_RUN/NOT_APPLICABLE.
    Not duplicated here."""
    import importlib
    importlib.import_module("tests.integration.test_underpowered_terminal_state")


def test_underpowered_misreported_as_negative_is_covered():
    """Already covered: tests/hotspot_stats/test_kcna2_regression.py reproduces
    the ORIGINATING defect exactly — the KCNA2 run 20260817T155124Z_df8ec59f
    that reported UNDERPOWERED (p_floor=7.1247e-04 > c_1=1.1062e-04 under the
    label-permutation null) as a scientific negative. Not duplicated here."""
    import importlib
    importlib.import_module("tests.hotspot_stats.test_kcna2_regression")


def test_radius_admissibility_significance_as_prerequisite_is_covered():
    """Already covered: tests/hotspot_stats/test_admissibility_v2.py proves a
    radius with zero significant centers stays admissible and remains a full
    Pareto candidate. Not duplicated here."""
    import importlib
    importlib.import_module("tests.hotspot_stats.test_admissibility_v2")


def test_degenerate_vacuous_selection_is_covered():
    """Already covered on BOTH selections: tests/hotspot_stats/test_admissibility_v2.py
    (r_hot) and tests/footprint/test_footprint_v2_selection_diagnostics.py (r_fp)
    each prove VACUOUS_PARETO_SELECTION fires on a single-member Pareto set / an
    all-degenerate objective set, and that degenerate-objective exclusion from
    distance-to-utopia is numerically neutral. Not duplicated here.

    The footprint module is checked by existence + content rather than import:
    it depends on tests/footprint/synthetic_footprint.py via a sys.path insert
    that tests/footprint/conftest.py performs for in-package collection, which
    a cross-directory importlib.import_module here would not have.
    """
    import importlib
    from pathlib import Path

    importlib.import_module("tests.hotspot_stats.test_admissibility_v2")

    footprint_test = (Path(__file__).resolve().parents[1] / "footprint"
                      / "test_footprint_v2_selection_diagnostics.py")
    assert footprint_test.is_file(), footprint_test
    text = footprint_test.read_text(encoding="utf-8")
    assert "VACUOUS_PARETO_SELECTION" in text
    assert "degenerate" in text.lower()


def test_deterministic_rendering_global_state_leakage_is_covered():
    """Already covered: tests/integration/test_figure_determinism.py reproduces
    the exact defect (Stage A's plt.rcdefaults() silently resetting
    svg.hashsalt, making every SVG rendered afterward non-reproducible) and
    proves the fix by rendering, poisoning, and re-rendering in one process.
    Not duplicated here."""
    import importlib
    importlib.import_module("tests.integration.test_figure_determinism")


def test_real_data_source_handoff_wiring_defects_are_covered():
    """Already covered: tests/data_structure/test_live_source.py (isoform-vs-
    fragment detection, AlphaFold version-template wiring, allow_network gating)
    and tests/data_structure/test_clinvar_download.py (resumable ranged download,
    checksum verification, retry accounting). Not duplicated here."""
    import importlib

    importlib.import_module("tests.data_structure.test_live_source")
    importlib.import_module("tests.data_structure.test_clinvar_download")


# --------------------------------------------------------------------------- #
# The one genuinely missing class: an invalid test fixture.                   #
# --------------------------------------------------------------------------- #
#
# Originating defect: an early draft of test_underpowered_terminal_state.py
# built a config overlay setting permutation.B_default=200 and passed it to
# run_pipeline(). The pipeline's own guard (assert_frozen_methodology, called
# unconditionally before any stage runs) correctly rejected it — the defect was
# never in production code. What was missing is a check that would have caught
# the mistake AT THE SOURCE, before ever executing a test, rather than only at
# run time after the fact.
#
# This backfills exactly that: every FROZEN key, individually overridden,
# proven to reject run_pipeline() before any run directory is even created —
# not just B_default, the one case this repository happened to hit.

from hotspot3d.utils.config import FROZEN_ASSERTIONS  # noqa: E402


def _wrong_value(frozen_value):
    """A value guaranteed to differ from the frozen one, of a matching shape."""
    if isinstance(frozen_value, bool):
        return not frozen_value
    if isinstance(frozen_value, (int, float)):
        return frozen_value + 1
    if frozen_value is None:
        return "not_none"
    if isinstance(frozen_value, str):
        return f"{frozen_value}_TAMPERED"
    return frozen_value  # unreached for the current FROZEN_ASSERTIONS shapes


@pytest.mark.parametrize("dotted_key,frozen_value", FROZEN_ASSERTIONS,
                         ids=[k for k, _ in FROZEN_ASSERTIONS])
def test_run_pipeline_rejects_every_frozen_key_before_creating_a_run_root(
        tmp_path, dotted_key, frozen_value):
    """An invalid fixture overlaying ANY frozen key must be rejected by
    run_pipeline() before a single byte of a run is written — proving the
    mechanism that caught the originating mistake covers the whole frozen list,
    not just the one key this repo happened to hit first.
    """
    import yaml

    from hotspot3d.orchestration.pipeline import run_pipeline
    from hotspot3d.utils.config import load_config
    from hotspot3d.utils.errors import BlockedError

    section, leaf = dotted_key.split(".", 1) if "." in dotted_key else (dotted_key, None)
    overlay_value = _wrong_value(frozen_value)
    overlay = {section: {leaf: overlay_value}} if leaf else {section: overlay_value}
    # Nested dotted keys (e.g. none currently in FROZEN_ASSERTIONS, but kept
    # general): only top-level.leaf pairs occur today, so this is exact for now
    # and safely general if a nested one is ever added.

    overlay_path = tmp_path / "invalid_fixture_overlay.yaml"
    overlay_path.write_text(yaml.safe_dump(overlay), encoding="utf-8")
    cfg = load_config("config/pipeline.yaml", gene_overlay=overlay_path)
    assert cfg.get(dotted_key) == overlay_value, \
        "the overlay did not even take effect — test setup is wrong, not the guard"

    results_root = tmp_path / "results"
    with pytest.raises(BlockedError, match="FROZEN methodology violation"):
        run_pipeline(gene="INVALIDFIXTURE", config=cfg, results_root=results_root,
                    synthetic=True, run_id="20250101T000000Z_regress1_regress2")

    assert not results_root.exists() or not any(results_root.rglob("*")), (
        "run_pipeline must reject a frozen-key overlay BEFORE creating any run "
        "artifact — not run partway and then fail")
