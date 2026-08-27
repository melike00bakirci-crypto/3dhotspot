"""Determinism is asserted by re-run, not by inspection (agent §7.3)."""
from __future__ import annotations

import pytest

from synthetic_footprint import synthetic_run, two_cluster_centers

from hotspot3d.footprint.stage import run_phase_c
from hotspot3d.robustness.stage import run_phase_d

pytestmark = [pytest.mark.integration, pytest.mark.slow]

DECISION_BEARING = [
    "09_ROBUSTNESS/perturbation_design.tsv",
    "09_ROBUSTNESS/perturbation_universe.tsv",
    "09_ROBUSTNESS/perturbation_results.tsv",
    "09_ROBUSTNESS/iteration_failures.tsv",
    "09_ROBUSTNESS/center_sensitivity.tsv",
    "09_ROBUSTNESS/robustness_summary.tsv",
    "09_ROBUSTNESS/perturbation_footprints.tar.gz",
]


def _run(root):
    ctx, upstream = synthetic_run(root, center_ids=two_cluster_centers()[:4])
    handoff_03, phase_c = run_phase_c(ctx, upstream=upstream)
    run_phase_d(ctx, phase_c=phase_c, handoff_03=handoff_03)
    return {rel: (ctx.full_results / rel).read_bytes() for rel in DECISION_BEARING}


def test_phase_d_is_byte_identical_on_re_run(tmp_path):
    """Same inputs and the same context-derived seeds -> identical outputs.

    The bundled archive is included deliberately: fixed member order, mtime and
    ownership make it reproducible too, which a naively written tar would not be.
    """
    first, second = _run(tmp_path / "a"), _run(tmp_path / "b")
    for rel in DECISION_BEARING:
        assert first[rel] == second[rel], rel
