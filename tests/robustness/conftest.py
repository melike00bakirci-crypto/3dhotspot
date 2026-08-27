"""Shared Stage D fixtures.

The synthetic upstream builder lives with the Stage C tests; it is imported rather
than duplicated so both phases are exercised against exactly the same fixture.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "footprint"))

from synthetic_footprint import synthetic_run, two_cluster_centers      # noqa: E402

from hotspot3d.footprint.stage import run_phase_c             # noqa: E402
from hotspot3d.robustness.params import RobustnessParams      # noqa: E402
from hotspot3d.robustness.stage import run_phase_d            # noqa: E402
from hotspot3d.utils.config import load_config                # noqa: E402
from hotspot3d.utils.seeds import SeedRegistry                # noqa: E402

CONFIG_PATH = Path(__file__).resolve().parents[2] / "config" / "pipeline.yaml"


@pytest.fixture(scope="session")
def config():
    return load_config(CONFIG_PATH)


@pytest.fixture(scope="session")
def rparams(config):
    return RobustnessParams.from_config(config)


@pytest.fixture(scope="session")
def seeds(config):
    return SeedRegistry(master_seed=int(config.get("seeding.MASTER_SEED")),
                        run_id="TEST_RUN")


@pytest.fixture(scope="session")
def phase_d_run(tmp_path_factory):
    """One complete Phase C -> FREEZE -> Phase D run on four centers."""
    tmp = tmp_path_factory.mktemp("phase_d")
    ctx, upstream = synthetic_run(tmp, center_ids=two_cluster_centers()[:4])
    handoff_03, phase_c = run_phase_c(ctx, upstream=upstream)
    handoff_04 = run_phase_d(ctx, phase_c=phase_c, handoff_03=handoff_03)
    return ctx, phase_c, handoff_03, handoff_04
