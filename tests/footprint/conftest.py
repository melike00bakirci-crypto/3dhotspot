"""Shared Stage C fixtures.

The full pipeline is run once per session and shared: every assertion below
inspects the *same* artifacts a real run would produce, rather than a hand-built
approximation of them.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from synthetic_footprint import synthetic_run, two_cluster_centers      # noqa: E402

from hotspot3d.footprint.params import FootprintParams        # noqa: E402
from hotspot3d.footprint.stage import run_phase_c             # noqa: E402
from hotspot3d.utils.config import load_config                # noqa: E402

CONFIG_PATH = Path(__file__).resolve().parents[2] / "config" / "pipeline.yaml"


@pytest.fixture(scope="session")
def config():
    return load_config(CONFIG_PATH)


@pytest.fixture(scope="session")
def params(config):
    return FootprintParams.from_config(config)


@pytest.fixture(scope="session")
def phase_c_run(tmp_path_factory):
    """One complete Phase C run on the four-center synthetic geometry."""
    tmp = tmp_path_factory.mktemp("phase_c")
    ctx, upstream = synthetic_run(tmp, center_ids=two_cluster_centers()[:4])
    handoff, result = run_phase_c(ctx, upstream=upstream)
    return ctx, handoff, result
