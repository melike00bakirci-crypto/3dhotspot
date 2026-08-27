"""Shared fixtures for the Stage B test suite."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from synthetic_hotspot import full_synthetic_run                      # noqa: E402


@pytest.fixture(scope="session", autouse=True)
def _isolate_logs(tmp_path_factory):
    """Keep the stage log out of the developer's working directory.

    Session-scoped, because the expensive stage runs live in module-scoped fixtures
    that are set up outside any function-scoped fixture and would otherwise write
    their logs into the repository.
    """
    root = tmp_path_factory.mktemp("stage_b_logs")
    previous = os.environ.get("HOTSPOT3D_LOG_ROOT")
    os.environ["HOTSPOT3D_LOG_ROOT"] = str(root)
    yield root
    if previous is None:
        os.environ.pop("HOTSPOT3D_LOG_ROOT", None)
    else:
        os.environ["HOTSPOT3D_LOG_ROOT"] = previous


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(20250101))


@pytest.fixture
def small_cohort():
    """A tiny labelled cohort with an obvious cluster — for unit-level checks."""
    coords = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0],
        [40.0, 0.0, 0.0], [41.0, 0.0, 0.0], [40.0, 1.0, 0.0], [80.0, 0.0, 0.0],
    ], dtype=np.float64)
    y = np.array([1, 1, 1, 1, 0, 0, 0, 0], dtype=np.int64)
    return coords, y


@pytest.fixture
def stage_b_run(tmp_path):
    """A complete synthetic Stage A surface plus a verified handoff_01."""
    return full_synthetic_run(tmp_path, B=200)
