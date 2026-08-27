"""Shared fixtures for the Stage E suite."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Put the repo root on sys.path so ``tests.fixtures.*`` imports as a namespace package.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tests.fixtures.synthetic_annotation import (  # noqa: E402
    DEFAULT_CONFIG,
    MockAnnotationSource,
    build_synthetic_run,
)

from hotspot3d.annotation.posthoc import PosthocGate  # noqa: E402
from hotspot3d.annotation.rubric import Rubric  # noqa: E402
from hotspot3d.utils.config import load_config  # noqa: E402


@pytest.fixture(scope="session")
def config():
    return load_config(DEFAULT_CONFIG)


@pytest.fixture(scope="session")
def rubric(config) -> Rubric:
    return Rubric.from_config(config)


@pytest.fixture(scope="session")
def gate(config) -> PosthocGate:
    return PosthocGate.from_config(config)


@pytest.fixture
def synthetic_run(tmp_path):
    return build_synthetic_run(tmp_path)


@pytest.fixture
def mock_source() -> MockAnnotationSource:
    return MockAnnotationSource()
