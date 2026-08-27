"""Shared Stage A test fixtures."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

from hotspot3d.data.synthetic import CASE_NAMES, make_synthetic_case
from hotspot3d.utils.config import load_config
from hotspot3d.utils.runctx import RunContext

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:                 # makes tests.fixtures.* importable
    sys.path.insert(0, str(REPO_ROOT))


@pytest.fixture(scope="session")
def config():
    return load_config(REPO_ROOT / "config" / "pipeline.yaml")


@pytest.fixture
def make_ctx(config, tmp_path):
    """Factory for a RunContext rooted in an isolated tmp directory."""
    counter = {"n": 0}

    def _make(gene: str = "SYNGENE", synthetic: bool = True) -> RunContext:
        counter["n"] += 1
        return RunContext.create(gene, config, results_root=tmp_path / "results",
                                 synthetic=synthetic, run_id=f"TST{counter['n']:03d}")

    return _make


@pytest.fixture(scope="session")
def cases():
    """All eight cases, built once."""
    return {name: make_synthetic_case(name) for name in CASE_NAMES}


@pytest.fixture(scope="session")
def stage_a_cases(cases):
    """The five Stage A behaviour cases.

    The three downstream-branch cases exercise Stage B/C conditions and one of
    them is 620 residues, so tests that run the whole stage per case use this
    subset and keep the suite fast; correctness for all eight is covered by the
    fixture-contract tests.
    """
    return {name: cases[name]
            for name in ("clustered", "no_cluster", "sparse", "conflict", "low_plddt")}


@pytest.fixture(scope="session")
def clustered():
    return make_synthetic_case("clustered")
