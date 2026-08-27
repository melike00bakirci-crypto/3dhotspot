"""Shared pytest fixtures.

No test in this suite may touch the network or real biological data. The
``no_network`` autouse fixture makes that structural rather than aspirational:
any attempt to open a socket fails the test.
"""
from __future__ import annotations

import socket
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


@pytest.fixture(autouse=True)
def no_network(monkeypatch):
    """Hard-block every socket connection for the whole suite."""
    def _blocked(*_args, **_kwargs):
        raise RuntimeError(
            "NETWORK ACCESS BLOCKED IN TESTS — no real biological data may be "
            "downloaded during implementation testing."
        )
    monkeypatch.setattr(socket.socket, "connect", _blocked)
    monkeypatch.setattr(socket, "create_connection", _blocked)


@pytest.fixture(scope="session")
def config():
    from hotspot3d.utils.config import load_config
    return load_config(REPO_ROOT / "config" / "pipeline.yaml")


@pytest.fixture
def run_ctx(tmp_path, config):
    """A throwaway RunContext rooted in tmp_path with a FIXED run_id.

    The fixed run_id keeps derived seeds identical across test invocations, which
    is what makes the determinism assertions meaningful.
    """
    from hotspot3d.utils.runctx import RunContext
    return RunContext.create(gene="TESTGENE", config=config,
                             results_root=tmp_path / "results", synthetic=True,
                             run_id="20250101T000000Z_deadbeef_cafebabe")


@pytest.fixture(scope="session")
def repo_root():
    return REPO_ROOT
