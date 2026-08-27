"""Integration-scope fixtures — TEST EXECUTION ONLY.

Overrides the repo-wide ``config`` fixture with one carrying the integration
perturbation budget. The reduction, and the reasoning behind it, live in
``tests/integration/_config.py``; nothing scientific changes.
"""
from __future__ import annotations

import pytest

from tests.integration._config import INTEGRATION_N_CAP, load_integration_config

__all__ = ["INTEGRATION_N_CAP"]


@pytest.fixture(scope="session")
def config(tmp_path_factory):
    """The production config with ``robustness.N_CAP`` reduced, via an overlay.

    Overrides the repo-wide ``config`` fixture for the integration package only.
    Every other parameter — including every FROZEN one — is read from the shipped
    ``config/pipeline.yaml`` unchanged.
    """
    return load_integration_config(tmp_path_factory.mktemp("integration_config"))
