"""The integration package's configuration overlay — TEST EXECUTION ONLY.

Nothing here touches scientific methodology. One operational budget knob is
reduced so the end-to-end chain is executable in a test suite; no formula,
threshold, seed, permutation count, radius, QC gate or expected output changes.

**Why the reduction exists.** The end-to-end tests run the whole chain, and the
``clustered`` fixture is deliberately decisive: it yields ``n_S = 32`` significant
centers. Stage D perturbs that set, so at the shipped ``robustness.N_CAP = 10000``
one full-chain run plans 10 000 footprint reconstructions — the stage records the
arithmetic itself in ``09_ROBUSTNESS/subset_design.json``::

    n_S = 32   K_MAX = 16   total_subset_space_sum_T_k = 2,448,023,842
    N_CAP = 10000   planned_iterations = 10000
    estimated_wall_seconds = 7950.4          (~2.2 h)

``test_end_to_end.py`` performs ten clustered runs and ``test_seed_determinism.py``
about fifteen more, so the production budget puts the integration package at more
than a day of wall time. This is the cost ``config/pipeline.yaml`` anticipates in
its ``max_wall_seconds`` note ("N_CAP iterations at ~1.5-15 s each is 4-40 h for a
real protein with n_S ~ 25").

**Why reducing it is legitimate.** ``robustness.N_CAP`` is a config knob with a
stated cost model, not a frozen constant: it is absent from ``FROZEN_ASSERTIONS``
in ``hotspot3d.utils.config``, and ``pipeline.yaml`` marks its neighbours
``FROZEN (A19)/(A20)`` while marking N_CAP itself as a knob. The methodology fixes
the *design rule* — ``K_MAX = min(n_S - 1, floor(n_S / 2))``, seeded combinatorial
unranking without replacement, budget allocated equally across remaining levels —
and every one of those still runs here, unmodified. Only how many subsets that
rule is asked to draw changes.

At ``N_CAP = 200`` with ``n_S = 32`` the design keeps k = 1 **exhaustive** (all 32
single-center removals, which the center-recurrence and center-sensitivity metrics
are computed from) and samples 11-12 subsets at each of k = 2..16. Both allocation
modes are exercised and every reported metric stays well-defined; the estimate
drops from ~2.2 h to ~3 min per run.

**Why it lives in its own module.** Two callers need it and neither may import the
other: the ``config`` fixture in ``conftest.py``, and ``test_seed_determinism.py``,
which builds its own config directly rather than taking the fixture. That second
caller is why this file exists — while it loaded ``config/pipeline.yaml`` itself it
silently ran at the production budget, and its ~15 full-chain runs dominated the
package.

The shipped value is asserted to still be 10 000 by
``test_lead_audit.test_production_config_ships_the_full_robustness_budget``, so the
reduction cannot leak into the production configuration unnoticed.
"""
from __future__ import annotations

from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
PIPELINE_CONFIG = REPO_ROOT / "config" / "pipeline.yaml"

#: Perturbation budget for integration runs. See the module docstring for why
#: this is a test-execution choice and not a methodological one.
INTEGRATION_N_CAP = 200

#: The overlay itself, as data. Written to disk rather than passed as an argument
#: so the deviation is declared the way a real one would be — the same mechanism
#: ``tests/hotspot_stats/synthetic_hotspot.make_context`` uses for a reduced ``B``.
INTEGRATION_OVERLAY = {"robustness": {"N_CAP": INTEGRATION_N_CAP}}


def write_integration_overlay(directory: Path) -> Path:
    """Write the overlay into ``directory`` and return its path."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "n_cap_overlay.yaml"
    path.write_text(yaml.safe_dump(INTEGRATION_OVERLAY), encoding="utf-8")
    return path


def load_integration_config(directory: Path):
    """The production config with ``robustness.N_CAP`` reduced, via an overlay.

    Every other parameter — including every FROZEN one — is read from the shipped
    ``config/pipeline.yaml`` unchanged.
    """
    from hotspot3d.utils.config import load_config

    return load_config(PIPELINE_CONFIG,
                       gene_overlay=write_integration_overlay(directory))
