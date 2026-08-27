"""Figure reproducibility across stage boundaries — Lead-owned cross-agent invariant.

matplotlib rcParams are process-global mutable state, so figure determinism is not
a property any one agent's module can hold on its own. ``svg.hashsalt`` decides
whether SVG element ids are derived from a fixed salt or from ``uuid.uuid4()``, and
whichever module touches it last wins for the rest of the process.

That is exactly how the defect this module pins arose. Stage A's ``_pyplot()``
calls ``plt.rcdefaults()`` for a clean baseline, which reset ``svg.hashsalt`` to
None; Stage B set it at import and never re-applied it. Stage A draws before every
other stage in a real run, so every vector figure the run wrote afterwards was
non-reproducible — Stage B's F1-F6, and, with no assertion covering them at all,
Stage C's F7-F8 and Stage D's F9-F10.

The existing Stage B determinism tests could not catch it: they run inside their
own package, and a targeted ``pytest tests/hotspot_stats`` executes Stage B before
any Stage A figure exists. The failure appeared only in a full suite, where
``tests/data_structure`` runs first — which reads as flakiness rather than as the
ordering dependency it is. **The ORDER is the test here**, so these assertions
provoke it directly instead of relying on collection order to do it.

Nothing here is scientific. Both parameters govern serialization, never a plotted
value.
"""
from __future__ import annotations

import hashlib
import re
import tempfile
from pathlib import Path

import numpy as np
import pytest

from hotspot3d.utils.figures import DETERMINISTIC_RCPARAMS

pytestmark = pytest.mark.integration

SRC = Path(__file__).resolve().parents[2] / "src" / "hotspot3d"

#: Every module that writes figures. Each must apply the shared settings, because
#: any of them may be the first to import matplotlib in a given process.
FIGURE_MODULES = ["data/figures.py", "hotspot/figures.py",
                  "footprint/figures.py", "robustness/figures.py"]


class _Stats:
    """The minimum F1 needs; the values are irrelevant, only that they are fixed."""

    def __init__(self):
        self.grid = np.linspace(1.0, 20.0, 40)
        self.k_obs = self.grid ** 2
        self.k_null_mean = self.grid ** 2 * 0.9
        self.k_env_lo = self.grid ** 2 * 0.8
        self.k_env_hi = self.grid ** 2 * 1.0
        self.l_obs = self.grid * 1.1
        self.l_env_lo = self.grid * 1.0
        self.l_env_hi = self.grid * 1.2
        self.n_points = 44
        self.t_obs = 3.14159
        self.p_global = 1e-4


def _render_f1() -> dict[str, str]:
    from hotspot3d.hotspot.figures import figure_f1_ripley

    out = Path(tempfile.mkdtemp())
    figure_f1_ripley(out, {"PLP": _Stats(), "BLB": _Stats()})
    return {p.suffix: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((out / "figures").iterdir())}


def test_stage_a_drawing_first_does_not_break_stage_b_figures():
    """The regression: Stage A draws, then Stage B must still be reproducible.

    Ordered deliberately — render, poison, render — so the assertion does not
    depend on which package pytest happens to collect first.
    """
    import matplotlib

    before_a, before_b = _render_f1(), _render_f1()
    assert before_a == before_b, "Stage B figures are not reproducible even alone"

    # Exactly what Stage A does before drawing any of its own figures.
    from hotspot3d.data.figures import _pyplot
    _pyplot()

    assert matplotlib.rcParams["svg.hashsalt"] == DETERMINISTIC_RCPARAMS["svg.hashsalt"], (
        "Stage A reset svg.hashsalt and did not restore it; every SVG written "
        "after Stage A in this process now has uuid4-derived element ids")
    assert matplotlib.rcParams["path.simplify"] == DETERMINISTIC_RCPARAMS["path.simplify"]

    after_a, after_b = _render_f1(), _render_f1()
    assert after_a == after_b, (
        "Stage B figures became non-reproducible after Stage A drew — the exact "
        "defect that failed 04_GLOBAL_CLUSTERING/figures/F1_ripleys_k.svg in the "
        "full suite while passing in every targeted run")
    assert after_a == before_a, (
        "Stage B figures changed across a Stage A call: same inputs, different bytes")


@pytest.mark.parametrize("module", FIGURE_MODULES)
def test_every_figure_module_applies_the_shared_settings(module):
    """A new figure module must not silently reintroduce the defect."""
    text = (SRC / module).read_text(encoding="utf-8")
    assert "apply_deterministic_rcparams" in text, (
        f"{module} writes figures without applying the shared determinism "
        f"settings from utils.figures — whichever module imports last would "
        f"decide whether this run's SVGs are reproducible")


@pytest.mark.parametrize("module", FIGURE_MODULES)
def test_no_figure_module_defines_the_settings_itself(module):
    """One definition, in utils — the same rule the Lead audit applies elsewhere."""
    text = (SRC / module).read_text(encoding="utf-8")
    for param in DETERMINISTIC_RCPARAMS:
        assert not re.search(rf"rcParams\[[\"']{re.escape(param)}[\"']\]\s*=", text), (
            f"{module} sets {param!r} directly; it belongs to "
            f"utils.figures.DETERMINISTIC_RCPARAMS so a reset has one place to "
            f"restore from")


@pytest.mark.parametrize("module", ["hotspot/figures.py", "footprint/figures.py",
                                    "robustness/figures.py"])
def test_svg_writers_suppress_the_timestamp(module):
    """An embedded <dc:date> makes two identical renders differ."""
    text = (SRC / module).read_text(encoding="utf-8")
    assert "SVG_METADATA" in text, (
        f"{module} writes SVG without the shared metadata suppression, so each "
        f"file carries the wall-clock time it was written")
