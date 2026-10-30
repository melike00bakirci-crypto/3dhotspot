"""The narrowed F10 exemption — DECISION-RFP-CONTAINMENT-0001.

`ReadOnlyHotspotRadius` exists so that `r_hot` cannot leak into Stage C by accident.
That decision did NOT remove it; it opened exactly one named door,
`derive_evidence_set`, through which only a per-residue boolean verdict passes. These
tests are the mechanical form of that promise: they fail if the door widens, if
another call site appears, or if the numeric radius starts escaping.

Modelled on the Stage A star audit in `sphare3d/data/audit.py`: the guarantee is
checked against the source tree, not against a config flag.
"""
from __future__ import annotations

import ast
import inspect
from pathlib import Path

import numpy as np
import pytest

from sphare3d.footprint import domain as domain_mod
from sphare3d.footprint.api import CenterSet, Universe
from sphare3d.footprint.inputs import EvidenceSet, ReadOnlyHotspotRadius
from sphare3d.utils.errors import BlockedError

pytestmark = pytest.mark.unit

SRC = Path(__file__).resolve().parents[2] / "src" / "sphare3d"
R_HOT = 9.5
#: A value no coordinate or distance in the fixture can coincide with, so that
#: "did r_hot leak" is decidable by equality rather than by eyeballing.
SENTINEL_R_HOT = 1234.5678


def _fixture():
    """Three P/LP residues: two inside r_hot of a center, one far outside."""
    ids = (10, 20, 30, 40)
    coords = np.array([[0.0, 0, 0],      # 10 — the center itself
                       [5.0, 0, 0],      # 20 — inside
                       [9.0, 0, 0],      # 30 — inside, furthest -> binding
                       [40.0, 0, 0]])    # 40 — outside
    universe = Universe(ids=ids, coords=coords)
    centers = CenterSet(ids=(10,), coords=coords[:1], universe_rows=(0,))
    return centers, universe, frozenset({10, 20, 30, 40})


# --------------------------------------------------------------------------- #
# the door is the only door                                                   #
# --------------------------------------------------------------------------- #

def test_derive_evidence_set_is_called_only_from_inputs_load():
    """The exemption may not leave the function that owns it.

    `load` derives two evidence sets — one at `r_hot`, one at the runner-up it
    reports against — so the invariant is not a call *count*. It is that no module
    other than `footprint/inputs.py`, and no function other than `load`, may invoke
    the method. A call anywhere else is a widening of the exemption and fails here.
    """
    call_sites = []
    for path in sorted(SRC.rglob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        enclosing = {}
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                for child in ast.walk(node):
                    enclosing.setdefault(id(child), node.name)
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "derive_evidence_set"):
                call_sites.append((str(path.relative_to(SRC)),
                                   enclosing.get(id(node)), node.lineno))

    assert call_sites, "the audit found no call at all — has the method been renamed?"
    offenders = [c for c in call_sites
                 if c[0] != "footprint/inputs.py" or c[1] != "load"]
    assert not offenders, (
        f"derive_evidence_set may only be called from footprint/inputs.py::load; "
        f"found {offenders}")


def test_numeric_dunders_remain_blocked():
    """The interlock is narrowed, not removed: arithmetic still raises."""
    r = ReadOnlyHotspotRadius(R_HOT)
    for op in (lambda: float(r), lambda: int(r),
               lambda: r < 1, lambda: r <= 1, lambda: r > 1, lambda: r >= 1,
               lambda: r + 1, lambda: r - 1, lambda: r * 2, lambda: r / 2):
        with pytest.raises(BlockedError):
            op()


def test_r_hot_value_never_escapes_the_evidence_set():
    """No field of the returned structure carries the radius."""
    centers, universe, plp = _fixture()
    ev = ReadOnlyHotspotRadius(SENTINEL_R_HOT).derive_evidence_set(
        centers=centers, universe=universe, cohort_plp=plp)

    def _walk(obj):
        if isinstance(obj, (int, float)) and not isinstance(obj, bool):
            yield float(obj)
        elif isinstance(obj, (list, tuple)):
            for item in obj:
                yield from _walk(item)
        elif hasattr(obj, "__dataclass_fields__"):
            for name in obj.__dataclass_fields__:
                yield from _walk(getattr(obj, name))

    assert SENTINEL_R_HOT not in set(_walk(ev)), (
        "the numeric r_hot appeared somewhere in EvidenceSet; the exemption "
        "promises only a per-residue boolean leaves the method")


def test_domain_module_takes_no_r_hot_anywhere():
    """footprint/domain.py keeps its original guarantee, verbatim."""
    offenders = []
    for name, fn in vars(domain_mod).items():
        if not callable(fn) or not hasattr(fn, "__code__"):
            continue
        if getattr(fn, "__module__", None) != domain_mod.__name__:
            continue
        params = inspect.signature(fn).parameters
        offenders += [f"{name}({p})" for p in params if "r_hot" in p.lower()]
    assert offenders == [], (
        f"r_hot reached footprint/domain.py: {offenders}. The containment radius is "
        f"computed from an EvidenceSet precisely so this module never sees it.")


# --------------------------------------------------------------------------- #
# the rule's lemmas                                                           #
# --------------------------------------------------------------------------- #

def _domain(params, r_hot=R_HOT):
    centers, universe, plp = _fixture()
    ev = ReadOnlyHotspotRadius(r_hot).derive_evidence_set(
        centers=centers, universe=universe, cohort_plp=plp)
    dom, rec = domain_mod.derive_containment_domain(
        centers.coords, universe.coords, ev, params)
    return dom, rec, ev


def test_L1_evidence_set_is_nonempty(params):
    _, rec, ev = _domain(params)
    assert rec["n_evidence_set"] >= 1
    assert len(ev.evidence) == 3          # residues 10, 20, 30
    assert rec["n_plp_outside_evidence"] == 1      # residue 40


def test_L2_r_fp_never_exceeds_r_hot(params):
    _, rec, _ = _domain(params)
    assert rec["r_contain_exact_A"] <= R_HOT


def test_L3_every_evidence_residue_is_inside(params):
    _, rec, ev = _domain(params)
    assert all(r.distance_A <= rec["r_contain_A"] + 1e-9 for r in ev.evidence)


def test_L4_one_step_smaller_excludes_the_binding_residue(params):
    """Minimality: the radius is the smallest that satisfies L3, not an argmax."""
    _, rec, ev = _domain(params)
    just_under = rec["r_contain_exact_A"] - 1e-9
    outside = [r for r in ev.evidence if r.distance_A > just_under]
    assert [r.residue_id for r in outside] == [rec["binding_residue"]]


def test_grid_rounding_is_upward_which_is_what_preserves_L3(params):
    """Rounding DOWN would drop the binding residue and break the whole rule."""
    _, rec, _ = _domain(params)
    assert rec["r_contain_A"] >= rec["r_contain_exact_A"]
    assert rec["grid_slack_A"] >= 0.0
    assert rec["grid_slack_A"] < float(params.step_fp_A)


def test_single_point_domain(params):
    dom, rec, _ = _domain(params)
    assert dom.grid == (rec["r_contain_A"],)
    assert dom.rho_min == dom.rho_max == rec["r_contain_A"]
    assert dom.binding_constraint == "containment"


def test_blb_cohort_cannot_change_the_radius(params):
    """B/LB enters no step — the signature cannot even accept it."""
    sig = inspect.signature(ReadOnlyHotspotRadius.derive_evidence_set).parameters
    assert "cohort_blb" not in sig
    assert set(sig) == {"self", "centers", "universe", "cohort_plp"}


def test_outlier_flag_is_recorded_not_left_to_the_reader(params):
    _, rec, _ = _domain(params)
    assert "binding_residue_is_outlier" in rec
    assert isinstance(rec["binding_residue_is_outlier"], bool)
    assert rec["outlier_fence_A"] == pytest.approx(
        rec["evidence_distance_q3_A"] + 1.5 * rec["evidence_distance_iqr_A"])
