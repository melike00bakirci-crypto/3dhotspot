"""II.4 — the permutation-resolution diagnostic (C0-C3, B_req, B_rec)."""
from __future__ import annotations

import numpy as np
import pytest

from hotspot3d.hotspot.resolution import K_TARGET, b_recommendation, evaluate

LADDER = [10_000, 100_000, 1_000_000, 10_000_000]
MAX_B = 10_000_000
S = 10
Q = 0.05


@pytest.mark.unit
def test_c1_fires_when_the_floor_exceeds_the_rank_one_critical_value():
    """Pre-flight: p_res = 1/(B+1) above q/m means no single center can ever be rejected."""
    out = evaluate("preflight", B=200, q=Q, m=5000, p_emp=None, rejected=None,
                   s=S, ladder=LADDER, max_b=MAX_B)
    assert out.p_res == pytest.approx(1 / 201)
    assert out.R == pytest.approx((1 / 201) * 5000 / Q)
    assert out.conditions["C1"] is True
    assert out.limited is True
    assert out.conditions["C2"] is None and out.conditions["C3"] is None


@pytest.mark.unit
def test_c1_silent_when_b_is_large_relative_to_the_family():
    out = evaluate("preflight", B=10_000, q=Q, m=200, p_emp=None, rejected=None,
                   s=S, ladder=LADDER, max_b=MAX_B)
    assert out.R < 1.0
    assert out.conditions["C1"] is False
    assert out.limited is False


@pytest.mark.unit
def test_c0_fires_when_the_floor_exceeds_q_itself():
    out = evaluate("preflight", B=10, q=Q, m=5, p_emp=None, rejected=None,
                   s=S, ladder=LADDER, max_b=MAX_B)
    assert out.p_res == pytest.approx(1 / 11)
    assert out.conditions["C0"] is True and out.limited is True


@pytest.mark.unit
def test_c2_fires_when_too_few_centers_sit_at_the_floor():
    B, m = 200, 5000
    p = np.full(m, 0.5)
    p[:3] = 1.0 / (B + 1)                       # 3 centers at the floor, R = 124
    rejected = np.zeros(m, dtype=bool)
    out = evaluate("posthoc", B=B, q=Q, m=m, p_emp=p, rejected=rejected,
                   s=S, ladder=LADDER, max_b=MAX_B)
    assert out.n_floor == 3
    assert out.conditions["C2"] is True
    assert out.limited is True


@pytest.mark.unit
def test_c3_fires_when_a_floor_center_is_not_rejected():
    B, m = 10_000, 100
    p = np.full(m, 0.9)
    p[:5] = 1.0 / (B + 1)
    rejected = np.zeros(m, dtype=bool)          # nothing rejected despite floor p-values
    out = evaluate("posthoc", B=B, q=Q, m=m, p_emp=p, rejected=rejected,
                   s=S, ladder=LADDER, max_b=MAX_B)
    assert out.n_floor == 5
    assert out.conditions["C3"] is True
    assert out.limited is True


@pytest.mark.unit
def test_no_condition_fires_when_the_floor_centers_are_all_rejected():
    B, m = 10_000, 100
    p = np.full(m, 0.9)
    p[:20] = 1.0 / (B + 1)
    rejected = np.zeros(m, dtype=bool)
    rejected[:20] = True
    out = evaluate("posthoc", B=B, q=Q, m=m, p_emp=p, rejected=rejected,
                   s=S, ladder=LADDER, max_b=MAX_B)
    assert out.conditions == {"C0": False, "C1": False, "C2": False, "C3": False}
    assert out.limited is False


@pytest.mark.unit
def test_b_req_and_ladder_recommendation():
    b_req, b_rec, exceeds = b_recommendation(m=1000, q=Q, s=S, ladder=LADDER, max_b=MAX_B)
    assert b_req == int(np.ceil(S * 1000 / (K_TARGET * Q))) - 1 == 199_999
    assert b_rec == 1_000_000                    # smallest ladder value >= B_req
    assert exceeds is False


@pytest.mark.unit
def test_b_req_above_the_cap_recommends_nothing_and_escalates():
    b_req, b_rec, exceeds = b_recommendation(m=10_000_000, q=Q, s=S, ladder=LADDER,
                                             max_b=MAX_B)
    assert exceeds is True and b_rec is None
    out = evaluate("posthoc", B=10_000, q=Q, m=10_000_000,
                   p_emp=np.array([0.5]), rejected=np.array([False]),
                   s=S, ladder=LADDER, max_b=MAX_B)
    assert out.escalate_above_ladder is True and out.B_rec is None


@pytest.mark.unit
def test_diagnostic_never_changes_b_or_the_fdr_method():
    out = evaluate("posthoc", B=200, q=Q, m=5000, p_emp=np.full(5, 1 / 201),
                   rejected=np.zeros(5, dtype=bool), s=S, ladder=LADDER, max_b=MAX_B)
    payload = out.as_json()
    assert payload["B"] == 200                   # unchanged, despite B_rec being larger
    assert payload["B_was_changed"] is False
    assert payload["fdr_method_was_changed"] is False
    assert payload["B_rec"] != 200
    assert payload["k_target"] == 1
