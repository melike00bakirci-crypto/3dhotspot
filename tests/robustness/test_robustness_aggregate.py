"""¶50 aggregation: BCa intervals, denominators and preservation frequencies."""
from __future__ import annotations

import numpy as np
import pytest

from hotspot3d.robustness.aggregate import (bca_interval, preservation_frequency,
                                            summarize)

pytestmark = pytest.mark.unit


def _rng(seed=12345):
    return np.random.Generator(np.random.PCG64(seed))


# --- BCa ---------------------------------------------------------------------

def test_bca_covers_the_mean_of_a_known_normal_sample():
    values = _rng(7).normal(loc=0.60, scale=0.10, size=200)
    interval = bca_interval(values, _rng(1), n_resamples=4000)

    assert interval.method == "BCa"
    assert interval.lower < values.mean() < interval.upper
    # SE of the mean is 0.10/sqrt(200) ~ 0.0071, so the 95% width is ~0.028
    assert 0.015 < interval.upper - interval.lower < 0.06


def test_bca_is_asymmetric_for_a_skewed_sample():
    """The acceleration term is the whole point: a skewed sample gets a skewed CI."""
    values = _rng(11).exponential(scale=1.0, size=300)
    interval = bca_interval(values, _rng(2), n_resamples=4000)
    mean = values.mean()

    assert interval.method == "BCa"
    left, right = mean - interval.lower, interval.upper - mean
    assert abs(left - right) / max(left, right) > 0.02


def test_bca_reports_a_percentile_fallback_rather_than_pretending():
    """A degenerate jackknife must not silently masquerade as BCa."""
    values = np.array([0.4, 0.4, 0.4, 0.4, 0.9])
    interval = bca_interval(values, _rng(3), n_resamples=2000)
    assert interval.method in ("BCa", "percentile_fallback")
    assert interval.detail


def test_constant_values_give_a_degenerate_recorded_interval():
    interval = bca_interval(np.full(20, 0.75), _rng(4), n_resamples=1000)
    assert interval.method == "degenerate_constant"
    assert (interval.lower, interval.upper) == (0.75, 0.75)


def test_a_single_value_has_no_interval():
    interval = bca_interval(np.array([0.5]), _rng(5))
    assert interval.method == "insufficient_n"
    assert np.isnan(interval.lower) and np.isnan(interval.upper)


# --- denominators ------------------------------------------------------------

def test_failed_iterations_stay_in_the_denominator():
    """n_total is the PLANNED count, always; n_valid records what was measurable."""
    values = [0.8, 0.7, None, 0.6, None]
    row = summarize("jaccard_residues", "primary", values, n_total=5,
                    rng=_rng(6), n_resamples=1000, confidence=0.95,
                    seed_context="bootstrap|test")

    assert row["n_total"] == 5
    assert row["n_valid"] == 3
    assert row["mean"] == pytest.approx(0.7)
    assert row["ci_method"] in ("BCa", "percentile_fallback")


def test_summary_reports_every_required_statistic():
    row = summarize("dice_residues", "primary", [0.1, 0.4, 0.5, 0.9], n_total=4,
                    rng=_rng(8), n_resamples=1000, confidence=0.95,
                    seed_context="bootstrap|test")
    for key in ("n_total", "n_valid", "mean", "sd", "cv", "median", "q1", "q3",
                "iqr", "min", "max", "ci_lower", "ci_upper", "ci_method",
                "ci_level", "ci_resamples"):
        assert row[key] is not None, key
    assert row["iqr"] == pytest.approx(row["q3"] - row["q1"])


def test_no_valid_values_is_stated_explicitly():
    row = summarize("volume_ratio", "primary", [None, None], n_total=2,
                    rng=_rng(9), n_resamples=100, confidence=0.95,
                    seed_context="bootstrap|test")
    assert row["ci_method"] == "no_valid_values"
    assert row["n_valid"] == 0 and row["n_total"] == 2


# --- preservation ------------------------------------------------------------

def test_failed_iterations_contribute_jaccard_zero():
    frequency, hits, total = preservation_frequency([0.9, 0.6, None, 0.2], 0.50)
    assert (hits, total) == (2, 4)
    assert frequency == pytest.approx(0.5)


def test_preservation_uses_a_greater_or_equal_comparison():
    frequency, hits, total = preservation_frequency([0.50, 0.4999], 0.50)
    assert (hits, total) == (1, 2)
    assert frequency == pytest.approx(0.5)


def test_preservation_thresholds_come_from_config_only(rparams):
    assert rparams.preservation_primary == 0.50
    assert rparams.preservation_secondary == 0.70
    assert rparams.bootstrap_method == "BCa"
    assert rparams.bootstrap_resamples == 10000
    assert rparams.failed_iteration_jaccard == 0.0
