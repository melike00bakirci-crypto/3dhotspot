"""¶50 aggregation — n, mean, SD, 95% BCa CI, CV, median and IQR.

Two rules that decide what the numbers mean:

* **Failed iterations stay in every denominator** and contribute ``J = 0``. Dropping
  them would silently condition the profile on the iterations that happened to
  work, which is exactly the direction that flatters the result.
* **The CI method is recorded, never silently substituted.** BCa needs a
  non-degenerate bootstrap distribution and a non-degenerate jackknife; when either
  is unavailable the interval falls back to the percentile method and says so in
  the output, because an interval whose construction is unknown is not usable.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

AGGREGATE_COLUMNS = [
    "metric", "reconstruction", "n_total", "n_valid", "mean", "sd", "cv",
    "median", "q1", "q3", "iqr", "min", "max", "ci_lower", "ci_upper",
    "ci_method", "ci_level", "ci_resamples", "seed_context",
]


@dataclass(frozen=True)
class Interval:
    lower: float
    upper: float
    method: str
    detail: str


def bca_interval(values: np.ndarray, rng: np.random.Generator,
                 n_resamples: int = 10000, confidence: float = 0.95) -> Interval:
    """Bias-corrected and accelerated bootstrap CI for the mean.

    ``z0`` corrects for median bias in the bootstrap distribution; ``a`` corrects
    for a skewness-driven change in the standard error, estimated by the jackknife.
    """
    from scipy.stats import norm

    values = np.asarray(values, dtype=np.float64)
    n = values.size
    if n < 2:
        return Interval(float("nan"), float("nan"), "insufficient_n",
                        f"n = {n}; a bootstrap interval is not defined")
    if float(np.ptp(values)) == 0.0:
        value = float(values[0])
        return Interval(value, value, "degenerate_constant",
                        "every iteration produced the same value; the interval is "
                        "the point itself")

    theta_hat = float(values.mean())
    idx = rng.integers(0, n, size=(int(n_resamples), n))
    theta_star = values[idx].mean(axis=1)

    alpha = 1.0 - confidence
    z_lo, z_hi = norm.ppf(alpha / 2.0), norm.ppf(1.0 - alpha / 2.0)

    proportion = float(np.count_nonzero(theta_star < theta_hat)) / n_resamples
    if proportion <= 0.0 or proportion >= 1.0:
        lower, upper = np.percentile(theta_star, [100 * alpha / 2,
                                                  100 * (1 - alpha / 2)])
        return Interval(float(lower), float(upper), "percentile_fallback",
                        "z0 is undefined: the bootstrap distribution lies entirely "
                        "on one side of the point estimate")
    z0 = float(norm.ppf(proportion))

    total = values.sum()
    jack = (total - values) / (n - 1)
    diff = jack.mean() - jack
    denom = 6.0 * float((diff ** 2).sum()) ** 1.5
    if denom == 0.0:
        lower, upper = np.percentile(theta_star, [100 * alpha / 2,
                                                  100 * (1 - alpha / 2)])
        return Interval(float(lower), float(upper), "percentile_fallback",
                        "the jackknife acceleration is undefined (zero variance)")
    a = float((diff ** 3).sum()) / denom

    def adjust(z: float) -> float:
        bottom = 1.0 - a * (z0 + z)
        if bottom == 0.0:
            return float("nan")
        return float(norm.cdf(z0 + (z0 + z) / bottom))

    a1, a2 = adjust(z_lo), adjust(z_hi)
    if not (np.isfinite(a1) and np.isfinite(a2)) or not 0.0 < a1 < a2 < 1.0:
        lower, upper = np.percentile(theta_star, [100 * alpha / 2,
                                                  100 * (1 - alpha / 2)])
        return Interval(float(lower), float(upper), "percentile_fallback",
                        "the BCa endpoint adjustment left the unit interval")
    lower, upper = np.percentile(theta_star, [100 * a1, 100 * a2])
    return Interval(float(lower), float(upper), "BCa",
                    f"z0 = {z0:.4f}, a = {a:.4f}")


def summarize(metric: str, reconstruction: str, values: list,
              n_total: int, rng: np.random.Generator, n_resamples: int,
              confidence: float, seed_context: str) -> dict:
    """One aggregate row. ``n_total`` is the FULL planned iteration count."""
    numeric = np.asarray([v for v in values if v is not None
                          and isinstance(v, (int, float, np.floating, np.integer))
                          and np.isfinite(float(v))], dtype=np.float64)
    row = {
        "metric": metric, "reconstruction": reconstruction,
        "n_total": int(n_total), "n_valid": int(numeric.size),
        "ci_level": confidence, "ci_resamples": int(n_resamples),
        "seed_context": seed_context,
    }
    if numeric.size == 0:
        row.update({k: None for k in ("mean", "sd", "cv", "median", "q1", "q3",
                                      "iqr", "min", "max", "ci_lower", "ci_upper")})
        row["ci_method"] = "no_valid_values"
        return row

    mean = float(numeric.mean())
    sd = float(numeric.std(ddof=1)) if numeric.size > 1 else 0.0
    q1, median, q3 = (float(v) for v in np.percentile(numeric, [25, 50, 75]))
    interval = bca_interval(numeric, rng, n_resamples, confidence)
    row.update({
        "mean": mean, "sd": sd,
        "cv": float(sd / abs(mean)) if mean != 0 else None,
        "median": median, "q1": q1, "q3": q3, "iqr": float(q3 - q1),
        "min": float(numeric.min()), "max": float(numeric.max()),
        "ci_lower": interval.lower, "ci_upper": interval.upper,
        "ci_method": interval.method,
    })
    return row


def preservation_frequency(jaccards: list, threshold: float) -> tuple[float, int, int]:
    """``P(J >= threshold)`` over ALL iterations; failures count as ``J = 0``.

    Only the two pre-registered thresholds are ever used, and they exist solely to
    compute the ¶51 preservation proportion — not to license a verdict.
    """
    values = [0.0 if v is None else float(v) for v in jaccards]
    if not values:
        return float("nan"), 0, 0
    hits = sum(1 for v in values if v >= threshold)
    return float(hits / len(values)), hits, len(values)


def distribution_summary(values: list) -> dict:
    """Median and IQR for the robustness profile.

    ``iqr_lo`` / ``iqr_hi`` are the interquartile *bounds* (Q1 and Q3), named as
    the Output Contract names them because the Lead's REVIEW_PACK and
    gene_summary read these keys verbatim; ``q1`` / ``q3`` / ``iqr`` are kept
    alongside so a reader is never left guessing whether IQR means the interval
    or its width.
    """
    numeric = np.asarray([v for v in values if v is not None
                          and np.isfinite(float(v))], dtype=np.float64)
    if numeric.size == 0:
        return {"n": 0, "mean": None, "sd": None, "median": None,
                "q1": None, "q3": None, "iqr_lo": None, "iqr_hi": None,
                "iqr": None, "min": None, "max": None}
    q1, median, q3 = (float(v) for v in np.percentile(numeric, [25, 50, 75]))
    return {
        "n": int(numeric.size), "mean": float(numeric.mean()),
        "sd": float(numeric.std(ddof=1)) if numeric.size > 1 else 0.0,
        "median": median, "q1": q1, "q3": q3,
        "iqr_lo": q1, "iqr_hi": q3, "iqr": float(q3 - q1),
        "min": float(numeric.min()), "max": float(numeric.max()),
    }
