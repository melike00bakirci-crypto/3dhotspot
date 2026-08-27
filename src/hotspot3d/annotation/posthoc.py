"""The post hoc mechanism-spatial gate and analysis — METHOD_SPEC II.13 / A24.

This is the **only** place in the entire pipeline where mechanism labels may enter a
statistical test, and it is terminal: nothing it produces feeds any upstream decision
(agent §3, §16).

The gate is deliberately unforgiving. It runs only when at least
``min_variants_per_category`` variants at ``min_evidence`` or better exist in each of
at least ``min_categories`` categories. Below that the test has no power, so it is
**not run** — recorded as a gate failure with its numbers, never as an omission. A
near-miss (9 where 10 is required) is a failure, not a pass; :func:`evaluate_gate`
has no tolerance parameter, because a tolerance is how a gate stops being one.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from ..utils.geometry import pairwise_distances
from .mechanisms import CuratedVariant, counts_at_or_above
from .rubric import ASSERTABLE, Evidence, parse_evidence

ANALYSIS_TYPE = "post_hoc"
SEED_CONTEXT = "posthoc|mechanism_spatial"


@dataclass(frozen=True)
class PosthocGate:
    """The frozen gate parameters, read from config. Absent config is BLOCKING."""

    min_variants_per_category: int
    min_evidence: Evidence
    min_categories: int
    separate_correction: bool
    feeds_upstream: bool

    @classmethod
    def from_config(cls, cfg) -> "PosthocGate":
        cfg.require_all(
            [
                "annotation.posthoc_gate.min_variants_per_category",
                "annotation.posthoc_gate.min_evidence",
                "annotation.posthoc_gate.min_categories",
                "annotation.posthoc_gate.separate_correction",
                "annotation.posthoc_gate.feeds_upstream",
            ]
        )
        gate = cls(
            min_variants_per_category=int(
                cfg.get("annotation.posthoc_gate.min_variants_per_category")),
            min_evidence=parse_evidence(cfg.get("annotation.posthoc_gate.min_evidence")),
            min_categories=int(cfg.get("annotation.posthoc_gate.min_categories")),
            separate_correction=bool(
                cfg.get("annotation.posthoc_gate.separate_correction")),
            feeds_upstream=bool(cfg.get("annotation.posthoc_gate.feeds_upstream")),
        )
        if gate.feeds_upstream:
            raise ValueError(
                "FROZEN methodology violation: annotation.posthoc_gate.feeds_upstream "
                "is TRUE. The post hoc analysis is terminal by construction."
            )
        return gate


@dataclass
class GateEvaluation:
    """The gate verdict *and* the numbers behind it — always recorded."""

    passed: bool
    counts: dict[str, int]
    qualifying_categories: list[str]
    n_qualifying_categories: int
    min_variants_per_category: int
    min_categories: int
    min_evidence: str
    n_eligible_variants: int
    reason: str

    def as_numbers(self) -> dict[str, Any]:
        return {
            "counts_by_category_at_or_above_min_evidence": self.counts,
            "qualifying_categories": self.qualifying_categories,
            "n_qualifying_categories": self.n_qualifying_categories,
            "n_eligible_variants": self.n_eligible_variants,
            "threshold_min_variants_per_category": self.min_variants_per_category,
            "threshold_min_categories": self.min_categories,
            "threshold_min_evidence": self.min_evidence,
            "passed": self.passed,
            "reason": self.reason,
        }


def evaluate_gate(curated: list[CuratedVariant], gate: PosthocGate) -> GateEvaluation:
    """Evaluate the gate. No tolerance, no near-miss allowance."""
    counts = counts_at_or_above(curated, gate.min_evidence, assertable_only=True)
    qualifying = sorted(
        name for name, n in counts.items() if n >= gate.min_variants_per_category
    )
    n_eligible = sum(
        1 for cv in curated
        if cv.mechanism in ASSERTABLE and cv.meets(gate.min_evidence)
    )
    passed = len(qualifying) >= gate.min_categories

    if passed:
        reason = (
            f"{len(qualifying)} categories ({', '.join(qualifying)}) each have "
            f">= {gate.min_variants_per_category} variants at "
            f">= {gate.min_evidence.value} evidence."
        )
    else:
        shortfall = ", ".join(
            f"{name}={n}" for name, n in sorted(counts.items())
        ) or "no assertable category"
        reason = (
            f"INSUFFICIENT EXPERIMENTALLY CHARACTERIZED VARIANTS — only "
            f"{len(qualifying)} of the required {gate.min_categories} categories reach "
            f"{gate.min_variants_per_category} variants at >= "
            f"{gate.min_evidence.value} evidence (counts: {shortfall}). The test has "
            f"no power at this sample size and was NOT run. This is a gate failure "
            f"with recorded numbers, not an omission, and a near-miss is not a pass."
        )

    return GateEvaluation(
        passed=passed,
        counts=counts,
        qualifying_categories=qualifying,
        n_qualifying_categories=len(qualifying),
        min_variants_per_category=gate.min_variants_per_category,
        min_categories=gate.min_categories,
        min_evidence=gate.min_evidence.value,
        n_eligible_variants=n_eligible,
        reason=reason,
    )


# --- the analysis (runs only on a passing gate) -----------------------------

def benjamini_hochberg(pvalues: list[float], q: float) -> tuple[list[float], list[bool]]:
    """BH step-up. Local to the post hoc family — corrected **separately** from Stage B.

    The separation is the requirement (A24): pooling these p-values with Stage B's
    would let a hypothesis-generating overlay alter the significance of the
    mechanism-blind discovery it is supposed to be independent of.
    """
    m = len(pvalues)
    if m == 0:
        return [], []
    order = sorted(range(m), key=lambda i: pvalues[i])
    adjusted = [1.0] * m
    running = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = m - rank + 1
        running = min(running, pvalues[idx] * m / i)
        adjusted[idx] = min(1.0, running)
    return adjusted, [adjusted[i] <= q for i in range(m)]


def run_posthoc_spatial(
    curated: list[CuratedVariant],
    coords: dict[int, tuple[float, float, float]],
    gate: PosthocGate,
    evaluation: GateEvaluation,
    *,
    rng: np.random.Generator,
    seed: int,
    n_permutations: int,
    q: float,
) -> dict[str, Any]:
    """Mechanism-specific spatial segregation test, labelled ``post_hoc``.

    Statistic: mean pairwise Cα distance among the variants of one mechanism category
    (spatial cohesion — small means the category is spatially concentrated). Null:
    mechanism labels are exchangeable *among the mechanism-labelled variants only*,
    so the test asks whether mechanisms segregate given where characterized variants
    are, not whether variants cluster at all — that question belongs to Stage B and
    was answered without any mechanism information.

    Hypothesis-generating only. Feeds nothing upstream.
    """
    if not evaluation.passed:
        raise ValueError(
            "run_posthoc_spatial called on a failed gate. The gate is not advisory."
        )

    eligible = [
        cv for cv in curated
        if cv.mechanism.value in evaluation.qualifying_categories
        and cv.meets(gate.min_evidence)
        and cv.residue_index in coords
    ]
    labels = np.array([cv.mechanism.value for cv in eligible])
    points = np.array([coords[cv.residue_index] for cv in eligible], dtype=float)
    categories = sorted(set(labels.tolist()))

    dist = pairwise_distances(points)
    observed = {c: _mean_intra_distance(dist, labels == c) for c in categories}

    n_extreme = {c: 0 for c in categories}
    valid = {c: 0 for c in categories}
    permuted = labels.copy()
    for _ in range(n_permutations):
        rng.shuffle(permuted)
        for c in categories:
            value = _mean_intra_distance(dist, permuted == c)
            if np.isnan(value):
                continue
            valid[c] += 1
            if value <= observed[c]:
                n_extreme[c] += 1

    names = [c for c in categories if not np.isnan(observed[c])]
    p_emp = [(1 + n_extreme[c]) / (1 + valid[c]) for c in names]
    q_bh, significant = benjamini_hochberg(p_emp, q)

    return {
        "analysis_type": ANALYSIS_TYPE,
        "hypothesis_generating": True,
        "feeds_upstream": False,
        "terminal": True,
        "statistic": "mean_intra_category_pairwise_ca_distance_angstrom",
        "alternative": "less (spatial concentration of a mechanism category)",
        "null_model": (
            "mechanism labels exchangeable among mechanism-labelled variants at "
            f">= {gate.min_evidence.value} evidence; positions held fixed"
        ),
        "n_permutations": n_permutations,
        "seed": seed,
        "seed_context": SEED_CONTEXT,
        "correction": {
            "method": "BH",
            "q": q,
            "family": "post_hoc_mechanism_spatial",
            "separate_from_stage_b": True,
        },
        "n_variants_tested": len(eligible),
        "categories": [
            {
                "mechanism": name,
                "n_variants": int((labels == name).sum()),
                "observed_mean_intra_distance_A": round(float(observed[name]), 6),
                "p_emp": round(float(p_emp[i]), 6),
                "q_bh": round(float(q_bh[i]), 6),
                "significant": bool(significant[i]),
            }
            for i, name in enumerate(names)
        ],
        "gate": evaluation.as_numbers(),
        "interpretation": (
            "Hypothesis-generating only. A significant result describes where "
            "experimentally characterized variants happen to sit; it is not evidence "
            "that any hotspot is real, and it did not and cannot influence hotspot "
            "discovery, r_hot, r_fp, footprint geometry or any robustness threshold."
        ),
    }


def _mean_intra_distance(dist: np.ndarray, mask: np.ndarray) -> float:
    """Mean pairwise distance within a group; NaN when fewer than two members."""
    idx = np.flatnonzero(mask)
    if idx.size < 2:
        return float("nan")
    sub = dist[np.ix_(idx, idx)]
    upper = sub[np.triu_indices(idx.size, k=1)]
    return float(upper.mean())


def gate_failure_record(evaluation: GateEvaluation, gate: PosthocGate) -> dict[str, Any]:
    """What is written when the gate fails: the numbers, not just a verdict."""
    return {
        "analysis_type": ANALYSIS_TYPE,
        "status": "NOT_RUN",
        "outcome_type": "SCIENTIFIC_NEGATIVE",
        "reason": evaluation.reason,
        "negative_result": {
            "condition": "INSUFFICIENT_CHARACTERIZED_VARIANTS_FOR_MECHANISM_ANALYSIS",
            "detail": (
                "Insufficient experimentally characterized variants for "
                "mechanism-specific spatial analysis. This is a valid, informative "
                "result — it says the literature does not yet support the test, not "
                "that the test failed."
            ),
        },
        "gate": evaluation.as_numbers(),
        "near_miss_is_not_a_pass": True,
        "feeds_upstream": False,
    }
