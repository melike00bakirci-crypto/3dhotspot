"""II.11 subset design — fixed, recorded, and never re-drawn (agent §3.8).

The design is decided by combinatorial arithmetic against ``N_CAP`` alone, written
out **in full before iteration 1**, and never revisited because early iterations
looked unstable. That is the difference between a pre-registered sampling plan and
design shopping.

Two properties are load-bearing:

* **Exactness.** Subsets are produced by combinatorial *unranking* (combinadic),
  not by drawing random subsets and rejecting duplicates. Distinct ranks give
  distinct subsets by construction, so a level is duplicate-free for free.
* **Order independence.** Every level's generator is seeded from its context
  string ``subset_design|k=<k>``, never from execution order, so the design is
  identical whether levels run sequentially, in parallel, or one at a time during
  a re-run of a single iteration.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import comb

import numpy as np

from ..utils.errors import BlockedError
from ..utils.seeds import SeedRegistry, ctx_subset_level


@dataclass(frozen=True)
class Level:
    k: int
    total_subsets: int              # T(k) = C(n_S, k), exact Python int
    n_evaluated: int
    mode: str                       # exhaustive | sampled
    sampling_fraction: float
    seed_context: str
    seed: int | None

    def as_row(self) -> dict:
        return {
            "k": self.k,
            "total_subsets_T_k": self.total_subsets,
            "n_evaluated": self.n_evaluated,
            "mode": self.mode,
            "sampling_fraction": self.sampling_fraction,
            "seed_context": self.seed_context,
            "seed": self.seed,
        }


@dataclass(frozen=True)
class SubsetDesign:
    n_S: int
    K_MAX: int
    N_CAP: int
    levels: tuple[Level, ...]
    subsets: tuple[tuple[int, ...], ...]      # removal sets, in iteration order

    @property
    def total_subset_space(self) -> int:
        return sum(level.total_subsets for level in self.levels)

    @property
    def total_evaluated(self) -> int:
        return sum(level.n_evaluated for level in self.levels)

    def as_dict(self) -> dict:
        return {
            "n_S": self.n_S,
            "K_MAX": self.K_MAX,
            "K_MAX_rule": "min(n_S - 1, floor(n_S / 2))",
            "N_CAP": self.N_CAP,
            "total_subset_space_sum_T_k": self.total_subset_space,
            "total_evaluated": self.total_evaluated,
            "budget_allocation": "equal_across_remaining_levels",
            "sampling": "seeded_combinatorial_unranking (without replacement)",
            "recorded_before_iteration_1": True,
        }


def k_max(n_s: int) -> int:
    """``K_MAX = min(n_S - 1, floor(n_S / 2))`` (A19, frozen)."""
    return min(n_s - 1, n_s // 2)


def build_design(n_s: int, n_cap: int, seeds: SeedRegistry,
                 center_ids: tuple[int, ...]) -> SubsetDesign:
    """Produce the complete design. Called once, before any iteration runs."""
    if n_s < 2:
        raise BlockedError(
            f"subset design requested for n_S = {n_s}; perturbation needs n_S >= 2. "
            f"The caller must emit ROBUSTNESS_NOT_EVALUABLE instead."
        )
    kmax = k_max(n_s)
    totals = {k: comb(n_s, k) for k in range(1, kmax + 1)}
    allocation, modes = _allocate(totals, kmax, n_cap)

    levels: list[Level] = []
    subsets: list[tuple[int, ...]] = []
    for k in range(1, kmax + 1):
        context = ctx_subset_level(k)
        n_eval = allocation[k]
        total = totals[k]
        seed = None
        if modes[k] == "exhaustive":
            ranks = range(total)
        elif n_eval == 0:
            ranks = ()
        else:
            seed = seeds.seed(context)
            ranks = sorted(_sample_ranks(seeds.rng(context), total, n_eval))
        for rank in ranks:
            subsets.append(tuple(int(center_ids[i])
                                 for i in unrank_subset(int(rank), n_s, k)))
        levels.append(Level(
            k=k, total_subsets=total, n_evaluated=n_eval, mode=modes[k],
            sampling_fraction=float(n_eval / total) if total else 0.0,
            seed_context=context, seed=seed,
        ))

    if len(subsets) != sum(level.n_evaluated for level in levels):
        raise BlockedError("subset design accounting mismatch")
    if len(set(subsets)) != len(subsets):
        raise BlockedError("subset design produced duplicate removal sets")
    return SubsetDesign(n_S=n_s, K_MAX=kmax, N_CAP=n_cap,
                        levels=tuple(levels), subsets=tuple(subsets))


def _allocate(totals: dict[int, int], kmax: int,
              n_cap: int) -> tuple[dict[int, int], dict[int, str]]:
    """Ascending exhaustion, then equal split of the remainder with redistribution."""
    allocation: dict[int, int] = {}
    modes: dict[int, str] = {}

    cumulative = 0
    exhaustive_upto = 0
    for k in range(1, kmax + 1):
        if cumulative + totals[k] <= n_cap:
            cumulative += totals[k]
            exhaustive_upto = k
        else:
            break
    for k in range(1, exhaustive_upto + 1):
        allocation[k] = totals[k]
        modes[k] = "exhaustive"

    remaining = [k for k in range(exhaustive_upto + 1, kmax + 1)]
    budget = n_cap - cumulative
    while remaining and budget > 0:
        base, extra = divmod(budget, len(remaining))
        tentative = {k: base + (1 if i < extra else 0)
                     for i, k in enumerate(remaining)}
        saturated = [k for k in remaining if tentative[k] >= totals[k]]
        if not saturated:
            allocation.update(tentative)
            for k in remaining:
                modes[k] = "sampled"
            budget, remaining = 0, []
            break
        # a saturated level becomes exhaustive and its surplus is redistributed
        for k in saturated:
            allocation[k] = totals[k]
            modes[k] = "exhaustive"
            budget -= totals[k]
        remaining = [k for k in remaining if k not in saturated]
    for k in remaining:
        allocation[k] = 0
        modes[k] = "sampled"
    for k in range(1, kmax + 1):
        allocation.setdefault(k, 0)
        modes.setdefault(k, "sampled")
    return allocation, modes


def unrank_subset(rank: int, n: int, k: int) -> tuple[int, ...]:
    """Lexicographic combinadic unranking of the ``rank``-th ``k``-subset of ``0..n-1``.

    Exact for arbitrarily large ``C(n, k)`` because Python integers are unbounded;
    a floating-point or int64 formulation would silently collide for large ``n_S``.
    """
    if not 0 <= rank < comb(n, k):
        raise ValueError(f"rank {rank} out of range for C({n}, {k})")
    out: list[int] = []
    x, b, start = rank, k, 0
    for _ in range(k):
        for element in range(start, n):
            count = comb(n - element - 1, b - 1)
            if x < count:
                out.append(element)
                start, b = element + 1, b - 1
                break
            x -= count
    return tuple(out)


def rank_subset(subset: tuple[int, ...], n: int) -> int:
    """Inverse of :func:`unrank_subset` — used only to verify exactness in tests."""
    k = len(subset)
    rank, previous, b = 0, -1, k
    for element in subset:
        for skipped in range(previous + 1, element):
            rank += comb(n - skipped - 1, b - 1)
        previous, b = element, b - 1
    return rank


def _sample_ranks(rng: np.random.Generator, total: int, n: int) -> list[int]:
    """``n`` DISTINCT uniform ranks in ``[0, total)``, without replacement.

    Uniform integers are drawn by rejection on the bit length rather than by a
    modulo reduction, which would be biased, and via Python integers rather than
    int64, which would overflow for large ``C(n_S, k)``.
    """
    if n > total:
        raise BlockedError(f"cannot draw {n} distinct ranks from {total}")
    if n == total:
        return list(range(total))
    bits = max(1, int(total - 1).bit_length())
    words = (bits + 31) // 32
    chosen: set[int] = set()
    guard = 0
    while len(chosen) < n:
        guard += 1
        if guard > 1000 * n + 10000:            # pragma: no cover - safety net
            raise BlockedError("rank sampling failed to converge")
        raw = rng.integers(0, 2 ** 32, size=words, dtype=np.uint64)
        value = 0
        for word in raw:
            value = (value << 32) | int(word)
        value &= (1 << bits) - 1
        if value < total:
            chosen.add(value)
    return sorted(chosen)


def estimate_wall_seconds(n_iterations: int, seconds_per_iteration: float) -> float:
    return float(n_iterations) * float(seconds_per_iteration)
