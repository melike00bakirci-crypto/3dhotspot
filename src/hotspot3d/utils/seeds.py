"""Reproducible seeding — METHOD_SPEC II.12 (FROZEN).

    seed(context) = int(sha256(f"{MASTER_SEED}|{run_id}|{context}")[:16], 16) mod 2**32
    generator     = numpy.random.Generator(PCG64(seed(context)))

A seed depends on its *context string* and never on execution order, so parallel
execution is bit-reproducible and any single iteration can be re-run in isolation.
Every derived seed is recorded in provenance.

Canonical context strings:
    global_null|PLP        global_null|BLB       scan|r=7.5
    final_detection        sensitivity|plddt70   sensitivity|star1
    iter|0123              bootstrap|jaccard     subset_design|k=3
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field

import numpy as np

_MOD = 2 ** 32


def derive_seed(master_seed: int, run_id: str, context: str) -> int:
    """FROZEN derivation. Do not alter — it defines reproducibility across runs."""
    digest = hashlib.sha256(f"{master_seed}|{run_id}|{context}".encode("utf-8")).hexdigest()
    return int(digest[:16], 16) % _MOD


@dataclass
class SeedRegistry:
    """Derives generators by context string and records every seed it issued."""

    master_seed: int
    run_id: str
    _issued: dict[str, int] = field(default_factory=dict)

    def seed(self, context: str) -> int:
        if context not in self._issued:
            self._issued[context] = derive_seed(self.master_seed, self.run_id, context)
        return self._issued[context]

    def rng(self, context: str) -> np.random.Generator:
        return np.random.Generator(np.random.PCG64(self.seed(context)))

    @property
    def issued(self) -> dict[str, int]:
        """Every (context -> seed) pair issued so far, for provenance."""
        return dict(sorted(self._issued.items()))

    def as_rows(self) -> list[dict]:
        """Rows for 12_REPRODUCIBILITY/random_seeds.tsv."""
        return [
            {"context": ctx, "seed": seed, "master_seed": self.master_seed,
             "run_id": self.run_id, "generator": "PCG64"}
            for ctx, seed in self.issued.items()
        ]


# --- canonical context-string builders (avoid free-form strings drifting) ----

def ctx_global_null(cohort: str) -> str:
    return f"global_null|{cohort}"


def ctx_scan(radius: float) -> str:
    return f"scan|r={radius:g}"


def ctx_final_detection() -> str:
    return "final_detection"


def ctx_sensitivity(name: str) -> str:
    return f"sensitivity|{name}"


def ctx_iteration(index: int) -> str:
    return f"iter|{index:04d}"


def ctx_bootstrap(metric: str) -> str:
    return f"bootstrap|{metric}"


def ctx_subset_level(k: int) -> str:
    return f"subset_design|k={k}"
