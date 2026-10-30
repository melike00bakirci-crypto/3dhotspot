"""Frozen Phase D parameters, read once from ``config/pipeline.yaml``.

The preservation thresholds, ``N_CAP``, the bootstrap settings and ``MASTER_SEED``
are all pre-registered constants. In particular the thresholds are read here and
nowhere else, so they cannot be chosen after the distribution has been seen.
"""
from __future__ import annotations

from dataclasses import dataclass

from ..utils.config import FrozenConfig
from ..utils.errors import BlockedError

REQUIRED_CONFIG_KEYS: list[str] = [
    "robustness.target",
    "robustness.perturbation_universe",
    "robustness.rerun_hotspot_pipeline",
    "robustness.modify_clinvar_records",
    "robustness.K_MAX_rule",
    "robustness.N_CAP",
    "robustness.sampling",
    "robustness.budget_allocation",
    "robustness.secondary_fixed_radius_reconstruction",
    "robustness.preservation_threshold_primary",
    "robustness.preservation_threshold_secondary",
    "robustness.bootstrap.method",
    "robustness.bootstrap.resamples",
    "robustness.bootstrap.confidence",
    "robustness.emit_categorical_verdict",
    "robustness.failed_iterations_jaccard",
    "robustness.min_n_S_for_evaluation",
    "loo_mcc.undefined_mcc_value",
    "seeding.MASTER_SEED",
    "execution.n_jobs",
    "execution.worker_nice",
]


@dataclass(frozen=True)
class RobustnessParams:
    n_cap: int
    preservation_primary: float
    preservation_secondary: float
    bootstrap_method: str
    bootstrap_resamples: int
    bootstrap_confidence: float
    failed_iteration_jaccard: float
    min_n_s: int
    undefined_mcc: float
    secondary_fixed_radius: bool
    max_wall_seconds: float | None
    n_workers: int
    worker_nice: int

    @classmethod
    def from_config(cls, cfg: FrozenConfig) -> "RobustnessParams":
        cfg.require_all(REQUIRED_CONFIG_KEYS)

        if cfg.get("robustness.target") != "final_footprint":
            raise BlockedError(
                "FROZEN (F11): the robustness target is the FINAL FOOTPRINT. This "
                "stage is geometric footprint robustness, never hotspot LOO."
            )
        if cfg.get("robustness.perturbation_universe") != "SIGNIFICANT_HOTSPOT_CENTERS":
            raise BlockedError(
                "FROZEN (A21): the perturbation universe is SIGNIFICANT_HOTSPOT_CENTERS "
                "— geometric center positions, never ClinVar records."
            )
        if not bool(cfg.get("robustness.secondary_fixed_radius_reconstruction")):
            raise BlockedError(
                "DECISION-STAGE-D-FIXED-RFP-0001 made fixed-radius reconstruction the "
                "SOLE per-iteration reconstruction: the re-deriving 'primary' path no "
                "longer exists, so there is nothing for this key to switch off. Setting "
                "it false would silently change nothing. If the intent is to restore "
                "per-iteration radius re-derivation, that needs a new decision record "
                "reversing DECISION-STAGE-D-FIXED-RFP-0001, not a config flag."
            )
        if bool(cfg.get("robustness.rerun_hotspot_pipeline")):
            raise BlockedError(
                "FROZEN (F11/A7): full-pipeline re-execution is a WITHDRAWN design, "
                "not an option. It must not return under any name."
            )
        if bool(cfg.get("robustness.modify_clinvar_records")):
            raise BlockedError(
                "FROZEN (A21): no ClinVar record is removed and no class label is "
                "altered at any point."
            )
        if bool(cfg.get("robustness.emit_categorical_verdict")):
            raise BlockedError(
                "FROZEN (II.11): no ROBUST / NOT_ROBUST verdict is emitted as a "
                "primary output. The scientific result is the continuous profile."
            )
        if cfg.get("robustness.bootstrap.method") != "BCa":
            raise BlockedError("FROZEN (A20): the bootstrap method is BCa.")
        if cfg.get("robustness.sampling") != "seeded_combinatorial_unranking":
            raise BlockedError("FROZEN (A19): subsets are drawn by seeded "
                               "combinatorial unranking, without replacement.")
        if cfg.get("robustness.budget_allocation") != "equal_across_remaining_levels":
            raise BlockedError("FROZEN (A19): the remaining budget is allocated "
                               "equally across every remaining level.")
        if float(cfg.get("robustness.failed_iterations_jaccard")) != 0.0:
            raise BlockedError("FROZEN (II.11): failed iterations contribute J = 0 "
                               "and stay in every denominator.")

        n_workers = int(cfg.get("execution.n_jobs"))
        if n_workers < 1:
            raise BlockedError(
                f"execution.n_jobs = {n_workers}; a worker count below 1 is not "
                f"meaningful. This is a machine-policy knob, never auto-scaled to "
                f"the node's CPU count — set it explicitly."
            )
        worker_nice = int(cfg.get("execution.worker_nice"))
        if worker_nice <= 0:
            raise BlockedError(
                f"execution.worker_nice = {worker_nice}; workers must run under a "
                f"POSITIVE nice value so they yield to other work on this shared "
                f"node. A value <= 0 would not yield and is rejected."
            )

        return cls(
            n_cap=int(cfg.get("robustness.N_CAP")),
            preservation_primary=float(
                cfg.get("robustness.preservation_threshold_primary")),
            preservation_secondary=float(
                cfg.get("robustness.preservation_threshold_secondary")),
            bootstrap_method=str(cfg.get("robustness.bootstrap.method")),
            bootstrap_resamples=int(cfg.get("robustness.bootstrap.resamples")),
            bootstrap_confidence=float(cfg.get("robustness.bootstrap.confidence")),
            failed_iteration_jaccard=float(
                cfg.get("robustness.failed_iterations_jaccard")),
            min_n_s=int(cfg.get("robustness.min_n_S_for_evaluation")),
            undefined_mcc=float(cfg.get("loo_mcc.undefined_mcc_value")),
            secondary_fixed_radius=bool(
                cfg.get("robustness.secondary_fixed_radius_reconstruction")),
            max_wall_seconds=_optional_float(cfg, "robustness.max_wall_seconds"),
            n_workers=n_workers,
            worker_nice=worker_nice,
        )

    def as_dict(self) -> dict:
        return {
            "N_CAP": self.n_cap,
            "preservation_threshold_primary": self.preservation_primary,
            "preservation_threshold_secondary": self.preservation_secondary,
            "bootstrap_method": self.bootstrap_method,
            "bootstrap_resamples": self.bootstrap_resamples,
            "bootstrap_confidence": self.bootstrap_confidence,
            "failed_iterations_jaccard": self.failed_iteration_jaccard,
            "min_n_S_for_evaluation": self.min_n_s,
            "undefined_mcc_value": self.undefined_mcc,
            "secondary_fixed_radius_reconstruction": self.secondary_fixed_radius,
            "max_wall_seconds": self.max_wall_seconds,
            "n_workers": self.n_workers,
            "worker_nice": self.worker_nice,
        }


def _optional_float(cfg: FrozenConfig, key: str) -> float | None:
    """Absence is meaningful here: no Lead-declared wall-time budget is set."""
    value = cfg.get_optional(key, None)
    return None if value is None else float(value)
