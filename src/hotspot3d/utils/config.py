"""Frozen configuration access.

Rules enforced here:
  * ``config/`` is Lead-owned; nothing in the pipeline writes it.
  * A missing parameter is BLOCKING, never defaulted (every agent §6).
  * The config is hashed once and the hash travels in every handoff.

**``sha256`` identifies the EFFECTIVE configuration.** It is the digest of the
canonical serialization of the merged parameter tree — base file plus any gene
overlay — not of the base file's bytes. That is what makes it usable as an
identity: two runs whose effective parameters are equal carry the same digest and
therefore the same ``run_id`` config component, and two runs that differ in any
parameter carry different ones, whichever file the difference came from.

The earlier composite form (``base[:32] + overlay[:32]``) failed that second
property in the place it mattered most: ``run_id`` takes ``config_sha256[0:8]``,
and those eight characters came entirely from the base file, so every overlay of
one base produced an identical run_id component. Different configurations were
not distinguishable by their run identity.

Byte-level provenance is not lost — ``source_sha256`` and ``overlay_sha256``
record the source files exactly as they were read, so a change confined to
comments (where this repository keeps its FROZEN annotations) is still detectable
even though it leaves the effective tree untouched.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from .errors import BlockedError
from .hashing import sha256_file, sha256_text

_MISSING = object()


def effective_config_yaml(data: dict) -> str:
    """The canonical text of an effective configuration.

    One serialization serves two purposes deliberately: it is what ``sha256``
    digests and it is what the run snapshots to
    ``12_REPRODUCIBILITY/config.yaml``. Because they are the same bytes, the
    snapshot is provably the configuration the digest identifies rather than a
    separately-produced description of it.

    ``sort_keys=True`` makes the text a function of the parameter tree alone, so
    the digest cannot move because two overlays merged their keys in a different
    order.
    """
    return yaml.safe_dump(data, sort_keys=True, default_flow_style=False,
                          allow_unicode=True)


@dataclass(frozen=True)
class FrozenConfig:
    """Read-only view over ``config/pipeline.yaml`` plus optional gene overlay."""

    data: dict
    sha256: str
    path: Path
    #: The gene overlay merged on top of ``path``, or None for a plain config.
    overlay: Path | None = None
    #: Digest of ``path``'s BYTES as read. Distinct from ``sha256``, which
    #: identifies the merged parameter tree; this one still moves when an edit
    #: touches only comments or formatting.
    source_sha256: str | None = None
    #: Digest of ``overlay``'s bytes as read, or None when there is no overlay.
    overlay_sha256: str | None = None

    # -- provenance ----------------------------------------------------------
    def effective_yaml(self) -> str:
        """The canonical text of the effective configuration this holds."""
        return effective_config_yaml(self.data)

    def recompute_sha256(self) -> str | None:
        """Re-read the sources, re-merge them, and re-derive ``sha256``.

        Returns None when a source is no longer readable, which is not by itself
        a freeze violation — the caller decides what an unavailable source means.
        Goes through ``load_config``, so the two can never drift apart.
        """
        if not Path(self.path).is_file():
            return None
        if self.overlay is not None and not Path(self.overlay).is_file():
            return None
        return load_config(self.path, self.overlay).sha256

    def recompute_source_digests(self) -> tuple[str | None, str | None]:
        """Re-hash the source files' bytes; None for a source that is gone."""
        base = sha256_file(self.path) if Path(self.path).is_file() else None
        if self.overlay is None:
            return base, None
        over = sha256_file(self.overlay) if Path(self.overlay).is_file() else None
        return base, over

    def describe_source(self) -> str:
        """Human-readable provenance for error messages."""
        if self.overlay is None:
            return str(self.path)
        return f"{self.path} + gene overlay {self.overlay}"

    # -- strict access -------------------------------------------------------
    def get(self, dotted: str) -> Any:
        """Fetch a parameter or raise BlockedError. There is no default."""
        node: Any = self.data
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                raise BlockedError(
                    f"Required config parameter '{dotted}' is missing or null in "
                    f"{self.path}. Missing parameters are BLOCKING and are never "
                    f"defaulted — escalate to the Lead."
                )
            node = node[part]
        if node is _MISSING:
            raise BlockedError(f"Required config parameter '{dotted}' is null.")
        return node

    def get_optional(self, dotted: str, default: Any = None) -> Any:
        """Only for parameters whose *absence is itself meaningful* (e.g. null filters)."""
        try:
            return self.get(dotted)
        except BlockedError:
            return default

    def require_all(self, dotted_keys: list[str]) -> None:
        """Validate a whole block up front so a stage fails fast, not mid-computation."""
        missing = []
        for key in dotted_keys:
            try:
                self.get(key)
            except BlockedError:
                missing.append(key)
        if missing:
            raise BlockedError(
                "BLOCKED — required config parameters absent: " + ", ".join(missing)
            )

    def assert_frozen_value(self, dotted: str, expected: Any) -> None:
        """Guard a FROZEN scientific constant against drift."""
        actual = self.get(dotted)
        if actual != expected:
            raise BlockedError(
                f"FROZEN methodology violation: config '{dotted}' is {actual!r} but the "
                f"approved Plan fixes it at {expected!r}. Changing it requires a "
                f"pre-registered methodological revision under a new RUN_ID."
            )


def load_config(path: str | Path = "config/pipeline.yaml",
                gene_overlay: str | Path | None = None) -> FrozenConfig:
    path = Path(path)
    if not path.is_file():
        raise BlockedError(f"BLOCKED — config not found: {path}")
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise BlockedError(f"BLOCKED — config did not parse to a mapping: {path}")
    source_digest = sha256_file(path)

    overlay_path: Path | None = None
    overlay_digest: str | None = None
    if gene_overlay is not None:
        overlay_path = Path(gene_overlay)
        if not overlay_path.is_file():
            raise BlockedError(f"BLOCKED — gene overlay not found: {overlay_path}")
        overlay = yaml.safe_load(overlay_path.read_text(encoding="utf-8")) or {}
        data = _deep_merge(data, overlay)
        overlay_digest = sha256_file(overlay_path)

    # The identity of a run's configuration is its EFFECTIVE parameter tree, so
    # that is what is digested. Hashing `path` here instead would give every
    # overlay of one base file the same identity — see the module docstring.
    return FrozenConfig(data=data, sha256=sha256_text(effective_config_yaml(data)),
                        path=path, overlay=overlay_path,
                        source_sha256=source_digest, overlay_sha256=overlay_digest)


def _deep_merge(base: dict, overlay: dict) -> dict:
    out = dict(base)
    for key, value in overlay.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = value
    return out


# --- FROZEN constants asserted at every run start ---------------------------
FROZEN_ASSERTIONS: list[tuple[str, Any]] = [
    ("representation.coordinate_atom", "CA"),
    ("representation.distance_metric", "euclidean"),
    ("clinvar.review_star_filter", None),
    ("clinvar.use_review_stars_for_inclusion", False),
    ("plddt.primary_filtering_enabled", False),
    ("plddt.sensitivity_threshold", 70.0),
    ("permutation.B_default", 100000),  # DECISION-B-DEFAULT-0001 (was 10000)
    ("fdr.method", "BH"),
    ("fdr.q", 0.05),
    ("loo_mcc.kappa", 2.0),
    ("loo_mcc.classification_threshold", None),
    ("loo_mcc.role", "radius_selection_only"),
    ("radius_domain.R_FLOOR", 5.0),
    ("radius_domain.R_CEIL", 25.0),
    ("radius_domain.step_hot_A", 0.5),
    ("radius_selection.w_k", 1.0),
    ("radius_selection.distance_metric", "L2"),
    ("radius_selection.normalization", "min_max"),
    # DECISION-FOOTPRINT-DOMAIN-0001: the r_fp domain is [5.0, min(0.25*D_max, 25.0)].
    ("footprint_domain.rho_min_A", 5.0),  # DECISION-FOOTPRINT-DOMAIN-0001 (was 3.0)
    ("footprint_domain.hard_ceiling_A", 25.0),  # DECISION-FOOTPRINT-DOMAIN-0001 (was 20.0)
    # The substantive half of that decision: rho_all measures centre packing, which is
    # not what r_fp is for. Pinned so it cannot be flipped back silently.
    ("footprint_domain.post_merge_factor_caps_rho_max", False),
    ("footprint_domain.step_fp_A", 0.5),
    ("footprint_domain.bound_below_by_r_hot", False),
    ("footprint_qc.QC_F1_max_coverage", 0.50),
    ("footprint_selection.w_k", 1.0),
    ("footprint_selection.distance_metric", "L2"),
    ("robustness.rerun_hotspot_pipeline", False),
    ("robustness.modify_clinvar_records", False),
    ("robustness.emit_categorical_verdict", False),
    ("robustness.preservation_threshold_primary", 0.50),
    ("seeding.MASTER_SEED", 20250101),
    ("annotation.may_influence_discovery", False),
    ("boundary_diagnostic.auto_widen_domain", False),
    ("radius_domain.automatic_domain_expansion", False),
]


def assert_frozen_methodology(cfg: FrozenConfig) -> None:
    """Fail the run if any FROZEN scientific constant has drifted.

    Called by the orchestrator before any stage executes, so drift cannot be
    introduced silently by editing config.
    """
    for dotted, expected in FROZEN_ASSERTIONS:
        actual = cfg.get_optional(dotted, _MISSING)
        if actual is _MISSING:
            raise BlockedError(f"FROZEN parameter '{dotted}' absent from config.")
        if actual != expected:
            raise BlockedError(
                f"FROZEN methodology violation: '{dotted}' is {actual!r}, "
                f"approved Plan fixes it at {expected!r}."
            )
