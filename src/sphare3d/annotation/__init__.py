"""Stage E — biological annotation (strictly downstream, terminal).

Public surface::

    run_stage_e(ctx, *, handoff_02, handoff_03, handoff_04, source=None) -> Handoff

Everything this package produces is an **overlay** on regions that were discovered
and validated without any biological input. Nothing here may reach hotspot discovery,
``r_hot``, ``r_fp``, footprint geometry or any robustness threshold.
"""
from .mechanisms import CuratedVariant, curate, mechanism_by_hotspot, mechanism_counts
from .posthoc import GateEvaluation, PosthocGate, evaluate_gate, run_posthoc_spatial
from .rubric import (
    Evidence,
    Mechanism,
    MechanismRecord,
    Resolution,
    Rubric,
    resolve_mechanism,
)
from .schema import assert_mandatory_context, assert_no_upstream_influence, assert_write_target
from .sources import AnnotationSource, NetworkAnnotationSource, SourceVersion, validate_source
from .stage import run_stage_e
from .verbatim import VerbatimAudit, load_verbatim_table, verbatim_copy

__all__ = [
    "run_stage_e",
    # rubric
    "Evidence", "Mechanism", "MechanismRecord", "Resolution", "Rubric",
    "resolve_mechanism",
    # curation
    "CuratedVariant", "curate", "mechanism_by_hotspot", "mechanism_counts",
    # post hoc
    "GateEvaluation", "PosthocGate", "evaluate_gate", "run_posthoc_spatial",
    # sources
    "AnnotationSource", "NetworkAnnotationSource", "SourceVersion", "validate_source",
    # guardrails
    "assert_mandatory_context", "assert_no_upstream_influence", "assert_write_target",
    "VerbatimAudit", "load_verbatim_table", "verbatim_copy",
]
