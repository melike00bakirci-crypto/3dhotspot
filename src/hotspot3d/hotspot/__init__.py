"""Stage B — primary spatial hotspot discovery (METHOD_SPEC II.3-II.7, Workflow v2).

A frozen, versioned, deterministic library. Other stages may **import** it; none may
edit it, and the robustness procedure of Stage D does not import it at all (F11: the
hotspot pipeline is never re-run inside robustness).

The primary per-center test is the **structure-aware positional null**
(:mod:`hotspot3d.hotspot.positional`, v2 §5.2). The label-permutation null
(:mod:`hotspot3d.hotspot.permutation`) is a SECONDARY reported analysis of a
different hypothesis. No run may be reported as a negative without a passing
:mod:`hotspot3d.hotspot.power` certificate (v2 §5.4).

Entry point::

    from hotspot3d.hotspot import run_stage_b
    handoff_02 = run_stage_b(ctx, upstream=handoff_01)
"""
from .detection import DetectionResult, final_detection, verify_object_nesting
from .fdr import FDRResult, bh_adjust, benjamini_hochberg, by_constant
from .loo import LooOutcome, kappa_is_boundary_neutral, leave_one_out
from .permutation import PermutationOutcome, permutation_pass, permute_labels
from .positional import PositionalOutcome, hypergeom_sf, positional_pass
from .power import PowerCertificate, UnderpoweredResult, certify
from .power import preflight as power_preflight
from .resolution import ResolutionDiagnostic, b_recommendation
from .resolution import evaluate as resolution_diagnostic
from .scan import RadiusResult, attach_neighbor_stability, scan_one_radius
from .selection import RadiusSelection, run_selection
from .sensitivity import SensitivityOutcome, run_restricted_detection
from .stage import is_underpowered, run_stage_b

__all__ = [
    "run_stage_b", "is_underpowered",
    "final_detection", "DetectionResult", "verify_object_nesting",
    "benjamini_hochberg", "bh_adjust", "by_constant", "FDRResult",
    "leave_one_out", "LooOutcome", "kappa_is_boundary_neutral",
    "permutation_pass", "permute_labels", "PermutationOutcome",
    "positional_pass", "PositionalOutcome", "hypergeom_sf",
    "certify", "power_preflight", "PowerCertificate", "UnderpoweredResult",
    "resolution_diagnostic", "ResolutionDiagnostic", "b_recommendation",
    "scan_one_radius", "attach_neighbor_stability", "RadiusResult",
    "run_selection", "RadiusSelection",
    "run_restricted_detection", "SensitivityOutcome",
]
