"""Lead integration audit — Phase 6.

Structural checks the Lead owns. These do not test the science (each agent tests
its own); they test that the four implementations *interlock* and that the frozen
methodology cannot drift into the code.
"""
from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[2] / "src" / "hotspot3d"
AGENT_PKGS = {
    "data-structure": ["data", "structure"],
    "hotspot-statistics": ["spatial", "hotspot"],
    "footprint-robustness": ["footprint", "robustness"],
    "biological-annotation": ["annotation"],
}
LEAD_PKGS = ["utils", "orchestration", "reporting"]


def _py_files(*packages: str) -> list[Path]:
    out: list[Path] = []
    for pkg in packages:
        out.extend(p for p in (SRC / pkg).rglob("*.py")
                   if "__pycache__" not in p.parts)
    return sorted(out)


def _all_agent_files() -> list[Path]:
    return _py_files(*[p for pkgs in AGENT_PKGS.values() for p in pkgs])


# --------------------------------------------------------------------------- #
# 1. no methodology drift — frozen constants live in config, not in code       #
# --------------------------------------------------------------------------- #

FROZEN_LITERALS = {
    # literal -> the config key that must be its only source
    "10000": "permutation.B_default / robustness.N_CAP / bootstrap.resamples",
    "0.05": "fdr.q (or a near-tie threshold, also config-driven)",
}
# Numeric literals that would indicate a hard-coded scientific constant.
SUSPICIOUS = [
    (r"\bq\s*=\s*0\.05\b", "FDR q hard-coded; must come from config fdr.q"),
    (r"\bB\s*=\s*10000\b", "permutation count hard-coded; config permutation.B_default"),
    (r"\bkappa\s*=\s*2(\.0)?\b", "LOO kappa hard-coded; config loo_mcc.kappa"),
    (r"coverage\s*>\s*0\.5(?!\d)", "coverage gate hard-coded; config *_max_coverage"),
    (r"\br_hot\s*=\s*\d", "r_hot assigned a literal radius"),
    (r"\br_fp\s*=\s*\d", "r_fp assigned a literal radius"),
]


def _code_lines(path: Path) -> list[tuple[int, str]]:
    """Executable source lines only — every string literal and comment removed.

    Documenting a frozen formula in a docstring, a report sentence or a comment is
    required by the methodology; only an *executable* hard-coded constant is a
    defect. Tokenizing is the reliable way to tell those apart: a regex cannot.
    """
    import io
    import tokenize

    text = path.read_text(encoding="utf-8")
    try:
        tokens = list(tokenize.generate_tokens(io.StringIO(text).readline))
    except (tokenize.TokenError, SyntaxError, IndentationError):
        return []

    lines: dict[int, list[str]] = {}
    for tok in tokens:
        if tok.type in (tokenize.STRING, tokenize.COMMENT, tokenize.NL,
                        tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT):
            continue
        if hasattr(tokenize, "FSTRING_START") and tok.type in (
                tokenize.FSTRING_START, tokenize.FSTRING_MIDDLE,
                tokenize.FSTRING_END):
            continue
        lines.setdefault(tok.start[0], []).append(tok.string)
    return [(no, " ".join(parts)) for no, parts in sorted(lines.items()) if parts]


def test_no_hardcoded_scientific_constants():
    """Configuration drives parameters rather than hidden constants."""
    offenders = []
    for path in _all_agent_files():
        for lineno, line in _code_lines(path):
            for pattern, why in SUSPICIOUS:
                if re.search(pattern, line) and not any(
                        tok in line for tok in ("config", "cfg", "default", "get(")):
                    offenders.append(
                        f"{path.relative_to(SRC)}:{lineno} {why} :: {line.strip()[:90]}")
    assert not offenders, "Hard-coded scientific constants found:\n" + "\n".join(offenders)


def test_frozen_methodology_assertions_hold(config):
    from hotspot3d.utils.config import assert_frozen_methodology
    assert_frozen_methodology(config)


def test_production_config_ships_the_full_robustness_budget(repo_root):
    """The SHIPPED config carries the production perturbation budget.

    The integration ``config`` fixture overlays a reduced ``robustness.N_CAP`` so
    the end-to-end chain is executable in a suite (see ``tests/integration/
    conftest.py``). That reduction is a test-execution choice and must never
    reach the file a real run reads, so the shipped value is asserted here
    against ``config/pipeline.yaml`` directly rather than through the fixture —
    the same guard ``test_production_default_permutation_count``
    provides for ``permutation.B_default``.

    The frozen assertions are re-checked on the shipped file too: the fixture
    version above proves the overlay does not violate them, which is a different
    claim from the shipped configuration being valid.
    """
    from hotspot3d.utils.config import assert_frozen_methodology, load_config

    shipped = load_config(repo_root / "config" / "pipeline.yaml")
    assert shipped.get("robustness.N_CAP") == 10000, (
        "config/pipeline.yaml no longer ships the production perturbation budget; "
        "a test-scoped reduction has leaked into the production configuration.")
    assert_frozen_methodology(shipped)


def test_production_config_freezes_the_compound_significance_ruling(repo_root):
    """The compound-label reading is frozen CONSERVATIVELY before the first real gene.

    It was the one parameter still marked PROVISIONAL, and it materially changes
    cohort membership on real data. The ruling is ``whole_field``: the primary
    P/LP-vs-B/LB comparison admits only records classifiable unambiguously, so a
    label carrying non-pathogenicity semantics never enters a binary class on the
    strength of one component.

    Asserted against the shipped file and against the behaviour it produces, since
    the value alone would not prove the classifier honours it.
    """
    from hotspot3d.data.classify import (REASON_NON_BINARY, SIG_BLB, SIG_PLP,
                                         SignificancePolicy)
    from hotspot3d.utils.config import load_config

    shipped = load_config(repo_root / "config" / "pipeline.yaml")
    assert shipped.get("clinvar.compound_significance_policy") == "whole_field", (
        "the conservative compound-significance ruling has been reverted")

    policy = SignificancePolicy.from_config(shipped)
    assert policy.component_matching is False

    # Unambiguous labels — bare and aggregate — still enter the primary cohort.
    # '/' is part of an aggregate label and is never a component separator.
    for label in ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"):
        assert policy.classify(label).significance_class == SIG_PLP, label
    for label in ("Benign", "Likely benign", "Benign/Likely benign"):
        assert policy.classify(label).significance_class == SIG_BLB, label

    # Compound labels carrying other semantics are excluded — but preserved, with
    # their components recoverable, never discarded.
    for label in ("Pathogenic, risk factor", "Likely benign; other",
                  "Benign, association", "Pathogenic, drug response",
                  "Likely pathogenic, affects"):
        call = policy.classify(label)
        assert call.significance_class not in (SIG_PLP, SIG_BLB), label
        assert call.eligible_primary is False, label
        assert call.exclusion_reason == REASON_NON_BINARY, label
        assert call.modifiers, f"{label}: components were dropped, not preserved"
        assert call.raw == label, f"{label}: the raw field was not preserved verbatim"

    # The surrounding policies this ruling must not disturb.
    assert shipped.get("clinvar.non_binary_disposition") == "preserve_outside_primary"
    assert shipped.get("clinvar.review_star_filter") is None
    assert shipped.get("clinvar.use_review_stars_for_inclusion") is False
    assert shipped.get("clinvar.residue_class_conflict.policy") == "exclude_from_primary"
    assert shipped.get("clinvar.residue_class_conflict.preserve_all_source_records") is True


def test_production_config_activates_the_stage_d_cost_preflight(repo_root):
    """The cost guard is ACTIVE, and it is advisory rather than a budget reducer.

    ``N_CAP`` stays at the production 10000; what changes is that an unusually
    expensive Stage D is escalated before iteration 1 instead of being discovered
    hours in. The threshold is the upper bound of the cost envelope the config
    itself states for this stage, not a new invented one.
    """
    from hotspot3d.robustness.params import RobustnessParams
    from hotspot3d.utils.config import load_config

    shipped = load_config(repo_root / "config" / "pipeline.yaml")
    assert shipped.get("robustness.N_CAP") == 10000, (
        "the cost guard must never be implemented by shrinking the budget")

    budget = shipped.get("robustness.max_wall_seconds")
    assert budget is not None, "the cost preflight is still dormant"
    assert budget == 144000, "40 h — the stated upper bound of this stage's envelope"

    params = RobustnessParams.from_config(shipped)
    assert params.max_wall_seconds == 144000 and params.n_cap == 10000

    # Advisory only: nothing about the design rule or the sampling is conditioned
    # on the budget. These stay exactly as pre-registered.
    assert shipped.get("robustness.K_MAX_rule") == "min(n_S - 1, floor(n_S / 2))"
    assert shipped.get("robustness.sampling") == "seeded_combinatorial_unranking"
    assert shipped.get("robustness.budget_allocation") == "equal_across_remaining_levels"


def test_overlaid_config_verifies_against_its_own_sources(repo_root, tmp_path):
    """A gene overlay is a declared deviation, not a freeze violation.

    ``hotspot3d --gene-config`` builds an overlaid config, and its digest is a
    COMPOSITE of both sources — yet it is 64 hex characters like any plain file
    hash, so the Stage A freeze check cannot infer the form from the digest. It
    must ask the config to recompute its own, which is what this pins.

    Before this was fixed the check hashed ``cfg.path`` unconditionally, so every
    overlaid run blocked at Stage A preflight with a message that printed the two
    truncated digests identically (a composite begins with the base file's hash).
    Stage A's acceptance of an overlay is covered end to end as well: the whole
    integration package now runs under one (``tests/integration/conftest.py``).
    """
    from hotspot3d.utils.config import load_config

    base = repo_root / "config" / "pipeline.yaml"
    overlay = tmp_path / "overlay.yaml"
    overlay.write_text("robustness:\n  N_CAP: 7\n", encoding="utf-8")

    plain = load_config(base)
    assert plain.overlay is None
    assert plain.recompute_sha256() == plain.sha256

    overlaid = load_config(base, gene_overlay=overlay)
    assert overlaid.get("robustness.N_CAP") == 7, "the overlay did not take effect"
    assert overlaid.sha256 != plain.sha256, "the overlay left the digest unchanged"
    assert len(overlaid.sha256) == 64, "a composite digest must stay digest-shaped"
    assert overlaid.recompute_sha256() == overlaid.sha256, (
        "an overlaid config failed to verify against its own sources")

    # ... and tampering with either source is still caught.
    overlay.write_text("robustness:\n  N_CAP: 8\n", encoding="utf-8")
    assert overlaid.recompute_sha256() != overlaid.sha256, (
        "editing the overlay after loading went undetected")


def test_deprecated_robustness_design_absent():
    """The withdrawn full-hotspot-pipeline-rerun robustness design must not return.

    Scans EXECUTABLE code only. Prose that states the prohibition — in a docstring,
    a comment, or a BlockedError message raised by a runtime guard — is not a
    violation; it is the prohibition being enforced.
    """
    banned_calls = [
        (r"\bripley\w*\s*\(", "calls a Ripley's K routine"),
        (r"\bpair_correlation\w*\s*\(", "calls pair correlation"),
        (r"\bradius_scan\w*\s*\(", "calls the r_hot radius scan"),
        (r"\bselect_r_hot\w*\s*\(", "re-selects r_hot"),
        (r"\bderive_(radius|r_hot)_domain\w*\s*\(", "re-derives the r_hot domain"),
        (r"\bbh_fdr\w*\s*\(", "recomputes BH-FDR"),
        (r"\blabel_permutation\w*\s*\(", "re-runs the per-center permutation test"),
    ]
    offenders = []
    for path in _py_files("robustness", "footprint"):
        for lineno, line in _code_lines(path):
            for pattern, why in banned_calls:
                if re.search(pattern, line):
                    offenders.append(f"{path.relative_to(SRC)}:{lineno} {why}")
    assert not offenders, "Deprecated robustness design detected:\n" + "\n".join(offenders)


def test_footprint_robustness_does_not_import_hotspot_statistics():
    """Phase D re-executes the FOOTPRINT methodology only.

    The hotspot-statistics packages must not be imported by footprint/robustness —
    this is the structural guarantee that the withdrawn full-pipeline-rerun design
    cannot come back, whatever the prose says.
    """
    forbidden = ("hotspot3d.spatial", "hotspot3d.hotspot")
    offenders = []
    for path in _py_files("footprint", "robustness"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            mod = None
            if isinstance(node, ast.ImportFrom) and node.module:
                mod = node.module
            elif isinstance(node, ast.Import):
                mod = node.names[0].name
            if mod and any(mod.startswith(f) for f in forbidden):
                offenders.append(f"{path.relative_to(SRC)}: imports {mod}")
    assert not offenders, (
        "footprint/robustness imports hotspot-statistics code:\n" + "\n".join(offenders))


def test_stage_b_does_not_import_footprint_or_annotation():
    """Information barrier: discovery cannot see geometry or biology."""
    forbidden = {"hotspot3d.footprint", "hotspot3d.robustness", "hotspot3d.annotation"}
    offenders = []
    for path in _py_files("spatial", "hotspot"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            mod = None
            if isinstance(node, ast.ImportFrom) and node.module:
                mod = node.module
            elif isinstance(node, ast.Import):
                mod = node.names[0].name
            if mod and any(mod.startswith(f) for f in forbidden):
                offenders.append(f"{path.relative_to(SRC)}: imports {mod}")
    assert not offenders, "Barrier violation:\n" + "\n".join(offenders)


def test_footprint_robustness_does_not_import_annotation():
    offenders = []
    for path in _py_files("footprint", "robustness"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            mod = None
            if isinstance(node, ast.ImportFrom) and node.module:
                mod = node.module
            elif isinstance(node, ast.Import):
                mod = node.names[0].name
            if mod and mod.startswith("hotspot3d.annotation"):
                offenders.append(f"{path.relative_to(SRC)}: imports {mod}")
    assert not offenders, "Barrier violation:\n" + "\n".join(offenders)


def test_annotation_never_imported_by_upstream():
    """Annotation is terminal — nothing upstream may depend on it."""
    offenders = []
    for path in _py_files("data", "structure", "spatial", "hotspot", "footprint",
                          "robustness"):
        text = path.read_text(encoding="utf-8")
        if "hotspot3d.annotation" in text:
            offenders.append(str(path.relative_to(SRC)))
    assert not offenders, f"Upstream depends on annotation: {offenders}"


def test_no_network_calls_in_barred_stages():
    """Stages B, C and D hold no network tools by design; their code must match."""
    banned = ("requests.", "urllib.request", "http.client", "socket.socket(",
              "urlopen(", "httpx.")
    offenders = []
    for path in _py_files("spatial", "hotspot", "footprint", "robustness"):
        text = path.read_text(encoding="utf-8")
        for token in banned:
            if token in text:
                offenders.append(f"{path.relative_to(SRC)}: {token}")
    assert not offenders, "Network access in a barred stage:\n" + "\n".join(offenders)


# --------------------------------------------------------------------------- #
# 2. no duplicated scientific implementation with conflicting formulas         #
# --------------------------------------------------------------------------- #

# Exact function names, anchored — a substring match would flag unrelated helpers
# (e.g. "indices" contains "dice").
SHARED_PRIMITIVES = {
    "mcc_from_confusion": r"^\s*def\s+mcc_from_confusion\s*\(",
    "classification_metrics": r"^\s*def\s+classification_metrics\s*\(",
    "jaccard": r"^\s*def\s+jaccard\s*\(",
    "dice": r"^\s*def\s+dice\s*\(",
    "pareto_front": r"^\s*def\s+pareto_front\s*\(",
    "min_max_normalize": r"^\s*def\s+min_max_normalize\s*\(",
    "distance_to_ideal": r"^\s*def\s+distance_to_ideal\s*\(",
    "pairwise_distances": r"^\s*def\s+pairwise_distances\s*\(",
    "mst_edges": r"^\s*def\s+mst_edges\s*\(",
    "hausdorff95": r"^\s*def\s+hausdorff95\s*\(",
    "boundary_diagnostic": r"^\s*def\s+boundary_diagnostic\s*\(",
}


@pytest.mark.parametrize("name,pattern", sorted(SHARED_PRIMITIVES.items()))
def test_shared_primitive_defined_once(name, pattern):
    """Each shared formula has exactly ONE definition, and it lives in utils."""
    definitions = []
    for path in _py_files(*LEAD_PKGS, *[p for pkgs in AGENT_PKGS.values() for p in pkgs]):
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if re.search(pattern, line):
                definitions.append(f"{path.relative_to(SRC)}:{lineno}")
    if not definitions:
        pytest.skip(f"{name} not yet implemented")
    assert len(definitions) == 1, (
        f"'{name}' is defined {len(definitions)} times — duplicated scientific "
        f"implementations can diverge:\n" + "\n".join(definitions))
    assert definitions[0].startswith("utils/"), (
        f"'{name}' must live in the Lead-owned utils package, found at {definitions[0]}")


def test_agents_use_shared_multiobjective_machinery():
    """Both selections must call the ONE frozen Pareto/normalization implementation."""
    for pkgs, label in ((["hotspot"], "r_hot"), (["footprint"], "r_fp")):
        files = _py_files(*pkgs)
        if not files:
            pytest.skip(f"{label} package not implemented yet")
        text = "\n".join(p.read_text(encoding="utf-8") for p in files)
        if "select(" not in text and "Candidate" not in text:
            pytest.skip(f"{label} selection not implemented yet")
        assert "multiobjective" in text, (
            f"{label} selection does not import utils.multiobjective — it must not "
            f"re-implement the frozen F15 machinery")


# --------------------------------------------------------------------------- #
# 3. canonical filenames and ownership                                         #
# --------------------------------------------------------------------------- #

def test_stage_ownership_partition_is_disjoint():
    from hotspot3d.utils.runctx import AGENT_WRITE_SCOPE, STAGE_OWNERS
    seen: dict[str, str] = {}
    for agent, stages in AGENT_WRITE_SCOPE.items():
        for stage in stages:
            assert stage not in seen, (
                f"{stage} claimed by both {seen[stage]} and {agent}")
            seen[stage] = agent
    assert set(seen) == set(STAGE_OWNERS), "write scope does not cover every stage"
    for stage, owner in STAGE_OWNERS.items():
        assert seen[stage] == owner, f"{stage}: owner mismatch {seen[stage]} vs {owner}"


def test_reads_flow_strictly_upstream():
    from hotspot3d.utils.runctx import AGENT_READ_SCOPE, AGENT_WRITE_SCOPE
    order = ["data-structure", "hotspot-statistics", "footprint-robustness",
             "biological-annotation"]
    for i, agent in enumerate(order):
        downstream = set()
        for later in order[i + 1:]:
            downstream |= AGENT_WRITE_SCOPE[later]
        overlap = AGENT_READ_SCOPE[agent] & downstream
        assert not overlap, f"{agent} may read downstream stages {sorted(overlap)}"


def test_manifest_expected_files_match_output_contract():
    """Every canonical filename named in the Output Contract is registered."""
    from hotspot3d.reporting.manifest import EXPECTED
    required = {
        "02_CLINVAR": {"variants_missense_all.tsv", "variants_residue_level.tsv",
                       "variants_excluded_from_primary.tsv",
                       "review_star_distribution.tsv"},
        "05_HOTSPOT_RADIUS": {"r_hot_scan.tsv"},
        "06_FINAL_HOTSPOTS": {"all_residue_center_tests.tsv",
                              "significant_hotspot_centers.tsv",
                              "hotspot_classified_variants.tsv",
                              "hotspot_covered_residues.tsv"},
        "07_FOOTPRINT_RADIUS": {"r_fp_scan.tsv"},
        "09_ROBUSTNESS": {"perturbation_results.tsv", "robustness_profile.json"},
        "10_ANNOTATION": {"hotspot_annotation.tsv", "functional_mechanism_variants.tsv"},
    }
    for stage, names in required.items():
        registered = {n for n, _t, _d in EXPECTED[stage]}
        missing = names - registered
        assert not missing, f"{stage} missing canonical files: {sorted(missing)}"


def test_stage_e_read_contract_reconciles_with_output_contract():
    """Schemas match across stage boundaries.

    Stage E declares every artifact it reads. Each must be a canonical file that an
    upstream stage actually emits, so an upstream rename fails here rather than
    surfacing as a BLOCKING runtime error mid-run.
    """
    from hotspot3d.annotation.inputs import UPSTREAM_PATHS
    from hotspot3d.reporting.manifest import EXPECTED

    unregistered = []
    for key, (stage, rel) in sorted(UPSTREAM_PATHS.items()):
        registered = {n for n, _t, _d in EXPECTED.get(stage, [])}
        name = rel.split("/")[-1]
        if name not in registered and not rel.startswith("structures/"):
            unregistered.append(f"{key}: {stage}/{rel}")
    assert not unregistered, (
        "Stage E reads artifacts that no upstream stage is registered to emit:\n"
        + "\n".join(unregistered))


def test_stage_e_reads_nothing_it_may_not_read():
    """Stage E's declared reads must all lie inside its permitted read scope."""
    from hotspot3d.annotation.inputs import UPSTREAM_PATHS
    from hotspot3d.utils.runctx import AGENT_READ_SCOPE

    allowed = AGENT_READ_SCOPE["biological-annotation"]
    violations = [f"{k}: {stage}" for k, (stage, _rel) in UPSTREAM_PATHS.items()
                  if stage not in allowed]
    assert not violations, f"Stage E declares out-of-scope reads: {violations}"


def test_handoff_required_keys_cover_the_chain():
    from hotspot3d.orchestration.contracts import HANDOFF_REQUIRED_KEYS
    assert set(HANDOFF_REQUIRED_KEYS) == {
        "handoff_01", "handoff_02", "handoff_03", "handoff_04", "handoff_05"}
    # the three residue objects and the r_hot/r_fp separation must be visible
    assert "n_significant_centers" in HANDOFF_REQUIRED_KEYS["handoff_02"]
    assert "n_centers_without_variant" in HANDOFF_REQUIRED_KEYS["handoff_02"]
    assert "footprint_radius" in HANDOFF_REQUIRED_KEYS["handoff_03"]
    assert "hotspot_radius" in HANDOFF_REQUIRED_KEYS["handoff_03"]
    assert "geometric_footprint_robustness_only" in HANDOFF_REQUIRED_KEYS["handoff_04"]


def test_handoff_05_cannot_influence_discovery():
    """Nothing in the terminal handoff may name a discovery parameter."""
    from hotspot3d.orchestration.contracts import HANDOFF_REQUIRED_KEYS
    banned = {"r_hot", "hotspot_radius", "r_fp", "footprint_radius", "q", "B",
              "fdr_method", "coverage_threshold"}
    assert not (set(HANDOFF_REQUIRED_KEYS["handoff_05"]) & banned)


# --------------------------------------------------------------------------- #
# 4. seeds and determinism primitives                                          #
# --------------------------------------------------------------------------- #

def test_seed_derivation_is_context_not_order_dependent():
    from hotspot3d.utils.seeds import SeedRegistry, derive_seed
    a = SeedRegistry(20250101, "RUN")
    b = SeedRegistry(20250101, "RUN")
    ctxs = ["scan|r=7.5", "final_detection", "iter|0003", "global_null|PLP"]
    forward = [a.seed(c) for c in ctxs]
    backward = [b.seed(c) for c in reversed(ctxs)][::-1]
    assert forward == backward, "seeds depend on issue order — parallel runs would diverge"
    assert derive_seed(20250101, "RUN", "x") != derive_seed(20250101, "RUN2", "x")
    assert all(0 <= s < 2 ** 32 for s in forward)


def test_warning_codes_are_frozen_enum():
    from hotspot3d.utils.status import WARNING_CODES, Warning_
    for code in ("PERMUTATION_RESOLUTION_LIMITED", "FALLBACK_RADIUS_DOMAIN",
                 "EXCESSIVE_FOOTPRINT_COVERAGE",
                 "FOOTPRINT_COVERAGE_CONSTRAINT_FAILURE", "NO_SIGNIFICANT_HOTSPOTS"):
        assert code in WARNING_CODES
    with pytest.raises(ValueError):
        Warning_(warning_code="INVENTED_CODE", stage="x", agent_owner="y", message="z")


def test_tsv_conventions():
    from hotspot3d.utils.io import fmt_value
    assert fmt_value(None) == "NA"
    assert fmt_value("") == "NA"
    assert fmt_value(True) == "TRUE"
    assert fmt_value(False) == "FALSE"
    assert fmt_value(float("nan")) == "NA"
    assert fmt_value(float("inf")) == "NA"
    assert fmt_value(1 / 3) == "0.333333"


def test_negative_and_failure_are_distinct_types():
    from hotspot3d.utils.errors import BlockedError, NegativeResult
    assert not issubclass(NegativeResult, BlockedError)
    assert not issubclass(BlockedError, NegativeResult)
