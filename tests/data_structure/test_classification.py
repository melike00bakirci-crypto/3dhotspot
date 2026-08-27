"""Significance normalization, HGVS parsing and transcript/reference-AA handling."""
from __future__ import annotations

import pytest

from hotspot3d.data.classify import (
    CLASS_BLB,
    CLASS_CONFLICT,
    CLASS_PLP,
    COMPOUND_POLICIES,
    POLICY_COMPONENT_WISE,
    POLICY_WHOLE_FIELD,
    REASON_NON_BINARY,
    REASON_NON_CANONICAL_TX,
    REASON_NON_MISSENSE_VARIANT_TYPE,
    REASON_REF_AA_MISMATCH,
    REASON_UNPARSEABLE_HGVS,
    REASON_VUS,
    SIG_BLB,
    SIG_CONFLICTING,
    SIG_OTHER,
    SIG_PLP,
    SIG_RECORD_CONFLICT,
    SIG_VUS,
    SignificancePolicy,
    compound_policy_cost,
    residue_class,
    split_significance,
)
from hotspot3d.data.clinvar import normalize_records, parse_variant_summary, transcript_matches
from hotspot3d.data.cohort import collapse_to_residues, evaluate_records
from hotspot3d.data.hgvs import parse_clinvar_name, parse_protein_hgvs
from hotspot3d.data.stage import _non_missense_true_cause

pytestmark = pytest.mark.unit


@pytest.fixture
def policy(config):
    return SignificancePolicy.from_config(config)


# --- significance normalization --------------------------------------------

@pytest.mark.parametrize("text,expected", [
    ("Pathogenic", SIG_PLP),
    ("Likely pathogenic", SIG_PLP),
    ("Pathogenic/Likely pathogenic", SIG_PLP),
    ("Benign", SIG_BLB),
    ("Likely benign", SIG_BLB),
    ("Benign/Likely benign", SIG_BLB),
    ("Uncertain significance", SIG_VUS),
    ("Conflicting classifications of pathogenicity", SIG_CONFLICTING),
    ("Conflicting interpretations of pathogenicity", SIG_CONFLICTING),
    ("drug response", SIG_OTHER),
    ("risk factor", SIG_OTHER),
    ("not provided", SIG_OTHER),
    ("association", SIG_OTHER),
])
def test_primary_labels_and_non_binary_categories(policy, text, expected):
    assert policy.classify(text).significance_class == expected


def test_aggregate_labels_are_matched_whole(policy):
    """'/' is part of the aggregate label, never a component separator."""
    assert split_significance("Pathogenic/Likely pathogenic") == \
        ["Pathogenic/Likely pathogenic"]
    call = policy.classify("Pathogenic/Likely pathogenic")
    assert call.significance_class == SIG_PLP
    assert call.exact_match is True


def test_compound_significance_is_matched_component_wise_and_flagged():
    """The component_wise READING, tested on its own terms.

    Built explicitly rather than from the shipped config — symmetric with
    ``test_whole_field_policy_keeps_compound_labels_out_of_the_binary_classes``
    below. The shipped ruling is now ``whole_field``, and this code path must stay
    covered regardless of which reading production selects.
    """
    loose = SignificancePolicy(
        pathogenic_labels=("Pathogenic", "Likely pathogenic"),
        benign_labels=("Benign", "Likely benign"),
        compound_policy=POLICY_COMPONENT_WISE,
    )
    call = loose.classify("Pathogenic, risk factor")
    assert call.significance_class == SIG_PLP
    assert call.exact_match is False, "the stricter whole-field reading must stay visible"
    assert call.modifiers == ("risk factor",)


# --- the compound-label reading is a config ruling, not a code decision -----

def test_config_freezes_the_compound_policy(config):
    assert config.get("clinvar.compound_significance_policy") in COMPOUND_POLICIES


def test_whole_field_policy_keeps_compound_labels_out_of_the_binary_classes():
    strict = SignificancePolicy(
        pathogenic_labels=("Pathogenic", "Likely pathogenic"),
        benign_labels=("Benign", "Likely benign"),
        compound_policy=POLICY_WHOLE_FIELD,
    )
    assert strict.classify("Pathogenic").significance_class == SIG_PLP
    assert strict.classify("Pathogenic").exact_match is True

    compound = strict.classify("Pathogenic, risk factor")
    assert compound.significance_class == SIG_OTHER
    assert compound.eligible_primary is False
    assert compound.exclusion_reason == REASON_NON_BINARY
    assert strict.classify("Uncertain significance").significance_class == SIG_VUS


def test_both_policies_agree_on_bare_and_aggregate_labels():
    kwargs = {"pathogenic_labels": ("Pathogenic", "Pathogenic/Likely pathogenic"),
              "benign_labels": ("Benign",)}
    loose = SignificancePolicy(**kwargs, compound_policy=POLICY_COMPONENT_WISE)
    strict = SignificancePolicy(**kwargs, compound_policy=POLICY_WHOLE_FIELD)
    for text in ("Pathogenic", "Pathogenic/Likely pathogenic", "Benign",
                 "Uncertain significance", "drug response"):
        assert loose.classify(text).significance_class == \
            strict.classify(text).significance_class, text


def test_an_unknown_compound_policy_is_rejected_not_guessed():
    with pytest.raises(ValueError, match="never guessed"):
        SignificancePolicy(pathogenic_labels=("Pathogenic",), benign_labels=("Benign",),
                           compound_policy="majority_vote")


@pytest.mark.parametrize("text", [
    "Conflicting interpretations of pathogenicity, risk factor",
    "Conflicting classifications of pathogenicity, other",
    "Uncertain significance, risk factor",
    "Uncertain significance; other",
    "not provided, risk factor",
])
@pytest.mark.parametrize("compound_policy", COMPOUND_POLICIES)
def test_component_matching_never_launders_a_non_binary_record(text, compound_policy):
    """The safety property that makes component-wise matching sound.

    A compound label whose *significance component* is itself non-binary must be
    excluded under BOTH readings. If component-wise matching could promote such a
    record into P/LP or B/LB it would launder conflicting and uncertain evidence
    into the primary comparison, which F12 forbids outright.
    """
    policy = SignificancePolicy(
        pathogenic_labels=("Pathogenic", "Likely pathogenic",
                           "Pathogenic/Likely pathogenic"),
        benign_labels=("Benign", "Likely benign", "Benign/Likely benign"),
        compound_policy=compound_policy,
    )
    call = policy.classify(text)
    assert call.eligible_primary is False
    assert call.significance_class not in (SIG_PLP, SIG_BLB)
    assert call.exclusion_reason is not None


def test_slash_is_never_a_component_separator():
    """'/' belongs to the aggregate labels; splitting on it would invent classes."""
    assert split_significance("Pathogenic/Likely pathogenic") == \
        ["Pathogenic/Likely pathogenic"]
    assert split_significance("Benign/Likely benign, other") == \
        ["Benign/Likely benign", "other"]


@pytest.mark.parametrize("compound_policy", COMPOUND_POLICIES)
def test_a_record_asserting_both_directions_is_never_binary(compound_policy):
    """Both readings exclude it; they differ only in the reason they record.

    component_wise sees two significance components pointing opposite ways and
    calls it RECORD_SIGNIFICANCE_CONFLICT. whole_field — the shipped ruling — never
    splits the field, so "Pathogenic; Benign" matches no bare label and is excluded
    as non-binary instead. The safety property that matters is identical under both
    and is what this pins: such a record never reaches the primary comparison.
    """
    both = SignificancePolicy(
        pathogenic_labels=("Pathogenic", "Likely pathogenic"),
        benign_labels=("Benign", "Likely benign"),
        compound_policy=compound_policy,
    )
    call = both.classify("Pathogenic; Benign")
    assert call.eligible_primary is False
    assert call.significance_class not in (SIG_PLP, SIG_BLB)
    assert call.exclusion_reason is not None
    expected = (SIG_RECORD_CONFLICT if compound_policy == POLICY_COMPONENT_WISE
                else SIG_OTHER)
    assert call.significance_class == expected


@pytest.mark.parametrize("text,reason", [
    ("Uncertain significance", REASON_VUS),
    ("Conflicting classifications of pathogenicity", "conflicting"),
    ("drug response", REASON_NON_BINARY),
])
def test_every_non_binary_category_carries_a_frozen_reason(policy, config, text, reason):
    call = policy.classify(text)
    assert call.exclusion_reason == reason
    assert reason in config.get("clinvar.exclusion_reason_enum")


def test_case_and_whitespace_are_normalized(policy):
    assert policy.classify("  likely PATHOGENIC ").significance_class == SIG_PLP


# --- residue class (F3) -----------------------------------------------------

@pytest.mark.parametrize("n_plp,n_blb,expected", [
    (3, 0, CLASS_PLP), (0, 2, CLASS_BLB), (1, 1, CLASS_CONFLICT), (5, 1, CLASS_CONFLICT),
])
def test_residue_class_never_resolves_a_conflict(n_plp, n_blb, expected):
    assert residue_class(n_plp, n_blb) == expected


def test_residue_class_rejects_an_empty_residue():
    with pytest.raises(ValueError):
        residue_class(0, 0)


# --- protein HGVS -----------------------------------------------------------

@pytest.mark.parametrize("text,ref,pos,alt", [
    ("p.Val600Glu", "V", 600, "E"),
    ("p.V600E", "V", 600, "E"),
    ("p.(Val600Glu)", "V", 600, "E"),
    ("  p.Gly12Asp ", "G", 12, "D"),
])
def test_missense_parsing(text, ref, pos, alt):
    change = parse_protein_hgvs(text)
    assert change.is_missense
    assert (change.ref_aa, change.position, change.alt_aa) == (ref, pos, alt)


@pytest.mark.parametrize("text,reason", [
    ("p.Val600=", "synonymous"),
    ("p.Val600Val", "synonymous"),
    ("p.Gln100Ter", "nonsense"),
    ("p.Gln100*", "nonsense"),
    ("p.Val600fs", "not_a_substitution:fs"),
    ("p.Lys100del", "not_a_substitution:del"),
    ("p.Ala10_Gly12dup", "not_a_substitution:dup"),
    ("c.1799T>A", "not_protein_level_hgvs"),
    ("", "empty_protein_hgvs"),
    (None, "empty_protein_hgvs"),
    ("p.Xyz12Abc", "non_standard_reference_aa:Xyz"),
    ("p.Val600Xyz", "non_standard_alternate_aa:Xyz"),
    ("p.600Glu", "hgvs_p_pattern_unmatched"),
    ("p.Sec12Ala", "non_standard_reference_aa:Sec"),
])
def test_non_missense_changes_are_rejected_with_a_reason(text, reason):
    change = parse_protein_hgvs(text)
    assert change.is_missense is False
    assert change.reason == reason


def test_clinvar_name_decomposition():
    parsed = parse_clinvar_name("NM_004333.6(BRAF):c.1799T>A (p.Val600Glu)")
    assert parsed.transcript == "NM_004333.6"
    assert parsed.gene == "BRAF"
    assert parsed.hgvs_c == "c.1799T>A"
    assert parsed.hgvs_p == "p.Val600Glu"


def test_transcript_comparison_ignores_version_only():
    assert transcript_matches("NM_004333.6", "NM_004333.5") is True
    assert transcript_matches("NM_004333.6", "NM_004334.6") is False
    assert transcript_matches(None, "NM_004333.6") is False


# --- record evaluation ------------------------------------------------------

def _evaluate(rows, policy, sequence, transcript="NM_900001.4"):
    return evaluate_records(normalize_records(rows, gene="SYNGENE"), policy=policy,
                            canonical_sequence=sequence, mane_transcript=transcript)


def test_reference_aa_mismatch_is_excluded_with_a_reason(policy):
    sequence = "MVAKLDEGGH"
    rows = [{"VariationID": "1", "GeneSymbol": "SYNGENE",
             "Name": "NM_900001.4(SYNGENE):c.7A>G (p.Trp3Gly)",
             "ClinicalSignificance": "Pathogenic",
             "ReviewStatus": "criteria provided, single submitter"}]
    ev = _evaluate(rows, policy, sequence)[0]
    assert ev.eligible_primary is False
    assert ev.exclusion_reason == REASON_REF_AA_MISMATCH
    assert ev.ref_aa_match is False
    assert ev.aa_ref_uniprot == "A"


def test_position_beyond_the_canonical_sequence_is_excluded(policy):
    rows = [{"VariationID": "1", "GeneSymbol": "SYNGENE",
             "Name": "NM_900001.4(SYNGENE):c.1A>G (p.Val900Glu)",
             "ClinicalSignificance": "Pathogenic", "ReviewStatus": "practice guideline"}]
    ev = _evaluate(rows, policy, "MVAKLDEGGH")[0]
    assert ev.exclusion_reason == REASON_REF_AA_MISMATCH
    assert "position_out_of_range" in ev.exclusion_detail


def test_non_mane_record_is_remapped_only_when_position_and_ref_aa_agree(policy):
    sequence = "MVAKLDEGGH"
    consistent = {"VariationID": "1", "GeneSymbol": "SYNGENE",
                  "Name": "NM_999999.1(SYNGENE):c.7A>G (p.Ala3Gly)",
                  "ClinicalSignificance": "Pathogenic",
                  "ReviewStatus": "criteria provided, single submitter"}
    inconsistent = {"VariationID": "2", "GeneSymbol": "SYNGENE",
                    "Name": "NM_999999.1(SYNGENE):c.7A>G (p.Trp3Gly)",
                    "ClinicalSignificance": "Pathogenic",
                    "ReviewStatus": "criteria provided, single submitter"}
    good, bad = _evaluate([consistent, inconsistent], policy, sequence)

    assert good.eligible_primary is True
    assert good.on_mane_transcript is False
    assert good.remapped_to_canonical is True

    assert bad.eligible_primary is False
    assert bad.exclusion_reason == REASON_NON_CANONICAL_TX
    assert bad.remapped_to_canonical is False


def test_frameshift_change_is_non_missense_variant_type_not_unparseable(policy):
    """[REVISED — Workflow v2 §1, Lead ruling 2026-08-18] a frameshift change
    parses CLEANLY into a non-missense consequence; it is not a parse failure
    and must not be reported as one. This test was previously pinned to the
    exact defect Workflow v2 §1 fixes (it asserted REASON_UNPARSEABLE_HGVS)."""
    rows = [{"VariationID": "1", "GeneSymbol": "SYNGENE",
             "Name": "NM_900001.4(SYNGENE):c.7del (p.Ala3fs)",
             "ClinicalSignificance": "Pathogenic", "ReviewStatus": "practice guideline"}]
    ev = _evaluate(rows, policy, "MVAKLDEGGH")[0]
    assert ev.exclusion_reason == REASON_NON_MISSENSE_VARIANT_TYPE
    assert ev.exclusion_detail == "not_a_substitution:fs"
    assert ev.residue_index is None


def test_genuinely_unreadable_protein_hgvs_stays_unparseable(policy):
    """A record whose HGVS string could not be interpreted at all — as opposed
    to one that parsed cleanly into a non-missense consequence — is still
    ``unparseable_hgvs``."""
    rows = [{"VariationID": "1", "GeneSymbol": "SYNGENE",
             "Name": "NM_900001.4(SYNGENE):c.7A>G (p.600Glu)",
             "ClinicalSignificance": "Pathogenic", "ReviewStatus": "practice guideline"}]
    ev = _evaluate(rows, policy, "MVAKLDEGGH")[0]
    assert ev.exclusion_reason == REASON_UNPARSEABLE_HGVS
    assert ev.exclusion_detail == "hgvs_p_pattern_unmatched"
    assert ev.residue_index is None


def test_vus_never_reaches_a_binary_class(policy):
    sequence = "MVAKLDEGGH"
    rows = [
        {"VariationID": "1", "GeneSymbol": "SYNGENE",
         "Name": "NM_900001.4(SYNGENE):c.7A>G (p.Ala3Gly)",
         "ClinicalSignificance": "Uncertain significance",
         "ReviewStatus": "reviewed by expert panel"},
        {"VariationID": "2", "GeneSymbol": "SYNGENE",
         "Name": "NM_900001.4(SYNGENE):c.10A>G (p.Lys4Glu)",
         "ClinicalSignificance": "Conflicting classifications of pathogenicity",
         "ReviewStatus": "criteria provided, conflicting classifications"},
    ]
    evals = _evaluate(rows, policy, sequence)
    assert all(not e.eligible_primary for e in evals)
    assert collapse_to_residues(evals, canonical_sequence=sequence) == []


# --- compound-significance policy cost (Workflow v2 §1) ---------------------

def test_compound_policy_cost_counts_records_the_frozen_policy_excludes(config):
    plp = config.get("clinvar.pathogenic_labels")
    blb = config.get("clinvar.benign_labels")
    rows = [
        {"VariationID": "1", "ClinicalSignificance": "Pathogenic, low penetrance"},
        {"VariationID": "2", "ClinicalSignificance": "Likely benign; other"},
        {"VariationID": "3", "ClinicalSignificance": "Pathogenic/Likely pathogenic"},
        {"VariationID": "4", "ClinicalSignificance": "Uncertain significance"},
        {"VariationID": "5", "ClinicalSignificance": "Pathogenic; Benign"},
    ]
    records = normalize_records(rows, gene="SYNGENE")

    result = compound_policy_cost(records, pathogenic_labels=plp, benign_labels=blb,
                                  configured_policy=POLICY_WHOLE_FIELD)

    assert result["configured_policy"] == POLICY_WHOLE_FIELD
    assert result["counterfactual_policy"] == POLICY_COMPONENT_WISE
    assert result["n_records_evaluated"] == 5
    # "Pathogenic, low penetrance" and "Likely benign; other" are compound labels
    # whole_field excludes (no exact bare-label match) that component_wise would
    # have read into PLP/BLB via their first component.
    assert result["n_excluded_by_configured_but_kept_by_counterfactual"] == 2
    assert set(result["example_clinical_significance_values_affected"]) == {
        "Pathogenic, low penetrance", "Likely benign; other",
    }
    # Aggregate labels ("Pathogenic/Likely pathogenic"), VUS and a genuine
    # both-directions conflict ("Pathogenic; Benign") are excluded from the
    # binary classes under BOTH policies, so they never contribute to the cost.
    assert result["n_kept_by_configured_but_excluded_by_counterfactual"] == 0


def test_compound_policy_cost_is_zero_with_no_compound_labels(config):
    plp = config.get("clinvar.pathogenic_labels")
    blb = config.get("clinvar.benign_labels")
    rows = [
        {"VariationID": "1", "ClinicalSignificance": "Pathogenic"},
        {"VariationID": "2", "ClinicalSignificance": "Benign"},
        {"VariationID": "3", "ClinicalSignificance": "Uncertain significance"},
    ]
    records = normalize_records(rows, gene="SYNGENE")

    result = compound_policy_cost(records, pathogenic_labels=plp, benign_labels=blb,
                                  configured_policy=POLICY_WHOLE_FIELD)
    assert result["n_excluded_by_configured_but_kept_by_counterfactual"] == 0
    assert result["example_clinical_significance_values_affected"] == []


def test_compound_policy_cost_is_symmetric_under_the_other_configured_policy(config):
    """The comparator is always "the other" policy, whichever is configured."""
    plp = config.get("clinvar.pathogenic_labels")
    blb = config.get("clinvar.benign_labels")
    rows = [{"VariationID": "1", "ClinicalSignificance": "Pathogenic, low penetrance"}]
    records = normalize_records(rows, gene="SYNGENE")

    result = compound_policy_cost(records, pathogenic_labels=plp, benign_labels=blb,
                                  configured_policy=POLICY_COMPONENT_WISE)
    assert result["configured_policy"] == POLICY_COMPONENT_WISE
    assert result["counterfactual_policy"] == POLICY_WHOLE_FIELD
    # component_wise KEEPS this record; whole_field would have excluded it —
    # the opposite direction of the frozen (whole_field) configuration's cost.
    assert result["n_kept_by_configured_but_excluded_by_counterfactual"] == 1
    assert result["n_excluded_by_configured_but_kept_by_counterfactual"] == 0


# --- exclusion-ledger true-cause accounting (Workflow v2 §1) ----------------

@pytest.mark.parametrize("detail,expected_bucket", [
    ("synonymous", "synonymous"),
    ("nonsense", "nonsense"),
    ("not_a_substitution:fs", "frameshift"),
    ("not_a_substitution:del", "indel_or_other_non_substitution"),
    ("not_a_substitution:dup", "indel_or_other_non_substitution"),
    ("not_a_substitution:ins", "indel_or_other_non_substitution"),
    ("not_a_substitution:delins", "indel_or_other_non_substitution"),
    ("reference_is_termination", "stop_loss_or_extension"),
    ("non_standard_reference_aa:Sec", "non_standard_amino_acid"),
    ("non_standard_alternate_aa:Pyl", "non_standard_amino_acid"),
    ("non_positive_position", "invalid_position"),
    ("hgvs_p_pattern_unmatched", "true_parse_failure"),
    ("not_protein_level_hgvs", "true_parse_failure"),
    ("empty_protein_hgvs", "true_parse_failure"),
])
def test_non_missense_true_cause_never_lumps_variant_types_as_parse_failures(
        detail, expected_bucket):
    bucket = _non_missense_true_cause(detail)
    assert bucket == expected_bucket
    # The four record types Workflow v2 §1 names explicitly must never land in
    # the parse-failure bucket.
    if detail in ("synonymous", "nonsense", "not_a_substitution:fs",
                 "not_a_substitution:del", "not_a_substitution:dup",
                 "not_a_substitution:ins", "not_a_substitution:delins"):
        assert bucket != "true_parse_failure"


def test_stage_a_exclusion_ledger_distinguishes_true_causes_end_to_end(policy):
    """Synonymous, nonsense, frameshift and indel records are excluded with
    DIFFERENT recorded true causes, never one undifferentiated parse-failure
    bucket — n_retrieved = n_eligible + n_excluded still holds exactly.
    """
    sequence = "M" + "A" * 20
    rows = [
        {"VariationID": "1", "GeneSymbol": "SYNGENE",
         "Name": "NM_900001.4(SYNGENE):c.15C>T (p.Ala5=)",
         "ClinicalSignificance": "Pathogenic", "ReviewStatus": "practice guideline"},
        {"VariationID": "2", "GeneSymbol": "SYNGENE",
         "Name": "NM_900001.4(SYNGENE):c.19C>T (p.Ala7Ter)",
         "ClinicalSignificance": "Pathogenic", "ReviewStatus": "practice guideline"},
        {"VariationID": "3", "GeneSymbol": "SYNGENE",
         "Name": "NM_900001.4(SYNGENE):c.25del (p.Ala9fs)",
         "ClinicalSignificance": "Pathogenic", "ReviewStatus": "practice guideline"},
        {"VariationID": "4", "GeneSymbol": "SYNGENE",
         "Name": "NM_900001.4(SYNGENE):c.31_33del (p.Ala11del)",
         "ClinicalSignificance": "Pathogenic", "ReviewStatus": "practice guideline"},
    ]
    evals = _evaluate(rows, policy, sequence)
    # All four parsed CLEANLY into a non-missense consequence — none is a parse
    # failure, so all four carry REASON_NON_MISSENSE_VARIANT_TYPE, not
    # REASON_UNPARSEABLE_HGVS (Workflow v2 §1, Lead ruling 2026-08-18).
    assert all(e.exclusion_reason == REASON_NON_MISSENSE_VARIANT_TYPE for e in evals)

    buckets = {ev.record.variation_id: _non_missense_true_cause(ev.exclusion_detail)
               for ev in evals}
    assert buckets == {
        "1": "synonymous", "2": "nonsense", "3": "frameshift",
        "4": "indel_or_other_non_substitution",
    }
    # every record landed in exactly one place: none eligible, all excluded.
    assert sum(1 for e in evals if e.eligible_primary) == 0
    assert len(evals) == 4


# --- retrieval-scope Type filter (Workflow v2 §1, Finding B) ----------------

def test_deletion_duplication_insertion_types_are_retrieved_and_excluded_with_true_cause(
        policy):
    """ClinVar ``Type=Deletion/Duplication/Insertion`` records resolving to a
    frameshift or in-frame indel must be RETRIEVED (present in n_retrieved) and
    then excluded with their true-cause reason — never silently absent before
    any accounting layer sees them, and never lumped as unparseable_hgvs.
    """
    sequence = "M" + "A" * 20
    tsv = (
        "VariationID\tGeneSymbol\tName\tType\tClinicalSignificance\tReviewStatus\tAssembly\n"
        "101\tSYNGENE\tNM_900001.4(SYNGENE):c.25del (p.Ala9fs)\tDeletion\t"
        "Pathogenic\tpractice guideline\tGRCh38\n"
        "102\tSYNGENE\tNM_900001.4(SYNGENE):c.31_33dup (p.Ala11dup)\tDuplication\t"
        "Pathogenic\tpractice guideline\tGRCh38\n"
        "103\tSYNGENE\tNM_900001.4(SYNGENE):c.40_42insGGT (p.Ala14_Val15insGly)\tInsertion\t"
        "Pathogenic\tpractice guideline\tGRCh38\n"
        "104\tSYNGENE\tNM_900001.4(SYNGENE):c.1_100del (p.?)\tMicrosatellite\t"
        "Pathogenic\tpractice guideline\tGRCh38\n"
    )

    retrieved = parse_variant_summary(tsv, gene="SYNGENE")
    retrieved_ids = {r["VariationID"] for r in retrieved}
    # Deletion/Duplication/Insertion now reach n_retrieved (the fix)...
    assert retrieved_ids == {"101", "102", "103"}
    # ...while Microsatellite (not a point-level protein-coding change) stays
    # excluded from retrieval entirely, exactly as before.
    assert "104" not in retrieved_ids

    evals = evaluate_records(
        normalize_records(retrieved, gene="SYNGENE"), policy=policy,
        canonical_sequence=sequence, mane_transcript="NM_900001.4")
    assert len(evals) == 3
    reasons = {ev.record.variation_id: ev.exclusion_reason for ev in evals}
    assert reasons == {
        "101": REASON_NON_MISSENSE_VARIANT_TYPE,
        "102": REASON_NON_MISSENSE_VARIANT_TYPE,
        "103": REASON_NON_MISSENSE_VARIANT_TYPE,
    }
    true_causes = {ev.record.variation_id: _non_missense_true_cause(ev.exclusion_detail)
                  for ev in evals}
    assert true_causes == {
        "101": "frameshift",
        "102": "indel_or_other_non_substitution",
        "103": "indel_or_other_non_substitution",
    }
    assert all(not ev.eligible_primary for ev in evals)
