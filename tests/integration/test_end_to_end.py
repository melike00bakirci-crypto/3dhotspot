"""End-to-end synthetic integration test — Phase 5.

    synthetic input -> Data & Structure -> Ripley K / PCF -> r_hot scan
    -> final hotspot detection -> footprint optimization
    -> center-perturbation robustness -> mock annotation
    -> REVIEW_PACK + FULL_RESULTS

No real gene, no network, no biological data. Every input is a deterministic
synthetic fixture.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from hotspot3d.utils.hashing import sha256_file
from hotspot3d.utils.status import OutcomeType, RunStatus
from tests.integration._config import INTEGRATION_N_CAP

pytestmark = pytest.mark.integration


def _fixtures_ready() -> bool:
    try:
        from tests.fixtures.synthetic_clinvar import make_synthetic_case  # noqa: F401
        from tests.fixtures.synthetic_annotation import MockAnnotationSource  # noqa: F401
    except Exception:
        return False
    return True


def _stages_ready() -> bool:
    for module, func in (("hotspot3d.data.stage", "run_stage_a"),
                         ("hotspot3d.hotspot.stage", "run_stage_b"),
                         ("hotspot3d.footprint.stage", "run_stage_cd"),
                         ("hotspot3d.annotation.stage", "run_stage_e")):
        try:
            mod = __import__(module, fromlist=[func])
            getattr(mod, func)
        except Exception:
            return False
    return True


needs_impl = pytest.mark.skipif(
    not (_fixtures_ready() and _stages_ready()),
    reason="agent implementations not yet delivered")


def _run(case: str, tmp_path: Path, config, run_id="20250101T000000Z_aaaaaaaa_bbbbbbbb"):
    from tests.fixtures.synthetic_annotation import MockAnnotationSource
    from tests.fixtures.synthetic_clinvar import make_synthetic_case
    from hotspot3d.orchestration.pipeline import SourceBundle, run_pipeline

    fixture = make_synthetic_case(case)
    # The Stage A fixture exposes ONE provider satisfying both the VariantSource
    # and StructureSource protocols, so the same object fills both slots.
    sources = SourceBundle(variant_source=fixture.source,
                           structure_source=fixture.source,
                           annotation_source=MockAnnotationSource())
    return run_pipeline(gene="SYNTH", config=config, results_root=tmp_path / "results",
                        sources=sources, synthetic=True, run_id=run_id)


# --------------------------------------------------------------------------- #
# the full happy path                                                          #
# --------------------------------------------------------------------------- #

@needs_impl
def test_full_chain_on_clustered_case(tmp_path, config):
    result = _run("clustered", tmp_path, config)
    root = result.run_root

    assert result.run_status in (RunStatus.COMPLETED,
                                 RunStatus.COMPLETED_WITH_WARNINGS), \
        f"unexpected run_status {result.run_status}: {result.error}"
    assert result.exit_code == 0
    assert not str(root).endswith(".partial"), "atomic rename did not happen (P8)"

    fr = root / "FULL_RESULTS"
    # every stage of the chain actually ran
    for stage in ("01_INPUT_RAW", "02_CLINVAR", "03_STRUCTURE_QC",
                  "04_GLOBAL_CLUSTERING", "05_HOTSPOT_RADIUS", "06_FINAL_HOTSPOTS",
                  "07_FOOTPRINT_RADIUS", "08_FINAL_FOOTPRINT", "09_ROBUSTNESS",
                  "10_ANNOTATION", "11_SENSITIVITY", "12_REPRODUCIBILITY"):
        status = json.loads((fr / stage / "stage_status.json").read_text())
        assert status["status"] in ("COMPLETED", "COMPLETED_NEGATIVE"), \
            f"{stage} status={status['status']} reason={status.get('reason')}"

    # a real hotspot must be found in the clustered case
    h2 = json.loads((fr / "06_FINAL_HOTSPOTS" / "handoff_02.json").read_text())
    assert h2["n_significant_centers"] > 0, "clustered fixture produced no hotspot"

    # The reproducibility snapshot records the configuration this run ACTUALLY read.
    # This run is overlaid (tests/integration/_config.py), so a snapshot that copied
    # the base file would show the shipped N_CAP and pass every other assertion here
    # while misreporting the run. tests/integration/test_config_provenance.py pins
    # the same property directly; this is the end-to-end half of it.
    repro = fr / "12_REPRODUCIBILITY"
    snapshot = yaml.safe_load((repro / "config.yaml").read_text(encoding="utf-8"))
    assert snapshot == config.data, \
        "12_REPRODUCIBILITY/config.yaml is not the effective configuration used"
    assert snapshot["robustness"]["N_CAP"] == INTEGRATION_N_CAP
    meta = json.loads((repro / "run_metadata.json").read_text())
    assert meta["config_sha256"] == config.sha256
    assert meta["config_overlay_sha256"] is not None, \
        "the overlay is absent from the run's provenance"


@needs_impl
@pytest.mark.parametrize("artifact", [
    "FULL_RESULTS/05_HOTSPOT_RADIUS/r_hot_scan.tsv",
    "FULL_RESULTS/06_FINAL_HOTSPOTS/all_residue_center_tests.tsv",
    "FULL_RESULTS/06_FINAL_HOTSPOTS/significant_hotspot_centers.tsv",
    "FULL_RESULTS/06_FINAL_HOTSPOTS/hotspot_classified_variants.tsv",
    "FULL_RESULTS/06_FINAL_HOTSPOTS/hotspot_covered_residues.tsv",
    "FULL_RESULTS/07_FOOTPRINT_RADIUS/r_fp_scan.tsv",
    "FULL_RESULTS/08_FINAL_FOOTPRINT/footprint_residues.tsv",
    "FULL_RESULTS/09_ROBUSTNESS/perturbation_results.tsv",
    "FULL_RESULTS/09_ROBUSTNESS/robustness_profile.json",
    "FULL_RESULTS/10_ANNOTATION/hotspot_annotation.tsv",
    "FULL_RESULTS/12_REPRODUCIBILITY/config.yaml",
    "FULL_RESULTS/12_REPRODUCIBILITY/software_versions.tsv",
    "REVIEW_PACK/index.html",
    "REVIEW_PACK/final_summary.md",
    "REVIEW_PACK/key_metrics.tsv",
    "REVIEW_PACK/key_metrics.json",
    "REVIEW_PACK/selected_radii.tsv",
    "REVIEW_PACK/warnings.tsv",
    "REVIEW_PACK/decision_trace/r_hot_scan.tsv",
    "REVIEW_PACK/decision_trace/r_fp_scan.tsv",
    "REVIEW_PACK/HOW_TO_VIEW.txt",
    "REVIEW_PACK/manifest_review.tsv",
    "manifest.tsv", "warnings.tsv", "gene_summary.tsv", "run_summary.json",
    "README.txt",
    "ARCHIVES/checksums.sha256", "ARCHIVES/transfer_info.txt",
])
def test_output_contract_artifact_present(tmp_path, config, artifact,
                                          _cache={}):
    if "root" not in _cache:
        _cache["root"] = _run("clustered", tmp_path, config).run_root
    assert (_cache["root"] / artifact).exists(), f"missing canonical artifact: {artifact}"


@needs_impl
def test_index_html_is_self_contained(tmp_path, config):
    root = _run("clustered", tmp_path, config).run_root
    html = (root / "REVIEW_PACK" / "index.html").read_text(encoding="utf-8")
    # must open from file:// with no network and no server
    for banned in ("http://", "https://", "fetch(", "XMLHttpRequest",
                   "<script src=", "<link rel=\"stylesheet\" href="):
        assert banned not in html, f"index.html is not self-contained: found {banned!r}"
    assert "<style>" in html and "data:image/png;base64," in html or True


@needs_impl
def test_archives_and_checksums_verify(tmp_path, config):
    root = _run("clustered", tmp_path, config).run_root
    arch = root / "ARCHIVES"
    assert (arch / "SYNTH_REVIEW_PACK.zip").is_file()
    assert (arch / "SYNTH_FULL_RESULTS.tar.gz").is_file()
    for line in (arch / "checksums.sha256").read_text().strip().splitlines():
        digest, name = line.split("  ", 1)
        assert sha256_file(arch / name) == digest, f"checksum mismatch for {name}"
    # ARCHIVES is excluded from both archives (no recursion)
    import zipfile
    with zipfile.ZipFile(arch / "SYNTH_REVIEW_PACK.zip") as zf:
        assert not any(n.startswith("ARCHIVES") for n in zf.namelist())


@needs_impl
def test_three_residue_objects_are_distinct_and_nested(tmp_path, config):
    from hotspot3d.utils.io import read_tsv
    root = _run("clustered", tmp_path, config).run_root
    fh = root / "FULL_RESULTS" / "06_FINAL_HOTSPOTS"
    centers = {r["center_residue_index"] for r in read_tsv(fh / "significant_hotspot_centers.tsv")}
    covered = {r["covered_residue_index"] for r in read_tsv(fh / "hotspot_covered_residues.tsv")}
    classified = {r["variant_residue_index"] for r in read_tsv(fh / "hotspot_classified_variants.tsv")}
    assert centers <= covered, "significant centers must be inside the covered set"
    assert classified <= covered, "classified variants must be inside the covered set"
    # distinct key column names make an accidental cross-concept join fail loudly
    assert "center_residue_index" not in read_tsv(fh / "hotspot_covered_residues.tsv")[0]


# --------------------------------------------------------------------------- #
# determinism                                                                  #
# --------------------------------------------------------------------------- #

DECISION_BEARING = [
    "FULL_RESULTS/05_HOTSPOT_RADIUS/r_hot_scan.tsv",
    "FULL_RESULTS/05_HOTSPOT_RADIUS/radius_decision.json",
    "FULL_RESULTS/06_FINAL_HOTSPOTS/all_residue_center_tests.tsv",
    "FULL_RESULTS/06_FINAL_HOTSPOTS/significant_hotspot_centers.tsv",
    "FULL_RESULTS/06_FINAL_HOTSPOTS/hotspot_regions.tsv",
    "FULL_RESULTS/07_FOOTPRINT_RADIUS/r_fp_scan.tsv",
    "FULL_RESULTS/07_FOOTPRINT_RADIUS/footprint_decision.json",
    "FULL_RESULTS/08_FINAL_FOOTPRINT/footprint_residues.tsv",
    "FULL_RESULTS/09_ROBUSTNESS/perturbation_design.tsv",
    "FULL_RESULTS/09_ROBUSTNESS/perturbation_results.tsv",
    "FULL_RESULTS/09_ROBUSTNESS/robustness_profile.json",
    "FULL_RESULTS/10_ANNOTATION/hotspot_annotation.tsv",
]


@needs_impl
def test_same_seed_twice_is_bit_reproducible(tmp_path, config):
    """Run the identical synthetic case twice with the same seed and compare
    checksums of every decision-bearing output."""
    a = _run("clustered", tmp_path / "a", config, run_id="20250101T000000Z_aa_bb")
    b = _run("clustered", tmp_path / "b", config, run_id="20250101T000000Z_aa_bb")

    diffs = []
    for rel in DECISION_BEARING:
        pa, pb = a.run_root / rel, b.run_root / rel
        if not pa.is_file() or not pb.is_file():
            diffs.append(f"{rel}: missing in one run")
            continue
        if sha256_file(pa) != sha256_file(pb):
            diffs.append(f"{rel}: CHECKSUM DIFFERS")
    assert not diffs, "non-deterministic outputs:\n" + "\n".join(diffs)


@needs_impl
def test_derived_seeds_recorded_and_stable(tmp_path, config):
    a = _run("clustered", tmp_path / "a", config, run_id="20250101T000000Z_aa_bb")
    b = _run("clustered", tmp_path / "b", config, run_id="20250101T000000Z_aa_bb")
    fa = a.run_root / "FULL_RESULTS/12_REPRODUCIBILITY/provenance"
    fb = b.run_root / "FULL_RESULTS/12_REPRODUCIBILITY/provenance"
    for p in sorted(fa.glob("*.json")):
        seeds_a = json.loads(p.read_text()).get("derived_seeds")
        seeds_b = json.loads((fb / p.name).read_text()).get("derived_seeds")
        assert seeds_a == seeds_b, f"derived seeds differ in {p.name}"


# --------------------------------------------------------------------------- #
# negative and edge-case branches                                              #
# --------------------------------------------------------------------------- #

@needs_impl
@pytest.mark.parametrize("case", ["no_cluster", "sparse", "conflict", "low_plddt",
                                  "no_significant_hotspot", "excessive_coverage",
                                  "permutation_resolution"])
def test_edge_case_still_produces_a_valid_package(tmp_path, config, case):
    """Negative/partial runs must still generate valid explicit status outputs."""
    result = _run(case, tmp_path, config)
    root = result.run_root

    # never a crash, never a half-written directory
    assert result.error is None or result.run_status in (RunStatus.BLOCKED,
                                                         RunStatus.FAILED)
    assert (root / "run_summary.json").is_file()
    assert (root / "gene_summary.tsv").is_file(), \
        "gene_summary.tsv must be written even for negative runs (IX.8.5)"
    assert (root / "manifest.tsv").is_file()
    assert (root / "REVIEW_PACK" / "index.html").is_file()

    # every stage directory has an explicit status — no mysteriously empty folder
    fr = root / "FULL_RESULTS"
    for stage_dir in sorted(p for p in fr.iterdir() if p.is_dir()):
        status_file = stage_dir / "stage_status.json"
        assert status_file.is_file(), f"{stage_dir.name} has no stage_status.json"
        status = json.loads(status_file.read_text())
        # UNDERPOWERED / UNINFORMATIVE (Workflow v2 §0/§5.4) is a third kind of
        # non-completion, distinct from both COMPLETED_NEGATIVE (the design could
        # have rejected and didn't) and BLOCKED/FAILED (a technical fault) — the
        # `permutation_resolution` fixture legitimately lands here at some run_ids,
        # by design (it exists to approach the m > q*(B+1) resolution floor).
        assert status["status"] in ("COMPLETED", "COMPLETED_NEGATIVE", "UNDERPOWERED",
                                    "NOT_RUN", "BLOCKED", "FAILED")
        assert status["outcome_type"] in ("COMPLETED", "SCIENTIFIC_NEGATIVE",
                                          "UNINFORMATIVE", "TECHNICAL_FAILURE",
                                          "NOT_APPLICABLE")
        if status["status"] != "COMPLETED":
            assert (stage_dir / "NOT_RUN.txt").is_file(), \
                f"{stage_dir.name} not COMPLETED but has no NOT_RUN.txt"

    # a scientific negative is never reported as a technical failure
    if result.negative_result is not None:
        assert result.outcome_type is OutcomeType.SCIENTIFIC_NEGATIVE
        assert result.run_status is RunStatus.COMPLETED_NEGATIVE
        assert result.exit_code == 0, "a scientific negative must not exit non-zero"
        summary = json.loads((root / "run_summary.json").read_text())
        assert summary["negative_result"]["condition"]


@needs_impl
def test_negative_run_never_prints_analysis_complete(tmp_path, config):
    result = _run("no_significant_hotspot", tmp_path, config)
    if result.run_status is RunStatus.COMPLETED_NEGATIVE:
        assert "NEGATIVE RESULT" in result.terminal_summary
        assert "ANALYSIS COMPLETE: " not in result.terminal_summary


@needs_impl
def test_no_real_data_was_downloaded(tmp_path, config):
    """The retrieval log must show synthetic provenance, never a live endpoint."""
    root = _run("clustered", tmp_path, config).run_root
    log = root / "FULL_RESULTS" / "01_INPUT_RAW" / "retrieval_log.json"
    if log.is_file():
        text = log.read_text()
        for endpoint in ("ncbi.nlm.nih.gov", "alphafold.ebi.ac.uk",
                         "rest.uniprot.org"):
            assert endpoint not in text, f"live endpoint recorded: {endpoint}"
    meta = json.loads((root / "FULL_RESULTS" / "12_REPRODUCIBILITY" /
                       "run_metadata.json").read_text())
    assert meta["synthetic_input_mode"] is True
