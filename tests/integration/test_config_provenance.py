"""Configuration identity and reproducibility-snapshot provenance.

Two defects are pinned here. Both were implementation faults in how a run RECORDS
and IDENTIFIES its configuration; neither touches scientific methodology, and no
frozen constant, formula, threshold or seed is involved.

1. **The snapshot was not the configuration the run used.**
   ``_snapshot_reproducibility`` copied ``config/pipeline.yaml`` verbatim, so an
   overlaid run wrote a ``12_REPRODUCIBILITY/config.yaml`` — labelled "frozen
   config snapshot" in the manifest — showing values the run never read.

2. **The digest did not identify the effective configuration.**
   An overlaid digest was ``sha256(base)[:32] + sha256(overlay)[:32]``, and
   ``run_id`` takes ``config_sha256[0:8]`` — eight characters that came entirely
   from the base file. Every overlay of one base produced the same run_id
   component, so two materially different configurations were indistinguishable by
   run identity.

The tests are ordered A-D after the requirements they discharge.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pytest
import yaml

from hotspot3d.utils.config import assert_frozen_methodology, load_config
from hotspot3d.utils.hashing import sha256_file
from hotspot3d.utils.runctx import RunContext, make_run_id
from tests.integration._config import INTEGRATION_N_CAP, PIPELINE_CONFIG

pytestmark = pytest.mark.integration

#: A value no shipped parameter holds, so "the overlay took effect" is never
#: satisfied by coincidence.
SENTINEL_N_CAP = 137

#: A fixed instant and a fixed code version, so a run_id comparison varies in the
#: configuration component alone.
_AT = datetime(2025, 1, 1, tzinfo=timezone.utc)
_CODE = "c0defeed" * 8


def _overlay(directory: Path, n_cap: int) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"overlay_{n_cap}.yaml"
    path.write_text(yaml.safe_dump({"robustness": {"N_CAP": n_cap}}), encoding="utf-8")
    return path


def _snapshot(tmp_path: Path, cfg) -> Path:
    """Run only the reproducibility snapshot and return its directory.

    Calls the private ``_snapshot_reproducibility`` directly. Requirement A is
    about what that function writes, and driving it through ``run_pipeline`` would
    mean a multi-minute full-chain run to assert something the function decides in
    isolation. The end-to-end coverage still exists:
    ``test_end_to_end.test_full_chain_on_clustered_case`` makes the same
    comparison against a real completed run.
    """
    from hotspot3d.orchestration.pipeline import _snapshot_reproducibility

    ctx = RunContext.create(gene="SYNTH", config=cfg, results_root=tmp_path / "results",
                            synthetic=True, run_id="20250101T000000Z_prov0001_prov0002")
    _snapshot_reproducibility(ctx)
    return ctx.full_results / "12_REPRODUCIBILITY"


# --------------------------------------------------------------------------- #
# A. the snapshot IS the effective runtime configuration                      #
# --------------------------------------------------------------------------- #

def test_a_snapshot_of_an_overlaid_run_is_the_effective_config(tmp_path):
    """config.yaml must be what the run read, overlay merged in."""
    cfg = load_config(PIPELINE_CONFIG,
                      gene_overlay=_overlay(tmp_path / "ov", SENTINEL_N_CAP))
    repro = _snapshot(tmp_path, cfg)

    snapshot = yaml.safe_load((repro / "config.yaml").read_text(encoding="utf-8"))
    assert snapshot == cfg.data, (
        "12_REPRODUCIBILITY/config.yaml is not the effective configuration the "
        "run used")
    assert snapshot["robustness"]["N_CAP"] == SENTINEL_N_CAP, (
        "the snapshot shows the base value, so the overlay is invisible in it — "
        "this is exactly the defect being pinned")

    # The snapshot is the very text the digest identifies, not a re-rendering.
    assert (repro / "config.yaml").read_text(encoding="utf-8") == cfg.effective_yaml()

    # ... and the annotated sources survive beside it, because the effective dump
    # cannot carry the YAML comments that state each parameter's FROZEN provenance.
    assert (repro / "config_base.yaml").read_bytes() == PIPELINE_CONFIG.read_bytes()
    assert (repro / "config_overlay.yaml").read_bytes() == Path(cfg.overlay).read_bytes()

    meta = json.loads((repro / "run_metadata.json").read_text(encoding="utf-8"))
    assert meta["config_sha256"] == cfg.sha256
    assert meta["config_base_sha256"] == sha256_file(PIPELINE_CONFIG)
    assert meta["config_overlay_sha256"] == sha256_file(cfg.overlay)


def test_a_snapshot_without_an_overlay_is_still_the_effective_config(tmp_path):
    """The plain case must not regress: same guarantee, and no stray overlay file."""
    cfg = load_config(PIPELINE_CONFIG)
    repro = _snapshot(tmp_path, cfg)

    snapshot = yaml.safe_load((repro / "config.yaml").read_text(encoding="utf-8"))
    assert snapshot == cfg.data
    assert snapshot["robustness"]["N_CAP"] == 10000
    assert (repro / "config_base.yaml").read_bytes() == PIPELINE_CONFIG.read_bytes()
    assert not (repro / "config_overlay.yaml").exists(), (
        "an overlay artifact appeared for a run that had no overlay")

    meta = json.loads((repro / "run_metadata.json").read_text(encoding="utf-8"))
    assert meta["config_overlay_path"] is None
    assert meta["config_overlay_sha256"] is None


# --------------------------------------------------------------------------- #
# B. equal effective configurations are equally identified                    #
# --------------------------------------------------------------------------- #

def test_b_identical_effective_configs_share_digest_and_run_id_component(tmp_path):
    """Identity depends on the effective parameters, not on how they were assembled.

    The strong form is asserted: an overlay that sets a parameter to the value the
    base file already holds changes nothing effective, and must therefore be
    indistinguishable from no overlay at all.
    """
    plain = load_config(PIPELINE_CONFIG)
    again = load_config(PIPELINE_CONFIG)
    assert again.sha256 == plain.sha256
    assert again.data == plain.data

    shipped_n_cap = plain.get("robustness.N_CAP")
    noop = load_config(PIPELINE_CONFIG,
                       gene_overlay=_overlay(tmp_path / "noop", shipped_n_cap))
    assert noop.data == plain.data, "the no-op overlay changed the effective tree"
    assert noop.sha256 == plain.sha256, (
        "two identical effective configurations were given different identities")

    assert make_run_id(noop.sha256, _CODE, _AT) == make_run_id(plain.sha256, _CODE, _AT)


# --------------------------------------------------------------------------- #
# C. different effective configurations are distinguishable                   #
# --------------------------------------------------------------------------- #

def test_c_different_overlays_give_different_digests_and_run_id_components(tmp_path):
    """The defect: run_id's config component came only from the base file."""
    base_digest = load_config(PIPELINE_CONFIG).sha256
    a = load_config(PIPELINE_CONFIG, gene_overlay=_overlay(tmp_path / "a", 200))
    b = load_config(PIPELINE_CONFIG, gene_overlay=_overlay(tmp_path / "b", 500))

    assert a.sha256 != b.sha256, "two different configurations share one digest"
    assert a.sha256 != base_digest and b.sha256 != base_digest

    ids = {make_run_id(d, _CODE, _AT) for d in (base_digest, a.sha256, b.sha256)}
    assert len(ids) == 3, (
        f"run_ids do not distinguish the three configurations: {sorted(ids)}")

    # The overlays differ only in a value nested under `robustness`, which is the
    # case the old composite digest could not see at all.
    assert a.get("robustness.N_CAP") == 200
    assert b.get("robustness.N_CAP") == 500


def test_c_digest_tracks_any_parameter_not_just_the_budget(tmp_path):
    """Nothing about the fix is specific to N_CAP; the digest covers the whole tree."""
    plain = load_config(PIPELINE_CONFIG)
    path = (tmp_path / "deep")
    path.mkdir()
    overlay = path / "deep.yaml"
    overlay.write_text(yaml.safe_dump(
        {"footprint_geometry": {"bbox_padding_A": 2.5}}), encoding="utf-8")

    changed = load_config(PIPELINE_CONFIG, gene_overlay=overlay)
    assert changed.get("footprint_geometry.bbox_padding_A") == 2.5
    assert changed.sha256 != plain.sha256


# --------------------------------------------------------------------------- #
# D. the production configuration is untouched                                #
# --------------------------------------------------------------------------- #

def test_d_production_config_is_unchanged_and_unwritten(tmp_path):
    """The shipped file keeps its production values and no run may write to it.

    ``config/`` is Lead-owned; the overlay mechanism exists precisely so a
    deviation never edits it. This asserts both the values and the immutability.
    """
    before = sha256_file(PIPELINE_CONFIG)

    shipped = load_config(PIPELINE_CONFIG)
    assert shipped.get("robustness.N_CAP") == 10000, (
        "config/pipeline.yaml no longer ships the production perturbation budget")
    # DECISION-B-DEFAULT-0001 raised this from 10000; the shipped value is the
    # production one and this assertion tracks it.
    assert shipped.get("permutation.B_default") == 100000
    assert_frozen_methodology(shipped)

    # A full snapshot, taken through an overlay, must leave the base file alone.
    _snapshot(tmp_path, load_config(
        PIPELINE_CONFIG, gene_overlay=_overlay(tmp_path / "ov", SENTINEL_N_CAP)))

    assert sha256_file(PIPELINE_CONFIG) == before, (
        "the production config was modified by a run")
    assert load_config(PIPELINE_CONFIG).get("robustness.N_CAP") == 10000
    assert INTEGRATION_N_CAP != 10000, (
        "the integration overlay no longer reduces anything, so the tests that "
        "rely on it are silently running the production budget")
