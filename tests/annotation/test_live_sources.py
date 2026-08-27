"""LiveAnnotationSource: graceful degradation and the literature_mechanisms
NOT_RUN contract (this task's two hard requirements).

No real network or Docker call is ever made here — the module-level
``no_network`` fixture (tests/conftest.py) hard-blocks sockets, and Docker
calls are mocked via subprocess. Every test proves a SPECIFIC failure mode
degrades to an explicit gap rather than a crash or a silent empty result
indistinguishable from a real negative.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from hotspot3d.annotation.live_sources import (
    LiveAnnotationSource,
    _parse_dssp,
)

pytestmark = pytest.mark.unit


def _source(tmp_path, allow_network=True) -> LiveAnnotationSource:
    clinvar_raw = tmp_path / "clinvar_raw.tsv"
    clinvar_raw.write_bytes(b"")          # never actually parsed in these tests
    variants = tmp_path / "variants_residue_level.tsv"
    variants.write_text("schema_version\tresidue_index\tclinvar_ids\n",
                        encoding="utf-8")
    return LiveAnnotationSource(
        allow_network=allow_network, gene="TESTGENE", uniprot_acc="P00000",
        clinvar_raw_path=clinvar_raw, variants_residue_level_path=variants,
    )


# --- literature_mechanisms is explicitly out of scope -----------------------

def test_literature_mechanisms_returns_empty_and_marks_not_run(tmp_path):
    source = _source(tmp_path)
    assert source.literature_mechanisms("TESTGENE") == []
    assert source.literature_mechanisms_not_run is True


# --- the network gate ---------------------------------------------------

def test_network_gated_methods_degrade_with_a_gap_when_network_disallowed(tmp_path):
    source = _source(tmp_path, allow_network=False)
    assert source.domains("P00000") == []
    assert source.functional_sites("P00000") == []
    assert source.conservation("P00000", [1, 2, 3]) == {}
    reasons = [g["source_attempted"] for g in source.source_gaps]
    assert reasons.count("domains") == 1
    assert reasons.count("functional_sites") == 1
    assert reasons.count("conservation") == 1
    for gap in source.source_gaps:
        assert gap["severity"] == "MAJOR"


def test_disease_associations_and_secondary_structure_do_not_need_network(tmp_path):
    """These two never call _require_network -- disease_associations reads
    already-downloaded ClinVar bytes, secondary_structure runs DSSP locally."""
    source = _source(tmp_path, allow_network=False)
    # empty clinvar_raw.tsv -> parse_variant_summary returns [] -> no crash
    assert source.disease_associations("TESTGENE") == []
    assert not any(g["source_attempted"] == "disease_associations"
                  for g in source.source_gaps)


# --- graceful degradation: DSSP -----------------------------------------

def test_secondary_structure_degrades_when_docker_is_unavailable(tmp_path):
    source = _source(tmp_path)
    structure = tmp_path / "structure.cif"
    structure.write_text("fake", encoding="utf-8")

    with patch("subprocess.run", side_effect=FileNotFoundError("docker: not found")):
        result = source.secondary_structure(structure)

    assert result == {}
    gap = next(g for g in source.source_gaps
              if g["source_attempted"] == "DSSP")
    assert gap["severity"] == "MAJOR"
    assert "docker" in gap["reason"].lower() or "not found" in gap["reason"].lower()


def test_secondary_structure_degrades_on_nonzero_dssp_exit(tmp_path):
    source = _source(tmp_path)
    structure = tmp_path / "structure.cif"
    structure.write_text("fake", encoding="utf-8")

    failed = MagicMock(returncode=1, stderr="mkdssp: unrecognized file format")
    with patch("subprocess.run", return_value=failed):
        result = source.secondary_structure(structure)

    assert result == {}
    gap = next(g for g in source.source_gaps if g["source_attempted"] == "DSSP")
    assert "unrecognized file format" in gap["reason"]


def test_secondary_structure_succeeds_when_docker_produces_real_output(tmp_path):
    """Proves the happy path end-to-end against a real (tiny) DSSP-format file,
    without actually invoking Docker: subprocess.run is mocked to write the
    fixture where the real container would have written its output."""
    source = _source(tmp_path)
    structure = tmp_path / "structure.cif"
    structure.write_text("fake", encoding="utf-8")

    dssp_text = (
        "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    "
        "N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA\n"
        "    1    1 A M              0   0  238      0, 0.0     0, 0.0     "
        "0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 114.6  -43.4   -4.6  -22.0\n"
        "    2    2 A T  H  4 S+     0   0  144      1,-0.1     2,-0.3     "
        "0, 0.0     0, 0.0   0.984 360.0-172.5  64.0 106.9  -46.9   -5.6  -20.8\n"
    )

    def fake_run(cmd, **kwargs):
        # locate the /out mount's host side (the -v flag ending in ":/out")
        host_out_dir = [a.split(":")[0] for a in cmd if a.endswith(":/out")][0]
        Path(host_out_dir, "out.dssp").write_text(dssp_text, encoding="utf-8")
        return MagicMock(returncode=0, stderr="")

    with patch("subprocess.run", side_effect=fake_run):
        result = source.secondary_structure(structure)

    assert result == {1: "-", 2: "H"}
    assert not any(g["source_attempted"] == "DSSP" for g in source.source_gaps)


# --- _parse_dssp unit coverage -------------------------------------------

def test_parse_dssp_reads_fixed_width_columns_correctly(tmp_path):
    path = tmp_path / "test.dssp"
    path.write_text(
        "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC\n"
        "    1    5 A G  E     -A   10   0A 100\n"
        "    2    6 A K              0   0  200\n",
        encoding="utf-8",
    )
    result = _parse_dssp(path)
    assert result == {5: "E", 6: "-"}


def test_parse_dssp_returns_empty_dict_for_malformed_file(tmp_path):
    path = tmp_path / "empty.dssp"
    path.write_text("not a real dssp file\n", encoding="utf-8")
    assert _parse_dssp(path) == {}


# --- disease_associations parses real-shaped ClinVar rows ------------------

def test_disease_associations_maps_phenotypes_to_the_right_residue(tmp_path, monkeypatch):
    source = _source(tmp_path)

    fake_rows = [
        {"VariationID": "1", "PhenotypeList": "Epilepsy|not provided",
         "RCVaccession": "RCV1|RCV2"},
    ]
    monkeypatch.setattr(
        "hotspot3d.annotation.live_sources.parse_variant_summary",
        lambda *a, **k: fake_rows)
    source._genomic_index.variants_residue_level_path.write_text(
        "schema_version\tresidue_index\tclinvar_ids\n1.0.0\t42\t1\n", encoding="utf-8")

    rows = source.disease_associations("TESTGENE")
    assert len(rows) == 1
    assert rows[0]["disease"] == "Epilepsy"
    assert rows[0]["residue_start"] == rows[0]["residue_end"] == 42
    assert rows[0]["reference"] == "RCV1;RCV2"
    # "not provided" is explicitly excluded, never reported as a real association
    assert not any(r["disease"] == "not provided" for r in rows)
