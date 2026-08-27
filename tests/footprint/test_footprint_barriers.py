"""Information barriers: no annotation reads, no network, no upstream writes."""
from __future__ import annotations

import pytest

from hotspot3d.footprint import inputs as inputs_mod
from hotspot3d.utils.errors import LeakageError
from hotspot3d.utils.runctx import AGENT_READ_SCOPE, AGENT_WRITE_SCOPE

pytestmark = pytest.mark.unit

AGENT = "footprint-robustness"


def test_annotation_is_outside_the_read_scope():
    assert "10_ANNOTATION" not in AGENT_READ_SCOPE[AGENT]
    assert AGENT_WRITE_SCOPE[AGENT] == {"07_FOOTPRINT_RADIUS", "08_FINAL_FOOTPRINT",
                                        "09_ROBUSTNESS"}


def test_opening_an_annotation_path_is_a_leakage_event(phase_c_run):
    ctx, _, _ = phase_c_run
    with pytest.raises(LeakageError):
        inputs_mod._guard_path(ctx.full_results / "10_ANNOTATION" / "anything.tsv")
    with pytest.raises(Exception):
        ctx.assert_may_read(AGENT, "10_ANNOTATION")


def test_no_upstream_stage_directory_is_written(phase_c_run):
    ctx, _, _ = phase_c_run
    for stage in ("01_INPUT_RAW", "02_CLINVAR", "03_STRUCTURE_QC",
                  "04_GLOBAL_CLUSTERING", "05_HOTSPOT_RADIUS", "06_FINAL_HOTSPOTS"):
        with pytest.raises(Exception):
            ctx.assert_may_write(AGENT, stage)


def test_the_package_holds_no_network_client():
    """The absence of network tooling is the mechanical form of the biology barrier."""
    import pathlib

    roots = [pathlib.Path(inputs_mod.__file__).parent,
             pathlib.Path(inputs_mod.__file__).parents[1] / "robustness"]
    banned = ("import requests", "from requests", "urllib.request", "httpx",
              "WebFetch", "WebSearch", "socket.socket")
    for root in roots:
        for path in sorted(root.rglob("*.py")):
            text = path.read_text(encoding="utf-8")
            for token in banned:
                assert token not in text, f"{path.name} references {token}"


def test_no_biological_or_mechanism_vocabulary_reaches_the_decision():
    import pathlib

    roots = [pathlib.Path(inputs_mod.__file__).parent,
             pathlib.Path(inputs_mod.__file__).parents[1] / "robustness"]
    banned = ("GOF", "LOF", "dbNSFP", "phyloP", "DSSP", "UniProt_domain")
    for root in roots:
        for path in sorted(root.rglob("*.py")):
            for line in path.read_text(encoding="utf-8").splitlines():
                if line.strip().startswith("#") or '"""' in line:
                    continue
                for token in banned:
                    assert token not in line, f"{path.name}: {line.strip()}"
