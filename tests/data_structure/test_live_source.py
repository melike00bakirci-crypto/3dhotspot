"""The live Stage A bundle: composition, delegation and the network gate.

No test here touches the network. ``tests/conftest.py`` hard-blocks sockets for
the whole suite; the ``allow_network=True`` cases replace ``requests.get`` with a
recorder and assert that the live path is *taken* — that the request is shaped
and dispatched — never that a real response arrives.
"""
from __future__ import annotations

import hashlib
import json
import types

import pytest

from hotspot3d.data.sources import (
    ClinVarClient,
    GeneResolution,
    LiveStageASource,
    StageASource,
    UniProtClient,
    assert_release_metadata,
    assert_source_complete,
)
from hotspot3d.data.stage import _resolve_source
from hotspot3d.structure.cif import render_mmcif
from hotspot3d.structure.sources import (
    AFDB_MODEL_TEMPLATE,
    AlphaFoldClient,
    Atom,
    StructureModel,
    StructureResidue,
    assert_model_admissible,
)
from hotspot3d.utils.errors import BlockedError, EscalationRequired

pytestmark = pytest.mark.unit

GENE = "TESTGENE"
ACC = "P00001"


# --- canned payloads (shapes only; no real biological data) -----------------

CLINVAR_TSV = (
    "#GeneSymbol\tAssembly\tType\tName\tClinicalSignificance\tReviewStatus\t"
    "VariationID\n"
    f"{GENE}\tGRCh38\tsingle nucleotide variant\tNM_000000.1(TESTGENE):c.1A>G "
    "(p.Met1Val)\tPathogenic\tcriteria provided, single submitter\t1\n"
    f"OTHERGENE\tGRCh38\tsingle nucleotide variant\tNM_9.9(OTHERGENE):c.2A>G "
    "(p.Lys2Arg)\tBenign\tno assertion criteria provided\t2\n"
)
CLINVAR_BYTES = CLINVAR_TSV.encode("utf-8")
#: What NCBI publishes beside the payload; the client verifies against it.
CLINVAR_MD5 = hashlib.md5(CLINVAR_BYTES, usedforsecurity=False).hexdigest()

UNIPROT_JSON = {
    "results": [{
        "primaryAccession": ACC,
        "uniProtkbId": "TEST_HUMAN",
        "entryAudit": {"entryVersion": "7", "lastSequenceUpdateDate": "2024-01-17"},
        "sequence": {"value": "MKVG"},
        "uniProtKBCrossReferences": [{"database": "MANE-Select", "id": "NM_000000.1"}],
        "proteinDescription": {"recommendedName": {"fullName": {"value": "Test protein"}}},
    }]
}

#: The shape AFDB really serves: the canonical protein AND a second isoform,
#: each a single-fragment (-F1) model. The composed live path must return one
#: fragment here — this is the payload shape the KCNA2 run blocked on.
AFDB_PREDICTION = [
    {
        "entryId": f"AF-{ACC}-F1",
        "uniprotAccession": ACC,
        "latestVersion": 4,
        "modelCreatedDate": "2022-06-01",
        "cifUrl": f"https://alphafold.ebi.ac.uk/files/AF-{ACC}-F1-model_v4.cif",
        "uniprotStart": 1,
        "uniprotEnd": 4,
    },
    {
        "entryId": f"AF-{ACC}-2-F1",
        "uniprotAccession": f"{ACC}-2",
        "latestVersion": 4,
        "modelCreatedDate": "2022-06-01",
        "cifUrl": f"https://alphafold.ebi.ac.uk/files/AF-{ACC}-2-F1-model_v4.cif",
        "uniprotStart": 1,
        "uniprotEnd": 3,
    },
]


def _model_cif() -> str:
    """A two-residue mmCIF, rendered by the repository's own writer."""
    residues = [
        StructureResidue(residue_index=i, aa3=aa3, chain_id="A",
                         atoms={"CA": Atom("CA", float(i), 0.0, 0.0, b_factor=90.0)})
        for i, aa3 in ((1, "MET"), (2, "LYS"))
    ]
    return render_mmcif(StructureModel(accession=ACC, entry_id=f"AF-{ACC}-F1",
                                       model_version="4", db_version="2022-06-01",
                                       residues=residues, n_fragments=1))


class _Response:
    def __init__(self, content: bytes, payload=None, headers: dict | None = None):
        self.content = content
        self._payload = payload
        self.headers = {"Last-Modified": "Mon, 10 Aug 2026 00:00:00 GMT",
                        **(headers or {})}
        self.status_code = 200

    def raise_for_status(self):
        return None

    def json(self):
        return self._payload

    def iter_content(self, chunk_size: int = 1):
        for start in range(0, len(self.content), chunk_size):
            yield self.content[start:start + chunk_size]


@pytest.fixture
def dispatched(monkeypatch):
    """Replace ``requests.get``/``head`` with a recorder serving canned payloads.

    The socket block stays in force; this fixture never reaches it, which is the
    point — the assertion is on the request that was shaped and dispatched.

    The ClinVar payload arrives through the resumable client, so the checksum NCBI
    publishes is served too: nothing is parsed here that was not verified, exactly
    as in a real run. This server declares no ``Accept-Ranges``, which exercises
    the single-GET fallback; the ranged transfer itself is covered in
    ``test_clinvar_download.py``.
    """
    import requests

    calls: list[dict] = []

    def _get(url, **kwargs):
        calls.append({"url": url, "method": "GET", **kwargs})
        if url.endswith(".md5"):
            return _Response(f"{CLINVAR_MD5}  variant_summary.txt.gz\n".encode())
        if "variant_summary" in url:
            return _Response(CLINVAR_BYTES)
        if "uniprot" in url:
            return _Response(b"{}", UNIPROT_JSON)
        if "/api/prediction/" in url:
            return _Response(json.dumps(AFDB_PREDICTION).encode("utf-8"))
        if url.endswith(".cif"):
            return _Response(_model_cif().encode("utf-8"))
        raise AssertionError(f"unexpected URL dispatched: {url}")

    def _head(url, **kwargs):
        calls.append({"url": url, "method": "HEAD"})
        return _Response(b"", headers={"Content-Length": str(len(CLINVAR_BYTES))})

    monkeypatch.setattr(requests, "get", _get)
    monkeypatch.setattr(requests, "head", _head)
    return calls


# --- composition ------------------------------------------------------------

def test_live_source_composes_the_three_clients():
    source = LiveStageASource(allow_network=False)
    assert isinstance(source.variants, ClinVarClient)
    assert isinstance(source.sequences, UniProtClient)
    assert isinstance(source.structures, AlphaFoldClient)


def test_live_source_satisfies_the_stage_a_source_contract():
    source = LiveStageASource(allow_network=False)
    assert_source_complete(source)                 # raises BlockedError if not
    assert isinstance(source, StageASource)


def test_each_method_delegates_to_its_own_client(monkeypatch):
    source = LiveStageASource(allow_network=False)
    monkeypatch.setattr(source.variants, "fetch_missense", lambda gene: [{"g": gene}])
    monkeypatch.setattr(source.variants, "release_metadata", lambda: {"who": "clinvar"})
    monkeypatch.setattr(source.sequences, "resolve_gene", lambda gene: f"res:{gene}")
    monkeypatch.setattr(source.sequences, "sequence_metadata", lambda: {"who": "uniprot"})
    monkeypatch.setattr(source.structures, "fetch_model", lambda acc: f"model:{acc}")
    monkeypatch.setattr(source.structures, "model_metadata", lambda: {"who": "alphafold"})

    assert source.fetch_missense(GENE) == [{"g": GENE}]
    assert source.release_metadata() == {"who": "clinvar"}
    assert source.resolve_gene(GENE) == f"res:{GENE}"
    assert source.sequence_metadata() == {"who": "uniprot"}
    assert source.fetch_model(ACC) == f"model:{ACC}"
    assert source.model_metadata() == {"who": "alphafold"}


def test_raw_payloads_are_offered_verbatim_per_kind(monkeypatch):
    source = LiveStageASource(allow_network=False)
    monkeypatch.setattr(source.variants, "raw_payload",
                        lambda kind: b"gz" if kind == "clinvar" else None)
    monkeypatch.setattr(source.sequences, "raw_payload",
                        lambda kind: b">fa" if kind == "uniprot" else None)
    assert source.raw_payload("clinvar") == b"gz"
    assert source.raw_payload("uniprot") == b">fa"
    assert source.raw_payload("structure") is None       # carried by raw_cif_text
    assert source.raw_payload("not_a_kind") is None


def test_no_star_filter_survives_the_composition():
    source = LiveStageASource(allow_network=False)
    request = source.variants.build_request(GENE)
    assert request["review_star_filter"] is None
    assert source.release_metadata()["use_review_stars_for_inclusion"] is False


# --- the network gate: closed -----------------------------------------------

def test_allow_network_false_reaches_every_client():
    source = LiveStageASource(allow_network=False)
    assert source.allow_network is False
    assert source.variants.allow_network is False
    assert source.sequences.allow_network is False
    assert source.structures.allow_network is False


def test_allow_network_false_blocks_every_external_retrieval():
    source = LiveStageASource(allow_network=False)
    for call in (lambda: source.fetch_missense(GENE),
                 lambda: source.resolve_gene(GENE),
                 lambda: source.fetch_model(ACC)):
        with pytest.raises(BlockedError, match="allow_network=FALSE"):
            call()


def test_a_blocked_retrieval_is_never_substituted():
    source = LiveStageASource(allow_network=False)
    with pytest.raises(BlockedError, match="never substituted|never replaced|never "
                                           "substitutes"):
        source.fetch_model(ACC)
    assert source.raw_payload("clinvar") is None
    assert source.raw_payload("structure") is None


# --- the network gate: open --------------------------------------------------

def test_allow_network_true_reaches_every_client():
    source = LiveStageASource(allow_network=True)
    assert (source.variants.allow_network, source.sequences.allow_network,
            source.structures.allow_network) == (True, True, True)


def test_allow_network_true_dispatches_the_clinvar_request(dispatched):
    source = LiveStageASource(allow_network=True)
    rows = source.fetch_missense(GENE)

    payload_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/" \
                  "variant_summary.txt.gz"
    # The checksum is fetched BEFORE the payload: a reassembled file is never
    # parsed unverified, so there is no point transferring it unverifiable.
    assert [(c["method"], c["url"]) for c in dispatched] == [
        ("GET", f"{payload_url}.md5"), ("HEAD", payload_url), ("GET", payload_url)]
    assert [r["GeneSymbol"] for r in rows] == [GENE]     # gene/assembly/type only
    assert rows[0]["ReviewStatus"] == "criteria provided, single submitter"
    meta = assert_release_metadata(source.release_metadata(), "ClinVar")
    assert meta["review_star_filter"] is None
    assert meta["transfer"]["integrity_verified"] is True
    assert source.raw_payload("clinvar") == CLINVAR_BYTES


def test_allow_network_true_dispatches_the_uniprot_request(dispatched):
    source = LiveStageASource(allow_network=True)
    resolution = source.resolve_gene(GENE)

    assert dispatched[0]["url"] == "https://rest.uniprot.org/uniprotkb/search"
    assert dispatched[0]["params"]["query"] == (
        f"gene_exact:{GENE} AND organism_id:9606 AND reviewed:true")
    assert isinstance(resolution, GeneResolution)
    assert (resolution.uniprot_acc, resolution.mane_transcript) == (ACC, "NM_000000.1")
    assert_release_metadata(source.sequence_metadata(), "UniProt")
    assert source.raw_payload("uniprot").startswith(f">sp|{ACC}|".encode())


def test_the_uniprot_request_asks_for_the_entry_version_and_dates():
    """UniProtKB returns only the requested fields, and §13 requires the version.

    Without these the entry version is recorded as "NA" and the payload is
    undated, which assert_release_metadata BLOCKS before any cohort is built.
    """
    fields = UniProtClient().build_request(GENE)["params"]["fields"].split(",")
    assert {"version", "date_sequence_modified"} <= set(fields)


def test_allow_network_true_dispatches_the_alphafold_requests(dispatched):
    source = LiveStageASource(allow_network=True)
    model = source.fetch_model(ACC)

    assert dispatched[0]["url"] == f"https://alphafold.ebi.ac.uk/api/prediction/{ACC}"
    assert isinstance(model, StructureModel)

    # The CIF URL is built from the PINNED template, never from the entryId or the
    # cifUrl in the payload. This fixture declares latestVersion 4 and a v4 cifUrl
    # while the client is pinned to the shipped v6 template, so the two are
    # distinguishable: if AFDB's answer were authoritative this would request v4.
    assert dispatched[1]["url"] == (
        f"https://alphafold.ebi.ac.uk/files/AF-{ACC}-F1-model_v6.cif")
    # The payload carries a second ISOFORM entry. End to end, through the
    # composed source, that must stay one fragment and must never redirect the
    # retrieval to the isoform's model — this is the KCNA2 false block.
    assert model.n_fragments == 1
    assert model.entry_id == f"AF-{ACC}-F1"
    assert f"AF-{ACC}-2-F1" not in dispatched[1]["url"]
    assert model.model_version == "6", "the RETRIEVED version must be recorded"
    meta = source.model_metadata()
    assert meta["afdb_isoform_accessions_excluded"] == [f"{ACC}-2"]
    assert meta["afdb_latest_version"] == "4", (
        "AFDB's current latest must still be recorded beside the pinned version, "
        "so a pin that has fallen behind is visible in provenance")
    assert meta["requested_model_id"] == f"AF-{ACC}-F1-model_v6"
    assert_model_admissible(model)                       # II.13 rules still apply
    assert model.raw_cif_text is not None                # verbatim for 01_INPUT_RAW
    assert source.model_metadata()["allow_network"] is True


def test_the_model_version_comes_from_the_template_it_was_given(dispatched):
    source = LiveStageASource(allow_network=True, model_template="AF-{acc}-F1-model_vX")
    assert source.structures.model_template == "AF-{acc}-F1-model_vX"
    assert source.structures.build_cif_request(ACC)["model_id"] == f"AF-{ACC}-F1-model_vX"
    assert LiveStageASource(allow_network=True).structures.model_template == (
        AFDB_MODEL_TEMPLATE)


# --- config -> clients, end to end -------------------------------------------

def _live_source_from_config(ctx):
    return _resolve_source(types.SimpleNamespace(ctx=ctx), None)


def test_config_allow_network_flows_to_every_client(make_ctx, config):
    source = _live_source_from_config(make_ctx(synthetic=False))
    expected = bool(config.get("execution.allow_network"))
    assert (source.allow_network, source.variants.allow_network,
            source.sequences.allow_network,
            source.structures.allow_network) == (expected,) * 4


def test_config_allow_network_false_closes_the_gate_end_to_end(make_ctx, monkeypatch):
    ctx = make_ctx(synthetic=False)
    monkeypatch.setitem(ctx.config.data["execution"], "allow_network", False)
    source = _live_source_from_config(ctx)

    assert source.variants.allow_network is False
    with pytest.raises(BlockedError, match="allow_network=FALSE"):
        source.fetch_missense(GENE)


def test_config_supplies_the_frozen_model_template(make_ctx, config):
    source = _live_source_from_config(make_ctx(synthetic=False))
    assert source.structures.model_template == config.get("structure.model_template")


# --- AFDB isoforms are not fragments (II.13 detection, not II.13 policy) -----

def _prediction(entry_id: str, acc: str, end: int, version: int = 6) -> dict:
    return {"entryId": entry_id, "uniprotAccession": acc, "latestVersion": version,
            "modelCreatedDate": "2025-08-01T00:00:00Z", "uniprotStart": 1,
            "uniprotEnd": end,
            "cifUrl": f"https://alphafold.ebi.ac.uk/files/{entry_id}-model_v{version}.cif",
            "modelEntityId": entry_id}


#: The real payload shape AFDB serves for P16389: the canonical protein plus a
#: shorter isoform, each a SINGLE-fragment model of a different molecule.
TWO_ISOFORM_PAYLOAD = [_prediction("AF-P16389-F1", "P16389", 499),
                       _prediction("AF-P16389-2-F1", "P16389-2", 356)]

#: A genuinely fragmented entry: one accession, several -F<n> fragments.
TWO_FRAGMENT_PAYLOAD = [_prediction("AF-Q8WZ42-F1", "Q8WZ42", 1400),
                        _prediction("AF-Q8WZ42-F2", "Q8WZ42", 2800)]


def test_isoform_entries_are_not_counted_as_fragments():
    """The regression: two isoforms are not a multi-fragment entry.

    AFDB returns one entry per isoform, so counting payload entries reported two
    fragments for a 499-residue protein — which AFDB cannot fragment, since it
    fragments only above ~2700 residues. The II.13 block fired on a condition that
    did not exist, and would have done so for any protein with a second isoform.
    """
    meta = AlphaFoldClient.parse_prediction_payload(TWO_ISOFORM_PAYLOAD, "P16389")
    assert meta["n_fragments"] == 1, "isoforms were counted as fragments again"
    assert meta["entry_id"] == "AF-P16389-F1"
    assert meta["uniprot_sequence_length"] == 499, "the isoform's length was taken"
    # The isoform is recorded rather than silently dropped.
    assert meta["afdb_isoform_accessions_excluded"] == ["P16389-2"]
    assert meta["afdb_accessions_returned"] == ["P16389", "P16389-2"]


def test_a_genuinely_fragmented_entry_still_blocks():
    """II.13 is untouched: real fragments of ONE accession still stop the run."""
    meta = AlphaFoldClient.parse_prediction_payload(TWO_FRAGMENT_PAYLOAD, "Q8WZ42")
    assert meta["n_fragments"] == 2
    assert meta["afdb_isoform_accessions_excluded"] == []

    model = StructureModel(accession="Q8WZ42", entry_id="AF-Q8WZ42-F1",
                           model_version="6", db_version="2025-08-01",
                           residues=[], n_fragments=2)
    with pytest.raises(BlockedError, match="MULTI_FRAGMENT_AFDB_ENTRY"):
        assert_model_admissible(model)


def test_the_canonical_entry_is_chosen_regardless_of_payload_order():
    """Fragment selection must not depend on the order AFDB happened to serve."""
    reversed_payload = list(reversed(TWO_ISOFORM_PAYLOAD))
    meta = AlphaFoldClient.parse_prediction_payload(reversed_payload, "P16389")
    assert meta["entry_id"] == "AF-P16389-F1"
    assert meta["n_fragments"] == 1

    frags = AlphaFoldClient.parse_prediction_payload(
        list(reversed(TWO_FRAGMENT_PAYLOAD)), "Q8WZ42")
    assert frags["entry_id"] == "AF-Q8WZ42-F1", "F2 was selected over F1"
    assert frags["afdb_fragment_entry_ids"] == ["AF-Q8WZ42-F1", "AF-Q8WZ42-F2"]


def test_a_payload_without_the_canonical_accession_escalates():
    """Only isoforms and no canonical entry is an escalation, never a guess."""
    with pytest.raises(EscalationRequired, match="none for that accession"):
        AlphaFoldClient.parse_prediction_payload(
            [_prediction("AF-P16389-2-F1", "P16389-2", 356)], "P16389")


def test_an_isoform_is_recognised_from_the_entry_id_alone():
    """The isoform test must not depend on uniprotAccession being spelled out.

    The entry id carries the same fact (AF-P16389-2-F1 -> P16389-2). Treating an
    entry that names no accession as canonical would count this isoform as a
    second fragment again, on a protein AFDB cannot fragment.
    """
    unlabelled = [{"entryId": "AF-P16389-F1", "latestVersion": 6, "uniprotEnd": 499},
                  {"entryId": "AF-P16389-2-F1", "latestVersion": 6, "uniprotEnd": 356}]
    meta = AlphaFoldClient.parse_prediction_payload(unlabelled, "P16389")
    assert meta["n_fragments"] == 1
    assert meta["entry_id"] == "AF-P16389-F1"
    assert meta["afdb_isoform_accessions_excluded"] == ["P16389-2"]

    # ...while real fragments of one accession are still counted from -F<n>.
    fragments = [{"entryId": "AF-Q8WZ42-F1"}, {"entryId": "AF-Q8WZ42-F2"}]
    assert AlphaFoldClient.parse_prediction_payload(
        fragments, "Q8WZ42")["n_fragments"] == 2

    # An entry that names no accession anywhere is still taken as the requested
    # one: exclusion is never inferred from the absence of all evidence.
    assert AlphaFoldClient.parse_prediction_payload(
        [{"entryId": "not-an-afdb-id"}], "P16389")["n_fragments"] == 1
