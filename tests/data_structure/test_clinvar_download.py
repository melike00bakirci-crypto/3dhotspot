"""Resumable, verified retrieval of the weekly ClinVar payload.

No test here touches the network: ``tests/conftest.py`` hard-blocks sockets and
every response below is scripted by :class:`_FakeNcbi`. What is asserted is the
behaviour of the transfer — the ranges requested, where a failure resumes from,
what is verified before anything is parsed, and what is recorded afterwards.

The payloads are a few kilobytes with a small ``chunk_bytes``; that exercises the
identical loop the 421 MB file takes, without pretending to move 421 MB.
"""
from __future__ import annotations

import gzip
import hashlib
import json

import pytest
import requests

from hotspot3d.data.sources import (
    CLINVAR_MD5_URL,
    CLINVAR_VARIANT_SUMMARY_URL,
    ClinVarClient,
    assert_release_metadata,
)
from hotspot3d.data.stage import run_stage_a
from hotspot3d.utils.errors import BlockedError
from hotspot3d.utils.io import read_json

pytestmark = pytest.mark.unit

URL = CLINVAR_VARIANT_SUMMARY_URL
PAYLOAD = bytes(range(256)) * 20                       # 5120 bytes, position-sensitive
LAST_MODIFIED = "Mon, 10 Aug 2026 00:00:00 GMT"


def _md5(data: bytes) -> str:
    return hashlib.md5(data, usedforsecurity=False).hexdigest()


class _FakeResponse:
    """One scripted HTTP response. ``fail_after`` resets the body mid-transfer."""

    def __init__(self, status: int, body: bytes = b"", headers: dict | None = None,
                 fail_after: int | None = None):
        self.status_code = status
        self.content = body
        self.headers = headers or {}
        self._fail_after = fail_after
        self.closed = False

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size: int = 1):
        sent = 0
        for start in range(0, len(self.content), chunk_size):
            part = self.content[start:start + chunk_size]
            if self._fail_after is not None and sent + len(part) > self._fail_after:
                head = part[:max(0, self._fail_after - sent)]
                if head:
                    yield head
                raise requests.exceptions.ChunkedEncodingError(
                    "Connection reset by peer")
            yield part
            sent += len(part)

    def close(self):
        self.closed = True


class _FakeNcbi:
    """A scripted stand-in for the NCBI host. Opens no socket, ever.

    ``range_script`` is consumed one entry per ranged request: ``None`` serves the
    range in full, an int serves that many bytes and then resets the connection,
    and ``"raise"`` fails before any byte arrives.
    """

    def __init__(self, payload: bytes = PAYLOAD, *, serves_ranges: bool = True,
                 ignores_range_header: bool = False, published_md5: str | None = None,
                 declared_length: int | None = None, head_fails: bool = False,
                 md5_body: bytes | None = None, md5_failures: int = 0,
                 range_script: list | None = None):
        self.payload = payload
        self.serves_ranges = serves_ranges
        self.ignores_range_header = ignores_range_header
        self.published_md5 = published_md5 or _md5(payload)
        self.declared_length = (len(payload) if declared_length is None
                                else declared_length)
        self.head_fails = head_fails
        self.md5_body = md5_body
        self.md5_failures = md5_failures
        self.range_script = list(range_script or [])
        self.calls: list[tuple[str, str, dict]] = []

    # -- what the client sees ------------------------------------------------
    def head(self, url, **kwargs):
        self.calls.append(("HEAD", url, {}))
        if self.head_fails:
            raise requests.ConnectionError("HEAD refused")
        headers = {"Content-Length": str(self.declared_length),
                   "Last-Modified": LAST_MODIFIED}
        if self.serves_ranges:
            headers["Accept-Ranges"] = "bytes"
        return _FakeResponse(200, b"", headers)

    def get(self, url, **kwargs):
        headers = dict(kwargs.get("headers") or {})
        self.calls.append(("GET", url, headers))
        if url.endswith(".md5"):
            return self._md5_response()
        if "Range" in headers and not self.ignores_range_header:
            return self._range_response(headers["Range"])
        return _FakeResponse(200, self.payload,
                             {"Last-Modified": LAST_MODIFIED,
                              "Content-Length": str(self.declared_length)})

    # -- internals -----------------------------------------------------------
    def _md5_response(self):
        if self.md5_failures > 0:
            self.md5_failures -= 1
            raise requests.ConnectionError("checksum fetch reset")
        body = (self.md5_body if self.md5_body is not None
                else f"{self.published_md5}  variant_summary.txt.gz\n".encode())
        return _FakeResponse(200, body)

    def _range_response(self, header: str):
        span = header.split("=", 1)[1]
        start, end = (int(part) for part in span.split("-"))
        if start >= len(self.payload):
            return _FakeResponse(416, b"")
        body = self.payload[start:end + 1]
        instruction = self.range_script.pop(0) if self.range_script else None
        if instruction == "raise":
            raise requests.ConnectionError("Connection reset by peer")
        return _FakeResponse(
            206, body,
            {"Content-Range": f"bytes {start}-{start + len(body) - 1}/"
                              f"{self.declared_length}",
             "Last-Modified": LAST_MODIFIED},
            fail_after=instruction if isinstance(instruction, int) else None)


@pytest.fixture
def ncbi(monkeypatch):
    """Install a scripted NCBI; returns a factory so each test scripts its own."""
    def _install(**kwargs) -> _FakeNcbi:
        server = _FakeNcbi(**kwargs)
        monkeypatch.setattr(requests, "get", server.get)
        monkeypatch.setattr(requests, "head", server.head)
        return server

    return _install


def _client(tmp_path, **kwargs) -> ClinVarClient:
    """A client whose scratch files land somewhere the test can inspect."""
    kwargs.setdefault("chunk_bytes", 1024)
    kwargs.setdefault("backoff_base_s", 0)              # no real sleeping in tests
    return ClinVarClient(allow_network=True, scratch_dir=str(tmp_path), **kwargs)


def _ranges(server: _FakeNcbi) -> list[str]:
    return [headers["Range"] for _, _, headers in server.calls if "Range" in headers]


# --- the happy path ---------------------------------------------------------

def test_ranged_download_assembles_and_verifies_the_payload(ncbi, tmp_path):
    server = ncbi()
    client = _client(tmp_path)

    assert client._download(URL) == PAYLOAD

    transfer = client.release_metadata()["transfer"]
    assert transfer["mode"] == "ranged_chunks"
    assert transfer["n_chunks"] == 5                    # 5120 bytes / 1024
    assert _ranges(server) == ["bytes=0-1023", "bytes=1024-2047", "bytes=2048-3071",
                               "bytes=3072-4095", "bytes=4096-5119"]
    assert transfer["integrity_verified"] is True
    assert transfer["md5_computed"] == transfer["md5_published"] == _md5(PAYLOAD)
    assert transfer["sha256"] == hashlib.sha256(PAYLOAD).hexdigest()
    assert transfer["n_bytes_assembled"] == transfer["content_length"] == len(PAYLOAD)
    assert client.release_date == LAST_MODIFIED


def test_the_published_checksum_is_fetched_per_run_and_never_assumed(ncbi, tmp_path):
    server = ncbi()
    client = _client(tmp_path)
    client._download(URL)

    assert ("GET", CLINVAR_MD5_URL, {}) in server.calls
    assert client.release_metadata()["transfer"]["md5_url"] == CLINVAR_MD5_URL
    assert CLINVAR_MD5_URL == f"{URL}.md5"
    # A second weekly release has a different digest; nothing here is pinned.
    other = ncbi(payload=PAYLOAD[:2048])
    fresh = _client(tmp_path)
    assert fresh._download(URL) == PAYLOAD[:2048]
    assert (fresh.release_metadata()["transfer"]["md5_published"]
            == other.published_md5 != _md5(PAYLOAD))


def test_the_final_chunk_is_clamped_to_the_end_of_the_file(ncbi, tmp_path):
    # 2500 bytes is not a multiple of the chunk size: the last range must stop at
    # the final byte rather than asking for bytes the file does not have.
    odd = PAYLOAD[:2500]
    server = ncbi(payload=odd)
    client = _client(tmp_path)

    assert client._download(URL) == odd
    assert _ranges(server) == ["bytes=0-1023", "bytes=1024-2047", "bytes=2048-2499"]
    assert client.release_metadata()["transfer"]["n_bytes_assembled"] == 2500


def test_the_scratch_file_is_deleted_once_the_payload_is_handed_on(ncbi, tmp_path):
    ncbi()
    _client(tmp_path)._download(URL)
    assert list(tmp_path.iterdir()) == []


# --- interrupted transfers --------------------------------------------------

def test_a_mid_transfer_failure_resumes_from_the_byte_that_landed(ncbi, tmp_path):
    #: first ranged response delivers 300 of 1024 bytes, then resets
    server = ncbi(range_script=[300])
    client = _client(tmp_path)

    assert client._download(URL) == PAYLOAD, "the resumed payload must be exact"

    assert _ranges(server)[:2] == ["bytes=0-1023", "bytes=300-1323"], (
        "the retry must resume from the byte that landed, not restart at 0")
    transfer = client.release_metadata()["transfer"]
    assert transfer["retries"] == 1
    assert transfer["failed_attempts"][0]["phase"] == "chunk_at_byte_0"
    assert "Connection reset" in transfer["failed_attempts"][0]["error"]
    assert transfer["integrity_verified"] is True


def test_repeated_failures_exhaust_the_budget_and_block(ncbi, tmp_path):
    server = ncbi(range_script=["raise"] * 6)
    client = _client(tmp_path, max_attempts=3)

    with pytest.raises(BlockedError, match="failed 3 consecutive times"):
        client._download(URL)

    transfer = client.release_metadata()["transfer"]
    assert transfer["retries"] == 3, "the recorded retry count must be the real one"
    assert [f["attempt"] for f in transfer["failed_attempts"]] == [1, 2, 3]
    assert len(_ranges(server)) == 3
    assert transfer.get("integrity_verified") is not True


def test_a_partial_download_is_never_kept_between_attempts(ncbi, tmp_path):
    ncbi(range_script=["raise"] * 6)
    client = _client(tmp_path, max_attempts=2)
    with pytest.raises(BlockedError):
        client._download(URL)

    assert list(tmp_path.iterdir()) == [], "a partial payload was left on disk"
    assert client.raw_payload("clinvar") is None


# --- verification is what makes reassembly safe ------------------------------

def test_a_checksum_mismatch_blocks_and_nothing_is_parsed(ncbi, tmp_path):
    ncbi(published_md5="0" * 32)
    client = _client(tmp_path)

    with pytest.raises(BlockedError, match="checksum mismatch"):
        client._download(URL)

    assert client.release_metadata()["transfer"]["integrity_verified"] is False
    assert client._payload is None, "a payload that failed verification was retained"
    assert list(tmp_path.iterdir()) == []


def test_a_short_payload_blocks_on_the_declared_length(ncbi, tmp_path):
    # No range support, so one GET arrives complete-looking but is 10 bytes short
    # of what the server declared. Length is checked before anything is parsed.
    ncbi(serves_ranges=False, declared_length=len(PAYLOAD) + 10)
    with pytest.raises(BlockedError, match="assembled 5120 bytes but the server "
                                           "declared 5130"):
        _client(tmp_path)._download(URL)


def test_an_unfetchable_checksum_blocks_before_the_payload_is_requested(ncbi, tmp_path):
    server = ncbi(md5_failures=99)
    client = _client(tmp_path, max_attempts=3)

    with pytest.raises(BlockedError, match="could not retrieve the published checksum"):
        client._download(URL)

    assert client.release_metadata()["transfer"]["retries"] == 3
    assert _ranges(server) == [], "the payload was fetched without a checksum to verify"


def test_an_unparseable_checksum_blocks(ncbi, tmp_path):
    ncbi(md5_body=b"<html>404 Not Found</html>\n")
    with pytest.raises(BlockedError, match="32-character MD5"):
        _client(tmp_path)._download(URL)


# --- servers that do not serve ranges ---------------------------------------

def test_a_server_without_range_support_falls_back_to_one_get(ncbi, tmp_path):
    server = ncbi(serves_ranges=False)
    client = _client(tmp_path)

    assert client._download(URL) == PAYLOAD
    assert client.release_metadata()["transfer"]["mode"] == "single_get"
    assert _ranges(server) == [], "ranges were requested from a server that refused them"


def test_a_server_that_ignores_the_range_header_is_handled_in_place(ncbi, tmp_path):
    server = ncbi(ignores_range_header=True)
    client = _client(tmp_path)

    assert client._download(URL) == PAYLOAD
    assert client.release_metadata()["transfer"]["mode"] == "single_get_range_ignored"
    payload_gets = [c for c in server.calls if c[0] == "GET" and not c[1].endswith(".md5")]
    assert len(payload_gets) == 1, "the body already arriving was requested twice"


def test_a_failed_probe_still_tries_ranges(ncbi, tmp_path):
    # A HEAD that fails leaves range support unknown; it is not evidence of refusal,
    # and the single-GET path is the one that is unreliable here.
    server = ncbi(head_fails=True)
    client = _client(tmp_path)

    assert client._download(URL) == PAYLOAD
    assert client.release_metadata()["transfer"]["mode"] == "ranged_chunks"
    assert _ranges(server)[0] == "bytes=0-1023"


def test_the_fallback_get_failing_blocks_rather_than_returning_short(monkeypatch,
                                                                    tmp_path):
    class _Broken(_FakeNcbi):
        def get(self, url, **kwargs):
            if url.endswith(".md5"):
                return super().get(url, **kwargs)
            self.calls.append(("GET", url, {}))
            raise requests.ConnectionError("reset")

    server = _Broken(serves_ranges=False)
    monkeypatch.setattr(requests, "get", server.get)
    monkeypatch.setattr(requests, "head", server.head)
    client = _client(tmp_path, max_attempts=2)

    with pytest.raises(BlockedError, match="does not serve ranged requests"):
        client._download(URL)

    assert client.release_metadata()["transfer"]["retries"] == 2
    assert list(tmp_path.iterdir()) == []


# --- the gate still comes first ---------------------------------------------

def test_allow_network_false_refuses_before_any_request(ncbi, tmp_path):
    server = ncbi()
    client = ClinVarClient(allow_network=False, scratch_dir=str(tmp_path))

    with pytest.raises(BlockedError, match="allow_network=FALSE"):
        client._download(URL)

    assert server.calls == [], "a refused run still talked to the network"
    assert list(tmp_path.iterdir()) == []


# --- integration with parsing and the raw-payload passthrough ---------------

def test_only_a_verified_payload_reaches_the_parser(ncbi, tmp_path):
    rows = ("#GeneSymbol\tAssembly\tType\tClinicalSignificance\tReviewStatus\n"
            "TESTGENE\tGRCh38\tsingle nucleotide variant\tPathogenic\t"
            "criteria provided, single submitter\n")
    gz = gzip.compress(rows.encode("utf-8"))
    ncbi(payload=gz)
    client = _client(tmp_path, chunk_bytes=16)

    records = client.fetch_missense("TESTGENE")

    assert [r["ClinicalSignificance"] for r in records] == ["Pathogenic"]
    assert client.raw_payload("clinvar") == gz, "01_INPUT_RAW must get the served bytes"
    meta = assert_release_metadata(client.release_metadata(), "ClinVar")
    assert meta["review_star_filter"] is None            # F12 untouched by any of this
    assert meta["transfer"]["integrity_verified"] is True
    assert list(tmp_path.iterdir()) == []


# --- retrieval_log.json tells the truth about retries -----------------------

def test_retrieval_log_records_the_retries_that_happened(clustered, make_ctx):
    """The log previously hardcoded failed_retrievals to [] whatever occurred."""
    class _Flaky:
        """The synthetic source, reporting a transfer that needed two retries."""

        def __init__(self, inner):
            self._inner = inner

        def __getattr__(self, name):
            return getattr(self._inner, name)

        def release_metadata(self) -> dict:
            return {**self._inner.release_metadata(), "retries": 2, "http_status": 200,
                    "failed_attempts": [
                        {"phase": "chunk_at_byte_0", "attempt": 1, "error": "reset"},
                        {"phase": "chunk_at_byte_300", "attempt": 1, "error": "reset"}]}

    ctx = make_ctx()
    run_stage_a(ctx, gene=clustered.gene, source=_Flaky(clustered.source))
    log = read_json(ctx.full_results / "01_INPUT_RAW" / "retrieval_log.json")

    assert [f["phase"] for f in log["failed_retrievals"]] == ["chunk_at_byte_0",
                                                              "chunk_at_byte_300"]
    assert {f["artifact"] for f in log["failed_retrievals"]} == {"clinvar_raw.tsv"}
    clinvar_entry = next(e for e in log["retrievals"] if e["artifact"] == "clinvar_raw.tsv")
    assert clinvar_entry["retries"] == 2
    assert json.dumps(log)                              # stays JSON-serializable
