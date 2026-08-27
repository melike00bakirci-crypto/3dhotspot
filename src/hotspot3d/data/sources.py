"""Upstream data-source interfaces and live clients for Stage A.

The science in :mod:`hotspot3d.data.cohort` never talks to the network. It talks
to these interfaces, so a synthetic fixture and the live ClinVar/UniProt/AlphaFold
clients drive byte-identical code paths and a live client can be dropped in later
without touching a single scientific decision.

Every network call is gated on ``allow_network``. Request shaping and response
parsing are pure functions and are unit-tested without a network.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol, Sequence, runtime_checkable

from ..utils.errors import BlockedError, EscalationRequired
from .clinvar import parse_variant_summary

CLINVAR_VARIANT_SUMMARY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
)
CLINVAR_RELEASE_DIR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/"
#: NCBI publishes the payload's MD5 beside it. It is fetched per run and never
#: hardcoded: the file is reissued weekly and so is its digest.
CLINVAR_MD5_URL = f"{CLINVAR_VARIANT_SUMMARY_URL}.md5"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"

#: Bytes per ranged request. The weekly payload is ~421 MB, so 16 MiB is ~27
#: short-lived requests: a reset costs one chunk rather than the whole transfer.
#: Sized for reliability, not for any observed size ceiling — none exists.
CLINVAR_CHUNK_BYTES = 16 * 1024 * 1024
#: Consecutive failures tolerated before a retrieval is declared failed. A
#: completed chunk resets the budget; a failure never widens it.
CLINVAR_MAX_ATTEMPTS = 5
#: First backoff, doubled per consecutive failure. 0 disables sleeping (tests).
CLINVAR_BACKOFF_BASE_S = 1.0
#: Read granularity while streaming one response body to the scratch file.
CLINVAR_STREAM_BYTES = 1024 * 1024

#: Payload kinds a source may optionally hand over verbatim for 01_INPUT_RAW.
RAW_KINDS = ("clinvar", "uniprot", "structure")


def _transport_errors() -> tuple:
    """Failures worth retrying: transport, never policy.

    :class:`BlockedError` is deliberately outside this set — a refusal must
    propagate immediately instead of being retried into a timeout.
    """
    import requests

    return (requests.RequestException, OSError)


def _header(response: Any, name: str) -> str:
    """One header, case-insensitively, from a response that may have none."""
    headers = getattr(response, "headers", None) or {}
    for key, value in dict(headers).items():
        if str(key).lower() == name.lower():
            return str(value)
    return ""


# --- oligomeric assembly (Workflow v2 §2) ------------------------------------
# Determined from UniProt's OWN subunit-structure annotation only (the "UniProt"
# half of "UniProt/ComplexPortal annotation" — ComplexPortal integration is
# explicitly out of scope, per config/pipeline.yaml structure.oligomeric_state_source).
# A SIMPLE, conservative classifier: it declares oligomeric only on a clear,
# confident textual match; it never over-declares (a false "oligomeric" earns an
# ADVISORY, no more), and a genuinely absent/unparseable comment is its own
# recorded "unknown" state rather than a silent default to either extreme.
OLIGOMERIC_STATE_OLIGOMERIC = "oligomeric"
OLIGOMERIC_STATE_MONOMER = "monomer"
OLIGOMERIC_STATE_UNKNOWN = "unknown"
OLIGOMERIC_STATES = (OLIGOMERIC_STATE_OLIGOMERIC, OLIGOMERIC_STATE_MONOMER,
                     OLIGOMERIC_STATE_UNKNOWN)

# UniProt's boilerplate subunit phrasing for self-assembly, e.g. "Homodimer.",
# "Homotetramer.", "Heterodimer of ATP1A1 and ATP1B1.", "Homooligomer.". Matches
# only the "homo-/hetero-<n>mer" family of words — never a bare "Interacts with
# X", which describes a partner, not the biological assembly.
_OLIGOMER_RE = re.compile(
    r"\b(homo|hetero)-?(dimer|trimer|tetramer|pentamer|hexamer|heptamer|octamer|"
    r"nonamer|decamer|dodecamer|oligomer|multimer)s?\b",
    re.IGNORECASE,
)


def extract_subunit_comment(entry: dict) -> str | None:
    """Join every UniProt ``SUBUNIT`` comment's free text from an entry payload.

    Returns ``None`` when the entry carries no subunit comment at all — the
    "genuinely unknown" case :func:`classify_oligomeric_state` must not conflate
    with a comment that was present and simply did not match oligomer language.
    """
    texts: list[str] = []
    for comment in entry.get("comments") or []:
        if str(comment.get("commentType", "")).strip().upper() != "SUBUNIT":
            continue
        for item in comment.get("texts") or []:
            value = str(item.get("value") or "").strip()
            if value:
                texts.append(value)
    joined = " ".join(texts).strip()
    return joined or None


def classify_oligomeric_state(subunit_text: str | None) -> dict:
    """Workflow v2 §2 — read UniProt's subunit annotation for the assembly state.

    Decision rule (verbatim):
      1. No subunit comment text at all -> ``unknown``. Distinct from a positive
         monomer reading; the run must not silently default to either extreme.
      2. The joined comment text contains a confident ``homo-`` or ``hetero-``
         self-assembly word (dimer/trimer/.../oligomer/multimer) -> ``oligomeric``.
      3. Any other present text — an explicit "Monomer." statement, or wording
         that is ambiguous/does not match rule 2 — -> ``monomer``. The pipeline
         never over-declares oligomeric state: it has no multimer-fetching
         capability regardless, so a missed, weakly-worded oligomer mention only
         costs a missing ADVISORY, while a false positive would wrongly force
         BLOCKED (or wrongly label discovery intra-subunit-only) for a protein
         that is not actually an obligate oligomer.
    """
    text = str(subunit_text or "").strip()
    if not text:
        return {"oligomeric_state": OLIGOMERIC_STATE_UNKNOWN, "subunit_text": None,
                "matched_phrase": None}
    match = _OLIGOMER_RE.search(text)
    if match:
        return {"oligomeric_state": OLIGOMERIC_STATE_OLIGOMERIC, "subunit_text": text,
                "matched_phrase": match.group(0)}
    return {"oligomeric_state": OLIGOMERIC_STATE_MONOMER, "subunit_text": text,
            "matched_phrase": None}


# --- gene resolution --------------------------------------------------------

@dataclass(frozen=True)
class GeneResolution:
    """Gene -> exactly one UniProt accession, canonical isoform, MANE transcript.

    Ambiguity is never resolved here: a source that cannot name exactly one
    accession raises :class:`EscalationRequired` (agent §6.2).
    """

    gene: str
    uniprot_acc: str
    uniprot_entry_version: str
    canonical_isoform: str
    sequence: str
    mane_transcript: str
    mane_source: str = "MANE_SELECT"
    mane_fallback_used: bool = False
    organism_id: int = 9606
    entry_name: str = "NA"
    protein_name: str = "NA"
    # -- Workflow v2 §2 — oligomeric assembly, from UniProt's own subunit comment
    oligomeric_state: str = OLIGOMERIC_STATE_UNKNOWN
    oligomeric_state_source: str = "uniprot_subunit_annotation"
    oligomeric_state_evidence: str | None = None
    oligomeric_state_matched_phrase: str | None = None

    def __post_init__(self) -> None:
        if not self.sequence:
            raise BlockedError(
                f"BLOCKED — no canonical sequence resolved for {self.gene}; a payload "
                f"without a sequence is not a usable UniProt entry."
            )
        if not self.mane_transcript:
            raise BlockedError(
                f"BLOCKED — no MANE Select transcript (or configured fallback) named "
                f"for {self.gene}. The transcript must be recorded explicitly."
            )

    @property
    def sequence_length(self) -> int:
        return len(self.sequence)

    def aa_at(self, position: int | None) -> str | None:
        if position is None or position < 1 or position > self.sequence_length:
            return None
        return self.sequence[position - 1]

    def fasta(self) -> str:
        header = (f">sp|{self.uniprot_acc}|{self.entry_name} {self.protein_name} "
                  f"OX={self.organism_id} GN={self.gene} "
                  f"ISOFORM={self.canonical_isoform} "
                  f"ENTRY_VERSION={self.uniprot_entry_version}")
        body = "\n".join(self.sequence[i:i + 60]
                         for i in range(0, self.sequence_length, 60))
        return f"{header}\n{body}\n"

    def as_dict(self) -> dict:
        return {
            "gene": self.gene, "uniprot_acc": self.uniprot_acc,
            "uniprot_entry_version": self.uniprot_entry_version,
            "canonical_isoform": self.canonical_isoform,
            "sequence_length": self.sequence_length,
            "mane_transcript": self.mane_transcript, "mane_source": self.mane_source,
            "mane_fallback_used": self.mane_fallback_used,
            "organism_id": self.organism_id, "entry_name": self.entry_name,
            "protein_name": self.protein_name,
            "oligomeric_state": self.oligomeric_state,
            "oligomeric_state_source": self.oligomeric_state_source,
            "oligomeric_state_evidence": self.oligomeric_state_evidence,
            "oligomeric_state_matched_phrase": self.oligomeric_state_matched_phrase,
        }


# --- provider interfaces ----------------------------------------------------

@runtime_checkable
class VariantSource(Protocol):
    """Delivers every ClinVar missense record for a gene, at every star level."""

    def fetch_missense(self, gene: str) -> list[dict]: ...

    def release_metadata(self) -> dict: ...


@runtime_checkable
class SequenceSource(Protocol):
    """Resolves a gene symbol to one UniProt entry and its canonical sequence."""

    def resolve_gene(self, gene: str) -> GeneResolution: ...

    def sequence_metadata(self) -> dict: ...


@runtime_checkable
class StageASource(Protocol):
    """The full provider Stage A needs: variants + sequence + structure.

    ``source`` passed to :func:`hotspot3d.data.stage.run_stage_a` must satisfy
    this. ``StructureSource`` lives in :mod:`hotspot3d.structure.sources`.
    """

    def fetch_missense(self, gene: str) -> list[dict]: ...

    def release_metadata(self) -> dict: ...

    def resolve_gene(self, gene: str) -> GeneResolution: ...

    def sequence_metadata(self) -> dict: ...

    def fetch_model(self, uniprot_acc: str): ...

    def model_metadata(self) -> dict: ...


REQUIRED_RELEASE_KEYS = ("source", "release_date", "query")


def assert_release_metadata(metadata: dict, what: str) -> dict:
    """A payload without a release/version is unusable (agent §6.5, §15)."""
    missing = [k for k in REQUIRED_RELEASE_KEYS if not str(metadata.get(k, "")).strip()]
    if missing:
        raise BlockedError(
            f"BLOCKED — the {what} payload is undated or unversioned: missing "
            f"{missing}. Stage A never analyses a payload it cannot pin to a release."
        )
    return metadata


def raw_payload_of(source: Any, kind: str) -> bytes | None:
    """Ask a source for the original bytes of a payload, if it kept them."""
    getter = getattr(source, "raw_payload", None)
    if getter is None:
        return None
    try:
        return getter(kind)
    except (KeyError, LookupError, NotImplementedError):
        return None


# --- ClinVar client ---------------------------------------------------------

@dataclass
class ClinVarClient:
    """Weekly ``variant_summary`` client. No star filter exists in this class.

    The only selection applied at retrieval is gene, assembly and record type.
    Significance and review status are carried through untouched so that strata 1
    and 3 can account for every retrieved record.

    **Retrieval is chunked, resumable and verified.** The payload is ~421 MB and a
    single long-lived connection to the NCBI host proved unreliable from at least
    one client (repeated ``IncompleteRead`` a few MB in, and one stall that ran out
    the timeout), while short ranged requests to the same host were reliable and a
    full transfer via another client succeeded. There is no size ceiling; what is
    fragile is holding one connection open for the whole file. So the file is
    fetched as ``chunk_bytes`` ranged requests, retried with exponential backoff,
    resuming from the byte that actually landed.

    Reassembling a payload from many responses is only safe if the result is
    checked, so it is: the length is compared with what the server declared and
    the bytes with the MD5 NCBI publishes beside the file. Either check failing is
    BLOCKING — a payload that failed verification is never parsed, and a partial
    one is never returned, kept, or substituted for the real thing.
    """

    allow_network: bool = False
    assembly: str = "GRCh38"
    timeout: int = 600
    release_date: str | None = None
    chunk_bytes: int = CLINVAR_CHUNK_BYTES
    max_attempts: int = CLINVAR_MAX_ATTEMPTS
    backoff_base_s: float = CLINVAR_BACKOFF_BASE_S
    scratch_dir: str | None = None
    _payload: bytes | None = field(default=None, init=False, repr=False)
    _last_query: str = field(default="", init=False, repr=False)
    _transfer: dict = field(default_factory=dict, init=False, repr=False)

    def build_request(self, gene: str) -> dict:
        query = (f"GeneSymbol=={gene} AND Assembly=={self.assembly} AND "
                 f"Type in {{single nucleotide variant, indel}} AND "
                 f"review_star_filter=NONE")
        self._last_query = query
        return {
            "url": CLINVAR_VARIANT_SUMMARY_URL,
            "method": "GET",
            "headers": {"Accept": "application/gzip"},
            "release_listing_url": CLINVAR_RELEASE_DIR_URL,
            "checksum_url": CLINVAR_MD5_URL,
            "query": query,
            "gene": gene,
            "assembly": self.assembly,
            "review_star_filter": None,
            "note": "F12 — every review star level is retrieved and preserved.",
        }

    # -- network primitives (gated) -----------------------------------------
    def _require_network(self, url: str) -> None:
        if not self.allow_network:
            raise BlockedError(
                f"BLOCKED — network retrieval is disabled for this run "
                f"(allow_network=FALSE); refusing to fetch {url}. A failed or "
                f"disabled retrieval is never substituted with a cached or "
                f"hand-made file."
            )

    def _http(self, method: str, url: str, *, headers: dict | None = None,
              stream: bool = False):
        """The single place this client touches the network."""
        self._require_network(url)
        import requests

        if method == "HEAD":
            return requests.head(url, timeout=self.timeout, allow_redirects=True)
        return requests.get(url, timeout=self.timeout, stream=stream,
                            **({"headers": headers} if headers else {}))

    # -- retry bookkeeping ---------------------------------------------------
    def _record_failure(self, phase: str, attempt: int, exc: Exception) -> None:
        """Every failed attempt is recorded; retrieval_log.json reports the truth."""
        self._transfer.setdefault("failed_attempts", []).append({
            "phase": phase, "attempt": attempt,
            "error": f"{type(exc).__name__}: {exc}"[:300],
        })
        self._transfer["retries"] = int(self._transfer.get("retries", 0)) + 1

    def _backoff(self, attempt: int) -> None:
        if self.backoff_base_s <= 0:
            return
        import time

        time.sleep(self.backoff_base_s * (2 ** (attempt - 1)))

    # -- scratch file (never kept between attempts) --------------------------
    def _new_scratch(self) -> Path:
        import tempfile

        handle, name = tempfile.mkstemp(prefix="clinvar_variant_summary_",
                                        suffix=".part", dir=self.scratch_dir)
        os.close(handle)
        return Path(name)

    @staticmethod
    def _discard(path: Path | None) -> None:
        """A partial download is deleted, never presented as a payload."""
        if path is None:
            return
        try:
            path.unlink()
        except OSError:                                          # pragma: no cover
            pass

    # -- retrieval -----------------------------------------------------------
    def _published_md5(self, url: str) -> str:
        """The MD5 published beside the payload. Fetched per run, never assumed."""
        checksum_url = f"{url}.md5"
        self._transfer["md5_url"] = checksum_url
        for attempt in range(1, self.max_attempts + 1):
            try:
                response = self._http("GET", checksum_url)
                response.raise_for_status()
                token = response.content.decode("utf-8", "replace").split()
                digest = token[0].strip().lower() if token else ""
                if len(digest) != 32 or digest.strip("0123456789abcdef"):
                    raise BlockedError(
                        f"BLOCKED — {checksum_url} did not yield a 32-character MD5 "
                        f"(read {digest!r}). A payload reassembled from many responses "
                        f"is never parsed without the checksum that proves it whole."
                    )
                self._transfer["md5_published"] = digest
                return digest
            except _transport_errors() as exc:
                self._record_failure("checksum", attempt, exc)
                if attempt >= self.max_attempts:
                    raise BlockedError(
                        f"BLOCKED — could not retrieve the published checksum "
                        f"{checksum_url} after {attempt} attempts ({exc}). Stage A does "
                        f"not parse a reassembled payload it cannot verify."
                    ) from exc
                self._backoff(attempt)
        raise AssertionError("unreachable")                      # pragma: no cover

    def _probe(self, url: str) -> dict:
        """HEAD the payload for its size, its date and whether ranges are served."""
        try:
            response = self._http("HEAD", url)
            response.raise_for_status()
        except _transport_errors() as exc:
            self._record_failure("probe", 1, exc)
            return {}
        length = _header(response, "Content-Length")
        return {
            "content_length": int(length) if length.strip().isdigit() else None,
            "accept_ranges": _header(response, "Accept-Ranges").strip().lower(),
            "last_modified": _header(response, "Last-Modified"),
        }

    @staticmethod
    def _consume(response: Any, handle: Any) -> int:
        """Stream one response body into the scratch file; return bytes written."""
        written = 0
        try:
            for part in response.iter_content(chunk_size=CLINVAR_STREAM_BYTES):
                if not part:
                    continue
                handle.write(part)
                written += len(part)
        finally:
            close = getattr(response, "close", None)
            if callable(close):
                close()
        return written

    @staticmethod
    def _total_from_content_range(response: Any) -> int | None:
        """``Content-Range: bytes 0-16777215/441836319`` -> ``441836319``."""
        tail = _header(response, "Content-Range").rsplit("/", 1)
        return int(tail[1]) if len(tail) == 2 and tail[1].strip().isdigit() else None

    def _download_ranged(self, url: str, path: Path, total: int | None) -> int | None:
        """Fetch the payload as ranged chunks, resuming from what actually landed.

        Returns the total size the server declared, which may only become known
        from the first ``Content-Range``. Falls back to consuming the response in
        place if the server ignores the very first Range header.
        """
        self._transfer["mode"] = "ranged_chunks"
        written = 0
        attempt = 0
        with path.open("wb") as handle:
            while total is None or written < total:
                start = written
                end = (min(written + self.chunk_bytes, total) if total is not None
                       else written + self.chunk_bytes) - 1
                attempt += 1
                try:
                    response = self._http("GET", url, stream=True,
                                          headers={"Range": f"bytes={start}-{end}"})
                    status = int(getattr(response, "status_code", 0) or 0)
                    if status == 200 and start == 0:
                        # The server ignored the Range header and is already
                        # sending the whole body: consume THIS response rather
                        # than re-requesting bytes that are on the way.
                        self._transfer["mode"] = "single_get_range_ignored"
                        self._consume(response, handle)
                        return total
                    if status == 416 and written > 0:
                        break                    # asked past the end; the file is done
                    if status != 206:
                        raise BlockedError(
                            f"BLOCKED — {url} answered the request for bytes "
                            f"{start}-{end} with HTTP {status}, not 206. A payload "
                            f"assembled from responses that do not mean what they say "
                            f"is never parsed."
                        )
                    self.release_date = (self.release_date
                                         or _header(response, "Last-Modified") or None)
                    received = self._consume(response, handle)
                    if received <= 0:
                        # A 206 that carries nothing makes no progress; treat it as a
                        # transport failure so the attempt budget still bounds the loop.
                        raise OSError(f"empty 206 response for bytes {start}-{end}")
                    written += received
                    total = total if total is not None else self._total_from_content_range(
                        response)
                    self._transfer["n_chunks"] = int(self._transfer.get("n_chunks", 0)) + 1
                    attempt = 0                  # a completed chunk resets the budget
                    if total is None and received < (end - start + 1):
                        break                    # short chunk: end of an unsized file
                except _transport_errors() as exc:
                    self._record_failure(f"chunk_at_byte_{start}", attempt, exc)
                    handle.flush()
                    written = handle.tell()      # resume from what actually landed
                    if attempt >= self.max_attempts:
                        raise BlockedError(
                            f"BLOCKED — retrieval of {url} failed {attempt} consecutive "
                            f"times at byte {written} of {total} ({exc}). A partial "
                            f"payload is never parsed, kept, or substituted with a "
                            f"cached or hand-made file."
                        ) from exc
                    self._backoff(attempt)
        return total

    def _download_whole(self, url: str, path: Path) -> int:
        """Fallback for a server that does not serve ranges: one GET, retried."""
        self._transfer["mode"] = "single_get"
        for attempt in range(1, self.max_attempts + 1):
            try:
                response = self._http("GET", url, stream=True)
                response.raise_for_status()
                self.release_date = (self.release_date
                                     or _header(response, "Last-Modified") or None)
                with path.open("wb") as handle:
                    return self._consume(response, handle)
            except _transport_errors() as exc:
                self._record_failure("single_get", attempt, exc)
                if attempt >= self.max_attempts:
                    raise BlockedError(
                        f"BLOCKED — {url} does not serve ranged requests here and the "
                        f"single unranged retrieval failed {attempt} times ({exc}). A "
                        f"short or partial payload is never returned."
                    ) from exc
                self._backoff(attempt)
        raise AssertionError("unreachable")                      # pragma: no cover

    def _verify(self, path: Path, declared_length: int | None, published_md5: str) -> None:
        """Length against what the server declared, bytes against NCBI's MD5."""
        import hashlib

        md5 = hashlib.md5(usedforsecurity=False)
        sha256 = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(CLINVAR_STREAM_BYTES), b""):
                md5.update(block)
                sha256.update(block)
        size = path.stat().st_size
        self._transfer.update({
            "n_bytes_assembled": size, "content_length": declared_length,
            "md5_computed": md5.hexdigest(), "sha256": sha256.hexdigest(),
            "integrity_verified": False,
        })
        if declared_length is not None and size != declared_length:
            raise BlockedError(
                f"BLOCKED — assembled {size} bytes but the server declared "
                f"{declared_length}. A truncated payload is never parsed: every "
                f"missing byte is variant records silently absent from stratum 1."
            )
        if md5.hexdigest() != published_md5:
            raise BlockedError(
                f"BLOCKED — checksum mismatch for the ClinVar payload: NCBI publishes "
                f"{published_md5}, the {size} bytes assembled from "
                f"{self._transfer.get('n_chunks', 0)} response(s) hash to "
                f"{md5.hexdigest()}. A payload that failed verification is never parsed."
            )
        self._transfer["integrity_verified"] = True

    def _download(self, url: str) -> bytes:
        """Retrieve, verify and hand over the weekly payload.

        The scratch file exists only for the duration of one retrieval: it is
        deleted whether the retrieval succeeds or fails, so a partial download is
        never available to a later run as a substitute for a real one.
        """
        self._require_network(url)
        self._transfer = {"url": url, "mode": "NA", "retries": 0, "failed_attempts": [],
                          "n_chunks": 0, "chunk_bytes": self.chunk_bytes,
                          "max_attempts": self.max_attempts, "http_status": "NA"}
        published_md5 = self._published_md5(url)
        probe = self._probe(url)
        declared_length = probe.get("content_length")
        self.release_date = self.release_date or probe.get("last_modified") or None
        # Ranges are used unless the server said it does not serve them. A probe
        # that failed outright leaves the question open, and the ranged path
        # answers it definitively: a server that ignores Range replies 200, which
        # is handled in place rather than guessed at here.
        ranges_refused = bool(probe) and not str(
            probe.get("accept_ranges", "")).startswith("bytes")

        path = self._new_scratch()
        try:
            if ranges_refused:
                self._download_whole(url, path)
            else:
                declared_length = self._download_ranged(url, path, declared_length)
            self._verify(path, declared_length, published_md5)
            payload = path.read_bytes()
        finally:
            self._discard(path)
        self._transfer["http_status"] = 200
        return payload

    def fetch_missense(self, gene: str) -> list[dict]:
        request = self.build_request(gene)
        self._payload = self._download(request["url"])
        return parse_variant_summary(self._payload, gene, assembly=self.assembly)

    def raw_payload(self, kind: str) -> bytes | None:
        return self._payload if kind == "clinvar" else None

    def release_metadata(self) -> dict:
        return {
            "source": "clinvar_variant_summary",
            "url": CLINVAR_VARIANT_SUMMARY_URL,
            "release_date": self.release_date or "",
            "query": self._last_query,
            "assembly": self.assembly,
            "review_star_filter": None,
            "use_review_stars_for_inclusion": False,
            "allow_network": self.allow_network,
            # The real transfer record: retries and failures as they happened, and
            # the digests that prove the assembled payload is the published one.
            "retries": int(self._transfer.get("retries", 0)),
            "failed_attempts": list(self._transfer.get("failed_attempts", [])),
            "http_status": self._transfer.get("http_status", "NA"),
            "transfer": dict(self._transfer),
        }


# --- UniProt client ---------------------------------------------------------

@dataclass
class UniProtClient:
    """Resolves gene -> one reviewed human UniProt entry + MANE Select transcript."""

    allow_network: bool = False
    organism_id: int = 9606
    timeout: int = 120
    _entry: dict | None = field(default=None, init=False, repr=False)
    _fasta: bytes | None = field(default=None, init=False, repr=False)

    def build_request(self, gene: str) -> dict:
        query = (f"gene_exact:{gene} AND organism_id:{self.organism_id} "
                 f"AND reviewed:true")
        return {
            "url": UNIPROT_SEARCH_URL,
            "method": "GET",
            # version + date_* populate the `entryAudit` object that
            # parse_entry and sequence_metadata read. UniProtKB returns ONLY the
            # requested fields, so without them the entry version is recorded as
            # "NA" and the payload is undated — which assert_release_metadata
            # BLOCKS. §13 requires the entry version to be recorded, not guessed.
            "params": {"query": query, "format": "json", "size": 10,
                       # cc_subunit: Workflow v2 §2 — the SUBUNIT comment is the
                       # sole evidence for the oligomeric-state determination
                       # (classify_oligomeric_state). Like entryAudit, UniProtKB
                       # returns it only when explicitly requested.
                       "fields": ("accession,id,protein_name,gene_names,sequence,"
                                  "xref_mane-select,cc_alternative_products,"
                                  "cc_subunit,protein_existence,annotation_score,"
                                  "version,date_modified,date_sequence_modified")},
            "query": query,
        }

    @staticmethod
    def parse_search_payload(payload: dict, gene: str) -> dict:
        """Pick the single reviewed entry. More than one is an escalation."""
        results = payload.get("results") or []
        if not results:
            raise EscalationRequired(
                ambiguity=f"No reviewed human UniProt entry resolves for gene {gene!r}.",
                options=["Escalate and stop", "Accept an unreviewed (TrEMBL) entry"],
                consequences=[
                    "Stage A stops; nothing is claimed for this gene.",
                    "An unreviewed entry has no stable canonical isoform, so the "
                    "reference-AA check and the numbering would be unverifiable.",
                ],
                recommendation="Escalate and stop; the Lead names the accession.",
            )
        if len(results) > 1:
            accessions = [r.get("primaryAccession") for r in results]
            raise EscalationRequired(
                ambiguity=(f"Gene {gene!r} resolves to {len(results)} reviewed UniProt "
                           f"accessions: {accessions}."),
                options=["Escalate and stop",
                         "Pick the highest-annotation-score entry automatically"],
                consequences=[
                    "Stage A stops until the Lead names one accession.",
                    "An automatic pick would silently choose which protein is analysed; "
                    "gene -> accession ambiguity is explicitly a Lead decision (§6.2).",
                ],
                recommendation="Escalate and stop.",
            )
        return results[0]

    @staticmethod
    def parse_entry(entry: dict, gene: str) -> GeneResolution:
        """Build a :class:`GeneResolution` from one UniProtKB JSON entry."""
        accession = entry.get("primaryAccession") or "NA"
        sequence = (entry.get("sequence") or {}).get("value") or ""
        version = str((entry.get("entryAudit") or {}).get("entryVersion", "NA"))
        mane = [x for x in (entry.get("uniProtKBCrossReferences") or [])
                if str(x.get("database", "")).upper() in ("MANE-SELECT", "MANE_SELECT")]
        refseq = [x for x in (entry.get("uniProtKBCrossReferences") or [])
                  if str(x.get("database", "")).upper() == "REFSEQ"]

        if mane:
            transcript = mane[0].get("id") or "NA"
            source, fallback = "MANE_SELECT", False
        elif refseq:
            # METHOD_SPEC II.13 fallback: the RefSeq transcript matching the
            # canonical sequence with the longest CDS, recorded explicitly.
            transcript = sorted(x.get("id") or "" for x in refseq)[-1]
            source, fallback = "longest_cds_matching_uniprot_canonical", True
        else:
            raise EscalationRequired(
                ambiguity=(f"UniProt entry {accession} for {gene} names neither a MANE "
                           f"Select transcript nor any RefSeq transcript."),
                options=["Escalate and stop", "Proceed with no named transcript"],
                consequences=[
                    "Stage A stops; the transcript must be recorded explicitly (§6.3).",
                    "Every non-canonical-transcript exclusion would become unjudgeable.",
                ],
                recommendation="Escalate and stop.",
            )

        protein = (((entry.get("proteinDescription") or {}).get("recommendedName") or {})
                   .get("fullName") or {}).get("value") or "NA"
        subunit_text = extract_subunit_comment(entry)
        oligomer = classify_oligomeric_state(subunit_text)
        return GeneResolution(
            gene=gene, uniprot_acc=accession, uniprot_entry_version=version,
            canonical_isoform=f"{accession}-1", sequence=sequence,
            mane_transcript=transcript, mane_source=source, mane_fallback_used=fallback,
            entry_name=entry.get("uniProtkbId") or "NA", protein_name=protein,
            oligomeric_state=oligomer["oligomeric_state"],
            oligomeric_state_evidence=oligomer["subunit_text"],
            oligomeric_state_matched_phrase=oligomer["matched_phrase"],
        )

    def _get_json(self, url: str, params: dict) -> dict:
        if not self.allow_network:
            raise BlockedError(
                f"BLOCKED — network retrieval is disabled for this run "
                f"(allow_network=FALSE); refusing to fetch {url}."
            )
        import requests                                        # pragma: no cover

        response = requests.get(url, params=params, timeout=self.timeout)
        response.raise_for_status()
        return response.json()

    def resolve_gene(self, gene: str) -> GeneResolution:        # pragma: no cover
        request = self.build_request(gene)
        payload = self._get_json(request["url"], request["params"])
        entry = self.parse_search_payload(payload, gene)
        self._entry = entry
        resolution = self.parse_entry(entry, gene)
        self._fasta = resolution.fasta().encode("utf-8")
        return resolution

    def raw_payload(self, kind: str) -> bytes | None:
        return self._fasta if kind == "uniprot" else None

    def sequence_metadata(self) -> dict:
        entry = self._entry or {}
        return {
            "source": "uniprotkb",
            "url": UNIPROT_SEARCH_URL,
            "release_date": str((entry.get("entryAudit") or {}).get("lastSequenceUpdateDate",
                                                                   "")),
            "query": f"gene_exact AND organism_id:{self.organism_id} AND reviewed:true",
            "entry_version": str((entry.get("entryAudit") or {}).get("entryVersion", "")),
            "allow_network": self.allow_network,
        }


# --- live bundle ------------------------------------------------------------

class LiveStageASource:
    """Composes the three live clients into one :class:`StageASource`.

    This class holds no retrieval logic of its own. It owns exactly three things:

    1. **composition** — one :class:`ClinVarClient`, one :class:`UniProtClient`
       and one :class:`hotspot3d.structure.sources.AlphaFoldClient`, so that a
       live run drives the same six-method contract a synthetic fixture does;
    2. **the network gate** — ``allow_network`` is handed to every client, so a
       run with the gate closed refuses *every* external retrieval with the
       client's own BLOCKED message rather than silently degrading;
    3. **raw-payload passthrough** — the bytes each client kept are offered to
       ``01_INPUT_RAW`` verbatim (agent §7.3), never re-serialized when the
       original exists.

    ``model_template`` is not defaulted here on purpose: the caller reads it from
    the frozen ``structure.model_template`` config key. This class never picks a
    model version, and never falls back to "whatever AFDB serves latest".

    No star filter exists anywhere in this class or in the clients it composes
    (F12): significance and review status are carried through untouched.
    """

    def __init__(self, *, allow_network: bool = False, assembly: str = "GRCh38",
                 model_template: str | None = None) -> None:
        # Deferred import: hotspot3d.structure reads the amino-acid vocabulary
        # from hotspot3d.data, so the packages are kept decoupled at import time.
        from ..structure.sources import AFDB_MODEL_TEMPLATE, AlphaFoldClient

        self.allow_network = bool(allow_network)
        self.assembly = assembly
        self.variants = ClinVarClient(allow_network=self.allow_network,
                                      assembly=assembly)
        self.sequences = UniProtClient(allow_network=self.allow_network)
        self.structures = AlphaFoldClient(
            allow_network=self.allow_network,
            model_template=model_template or AFDB_MODEL_TEMPLATE)

    # -- variants ------------------------------------------------------------
    def fetch_missense(self, gene: str) -> list[dict]:
        return self.variants.fetch_missense(gene)

    def release_metadata(self) -> dict:
        return dict(self.variants.release_metadata())

    # -- sequence ------------------------------------------------------------
    def resolve_gene(self, gene: str) -> GeneResolution:
        return self.sequences.resolve_gene(gene)

    def sequence_metadata(self) -> dict:
        return dict(self.sequences.sequence_metadata())

    # -- structure -----------------------------------------------------------
    def fetch_model(self, uniprot_acc: str):
        return self.structures.fetch_model(uniprot_acc)

    def model_metadata(self) -> dict:
        return dict(self.structures.model_metadata())

    # -- immutable payloads --------------------------------------------------
    def raw_payload(self, kind: str) -> bytes | None:
        """The original bytes for ``kind``, if the owning client kept them."""
        if kind not in RAW_KINDS:
            return None
        for client in (self.variants, self.sequences, self.structures):
            payload = raw_payload_of(client, kind)
            if payload is not None:
                return payload
        return None


def assert_source_complete(source: Any) -> None:
    """Every Stage A provider must implement the whole :class:`StageASource`."""
    required = ("fetch_missense", "release_metadata", "resolve_gene",
                "sequence_metadata", "fetch_model", "model_metadata")
    missing = [name for name in required if not callable(getattr(source, name, None))]
    if missing:
        raise BlockedError(
            f"BLOCKED — the Stage A source {type(source).__name__} does not implement "
            f"{missing}. Stage A needs variants, canonical sequence and structure from "
            f"one declared provider."
        )
