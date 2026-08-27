"""Live ``AnnotationSource`` for Stage E — the five config-pinned methods.

Modelled on :class:`hotspot3d.data.sources.LiveStageASource`: this class holds no
retrieval logic of its own beyond composing four independent, single-purpose
clients and gating them behind ``allow_network`` exactly as
:class:`hotspot3d.annotation.sources.NetworkAnnotationSource` requires.

Sources, each already named/pinned in ``config/pipeline.yaml`` §II.13 or the
agent charter (§7.4):

* **secondary_structure** — DSSP, run locally (no network) on the already
  -downloaded AlphaFold model via a pinned Docker image (biocontainers has no
  DSSP recipe; conda-forge does, so the image is built from
  ``condaforge/miniforge3`` + ``conda install -c conda-forge dssp=4.5.3``).
* **domains** / **functional_sites** — UniProtKB REST, the same accession
  Stage A already resolved, requesting features (``ft_domain``, ``ft_region``,
  ``ft_motif``, ``ft_binding``, ``ft_site``) and InterPro/Pfam cross-references.
* **conservation** — dbNSFP (GERP++ RS, phyloP100way) via myvariant.info,
  queried by genomic position. The genomic position for a residue is recovered
  from the SAME raw ClinVar payload Stage A already downloaded (never a new
  ClinVar fetch): a residue only has one if a ClinVar record maps to it, so
  conservation is necessarily incomplete over ``U_struct`` — the residues it
  cannot resolve are reported through the EXISTING gap mechanism
  (``annotation/stage.py:_cons_rows``), not fabricated.
* **disease_associations** — the SAME already-downloaded raw ClinVar payload's
  ``PhenotypeList``, mapped to residues via ``variants_residue_level.tsv``'s
  ``clinvar_ids``. No new network call.

``literature_mechanisms`` is deliberately NOT implemented here — out of scope
per this task. It returns ``[]`` and sets ``literature_mechanisms_not_run =
True``, a marker ``annotation/stage.py`` checks to record an EXPLICIT
"not run" warning, distinct from "searched and found nothing" (agent
§9.7/§17).

Graceful degradation (the offline/gap-reporting policy already established by
``annotation_gaps.tsv``): every method independently catches its own retrieval
failures — a missing Docker/DSSP install, a network error, an unreachable
host, a malformed response — and returns an EMPTY result for that method
rather than raising, so one source being down never fails the whole stage.
Every such failure is recorded in ``self.source_gaps`` (the same shape
:func:`annotation.stage._gap` produces) with a severity reflecting how much
of the method failed, and merged into ``annotation_gaps.tsv`` by the stage.
This is a stronger claim than an empty list returned by a healthy source
(which means "queried successfully, nothing found" — agent §9.7's own
informative negative): a source-level failure is always visible as an
explicit gap, never silently indistinguishable from a real negative.
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import Any

from ..data.clinvar import parse_variant_summary
from ..utils.io import read_tsv
from .rubric import MechanismRecord
from .sources import SourceVersion

DSSP_DOCKER_IMAGE = "hotspot3d-dssp:4.5.3"
DSSP_VERSION = "4.5.3"

UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"
#: The feature/xref fields this client actually consumes; UniProtKB returns
#: only what is requested.
UNIPROT_FIELDS = ("ft_domain,ft_region,ft_motif,ft_binding,ft_site,"
                  "cc_disease,xref_interpro,xref_pfam,version")

MYVARIANT_BATCH_URL = "https://myvariant.info/v1/variant"
MYVARIANT_METADATA_URL = "https://myvariant.info/v1/metadata"
MYVARIANT_FIELDS = "dbnsfp.gerp++.rs,dbnsfp.phylop.100way_vertebrate.score"
#: One POST per this many variants (myvariant.info's own documented batch cap
#: is 1000; kept well under it for reliability over a single request).
MYVARIANT_BATCH_SIZE = 200

#: UniProt feature `type` -> Stage E bucket. Anything not listed here is
#: skipped for that bucket, never guessed into one.
_DOMAIN_FEATURE_TYPES = {"Domain"}
_SITE_FEATURE_TYPES = {
    "Region": "functional_region", "Motif": "functional_region",
    "Coiled coil": "functional_region", "Repeat": "functional_region",
    "Zinc finger": "functional_region", "DNA binding": "functional_region",
    "Site": "functional_region", "Binding site": "ligand_binding",
}


def _gap(scope: str, identifier: str, missing: str, reason: str, source: str,
        severity: str) -> dict:
    return {"scope": scope, "identifier": identifier, "missing_item": missing,
            "reason": reason, "source_attempted": source, "severity": severity}


class _DsspClient:
    """Secondary structure via a pinned local DSSP, never a network call."""

    def __init__(self, docker_image: str = DSSP_DOCKER_IMAGE):
        self.docker_image = docker_image

    def secondary_structure(self, structure_path: Path) -> tuple[dict[int, str], list[dict]]:
        structure_path = Path(structure_path).resolve()
        gaps: list[dict] = []
        try:
            with tempfile.TemporaryDirectory() as tmp:
                out_path = Path(tmp) / "out.dssp"
                proc = subprocess.run(
                    ["docker", "run", "--rm",
                     "-v", f"{structure_path.parent}:/data:ro",
                     "-v", f"{tmp}:/out",
                     self.docker_image,
                     f"/data/{structure_path.name}", "/out/out.dssp"],
                    capture_output=True, text=True, timeout=600,
                )
                if proc.returncode != 0 or not out_path.is_file():
                    gaps.append(_gap(
                        "source", "secondary_structure", "secondary_structure",
                        f"DSSP failed (exit {proc.returncode}): "
                        f"{proc.stderr.strip()[:300] or 'no output produced'}",
                        "DSSP", "MAJOR"))
                    return {}, gaps
                return _parse_dssp(out_path), gaps
        except (OSError, subprocess.SubprocessError) as exc:
            gaps.append(_gap(
                "source", "secondary_structure", "secondary_structure",
                f"DSSP could not be run: {type(exc).__name__}: {exc}",
                "DSSP", "MAJOR"))
            return {}, gaps


def _parse_dssp(path: Path) -> dict[int, str]:
    """Classic DSSP fixed-width text format: column 17 is the SS code, the PDB
    residue number is columns 6-10. Header ends at the ``#  RESIDUE`` line."""
    ss_map: dict[int, str] = {}
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    start = next((i for i, ln in enumerate(lines) if ln.startswith("  #  RESIDUE")), None)
    if start is None:
        return ss_map
    for line in lines[start + 1:]:
        if len(line) < 17:
            continue
        resnum = line[5:10].strip()
        code = line[16]
        if not resnum or resnum == "!":       # chain break marker
            continue
        try:
            idx = int(resnum)
        except ValueError:
            continue
        ss_map[idx] = code if code != " " else "-"
    return ss_map


class _UniProtFeatureClient:
    """Domains and functional sites from the same UniProtKB entry Stage A
    already resolved the accession for. One request, cached."""

    def __init__(self, *, timeout: int = 60):
        self.timeout = timeout
        self._entry: dict | None = None
        self._acc: str | None = None
        self.version = "NA"

    def _fetch(self, acc: str) -> dict:
        if self._entry is not None and self._acc == acc:
            return self._entry
        import requests

        url = UNIPROT_ENTRY_URL.format(acc=acc)
        resp = requests.get(url, params={"fields": UNIPROT_FIELDS}, timeout=self.timeout)
        resp.raise_for_status()
        entry = resp.json()
        self._entry, self._acc = entry, acc
        audit = entry.get("entryAudit", {}) or {}
        self.version = str(audit.get("entryVersion") or entry.get("annotationScore") or "NA")
        return entry

    def domains(self, acc: str) -> tuple[list[dict], list[dict]]:
        gaps: list[dict] = []
        try:
            entry = self._fetch(acc)
        except Exception as exc:                                    # noqa: BLE001
            gaps.append(_gap("source", "domains", "domains",
                             f"UniProtKB feature retrieval failed: "
                             f"{type(exc).__name__}: {exc}", "UniProt", "MAJOR"))
            return [], gaps
        try:
            rows = []
            interpro = {x.get("id"): _prop(x, "EntryName")
                       for x in entry.get("uniProtKBCrossReferences", [])
                       if x.get("database") == "InterPro"}
            for feat in entry.get("features", []):
                if feat.get("type") not in _DOMAIN_FEATURE_TYPES:
                    continue
                start, end = _feature_span(feat)
                if start is None:
                    continue
                rows.append({
                    "domain_id": _first(interpro) or "NA",
                    "name": feat.get("description") or feat.get("type"),
                    "start": start, "end": end, "source": "UniProt",
                })
            return rows, gaps
        except Exception as exc:                                     # noqa: BLE001
            gaps.append(_gap("source", "domains", "domains",
                             f"UniProtKB domain features could not be parsed: "
                             f"{type(exc).__name__}: {exc}", "UniProt", "MAJOR"))
            return [], gaps

    def functional_sites(self, acc: str) -> tuple[list[dict], list[dict]]:
        gaps: list[dict] = []
        try:
            entry = self._fetch(acc)
        except Exception as exc:                                    # noqa: BLE001
            gaps.append(_gap("source", "functional_sites", "functional_sites",
                             f"UniProtKB feature retrieval failed: "
                             f"{type(exc).__name__}: {exc}", "UniProt", "MAJOR"))
            return [], gaps
        try:
            rows = []
            for i, feat in enumerate(entry.get("features", [])):
                site_type = _SITE_FEATURE_TYPES.get(feat.get("type"))
                if site_type is None:
                    continue
                start, end = _feature_span(feat)
                if start is None:
                    continue
                rows.append({
                    "site_id": f"UniProt:{acc}:{i}", "site_type": site_type,
                    "residue_start": start, "residue_end": end,
                    "description": feat.get("description") or feat.get("type"),
                    "source": "UniProt",
                })
            return rows, gaps
        except Exception as exc:                                     # noqa: BLE001
            gaps.append(_gap("source", "functional_sites", "functional_sites",
                             f"UniProtKB site features could not be parsed: "
                             f"{type(exc).__name__}: {exc}", "UniProt", "MAJOR"))
            return [], gaps


def _feature_span(feat: dict) -> tuple[int | None, int | None]:
    loc = feat.get("location") or {}
    start = (loc.get("start") or {}).get("value")
    end = (loc.get("end") or {}).get("value")
    return (int(start), int(end)) if start is not None and end is not None else (None, None)


def _prop(xref: dict, name: str) -> str | None:
    for p in xref.get("properties", []):
        if p.get("key") == name:
            return p.get("value")
    return None


def _first(mapping: dict) -> str | None:
    return next(iter(mapping), None)


class _ClinVarGenomicIndex:
    """Genomic position and phenotype per residue, from the raw ClinVar bytes
    Stage A already downloaded — never a new ClinVar fetch."""

    def __init__(self, *, clinvar_raw_path: Path, variants_residue_level_path: Path,
                gene: str, assembly: str = "GRCh38"):
        self.clinvar_raw_path = Path(clinvar_raw_path)
        self.variants_residue_level_path = Path(variants_residue_level_path)
        self.gene = gene
        self.assembly = assembly
        self._by_variation_id: dict[str, dict] | None = None
        self._residue_variation_ids: dict[int, list[str]] | None = None

    def release_version(self) -> str:
        """The SAME release Stage A already recorded for this raw payload —
        never a new claim about ClinVar's current state."""
        import json

        log_path = self.clinvar_raw_path.parent / "retrieval_log.json"
        try:
            log = json.loads(log_path.read_text(encoding="utf-8"))
            for r in log.get("retrievals", []):
                if r.get("artifact") == self.clinvar_raw_path.name:
                    release = r.get("release_date")
                    if release:
                        return str(release)
        except (OSError, ValueError, KeyError):
            pass
        return "NA"

    def _load(self) -> None:
        if self._by_variation_id is not None:
            return
        payload = self.clinvar_raw_path.read_bytes()
        rows = parse_variant_summary(payload, gene=self.gene, assembly=self.assembly)
        self._by_variation_id = {str(r.get("VariationID")): r for r in rows if r.get("VariationID")}

        residue_ids: dict[int, list[str]] = {}
        for row in read_tsv(self.variants_residue_level_path):
            raw = str(row.get("clinvar_ids") or "").strip()
            if not raw:
                continue
            idx = int(row["residue_index"])
            residue_ids[idx] = [v for v in raw.split(";") if v]
        self._residue_variation_ids = residue_ids

    def genomic_position(self, residue_index: int) -> dict | None:
        """Chromosome/PositionVCF/ref/alt for ONE representative ClinVar record
        at this residue, or None if the residue has no ClinVar record at all
        (most of U_struct — a structural fact, not a retrieval failure)."""
        self._load()
        for vid in self._residue_variation_ids.get(residue_index, []):
            row = self._by_variation_id.get(vid)
            if row and row.get("Chromosome") and row.get("PositionVCF"):
                return row
        return None

    def diseases_for_residue(self, residue_index: int) -> list[dict]:
        self._load()
        out = []
        for vid in self._residue_variation_ids.get(residue_index, []):
            row = self._by_variation_id.get(vid)
            if not row:
                continue
            phenotypes = [p for p in str(row.get("PhenotypeList", "")).split("|")
                         if p and p.strip().lower() != "not provided"]
            for disease in phenotypes:
                out.append({
                    "disease": disease, "residue_start": residue_index,
                    "residue_end": residue_index,
                    "reference": ";".join(
                        r for r in str(row.get("RCVaccession", "")).split("|") if r) or "NA",
                    "source": "ClinVar",
                })
        return out


class _ConservationClient:
    """dbNSFP (GERP++ RS, phyloP100way) via myvariant.info, by genomic position.

    Only residues with a resolvable ClinVar-derived genomic position are ever
    queried; every other residue is simply absent from the returned dict, which
    the existing ``_cons_rows`` gap path already reports honestly.
    """

    def __init__(self, genomic_index: _ClinVarGenomicIndex, *, timeout: int = 30):
        self.genomic_index = genomic_index
        self.timeout = timeout
        self.dbnsfp_version = "NA"
        self._version_fetched = False

    def _fetch_version(self) -> None:
        if self._version_fetched:
            return
        self._version_fetched = True
        try:
            import requests

            resp = requests.get(MYVARIANT_METADATA_URL, timeout=self.timeout)
            resp.raise_for_status()
            self.dbnsfp_version = str(
                resp.json().get("src", {}).get("dbnsfp", {}).get("version", "NA"))
        except Exception:                                            # noqa: BLE001
            pass          # non-fatal: individual queries still carry their own gap on failure

    def conservation(self, residues: list[int]) -> tuple[dict[int, dict], list[dict]]:
        import requests

        self._fetch_version()
        out: dict[int, dict] = {}
        gaps: list[dict] = []

        hgvs_by_residue: dict[int, str] = {}
        try:
            for idx in residues:
                pos = self.genomic_index.genomic_position(idx)
                if pos is None:
                    continue                   # no ClinVar record here; not a failure
                hgvs_by_residue[idx] = (
                    f"chr{pos['Chromosome']}:g.{pos['PositionVCF']}"
                    f"{pos['ReferenceAlleleVCF']}>{pos['AlternateAlleleVCF']}")
        except Exception as exc:                                      # noqa: BLE001
            gaps.append(_gap(
                "source", "conservation", "conservation",
                f"could not resolve genomic positions from the ClinVar-derived "
                f"index: {type(exc).__name__}: {exc}", "dbNSFP", "MAJOR"))
            return out, gaps

        attempted = len(hgvs_by_residue)
        if attempted == 0:
            return out, gaps

        by_hgvs = {v: k for k, v in hgvs_by_residue.items()}     # 1:1 within one gene
        all_hgvs = list(hgvs_by_residue.values())
        failed_batches = 0
        for chunk_start in range(0, len(all_hgvs), MYVARIANT_BATCH_SIZE):
            chunk = all_hgvs[chunk_start:chunk_start + MYVARIANT_BATCH_SIZE]
            try:
                resp = requests.post(
                    MYVARIANT_BATCH_URL,
                    data={"ids": ",".join(chunk), "assembly": "hg38",
                          "fields": MYVARIANT_FIELDS},
                    timeout=self.timeout)
                resp.raise_for_status()
                for hit in resp.json():
                    idx = by_hgvs.get(hit.get("query"))
                    if idx is None or hit.get("notfound"):
                        continue
                    dbnsfp = hit.get("dbnsfp", {}) or {}
                    gerp = dbnsfp.get("gerp++", {}).get("rs")
                    phylop = (dbnsfp.get("phylop", {}).get("100way_vertebrate", {})
                             or {}).get("score")
                    if gerp is None and phylop is None:
                        continue    # queried successfully, no dbNSFP entry at this position
                    out[idx] = {"GERP++_RS": gerp, "phyloP100way": phylop}
            except Exception:                                        # noqa: BLE001
                failed_batches += 1
        if failed_batches:
            resolved = len(out)
            gaps.append(_gap(
                "source", "conservation", "conservation",
                f"dbNSFP/myvariant.info batch retrieval failed for {failed_batches} "
                f"chunk(s); {resolved} of {attempted} residues with a resolvable "
                f"genomic position were annotated",
                "dbNSFP", "ADVISORY" if resolved else "MAJOR"))
        return out, gaps


class LiveAnnotationSource:
    """Composes the four live clients into the five config-pinned methods.

    ``literature_mechanisms`` is out of scope for this build (see module
    docstring); it always returns ``[]`` and sets
    ``literature_mechanisms_not_run = True``.
    """

    def __init__(self, *, allow_network: bool, gene: str, uniprot_acc: str,
                clinvar_raw_path: Path, variants_residue_level_path: Path,
                assembly: str = "GRCh38", dssp_docker_image: str = DSSP_DOCKER_IMAGE):
        self.allow_network = bool(allow_network)
        self.gene = gene
        self.uniprot_acc = uniprot_acc
        self.literature_mechanisms_not_run = True
        self.source_gaps: list[dict] = []

        self._dssp = _DsspClient(dssp_docker_image)
        self._uniprot = _UniProtFeatureClient()
        self._genomic_index = _ClinVarGenomicIndex(
            clinvar_raw_path=clinvar_raw_path,
            variants_residue_level_path=variants_residue_level_path,
            gene=gene, assembly=assembly)
        self._conservation = _ConservationClient(self._genomic_index)

    def _require_network(self, what: str) -> bool:
        if self.allow_network:
            return True
        self.source_gaps.append(_gap(
            "source", what, what,
            "execution.allow_network is FALSE; Stage E holds network privileges "
            "only when the run explicitly enables them",
            what, "MAJOR"))
        return False

    def versions(self) -> dict[str, SourceVersion]:
        """Eagerly resolves UniProt's/dbNSFP's version if not already fetched
        (stage.py validates every source's versions BEFORE calling any of the
        five retrieval methods, so this cannot rely on their normal lazy
        population as a side effect of domains()/conservation() running
        first). Network-gated and individually best-effort: a version that
        cannot be resolved stays "NA" and is caught by validate_source's own
        ANNOTATION_SOURCE_UNVERSIONED check, not silently accepted here.
        """
        if self.allow_network and self._uniprot.version == "NA":
            try:
                self._uniprot._fetch(self.uniprot_acc)
            except Exception:                                        # noqa: BLE001
                pass
        if self.allow_network and self._conservation.dbnsfp_version == "NA":
            self._conservation._fetch_version()
        return {
            "DSSP": SourceVersion(name="DSSP", version=DSSP_VERSION,
                                  url="https://github.com/PDB-REDO/dssp"),
            "UniProt": SourceVersion(name="UniProt", version=self._uniprot.version,
                                     accession=self.uniprot_acc),
            "dbNSFP": SourceVersion(name="dbNSFP", version=self._conservation.dbnsfp_version,
                                    url="https://myvariant.info"),
            "ClinVar": SourceVersion(name="ClinVar",
                                     version=self._genomic_index.release_version(),
                                     url="https://www.ncbi.nlm.nih.gov/clinvar/"),
        }

    # -- the five methods -----------------------------------------------
    def secondary_structure(self, structure: Any) -> dict[int, str]:
        ss_map, gaps = self._dssp.secondary_structure(structure)
        self.source_gaps.extend(gaps)
        return ss_map

    def domains(self, acc: str) -> list[dict]:
        if not self._require_network("domains"):
            return []
        rows, gaps = self._uniprot.domains(acc)
        self.source_gaps.extend(gaps)
        return rows

    def functional_sites(self, acc: str) -> list[dict]:
        if not self._require_network("functional_sites"):
            return []
        rows, gaps = self._uniprot.functional_sites(acc)
        self.source_gaps.extend(gaps)
        return rows

    def conservation(self, acc: str, residues: Any) -> dict[int, dict]:
        if not self._require_network("conservation"):
            return {}
        cons, gaps = self._conservation.conservation(list(residues))
        self.source_gaps.extend(gaps)
        return cons

    def disease_associations(self, gene: str) -> list[dict]:
        try:
            out = []
            self._genomic_index._load()
            for idx in self._genomic_index._residue_variation_ids:
                out.extend(self._genomic_index.diseases_for_residue(idx))
            return out
        except Exception as exc:                                     # noqa: BLE001
            self.source_gaps.append(_gap(
                "source", "disease_associations", "disease_associations",
                f"ClinVar-derived disease association retrieval failed: "
                f"{type(exc).__name__}: {exc}", "ClinVar", "MAJOR"))
            return []

    def literature_mechanisms(self, gene: str) -> list[MechanismRecord]:
        return []
