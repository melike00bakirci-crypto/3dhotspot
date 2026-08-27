"""Archives, checksums and transfer info — Output Contract IX.9 (Lead-owned).

P10 No silent transfer: the pipeline NEVER copies anything off the machine. It
emits transfer_info.txt with placeholder commands for the user to run.
ARCHIVES/ is excluded from both archives (no recursion).
"""
from __future__ import annotations

import tarfile
import zipfile
from pathlib import Path

from ..utils.hashing import sha256_file
from ..utils.io import write_text
from ..utils.runctx import RunContext


def build_archives(ctx: RunContext) -> dict:
    """Create REVIEW_PACK.zip, FULL_RESULTS.tar.gz, checksums and transfer_info."""
    arch = ctx.archives
    arch.mkdir(parents=True, exist_ok=True)
    produced: dict[str, Path] = {}

    # -- REVIEW_PACK.zip — Windows opens zip natively by double-click ---------
    if ctx.review_pack.is_dir():
        zip_path = arch / f"{ctx.gene}_REVIEW_PACK.zip"
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for p in sorted(ctx.review_pack.rglob("*")):
                if p.is_file():
                    # single top-level entry REVIEW_PACK/ so unzipping never scatters
                    zf.write(p, Path("REVIEW_PACK") / p.relative_to(ctx.review_pack))
        produced["review_pack_zip"] = zip_path

    # -- FULL_RESULTS.tar.gz — tar.gz for universality ----------------------
    if ctx.full_results.is_dir():
        tar_path = arch / f"{ctx.gene}_FULL_RESULTS.tar.gz"
        with tarfile.open(tar_path, "w:gz") as tf:
            for p in sorted(ctx.full_results.rglob("*")):
                if p.is_file():
                    tf.add(p, arcname=str(Path("FULL_RESULTS") /
                                          p.relative_to(ctx.full_results)))
        produced["full_results_tar"] = tar_path

    # -- checksums.sha256 in canonical `sha256sum -c` format ----------------
    lines = [f"{sha256_file(p)}  {p.name}"
             for p in sorted(produced.values()) if p.is_file()]
    checksums = arch / "checksums.sha256"
    checksums.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    produced["checksums"] = checksums

    # -- transfer_info.txt (CRLF; placeholders only, never real credentials) -
    write_text(arch / "transfer_info.txt", _transfer_text(ctx, produced), crlf=True)
    return {k: str(v) for k, v in produced.items()}


def _transfer_text(ctx: RunContext, produced: dict[str, Path]) -> str:
    rows = []
    for key in ("review_pack_zip", "full_results_tar"):
        p = produced.get(key)
        if p and p.is_file():
            rows.append(f"  {p.name}\n"
                        f"    size   : {p.stat().st_size:,} bytes\n"
                        f"    sha256 : {sha256_file(p)}\n"
                        f"    path   : {p.resolve()}")
    listing = "\n".join(rows) or "  (no archives were produced for this run)"
    return f"""TRANSFER INFORMATION — {ctx.gene} / {ctx.run_id}
================================================================

This pipeline NEVER transfers anything off this machine. The commands below are
placeholders for YOU to run; substitute your own username, host and destination.

ARCHIVES
{listing}

COPY TO A LOCAL MACHINE (choose one)

  scp <user>@<host>:{ctx.archives.resolve()}/{ctx.gene}_REVIEW_PACK.zip <local_dest>

  rsync -avP --partial <user>@<host>:{ctx.archives.resolve()}/ <local_dest>

VERIFY AFTER TRANSFER

  Linux / macOS :  sha256sum -c checksums.sha256
  Windows PS    :  Get-FileHash .\\{ctx.gene}_REVIEW_PACK.zip -Algorithm SHA256
  Windows cmd   :  certutil -hashfile {ctx.gene}_REVIEW_PACK.zip SHA256

THEN

  1. Unzip {ctx.gene}_REVIEW_PACK.zip
  2. Open REVIEW_PACK/index.html by double-clicking it.
     No server and no internet connection are required.
"""


def write_readme(ctx: RunContext, run_status: str, banner: str) -> Path:
    """README.txt at the run root — 20 lines, CRLF for naive Windows editors."""
    return write_text(ctx.run_root / "README.txt", f"""{banner}
================================================================

WHAT THIS IS
  Output of the 3D hotspot / footprint pipeline for {ctx.gene}.
  RUN_ID: {ctx.run_id}

START HERE
  REVIEW_PACK/index.html    open by double-click; works offline
  REVIEW_PACK/final_summary.md

IF YOU NEED THE DETAIL
  FULL_RESULTS/             complete audit trail, one directory per stage
  FULL_RESULTS/12_REPRODUCIBILITY/   config, seeds, versions, checksums
  manifest.tsv              every file, including expected-but-absent ones
  warnings.tsv              severity-sorted; read MAJOR and BLOCKING first

TO SHARE
  ARCHIVES/{ctx.gene}_REVIEW_PACK.zip        small, for reviewers
  ARCHIVES/{ctx.gene}_FULL_RESULTS.tar.gz    complete
  ARCHIVES/transfer_info.txt                 copy + verify commands

run_status = {run_status}
""", crlf=True)
