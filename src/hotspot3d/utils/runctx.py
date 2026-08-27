"""RunContext — the object every stage receives from the Lead.

Owns the canonical directory layout (Output Contract IX.1), the frozen config,
the seed registry, and the provenance recorder. Stages never construct paths by
hand; they ask the context, which is what keeps one-writer-per-path true.
"""
from __future__ import annotations

import os
import platform
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from .config import FrozenConfig
from .errors import BlockedError
from .hashing import code_version, env_fingerprint, sha256_file
from .io import write_json
from .seeds import SeedRegistry

# Canonical FULL_RESULTS stage directories and their owning agent (IX.12).
STAGE_OWNERS: dict[str, str] = {
    "00_RUN_SUMMARY": "lead",
    "01_INPUT_RAW": "data-structure",
    "02_CLINVAR": "data-structure",
    "03_STRUCTURE_QC": "data-structure",
    "04_GLOBAL_CLUSTERING": "hotspot-statistics",
    "05_HOTSPOT_RADIUS": "hotspot-statistics",
    "06_FINAL_HOTSPOTS": "hotspot-statistics",
    "07_FOOTPRINT_RADIUS": "footprint-robustness",
    "08_FINAL_FOOTPRINT": "footprint-robustness",
    "09_ROBUSTNESS": "footprint-robustness",
    "10_ANNOTATION": "biological-annotation",
    "11_SENSITIVITY": "hotspot-statistics",
    "12_REPRODUCIBILITY": "lead",
}

# Which stage directories each agent may WRITE. Enforced by assert_may_write.
AGENT_WRITE_SCOPE: dict[str, set[str]] = {
    "data-structure": {"01_INPUT_RAW", "02_CLINVAR", "03_STRUCTURE_QC"},
    "hotspot-statistics": {"04_GLOBAL_CLUSTERING", "05_HOTSPOT_RADIUS",
                           "06_FINAL_HOTSPOTS", "11_SENSITIVITY"},
    "footprint-robustness": {"07_FOOTPRINT_RADIUS", "08_FINAL_FOOTPRINT",
                             "09_ROBUSTNESS"},
    "biological-annotation": {"10_ANNOTATION"},
    "lead": {"00_RUN_SUMMARY", "12_REPRODUCIBILITY"},
}

# Which stage directories each agent may READ. Reads flow strictly upstream.
AGENT_READ_SCOPE: dict[str, set[str]] = {
    "data-structure": set(),
    "hotspot-statistics": {"01_INPUT_RAW", "02_CLINVAR", "03_STRUCTURE_QC"},
    "footprint-robustness": {"01_INPUT_RAW", "02_CLINVAR", "03_STRUCTURE_QC",
                             "04_GLOBAL_CLUSTERING", "05_HOTSPOT_RADIUS",
                             "06_FINAL_HOTSPOTS", "11_SENSITIVITY"},
    "biological-annotation": {"01_INPUT_RAW", "02_CLINVAR", "03_STRUCTURE_QC",
                              "04_GLOBAL_CLUSTERING", "05_HOTSPOT_RADIUS",
                              "06_FINAL_HOTSPOTS", "07_FOOTPRINT_RADIUS",
                              "08_FINAL_FOOTPRINT", "09_ROBUSTNESS", "11_SENSITIVITY"},
    "lead": set(STAGE_OWNERS),
}


def make_run_id(config_sha256: str, code_ver: str, when: datetime | None = None) -> str:
    """``<YYYYMMDDTHHMMSSZ>_<config_sha[0:8]>_<code_version[0:8]>`` (IX.0).

    Deterministic and sortable: encodes *when*, *under which frozen configuration*
    and *with which pipeline version*, so two runs that should be identical are
    visibly identical.
    """
    when = when or datetime.now(timezone.utc)
    return f"{when.strftime('%Y%m%dT%H%M%SZ')}_{config_sha256[:8]}_{code_ver[:8]}"


@dataclass
class RunContext:
    gene: str
    run_id: str
    run_root: Path                    # the .partial directory during execution
    config: FrozenConfig
    code_version: str
    seeds: SeedRegistry
    started_utc: str
    synthetic: bool = False
    _provenance: dict = field(default_factory=dict)

    # -- construction --------------------------------------------------------
    @classmethod
    def create(cls, gene: str, config: FrozenConfig, results_root: str | Path = "results",
               synthetic: bool = False, when: datetime | None = None,
               run_id: str | None = None) -> "RunContext":
        code_ver = code_version()
        rid = run_id or make_run_id(config.sha256, code_ver, when)
        root = Path(results_root) / gene / f"{gene}_{rid}.partial"
        root.mkdir(parents=True, exist_ok=True)
        ctx = cls(
            gene=gene, run_id=rid, run_root=root, config=config,
            code_version=code_ver,
            seeds=SeedRegistry(master_seed=int(config.get("seeding.MASTER_SEED")),
                               run_id=rid),
            started_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            synthetic=synthetic,
        )
        ctx.full_results.mkdir(parents=True, exist_ok=True)
        return ctx

    # -- canonical paths -----------------------------------------------------
    @property
    def full_results(self) -> Path:
        return self.run_root / "FULL_RESULTS"

    @property
    def review_pack(self) -> Path:
        return self.run_root / "REVIEW_PACK"

    @property
    def archives(self) -> Path:
        return self.run_root / "ARCHIVES"

    def stage_dir(self, stage: str) -> Path:
        if stage not in STAGE_OWNERS:
            raise BlockedError(f"unknown canonical stage directory: {stage!r}")
        path = self.full_results / stage
        path.mkdir(parents=True, exist_ok=True)
        (path / "figures").mkdir(exist_ok=True)
        return path

    def provenance_path(self, stage: str) -> Path:
        d = self.full_results / "12_REPRODUCIBILITY" / "provenance"
        d.mkdir(parents=True, exist_ok=True)
        return d / f"provenance_{stage.lower()}.json"

    def handoff_path(self, name: str) -> Path:
        """Handoffs live in the producing agent's terminal stage directory."""
        location = {
            "handoff_01": "03_STRUCTURE_QC", "handoff_02": "06_FINAL_HOTSPOTS",
            "handoff_03": "08_FINAL_FOOTPRINT", "handoff_04": "09_ROBUSTNESS",
            "handoff_05": "10_ANNOTATION",
        }[name]
        return self.stage_dir(location) / f"{name}.json"

    # -- ownership enforcement ----------------------------------------------
    def assert_may_write(self, agent: str, stage: str) -> None:
        if stage not in AGENT_WRITE_SCOPE.get(agent, set()):
            raise BlockedError(
                f"OWNERSHIP violation: agent {agent!r} may not write {stage!r}. "
                f"Exactly one writer per path."
            )

    def assert_may_read(self, agent: str, stage: str) -> None:
        from .errors import LeakageError
        if stage not in AGENT_READ_SCOPE.get(agent, set()):
            raise LeakageError(
                f"INFORMATION BARRIER violation: agent {agent!r} may not read {stage!r}. "
                f"Reads flow strictly upstream."
            )

    # -- provenance ----------------------------------------------------------
    def record_provenance(self, stage: str, agent: str, *, inputs: dict | None = None,
                          outputs: dict | None = None, parameters: dict | None = None,
                          commands: list[str] | None = None,
                          warnings: list[dict] | None = None,
                          extra: dict | None = None) -> Path:
        payload = {
            "stage": stage, "agent_owner": agent, "run_id": self.run_id,
            "gene": self.gene, "config_sha256": self.config.sha256,
            "code_version": self.code_version,
            "started_utc": self.started_utc,
            "ended_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "environment": env_fingerprint(),
            "software_versions": software_versions(),
            "synthetic_input_mode": self.synthetic,
            "network_used": False if self.synthetic else None,
            "inputs": inputs or {},
            "outputs": outputs or {},
            "parameters": parameters or {},
            "derived_seeds": self.seeds.issued or {"note": "no stochastic context used"},
            "commands": commands or [],
            "warnings": warnings or [],
        }
        if extra:
            payload.update(extra)
        return write_json(self.provenance_path(stage), payload)

    def finalize(self, ok: bool = True) -> Path:
        """Atomic rename ``.partial`` -> final (or ``.failed``) — P8."""
        suffix = "" if ok else ".failed"
        final = self.run_root.parent / (self.run_root.name.replace(".partial", "") + suffix)
        if final.exists():
            import shutil
            shutil.rmtree(final)
        self.run_root.rename(final)
        object.__setattr__(self, "run_root", final)
        return final


def software_versions() -> dict[str, str]:
    """Package versions for 12_REPRODUCIBILITY/software_versions.tsv."""
    import importlib.metadata as md
    out = {"python": sys.version.split()[0], "platform": platform.platform()}
    for pkg in ("numpy", "scipy", "pandas", "pyyaml", "scikit-image",
                "scikit-learn", "matplotlib", "networkx", "biopython"):
        try:
            out[pkg] = md.version(pkg)
        except Exception:
            out[pkg] = "NOT_INSTALLED"
    return out


def uv_pip_freeze() -> str:
    try:
        res = subprocess.run(["uv", "pip", "freeze"], capture_output=True, text=True,
                             timeout=60, cwd=os.getcwd())
        return res.stdout.strip() if res.returncode == 0 else "uv pip freeze failed"
    except Exception as exc:                                  # pragma: no cover
        return f"unavailable: {exc}"
