"""Canonical command-line entry point.

    hotspot3d run --gene GENE
    hotspot3d run --gene GENE --synthetic-case clustered
    hotspot3d validate-config
    hotspot3d show-config
    hotspot3d version

Gene-independent by construction: nothing in this module (or anywhere in the
pipeline) branches on a specific gene or protein.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from ..utils.config import assert_frozen_methodology, load_config
from ..utils.errors import BlockedError
from ..utils.hashing import code_version
from .pipeline import SourceBundle, run_pipeline

DESCRIPTION = (
    "Gene-independent 3D spatial hotspot and footprint pipeline for pathogenic "
    "missense variants. Methodology is FROZEN and pre-registered in config/pipeline.yaml."
)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="hotspot3d", description=DESCRIPTION)
    p.add_argument("--version", action="store_true", help="print version and code_version")
    sub = p.add_subparsers(dest="command")

    run = sub.add_parser("run", help="run the full pipeline for one gene")
    run.add_argument("--gene", required=True, help="HGNC gene symbol (e.g. TP53)")
    run.add_argument("--config", default="config/pipeline.yaml",
                     help="frozen configuration (default: config/pipeline.yaml)")
    run.add_argument("--gene-config", default=None,
                     help="optional per-gene overlay, config/genes/<GENE>.yaml")
    run.add_argument("--results-root", default="results", help="output root")
    run.add_argument("--run-id", default=None,
                     help="override RUN_ID (used for deterministic reproduction checks)")
    run.add_argument("--synthetic-case", default=None,
                     help="run against a named synthetic fixture instead of live data "
                          "(clustered | no_cluster | sparse | conflict | low_plddt | ...)")
    run.add_argument("--quiet", action="store_true", help="suppress the terminal summary")

    sub.add_parser("validate-config", help="assert the FROZEN methodology constants")
    show = sub.add_parser("show-config", help="print a resolved config value")
    show.add_argument("key", nargs="?", help="dotted key, e.g. fdr.q")
    sub.add_parser("version", help="print version and code_version")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.version or args.command == "version":
        print(f"hotspot3d 0.1.0\ncode_version {code_version()}")
        return 0

    if args.command == "validate-config":
        return _validate_config()

    if args.command == "show-config":
        cfg = load_config()
        if args.key:
            print(cfg.get(args.key))
        else:
            import yaml
            print(yaml.safe_dump(cfg.data, sort_keys=False))
        return 0

    if args.command != "run":
        build_parser().print_help()
        return 2

    return _run(args)


def _validate_config() -> int:
    try:
        cfg = load_config()
        assert_frozen_methodology(cfg)
    except BlockedError as exc:
        print(f"CONFIG INVALID\n{exc}", file=sys.stderr)
        return 2
    print(f"CONFIG VALID — frozen methodology constants intact\n"
          f"config_sha256 {cfg.sha256}")
    return 0


def _run(args) -> int:
    gene_overlay = args.gene_config
    if gene_overlay is None:
        candidate = Path("config/genes") / f"{args.gene}.yaml"
        gene_overlay = str(candidate) if candidate.is_file() else None

    try:
        cfg = load_config(args.config, gene_overlay)
        assert_frozen_methodology(cfg)
    except BlockedError as exc:
        print(f"RUN BLOCKED: {args.gene}\n{exc}", file=sys.stderr)
        return 2

    sources, synthetic = SourceBundle(), False
    if args.synthetic_case:
        sources, synthetic = _synthetic_sources(args.synthetic_case)

    result = run_pipeline(
        gene=args.gene, config=cfg, results_root=args.results_root,
        sources=sources, synthetic=synthetic, run_id=args.run_id,
    )

    if not args.quiet:
        print(result.terminal_summary)
    if result.escalations:
        print("\nESCALATIONS REQUIRING A LEAD DECISION", file=sys.stderr)
        for esc in result.escalations:
            print(f"  [{esc['stage']} / {esc['agent']}] {esc['ambiguity']}\n"
                  f"    options       : {esc['options']}\n"
                  f"    consequences  : {esc['consequences']}\n"
                  f"    recommendation: {esc['recommendation']}", file=sys.stderr)
    return result.exit_code


def _synthetic_sources(case_name: str) -> tuple[SourceBundle, bool]:
    """Load a deterministic synthetic fixture. No network, no real biological data."""
    try:
        from tests.fixtures.synthetic_clinvar import make_synthetic_case
        from tests.fixtures.synthetic_annotation import MockAnnotationSource
    except ImportError as exc:                                # pragma: no cover
        raise BlockedError(
            f"synthetic fixtures are not available ({exc}); "
            f"run from the repository root with the dev environment active"
        ) from exc
    case = make_synthetic_case(case_name)
    # The Stage A fixture exposes ONE provider satisfying both the VariantSource and
    # the StructureSource protocol (see tests/fixtures/synthetic_clinvar.py), so the
    # same object fills both slots. There is no `case.variant_source`.
    return SourceBundle(variant_source=case.source,
                        structure_source=case.source,
                        annotation_source=MockAnnotationSource()), True


if __name__ == "__main__":                                    # pragma: no cover
    raise SystemExit(main())
