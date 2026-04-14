from __future__ import annotations

import argparse
import sys
from pathlib import Path

from config_utils import load_config
from run_burden import run_burden
from run_gwas_gene_catalog import run_gwas_gene_catalog
from run_quant_assoc import run_quant_assoc
from run_skato import run_skato


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Genetic Support Tool: Python orchestrator + R statistics backend."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    burden_parser = subparsers.add_parser(
        "burden",
        help="Run the burden test module.",
    )
    burden_parser.add_argument(
        "--config",
        required=True,
        help="Path to YAML config file.",
    )

    skato_parser = subparsers.add_parser(
        "skato",
        help="Run the SKAT-O module.",
    )
    skato_parser.add_argument(
        "--config",
        required=True,
        help="Path to YAML config file.",
    )

    quant_parser = subparsers.add_parser(
        "quant-assoc",
        help="Run the quantitative phenotype association module.",
    )
    quant_parser.add_argument(
        "--config",
        required=True,
        help="Path to YAML config file.",
    )

    gene_catalog_parser = subparsers.add_parser(
        "gwas-gene-catalog",
        help="Annotate gene-level results with GWAS Catalog evidence.",
    )
    gene_catalog_parser.add_argument(
        "--config",
        required=True,
        help="Path to YAML config file.",
    )

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parents[2]

    if args.command == "burden":
        config = load_config(project_root / args.config)
        result_path = run_burden(config, project_root)
        print(f"Burden analysis completed: {result_path}")
        return 0

    if args.command == "skato":
        config = load_config(project_root / args.config)
        result_path = run_skato(config, project_root)
        print(f"SKAT-O analysis completed: {result_path}")
        return 0

    if args.command == "quant-assoc":
        config = load_config(project_root / args.config)
        result_path = run_quant_assoc(config, project_root)
        print(f"Quantitative association analysis completed: {result_path}")
        return 0

    if args.command == "gwas-gene-catalog":
        config = load_config(project_root / args.config)
        result_path = run_gwas_gene_catalog(config, project_root)
        print(f"GWAS gene catalog annotation completed: {result_path}")
        return 0

    parser.error(f"Unsupported command: {args.command}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
