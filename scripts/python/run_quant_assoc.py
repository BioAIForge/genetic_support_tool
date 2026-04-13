from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from config_utils import get_required
from path_utils import resolve_output_root


def _ensure_exists(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _resolve_plink_input(config: dict[str, Any], project_root: Path) -> tuple[list[str], str]:
    input_cfg = config.get("input", {})
    if "pfile_prefix" in input_cfg:
        prefix = project_root / input_cfg["pfile_prefix"]
        for suffix in [".pgen", ".pvar", ".psam"]:
            _ensure_exists(prefix.with_suffix(suffix), f"PLINK2 pfile component {suffix}")
        return ["--pfile", str(prefix)], "pfile"
    if "bfile_prefix" in input_cfg:
        prefix = project_root / input_cfg["bfile_prefix"]
        for suffix in [".bed", ".bim", ".fam"]:
            _ensure_exists(prefix.with_suffix(suffix), f"PLINK binary component {suffix}")
        return ["--bfile", str(prefix)], "bfile"
    if "pedmap_prefix" in input_cfg:
        prefix = project_root / input_cfg["pedmap_prefix"]
        for suffix in [".ped", ".map"]:
            _ensure_exists(prefix.with_suffix(suffix), f"PLINK PED/MAP component {suffix}")
        return ["--pedmap", str(prefix)], "pedmap"
    raise KeyError(
        "Quantitative association requires input.pfile_prefix, "
        "input.bfile_prefix, or input.pedmap_prefix in config."
    )


def _parse_plink_glm(output_prefix: Path, result_tsv: Path) -> Path:
    glm_candidates = list(output_prefix.parent.glob(f"{output_prefix.name}*.glm.linear"))
    if not glm_candidates:
        raise FileNotFoundError(
            "PLINK2 did not produce a .glm.linear result file. "
            "Check whether plink2 ran successfully and whether the phenotype is continuous."
        )
    glm_file = glm_candidates[0]
    rows = _read_tsv(glm_file)
    parsed_rows: list[dict[str, Any]] = []
    for row in rows:
        if row.get("TEST", "") not in {"", "ADD"}:
            continue
        parsed_rows.append(
            {
                "variant_id": row.get("ID", ""),
                "chr": row.get("#CHROM", row.get("CHROM", "")),
                "pos": row.get("POS", ""),
                "ref": row.get("REF", ""),
                "alt": row.get("ALT1", row.get("ALT", "")),
                "n_samples": row.get("OBS_CT", ""),
                "beta": row.get("BETA", ""),
                "se": row.get("SE", ""),
                "p_value": row.get("P", ""),
            }
        )
    _write_tsv(
        result_tsv,
        parsed_rows,
        ["variant_id", "chr", "pos", "ref", "alt", "n_samples", "beta", "se", "p_value"],
    )
    return glm_file


def run_quant_assoc(config: dict[str, Any], project_root: Path) -> Path:
    phenotype_file = project_root / get_required(config, "input", "phenotype_file")
    covariate_file = project_root / get_required(config, "input", "covariate_file")
    phenotype_name = get_required(config, "phenotype", "name")
    covariates = get_required(config, "quant_assoc", "covariates")
    engine = config.get("quant_assoc", {}).get("engine", "plink2")
    plink2_bin = config.get("quant_assoc", {}).get("plink2_bin", "plink2")
    variance_standardize = bool(config.get("quant_assoc", {}).get("covar_variance_standardize", True))
    output_root = resolve_output_root(project_root, get_required(config, "output", "root_dir"))
    output_dir = output_root / "quant_assoc"

    if engine != "plink2":
        raise ValueError("quant_assoc currently supports only engine=plink2")

    if shutil.which(plink2_bin) is None:
        raise FileNotFoundError(
            f"PLINK2 executable not found: {plink2_bin}. "
            "Install plink2 and ensure it is available on PATH, or set quant_assoc.plink2_bin."
        )

    plink_input_args, plink_input_mode = _resolve_plink_input(config, project_root)
    _ensure_exists(phenotype_file, "Phenotype file")
    _ensure_exists(covariate_file, "Covariate file")

    output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_dir / "plink2_quant_assoc"

    covar_name_args: list[str] = []
    if covariates:
        covar_name_args = ["--covar-name", ",".join(covariates)]

    covar_scale_args: list[str] = []
    if covariates and variance_standardize:
        covar_scale_args = ["--covar-variance-standardize", *covariates]

    command = [
        plink2_bin,
        *plink_input_args,
        "--pheno",
        str(phenotype_file),
        "--pheno-name",
        phenotype_name,
        "--covar",
        str(covariate_file),
        *covar_name_args,
        *covar_scale_args,
        "--glm",
        "hide-covar",
        "--allow-extra-chr",
        "--out",
        str(output_prefix),
    ]

    completed = subprocess.run(command, cwd=project_root, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            "PLINK2 quantitative association run failed. "
            "Check whether plink2 is installed and the input PLINK files are valid."
        )

    result_tsv = output_dir / "quant_assoc_result.tsv"
    raw_glm_file = _parse_plink_glm(output_prefix=output_prefix, result_tsv=result_tsv)

    metadata = {
        "method": "quant_assoc",
        "engine": engine,
        "plink2_bin": plink2_bin,
        "plink_input_mode": plink_input_mode,
        "phenotype_file": str(phenotype_file),
        "covariate_file": str(covariate_file),
        "phenotype_name": phenotype_name,
        "covariates": covariates,
        "covar_variance_standardize": variance_standardize,
        "raw_glm_file": str(raw_glm_file),
        "result_file": str(result_tsv),
    }
    metadata_json = output_dir / "run_metadata.json"
    metadata_json.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return result_tsv


def main() -> int:
    raise SystemExit("This module is intended to be imported by main.py, not run directly.")


if __name__ == "__main__":
    sys.exit(main())
