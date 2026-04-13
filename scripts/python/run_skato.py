from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from config_utils import get_required
from path_utils import resolve_output_root
from run_regenie import (
    build_regenie_metadata,
    ensure_exists as ensure_regenie_exists,
    run_regenie_set_tests,
    write_metadata,
)


def _ensure_exists(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def _run_command(command: list[str], *, cwd: Path, tool_name: str, install_hint: str) -> None:
    executable = command[0]
    if shutil.which(executable) is None:
        raise FileNotFoundError(
            f"{tool_name} executable not found: {executable}. {install_hint}"
        )

    completed = subprocess.run(
        command,
        cwd=cwd,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"{tool_name} backend failed.\n"
            f"Command: {' '.join(command)}\n"
            f"Exit code: {completed.returncode}\n"
            f"STDOUT:\n{completed.stdout}\n"
            f"STDERR:\n{completed.stderr}"
        )


def run_skato(config: dict[str, Any], project_root: Path) -> Path:
    phenotype_file = project_root / get_required(config, "input", "phenotype_file")
    phenotype_name = get_required(config, "phenotype", "name")
    phenotype_type = get_required(config, "phenotype", "type")
    set_id = get_required(config, "skato", "set_id")
    engine = config.get("skato", {}).get("engine", "skat")
    output_root = resolve_output_root(project_root, get_required(config, "output", "root_dir"))
    output_dir = output_root / "skato"

    _ensure_exists(phenotype_file, "Phenotype file")

    output_dir.mkdir(parents=True, exist_ok=True)

    result_tsv = output_dir / "skato_result.tsv"
    metadata_json = output_dir / "run_metadata.json"

    if engine == "regenie":
        genotype_format = config.get("regenie", {}).get("genotype_format", "pfile")
        if genotype_format not in {"pfile", "bfile"}:
            raise ValueError("regenie.genotype_format must be either 'pfile' or 'bfile'.")
        genotype_key = "pfile_prefix" if genotype_format == "pfile" else "bfile_prefix"
        genotype_prefix = project_root / get_required(config, "input", genotype_key)
        covariate_value = config.get("input", {}).get("covariate_file")
        covariate_file = project_root / covariate_value if covariate_value else None
        anno_file = project_root / get_required(config, "input", "anno_file")
        set_list = project_root / get_required(config, "input", "set_list")
        mask_def = project_root / get_required(config, "input", "mask_def")
        covariates = config.get("skato", {}).get("covariates", [])
        regenie_bin = config.get("skato", {}).get("regenie_bin", "regenie")
        aaf_bins = str(config.get("skato", {}).get("regenie_aaf_bins", "0.01"))
        build_mask = str(config.get("skato", {}).get("regenie_build_mask", "max"))
        ignore_pred = bool(config.get("skato", {}).get("regenie_ignore_pred", False))
        pred_value = config.get("input", {}).get("pred_file")
        pred_file = project_root / pred_value if pred_value else None
        bsize = int(config.get("skato", {}).get("regenie_bsize", 200))
        extra_args = config.get("skato", {}).get("regenie_extra_args", [])
        regenie_prefix = output_dir / "regenie_skato"

        ensure_regenie_exists(anno_file, "regenie annotation file")
        ensure_regenie_exists(set_list, "regenie set-list file")
        ensure_regenie_exists(mask_def, "regenie mask-def file")
        if covariate_file is not None:
            ensure_regenie_exists(covariate_file, "Covariate file")
        if pred_file is not None:
            ensure_regenie_exists(pred_file, "regenie prediction file")
        elif not ignore_pred:
            raise FileNotFoundError(
                "regenie prediction file not configured. Set input.pred_file or "
                "skato.regenie_ignore_pred=true for a small demo run."
            )
        if genotype_format == "pfile":
            ensure_regenie_exists(genotype_prefix.with_suffix(".pgen"), "regenie PGEN file")
            ensure_regenie_exists(genotype_prefix.with_suffix(".pvar"), "regenie PVAR file")
            ensure_regenie_exists(genotype_prefix.with_suffix(".psam"), "regenie PSAM file")
        else:
            ensure_regenie_exists(genotype_prefix.with_suffix(".bed"), "regenie BED file")
            ensure_regenie_exists(genotype_prefix.with_suffix(".bim"), "regenie BIM file")
            ensure_regenie_exists(genotype_prefix.with_suffix(".fam"), "regenie FAM file")

        regenie_result = run_regenie_set_tests(
            regenie_bin=regenie_bin,
            genotype_prefix=genotype_prefix,
            genotype_format=genotype_format,
            phenotype_file=phenotype_file,
            phenotype_name=phenotype_name,
            phenotype_type=phenotype_type,
            covariate_file=covariate_file,
            covariate_names=covariates,
            anno_file=anno_file,
            set_list=set_list,
            mask_def=mask_def,
            output_prefix=regenie_prefix,
            vc_tests=["skato"],
            aaf_bins=aaf_bins,
            build_mask=build_mask,
            ignore_pred=ignore_pred,
            pred_file=pred_file,
            bsize=bsize,
            extra_args=extra_args,
            working_directory=project_root,
        )
        metadata = build_regenie_metadata(
            method="skato",
            engine=engine,
            phenotype_name=phenotype_name,
            phenotype_type=phenotype_type,
            genotype_format=genotype_format,
            genotype_prefix=genotype_prefix,
            pred_file=pred_file,
            phenotype_file=phenotype_file,
            covariate_file=covariate_file,
            anno_file=anno_file,
            set_list=set_list,
            mask_def=mask_def,
            result_file=regenie_result,
        )
        metadata["result_file"] = str(regenie_result)
        write_metadata(metadata_json, metadata)
        return regenie_result

    genotype_matrix = project_root / get_required(config, "input", "genotype_matrix")
    covariate_file = project_root / get_required(config, "input", "covariate_file")
    covariates = get_required(config, "skato", "covariates")
    _ensure_exists(genotype_matrix, "Genotype matrix")
    _ensure_exists(covariate_file, "Covariate file")

    covariates_csv = ",".join(covariates)

    r_script = project_root / "scripts" / "r" / "skato.R"
    command = [
        "Rscript",
        str(r_script),
        "--geno", str(genotype_matrix),
        "--pheno", str(phenotype_file),
        "--covar", str(covariate_file),
        "--pheno-name", str(phenotype_name),
        "--pheno-type", str(phenotype_type),
        "--set-id", str(set_id),
        "--covariates", covariates_csv,
        "--result-out", str(result_tsv),
    ]

    _run_command(
        command,
        cwd=project_root,
        tool_name="SKAT-O R",
        install_hint="Install R and make sure 'Rscript' is available on PATH.",
    )

    metadata = {
        "method": "skato",
        "engine": engine,
        "phenotype_name": phenotype_name,
        "phenotype_type": phenotype_type,
        "genotype_matrix": str(genotype_matrix),
        "phenotype_file": str(phenotype_file),
        "covariate_file": str(covariate_file),
        "result_file": str(result_tsv),
    }
    metadata_json.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return result_tsv


def main() -> int:
    raise SystemExit("This module is intended to be imported by main.py, not run directly.")


if __name__ == "__main__":
    sys.exit(main())
