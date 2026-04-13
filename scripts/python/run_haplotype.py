from __future__ import annotations

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


def run_haplotype(config: dict[str, Any], project_root: Path) -> Path:
    genotype_matrix = project_root / get_required(config, "input", "genotype_matrix")
    phenotype_file = project_root / get_required(config, "input", "phenotype_file")
    phenotype_name = get_required(config, "phenotype", "name")
    phenotype_type = get_required(config, "phenotype", "type")
    hap_set_id = get_required(config, "haplotype", "set_id")
    covariates = config.get("haplotype", {}).get("covariates", [])
    covariate_file_value = config.get("input", {}).get("covariate_file")
    covariate_file = project_root / covariate_file_value if covariate_file_value else None
    output_root = resolve_output_root(project_root, get_required(config, "output", "root_dir"))
    output_dir = output_root / "haplotype"

    _ensure_exists(genotype_matrix, "Genotype matrix")
    _ensure_exists(phenotype_file, "Phenotype file")
    if covariate_file is not None:
        _ensure_exists(covariate_file, "Covariate file")

    output_dir.mkdir(parents=True, exist_ok=True)

    covariates_csv = ",".join(covariates)
    result_tsv = output_dir / "haplotype_result.tsv"
    freq_tsv = output_dir / "haplotype_frequency.tsv"
    metadata_json = output_dir / "run_metadata.json"

    r_script = project_root / "scripts" / "r" / "haplotype.R"
    command = [
        "Rscript",
        str(r_script),
        "--geno", str(genotype_matrix),
        "--pheno", str(phenotype_file),
        "--pheno-name", str(phenotype_name),
        "--pheno-type", str(phenotype_type),
        "--set-id", str(hap_set_id),
        "--covariates", covariates_csv,
        "--result-out", str(result_tsv),
        "--freq-out", str(freq_tsv),
    ]
    if covariate_file is not None:
        command.extend(["--covar", str(covariate_file)])

    _run_command(
        command,
        cwd=project_root,
        tool_name="Haplotype R",
        install_hint="Install R and make sure 'Rscript' is available on PATH.",
    )

    metadata = {
        "method": "haplotype",
        "phenotype_name": phenotype_name,
        "phenotype_type": phenotype_type,
        "genotype_matrix": str(genotype_matrix),
        "phenotype_file": str(phenotype_file),
        "covariate_file": str(covariate_file) if covariate_file else None,
        "result_file": str(result_tsv),
        "frequency_file": str(freq_tsv),
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
