from __future__ import annotations

from pathlib import Path
from typing import Any

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover - depends on local environment
    yaml = None


def load_config(config_path: str | Path) -> dict[str, Any]:
    if yaml is None:
        raise RuntimeError(
            "PyYAML is not installed. Install Python dependencies first, for example:\n"
            "  pip install -r requirements.txt"
        )

    path = Path(config_path)
    try:
        with path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid YAML config file: {path}\n{exc}") from exc

    if not isinstance(data, dict):
        raise ValueError("Config root must be a mapping.")
    validate_config(data, path)
    return data


def get_required(config: dict[str, Any], *keys: str) -> Any:
    current: Any = config
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            dotted = ".".join(keys)
            raise KeyError(f"Missing required config key: {dotted}")
        current = current[key]
    return current


def validate_config(config: dict[str, Any], config_path: Path) -> None:
    phenotype = config.get("phenotype", {})
    phenotype_type = phenotype.get("type")
    if phenotype_type is not None and phenotype_type not in {"continuous", "binary"}:
        raise ValueError(
            f"Invalid phenotype.type in {config_path}: {phenotype_type!r}. "
            "Supported values are 'continuous' and 'binary'."
        )

    output_cfg = config.get("output", {})
    root_dir = output_cfg.get("root_dir")
    if root_dir is not None and not isinstance(root_dir, str):
        raise ValueError(f"output.root_dir must be a string in {config_path}")

    regenie_cfg = config.get("regenie", {})
    genotype_format = regenie_cfg.get("genotype_format")
    if genotype_format is not None and genotype_format not in {"pfile", "bfile"}:
        raise ValueError(
            f"Invalid regenie.genotype_format in {config_path}: {genotype_format!r}. "
            "Supported values are 'pfile' and 'bfile'."
        )

    for section_name in ("burden", "skato", "haplotype", "quant_assoc"):
        section = config.get(section_name)
        if section is None:
            continue
        if not isinstance(section, dict):
            raise ValueError(f"Config section {section_name} must be a mapping in {config_path}")

        covariates = section.get("covariates")
        if covariates is not None:
            if not isinstance(covariates, list) or any(not isinstance(item, str) for item in covariates):
                raise ValueError(
                    f"{section_name}.covariates must be a list of strings in {config_path}"
                )

    input_cfg = config.get("input", {})
    if not isinstance(input_cfg, dict):
        raise ValueError(f"Config section input must be a mapping in {config_path}")

    quant_cfg = config.get("quant_assoc", {})
    if isinstance(quant_cfg, dict) and quant_cfg.get("engine", "plink2") == "plink2":
        quant_input_keys = {"pfile_prefix", "bfile_prefix", "pedmap_prefix"}
        present_keys = sorted(key for key in quant_input_keys if key in input_cfg)
        if len(present_keys) > 1:
            raise ValueError(
                f"Quant association input is ambiguous in {config_path}: "
                f"found multiple input prefixes {present_keys}. Use only one."
            )
