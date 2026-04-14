from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover - depends on local environment
    yaml = None


def _expand_env_vars(value: Any) -> Any:
    """递归展开字符串中的环境变量，如 ${VAR} 或 $VAR"""
    if isinstance(value, str):
        # 匹配 ${VAR} 或 $VAR 格式
        def replace_env_var(match):
            var_expr = match.group(0)
            # 提取变量名
            if var_expr.startswith("${"):
                var_name = var_expr[2:-1]
            else:
                var_name = var_expr[1:]
            return os.environ.get(var_name, var_expr)
        return re.sub(r'\$\{[^}]+\}|\$[a-zA-Z_][a-zA-Z0-9_]*', replace_env_var, value)
    elif isinstance(value, dict):
        return {k: _expand_env_vars(v) for k, v in value.items()}
    elif isinstance(value, list):
        return [_expand_env_vars(item) for item in value]
    return value


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

    # 展开环境变量
    data = _expand_env_vars(data)

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

    for section_name in ("burden", "skato", "haplotype", "quant_assoc", "gwas_gene_catalog"):
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

    burden_cfg = config.get("burden", {})
    if isinstance(burden_cfg, dict):
        burden_engine = burden_cfg.get("engine")
        if burden_engine is not None and burden_engine not in {"skat", "regenie"}:
            raise ValueError(
                f"Invalid burden.engine in {config_path}: {burden_engine!r}. "
                "Supported values are 'skat' and 'regenie'."
            )

    skato_cfg = config.get("skato", {})
    if isinstance(skato_cfg, dict):
        skato_engine = skato_cfg.get("engine")
        if skato_engine is not None and skato_engine not in {"skat", "regenie"}:
            raise ValueError(
                f"Invalid skato.engine in {config_path}: {skato_engine!r}. "
                "Supported values are 'skat' and 'regenie'."
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

    gwas_gene_catalog_cfg = config.get("gwas_gene_catalog", {})
    if isinstance(gwas_gene_catalog_cfg, dict):
        match_fields = gwas_gene_catalog_cfg.get("match_fields")
        if match_fields is not None:
            allowed_fields = {"mapped_genes", "reported_genes"}
            if not isinstance(match_fields, list) or any(not isinstance(item, str) for item in match_fields):
                raise ValueError(
                    f"gwas_gene_catalog.match_fields must be a list of strings in {config_path}"
                )
            invalid_fields = sorted(set(match_fields) - allowed_fields)
            if invalid_fields:
                raise ValueError(
                    f"Unsupported gwas_gene_catalog.match_fields in {config_path}: {invalid_fields}. "
                    f"Supported values are {sorted(allowed_fields)}."
                )

        source_mode = gwas_gene_catalog_cfg.get("source_mode")
        if source_mode is not None and source_mode not in {"tsv", "api"}:
            raise ValueError(
                f"Invalid gwas_gene_catalog.source_mode in {config_path}: {source_mode!r}. "
                "Supported values are 'tsv' and 'api'."
            )

        api_base_url = gwas_gene_catalog_cfg.get("gwas_catalog_api_base_url")
        if api_base_url is not None and not isinstance(api_base_url, str):
            raise ValueError(
                f"gwas_gene_catalog.gwas_catalog_api_base_url must be a string in {config_path}"
            )

        api_page_size = gwas_gene_catalog_cfg.get("api_page_size")
        if api_page_size is not None and (not isinstance(api_page_size, int) or api_page_size <= 0):
            raise ValueError(
                f"gwas_gene_catalog.api_page_size must be a positive integer in {config_path}"
            )

        api_max_pages = gwas_gene_catalog_cfg.get("api_max_pages")
        if api_max_pages is not None and (not isinstance(api_max_pages, int) or api_max_pages <= 0):
            raise ValueError(
                f"gwas_gene_catalog.api_max_pages must be a positive integer in {config_path}"
            )

        api_extended_geneset = gwas_gene_catalog_cfg.get("api_extended_geneset")
        if api_extended_geneset is not None and not isinstance(api_extended_geneset, bool):
            raise ValueError(
                f"gwas_gene_catalog.api_extended_geneset must be a boolean in {config_path}"
            )

        if source_mode == "api" and match_fields is not None and "reported_genes" in match_fields:
            has_reference_file = bool(gwas_gene_catalog_cfg.get("gwas_reference_file"))
            has_tsv_url = bool(gwas_gene_catalog_cfg.get("gwas_catalog_tsv_url"))
            if not has_reference_file and not has_tsv_url:
                raise ValueError(
                    f"API mode with gwas_gene_catalog.match_fields including 'reported_genes' "
                    f"requires gwas_reference_file or gwas_catalog_tsv_url in {config_path}, "
                    "because the official GWAS Catalog API only supports direct mapped_gene queries."
                )
