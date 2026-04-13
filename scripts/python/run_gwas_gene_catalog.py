from __future__ import annotations

import csv
import json
import re
import sys
import urllib.request
from pathlib import Path
from typing import Any

from config_utils import get_required
from path_utils import resolve_output_root


def _ensure_exists(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def _download_file(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=300) as response:
        data = response.read()
    target.write_bytes(data)


def _detect_delimiter(path: Path) -> str:
    if path.suffix.lower() == ".csv":
        return ","
    return "\t"


def _read_table(path: Path) -> list[dict[str, str]]:
    delimiter = _detect_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def _write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _get_first_present(row: dict[str, str], candidates: list[str]) -> str:
    for key in candidates:
        value = row.get(key)
        if value not in (None, ""):
            return value
    return ""


def _normalize_gene_token(value: str) -> str:
    return value.strip().upper()


def _split_gene_field(value: str) -> list[str]:
    if not value:
        return []
    normalized = value.strip()
    if not normalized or normalized.upper() in {"NR", "NA", "N/A", "NONE", "INTERGENIC"}:
        return []

    tokens = re.split(r"\s*[,;/]\s*|\s+-\s+|\s+x\s+", normalized)
    cleaned: list[str] = []
    for token in tokens:
        gene = _normalize_gene_token(token)
        if gene and gene not in {"NR", "NA", "N/A", "NONE", "INTERGENIC"}:
            cleaned.append(gene)
    return cleaned


def _parse_pvalue(value: str) -> float | None:
    if value is None:
        return None
    text = value.strip()
    if not text or text.upper() in {"NA", "N/A"}:
        return None
    try:
        return float(text)
    except ValueError:
        pass

    normalized = text.replace("×", "x").replace("*", "x")
    match = re.fullmatch(r"([0-9]*\.?[0-9]+)\s*x\s*10\^?(-?[0-9]+)", normalized, flags=re.IGNORECASE)
    if match:
        base = float(match.group(1))
        exponent = int(match.group(2))
        return base * (10 ** exponent)
    return None


def _format_pvalue(value: float | None, raw_value: str) -> str:
    if value is None:
        return raw_value
    return f"{value:.12g}"


def _load_catalog_source(config: dict[str, Any], project_root: Path, work_dir: Path) -> tuple[Path, str]:
    catalog_cfg = config.get("gwas_gene_catalog", {})
    if "gwas_reference_file" in catalog_cfg:
        ref_path = project_root / catalog_cfg["gwas_reference_file"]
        _ensure_exists(ref_path, "GWAS Catalog reference file")
        return ref_path, "local_file"
    if "gwas_catalog_tsv_url" in catalog_cfg:
        target = work_dir / "gwas_catalog_download.tsv"
        _download_file(catalog_cfg["gwas_catalog_tsv_url"], target)
        return target, "remote_download"
    raise KeyError(
        "GWAS gene catalog annotation requires either "
        "gwas_gene_catalog.gwas_reference_file or gwas_gene_catalog.gwas_catalog_tsv_url"
    )


def _normalize_catalog_rows(
    source_rows: list[dict[str, str]],
    *,
    trait_keyword: str | None,
) -> list[dict[str, Any]]:
    normalized_rows: list[dict[str, Any]] = []
    for row in source_rows:
        mapped_genes_raw = _get_first_present(
            row,
            ["MAPPED_GENE", "MAPPED_GENE(S)", "Mapped gene", "MAPPED GENES"],
        )
        reported_genes_raw = _get_first_present(
            row,
            ["REPORTED GENE(S)", "Reported gene(s)", "REPORTED GENES"],
        )
        trait_value = _get_first_present(
            row,
            ["MAPPED_TRAIT", "MAPPED TRAIT", "DISEASE/TRAIT", "Trait(s)", "Trait"],
        )
        reported_trait_value = _get_first_present(
            row,
            ["DISEASE/TRAIT", "Reported trait", "REPORTED TRAIT", "REPORTED TRAIT(S)"],
        )
        if trait_keyword:
            haystacks = " ".join([trait_value, reported_trait_value]).lower()
            if trait_keyword.lower() not in haystacks:
                continue

        pvalue_raw = _get_first_present(row, ["P-VALUE", "P VALUE", "P-VALUE (TEXT)", "P-value"])
        pvalue_numeric = _parse_pvalue(pvalue_raw)
        variant_value = _get_first_present(
            row,
            ["SNPS", "SNP_ID_CURRENT", "STRONGEST SNP-RISK ALLELE", "Variant and risk allele"],
        )
        accession_value = _get_first_present(row, ["STUDY ACCESSION", "Study accession", "ACCESSION"])
        first_author = _get_first_present(row, ["FIRST AUTHOR", "First author"])
        pubmed_id = _get_first_present(row, ["PUBMEDID", "PUBMED ID", "PubMed ID"])

        normalized_rows.append(
            {
                "mapped_genes_raw": mapped_genes_raw,
                "reported_genes_raw": reported_genes_raw,
                "mapped_gene_tokens": _split_gene_field(mapped_genes_raw),
                "reported_gene_tokens": _split_gene_field(reported_genes_raw),
                "trait": trait_value,
                "reported_trait": reported_trait_value,
                "variant": variant_value,
                "pvalue_raw": pvalue_raw,
                "pvalue_numeric": pvalue_numeric,
                "study_accession": accession_value,
                "first_author": first_author,
                "pubmed_id": pubmed_id,
            }
        )
    return normalized_rows


def _build_gene_index(rows: list[dict[str, Any]], match_fields: list[str]) -> dict[str, list[dict[str, Any]]]:
    index: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        tokens: list[str] = []
        if "mapped_genes" in match_fields:
            tokens.extend(row["mapped_gene_tokens"])
        if "reported_genes" in match_fields:
            tokens.extend(row["reported_gene_tokens"])
        for token in sorted(set(tokens)):
            index.setdefault(token, []).append(row)
    return index


def _unique_join(values: list[str]) -> str:
    seen: list[str] = []
    for value in values:
        cleaned = value.strip()
        if cleaned and cleaned not in seen:
            seen.append(cleaned)
    return "; ".join(seen)


def _summarize_matches(gene: str, matches: list[dict[str, Any]]) -> dict[str, Any]:
    if not matches:
        return {
            "gene": gene,
            "gwas_catalog_found": "FALSE",
            "gwas_catalog_result_count": 0,
            "gwas_catalog_traits": "",
            "gwas_catalog_reported_traits": "",
            "gwas_catalog_mapped_genes": "",
            "gwas_catalog_top_variant": "",
            "gwas_catalog_top_pvalue": "",
            "gwas_catalog_top_study": "",
            "gwas_catalog_study_accessions": "",
        }

    top_match = min(
        matches,
        key=lambda row: (row["pvalue_numeric"] is None, row["pvalue_numeric"] if row["pvalue_numeric"] is not None else float("inf")),
    )
    top_study = top_match["study_accession"]
    if top_match["first_author"]:
        top_study = f"{top_study} ({top_match['first_author']})" if top_study else top_match["first_author"]

    return {
        "gene": gene,
        "gwas_catalog_found": "TRUE",
        "gwas_catalog_result_count": len(matches),
        "gwas_catalog_traits": _unique_join([row["trait"] for row in matches]),
        "gwas_catalog_reported_traits": _unique_join([row["reported_trait"] for row in matches]),
        "gwas_catalog_mapped_genes": _unique_join([row["mapped_genes_raw"] for row in matches]),
        "gwas_catalog_top_variant": top_match["variant"],
        "gwas_catalog_top_pvalue": _format_pvalue(top_match["pvalue_numeric"], top_match["pvalue_raw"]),
        "gwas_catalog_top_study": top_study,
        "gwas_catalog_study_accessions": _unique_join([row["study_accession"] for row in matches]),
    }


def run_gwas_gene_catalog(config: dict[str, Any], project_root: Path) -> Path:
    catalog_cfg = config.get("gwas_gene_catalog", {})
    project_result = project_root / get_required(config, "gwas_gene_catalog", "project_result_file")
    gene_column = catalog_cfg.get("gene_column", "gene")
    match_fields = catalog_cfg.get("match_fields", ["mapped_genes", "reported_genes"])
    trait_keyword = catalog_cfg.get("trait_keyword")
    output_root = resolve_output_root(project_root, get_required(config, "output", "root_dir"))
    output_dir = output_root / "gwas_gene_catalog"

    _ensure_exists(project_result, "Project result file")
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = output_dir / "tmp_catalog"
    work_dir.mkdir(parents=True, exist_ok=True)

    project_rows = _read_table(project_result)
    if not project_rows:
        raise ValueError(f"Project result file has no data rows: {project_result}")
    if gene_column not in project_rows[0]:
        raise KeyError(f"Gene column not found in project result file: {gene_column}")

    source_path, source_mode = _load_catalog_source(config, project_root, work_dir)
    catalog_rows = _normalize_catalog_rows(_read_table(source_path), trait_keyword=trait_keyword)
    gene_index = _build_gene_index(catalog_rows, match_fields=match_fields)

    summary_columns = [
        "gwas_catalog_found",
        "gwas_catalog_result_count",
        "gwas_catalog_traits",
        "gwas_catalog_reported_traits",
        "gwas_catalog_mapped_genes",
        "gwas_catalog_top_variant",
        "gwas_catalog_top_pvalue",
        "gwas_catalog_top_study",
        "gwas_catalog_study_accessions",
    ]

    annotated_rows: list[dict[str, Any]] = []
    detail_rows: list[dict[str, Any]] = []
    for row_index, row in enumerate(project_rows, start=1):
        gene_value = row.get(gene_column, "")
        normalized_gene = _normalize_gene_token(gene_value)
        matches = gene_index.get(normalized_gene, []) if normalized_gene else []
        summary = _summarize_matches(gene_value, matches)

        annotated_row = dict(row)
        for column in summary_columns:
            annotated_row[column] = summary[column]
        annotated_rows.append(annotated_row)

        for match in matches:
            detail_rows.append(
                {
                    "project_row_index": row_index,
                    "gene": gene_value,
                    "gwas_catalog_trait": match["trait"],
                    "gwas_catalog_reported_trait": match["reported_trait"],
                    "gwas_catalog_mapped_genes": match["mapped_genes_raw"],
                    "gwas_catalog_reported_genes": match["reported_genes_raw"],
                    "gwas_catalog_variant": match["variant"],
                    "gwas_catalog_pvalue": _format_pvalue(match["pvalue_numeric"], match["pvalue_raw"]),
                    "gwas_catalog_study_accession": match["study_accession"],
                    "gwas_catalog_first_author": match["first_author"],
                    "gwas_catalog_pubmed_id": match["pubmed_id"],
                }
            )

    annotated_path = output_dir / "project_results_with_gwas_catalog.tsv"
    detail_path = output_dir / "gwas_catalog_gene_hits.tsv"
    metadata_path = output_dir / "run_metadata.json"

    project_fieldnames = list(project_rows[0].keys())
    annotated_fieldnames = project_fieldnames + [column for column in summary_columns if column not in project_fieldnames]
    _write_tsv(annotated_path, annotated_rows, annotated_fieldnames)

    detail_fieldnames = [
        "project_row_index",
        "gene",
        "gwas_catalog_trait",
        "gwas_catalog_reported_trait",
        "gwas_catalog_mapped_genes",
        "gwas_catalog_reported_genes",
        "gwas_catalog_variant",
        "gwas_catalog_pvalue",
        "gwas_catalog_study_accession",
        "gwas_catalog_first_author",
        "gwas_catalog_pubmed_id",
    ]
    _write_tsv(detail_path, detail_rows, detail_fieldnames)

    metadata = {
        "method": "gwas_gene_catalog",
        "project_result_file": str(project_result),
        "gene_column": gene_column,
        "match_fields": match_fields,
        "trait_keyword": trait_keyword,
        "source_mode": source_mode,
        "gwas_reference_source": str(source_path),
        "annotated_result_file": str(annotated_path),
        "detail_result_file": str(detail_path),
        "result_file": str(annotated_path),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    return annotated_path


def main() -> int:
    raise SystemExit("This module is intended to be imported by main.py, not run directly.")


if __name__ == "__main__":
    sys.exit(main())
