from __future__ import annotations

import csv
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from config_utils import get_required


def _ensure_exists(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def _download_file(url: str, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as response, destination.open("wb") as handle:
        handle.write(response.read())
    return destination


def _load_json(url: str) -> dict[str, Any]:
    with urllib.request.urlopen(url) as response:
        content = response.read().decode("utf-8")
    data = json.loads(content)
    if not isinstance(data, dict):
        raise ValueError(f"Expected JSON object from {url}")
    return data


def _detect_delimiter(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def _read_table(path: Path) -> list[dict[str, str]]:
    delimiter = _detect_delimiter(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"Input table has no header: {path}")
        return [{key: value or "" for key, value in row.items()} for row in reader]


def _write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def _get_first_present(row: dict[str, Any], *keys: str) -> str:
    for key in keys:
        value = row.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def _normalize_gene_token(token: str) -> str:
    return re.sub(r"\s+", "", token).upper()


def _split_gene_field(value: str) -> list[str]:
    if not value:
        return []
    tokens = re.split(r"\s*[,;/]\s*|\s+-\s+|\s+x\s+", value)
    return [token.strip() for token in tokens if token and token.strip()]


def _parse_pvalue(raw_value: str) -> float | None:
    value = raw_value.strip()
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        pass

    sci_match = re.fullmatch(r"([0-9.]+)\s*[xX*]\s*10\^?(-?[0-9]+)", value)
    if sci_match:
        return float(sci_match.group(1)) * (10 ** int(sci_match.group(2)))
    return None


def _format_pvalue(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6g}"


def _extract_strings(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value] if value.strip() else []
    if isinstance(value, (int, float)):
        return [str(value)]
    if isinstance(value, list):
        result: list[str] = []
        for item in value:
            result.extend(_extract_strings(item))
        return result
    if isinstance(value, dict):
        result: list[str] = []
        preferred_keys = (
            "trait",
            "label",
            "name",
            "value",
            "mappedGene",
            "reportedGene",
            "geneName",
            "rsId",
            "variantId",
            "accessionId",
            "fullName",
            "shortForm",
        )
        for key in preferred_keys:
            if key in value:
                result.extend(_extract_strings(value[key]))
        if result:
            return result
        for item in value.values():
            result.extend(_extract_strings(item))
        return result
    return []


def _normalize_catalog_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    normalized_rows: list[dict[str, Any]] = []
    for row in rows:
        mapped_genes = _split_gene_field(
            _get_first_present(row, "MAPPED_GENE", "MAPPED_GENE(S)", "mapped_gene", "mapped_genes")
        )
        reported_genes = _split_gene_field(
            _get_first_present(row, "REPORTED GENE(S)", "reported_gene", "reported_genes")
        )
        normalized_rows.append(
            {
                "mapped_genes": mapped_genes,
                "reported_genes": reported_genes,
                "trait": _get_first_present(row, "MAPPED_TRAIT", "DISEASE/TRAIT", "trait"),
                "reported_trait": _get_first_present(row, "DISEASE/TRAIT", "reported_trait"),
                "variant": _get_first_present(row, "SNPS", "SNP_ID_CURRENT", "variant"),
                "pvalue": _parse_pvalue(_get_first_present(row, "P-VALUE", "pvalue")),
                "study_accession": _get_first_present(row, "STUDY ACCESSION", "study_accession"),
                "first_author": _get_first_present(row, "FIRST AUTHOR", "first_author"),
                "pubmed_id": _get_first_present(row, "PUBMEDID", "pubmed_id"),
            }
        )
    return normalized_rows


def _extract_api_gene_list(item: Any, fallback_text: str = "") -> list[str]:
    values = _extract_strings(item)
    if not values and fallback_text:
        values = [fallback_text]
    result: list[str] = []
    for value in values:
        result.extend(_split_gene_field(value))
    return result


def _normalize_api_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    normalized_rows: list[dict[str, Any]] = []
    for row in rows:
        mapped_gene_field = row.get("mappedGenes") or row.get("mappedGene") or row.get("mapped_genes")
        reported_gene_field = (
            row.get("reportedGenes")
            or row.get("reportedGene")
            or row.get("reported_genes")
            or row.get("reportedTrait")  # harmless fallback when genes are absent
        )
        trait_field = row.get("diseaseTrait") or row.get("trait") or row.get("efoTraits") or row.get("efoTrait")
        variant_field = row.get("variantId") or row.get("rsId") or row.get("snps")
        study_field = row.get("studyAccession") or row.get("study") or row.get("study accession")
        author_field = row.get("firstAuthor") or row.get("author")
        pubmed_field = row.get("pubmedId") or row.get("pubmed_id")

        study_values = _extract_strings(study_field)
        study_accession = ""
        if study_values:
            study_accession = study_values[0]

        first_author_values = _extract_strings(author_field)
        first_author = first_author_values[0] if first_author_values else ""
        pubmed_values = _extract_strings(pubmed_field)
        pubmed_id = pubmed_values[0] if pubmed_values else ""

        pvalue = row.get("pvalue")
        if pvalue is None:
            pvalue = row.get("pValue")
        pvalue_value = _parse_pvalue(str(pvalue)) if pvalue is not None else None

        normalized_rows.append(
            {
                "mapped_genes": _extract_api_gene_list(mapped_gene_field),
                "reported_genes": _extract_api_gene_list(reported_gene_field),
                "trait": "; ".join(dict.fromkeys(_extract_strings(trait_field))),
                "reported_trait": "; ".join(
                    dict.fromkeys(_extract_strings(row.get("reportedTrait") or row.get("diseaseTrait")))
                ),
                "variant": "; ".join(dict.fromkeys(_extract_strings(variant_field))),
                "pvalue": pvalue_value,
                "study_accession": study_accession,
                "first_author": first_author,
                "pubmed_id": pubmed_id,
            }
        )
    return normalized_rows


def _build_gene_index(
    rows: list[dict[str, Any]],
    match_fields: list[str],
) -> dict[str, list[dict[str, Any]]]:
    gene_index: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        matched_tokens: set[str] = set()
        for field in match_fields:
            for gene in row.get(field, []):
                normalized_gene = _normalize_gene_token(gene)
                if normalized_gene:
                    matched_tokens.add(normalized_gene)
        for token in matched_tokens:
            gene_index.setdefault(token, []).append(row)
    return gene_index


def _unique_join(values: list[str]) -> str:
    seen: dict[str, None] = {}
    for value in values:
        cleaned = value.strip()
        if cleaned and cleaned not in seen:
            seen[cleaned] = None
    return "; ".join(seen.keys())


def _summarize_matches(gene: str, matches: list[dict[str, Any]]) -> dict[str, str]:
    if not matches:
        return {
            "gwas_catalog_found": "FALSE",
            "gwas_catalog_result_count": "0",
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
        key=lambda row: (row["pvalue"] is None, row["pvalue"] if row["pvalue"] is not None else float("inf")),
    )
    top_study_parts = [part for part in [top_match["study_accession"], top_match["first_author"]] if part]
    if not top_study_parts:
        top_study = ""
    elif len(top_study_parts) == 1:
        top_study = top_study_parts[0]
    else:
        top_study = f"{top_study_parts[0]} ({top_study_parts[1]})"
    return {
        "gwas_catalog_found": "TRUE",
        "gwas_catalog_result_count": str(len(matches)),
        "gwas_catalog_traits": _unique_join([row["trait"] for row in matches]),
        "gwas_catalog_reported_traits": _unique_join([row["reported_trait"] for row in matches]),
        "gwas_catalog_mapped_genes": _unique_join(
            [gene_name for row in matches for gene_name in row.get("mapped_genes", [])]
        ),
        "gwas_catalog_top_variant": top_match["variant"],
        "gwas_catalog_top_pvalue": _format_pvalue(top_match["pvalue"]),
        "gwas_catalog_top_study": top_study,
        "gwas_catalog_study_accessions": _unique_join([row["study_accession"] for row in matches]),
    }


def _load_catalog_source(config: dict[str, Any], project_root: Path) -> tuple[str, Path]:
    module_cfg = config.get("gwas_gene_catalog", {})
    source_mode = str(module_cfg.get("source_mode", "tsv")).strip().lower()

    if source_mode == "api":
        return "api", Path()

    reference_file = module_cfg.get("gwas_reference_file")
    if reference_file:
        source_path = project_root / str(reference_file)
        _ensure_exists(source_path, "GWAS Catalog reference file")
        return "local_file", source_path

    tsv_url = module_cfg.get("gwas_catalog_tsv_url")
    if tsv_url:
        cache_dir = project_root / "output" / "gwas_gene_catalog" / "_downloads"
        source_path = _download_file(str(tsv_url), cache_dir / "gwas_catalog_reference.tsv")
        return "remote_download", source_path

    raise KeyError(
        "Missing GWAS Catalog source. Provide one of "
        "'gwas_gene_catalog.gwas_reference_file', "
        "'gwas_gene_catalog.gwas_catalog_tsv_url', or set "
        "'gwas_gene_catalog.source_mode: api'."
    )


def _get_next_link(data: dict[str, Any]) -> str:
    links = data.get("_links")
    if not isinstance(links, dict):
        return ""
    next_link = links.get("next")
    if isinstance(next_link, dict):
        href = next_link.get("href")
        if isinstance(href, str):
            return href
    return ""


def _get_embedded_associations(data: dict[str, Any]) -> list[dict[str, Any]]:
    embedded = data.get("_embedded")
    if not isinstance(embedded, dict):
        return []
    for key in ("associations", "associationList", "singleNucleotidePolymorphisms"):
        value = embedded.get(key)
        if isinstance(value, list):
            return [item for item in value if isinstance(item, dict)]
    return []


def _fetch_api_matches_for_gene(
    gene: str,
    api_base_url: str,
    api_page_size: int,
    api_extended_geneset: bool,
    api_max_pages: int | None,
) -> list[dict[str, Any]]:
    encoded_gene = urllib.parse.quote(gene)
    next_url = (
        f"{api_base_url}?mapped_gene={encoded_gene}"
        f"&size={api_page_size}"
        f"&extended_geneset={'true' if api_extended_geneset else 'false'}"
    )
    api_rows: list[dict[str, Any]] = []
    visited_urls: set[str] = set()
    page_count = 0

    while next_url and next_url not in visited_urls:
        if api_max_pages is not None and page_count >= api_max_pages:
            break
        visited_urls.add(next_url)
        payload = _load_json(next_url)
        api_rows.extend(_get_embedded_associations(payload))
        next_url = _get_next_link(payload)
        page_count += 1

    return _normalize_api_rows(api_rows)


def run_gwas_gene_catalog(config: dict[str, Any], project_root: Path) -> Path:
    module_cfg = config.get("gwas_gene_catalog", {})
    project_result_file = project_root / str(get_required(module_cfg, "project_result_file"))
    gene_column = str(module_cfg.get("gene_column", "gene"))
    match_fields = list(module_cfg.get("match_fields", ["mapped_genes", "reported_genes"]))
    trait_keyword = str(module_cfg.get("trait_keyword", "")).strip().lower()
    api_base_url = str(
        module_cfg.get("gwas_catalog_api_base_url", "https://www.ebi.ac.uk/gwas/rest/api/v2/associations")
    ).strip()
    api_page_size = int(module_cfg.get("api_page_size", 100))
    api_extended_geneset = bool(module_cfg.get("api_extended_geneset", False))
    api_max_pages_raw = module_cfg.get("api_max_pages")
    api_max_pages = int(api_max_pages_raw) if api_max_pages_raw is not None else None

    _ensure_exists(project_result_file, "Project result file")

    output_root = project_root / str(config.get("output", {}).get("root_dir", "output"))
    output_dir = output_root / "gwas_gene_catalog"
    output_dir.mkdir(parents=True, exist_ok=True)

    project_rows = _read_table(project_result_file)
    if not project_rows:
        raise ValueError(f"Project result file is empty: {project_result_file}")
    if gene_column not in project_rows[0]:
        raise KeyError(f"Gene column '{gene_column}' not found in {project_result_file}")

    source_kind, source_path = _load_catalog_source(config, project_root)
    normalized_catalog_rows: list[dict[str, Any]] = []
    gene_index: dict[str, list[dict[str, Any]]] = {}

    if source_kind == "api":
        catalog_source_label = api_base_url
    else:
        catalog_rows = _read_table(source_path)
        normalized_catalog_rows = _normalize_catalog_rows(catalog_rows)
        gene_index = _build_gene_index(normalized_catalog_rows, match_fields)
        catalog_source_label = str(source_path)

    enriched_rows: list[dict[str, Any]] = []
    detail_rows: list[dict[str, Any]] = []

    for row_index, row in enumerate(project_rows, start=1):
        gene_value = str(row.get(gene_column, "")).strip()
        normalized_gene = _normalize_gene_token(gene_value)
        matches = gene_index.get(normalized_gene, []) if source_kind != "api" else _fetch_api_matches_for_gene(
            gene_value,
            api_base_url=api_base_url,
            api_page_size=api_page_size,
            api_extended_geneset=api_extended_geneset,
            api_max_pages=api_max_pages,
        )

        if trait_keyword:
            matches = [
                match
                for match in matches
                if trait_keyword in match.get("trait", "").lower()
                or trait_keyword in match.get("reported_trait", "").lower()
            ]

        summary = _summarize_matches(gene_value, matches)
        enriched_row = dict(row)
        enriched_row.update(summary)
        enriched_rows.append(enriched_row)

        for match in matches:
            detail_rows.append(
                {
                    "project_row_index": str(row_index),
                    "gene": gene_value,
                    "gwas_catalog_trait": match["trait"],
                    "gwas_catalog_reported_trait": match["reported_trait"],
                    "gwas_catalog_mapped_genes": _unique_join(match.get("mapped_genes", [])),
                    "gwas_catalog_reported_genes": _unique_join(match.get("reported_genes", [])),
                    "gwas_catalog_variant": match["variant"],
                    "gwas_catalog_pvalue": _format_pvalue(match["pvalue"]),
                    "gwas_catalog_study_accession": match["study_accession"],
                    "gwas_catalog_first_author": match["first_author"],
                    "gwas_catalog_pubmed_id": match["pubmed_id"],
                }
            )

    summary_columns = list(project_rows[0].keys()) + [
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
    detail_columns = [
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

    summary_path = output_dir / "project_results_with_gwas_catalog.tsv"
    detail_path = output_dir / "gwas_catalog_gene_hits.tsv"
    metadata_path = output_dir / "run_metadata.json"

    _write_tsv(summary_path, enriched_rows, summary_columns)
    _write_tsv(detail_path, detail_rows, detail_columns)

    metadata = {
        "method": "gwas_gene_catalog",
        "project_result_file": str(project_result_file),
        "gene_column": gene_column,
        "catalog_source_kind": source_kind,
        "catalog_source": catalog_source_label,
        "result_file": str(summary_path),
        "detail_file": str(detail_path),
        "n_input_rows": len(project_rows),
        "n_detail_rows": len(detail_rows),
        "match_fields": match_fields,
        "trait_keyword": trait_keyword,
    }
    if source_kind == "api":
        metadata["gwas_catalog_api_base_url"] = api_base_url
        metadata["api_page_size"] = api_page_size
        metadata["api_extended_geneset"] = api_extended_geneset
        metadata["api_max_pages"] = api_max_pages

    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary_path
