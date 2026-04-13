from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import urllib.request
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


def _download_file(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=300) as response:
        data = response.read()
    target.write_bytes(data)


def _load_gwas_reference(config: dict[str, Any], project_root: Path, work_dir: Path) -> tuple[Path, str]:
    gwas_cfg = config.get("gwas_overlap", {})
    if "gwas_reference_file" in gwas_cfg:
        ref_path = project_root / gwas_cfg["gwas_reference_file"]
        _ensure_exists(ref_path, "GWAS reference file")
        return ref_path, "local_file"
    if "gwas_catalog_tsv_url" in gwas_cfg:
        target = work_dir / "gwas_catalog_download.tsv"
        _download_file(gwas_cfg["gwas_catalog_tsv_url"], target)
        return target, "remote_download"
    raise KeyError(
        "GWAS overlap requires either gwas_overlap.gwas_reference_file "
        "or gwas_overlap.gwas_catalog_tsv_url"
    )


def _normalize_project_hits(input_path: Path, output_path: Path) -> list[dict[str, str]]:
    rows = _read_tsv(input_path)
    normalized: list[dict[str, str]] = []
    for row in rows:
        normalized.append(
            {
                "variant_id": row["variant_id"],
                "chr": row["chr"].replace("chr", ""),
                "pos": row["pos"],
            }
        )
    _write_tsv(output_path, normalized, ["variant_id", "chr", "pos"])
    return normalized


def _normalize_gwas_reference(
    input_path: Path,
    output_path: Path,
    trait_keyword: str | None,
) -> list[dict[str, str]]:
    rows = _read_tsv(input_path)
    normalized: list[dict[str, str]] = []
    for row in rows:
        if "chr" not in row or "pos" not in row:
            continue
        trait_value = row.get("trait", row.get("DISEASE/TRAIT", ""))
        if trait_keyword and trait_keyword.lower() not in trait_value.lower():
            continue
        normalized.append(
            {
                "variant_id": row.get("variant_id", row.get("SNP_ID_CURRENT", row.get("SNPS", ""))),
                "chr": row["chr"].replace("chr", ""),
                "pos": row["pos"],
                "trait": trait_value,
            }
        )
    _write_tsv(output_path, normalized, ["variant_id", "chr", "pos", "trait"])
    return normalized


def _write_bed(
    rows: list[dict[str, str]],
    path: Path,
    window_bp: int = 0,
    include_trait: bool = False,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        for row in rows:
            pos = int(row["pos"])
            start = max(0, pos - 1 - window_bp)
            end = pos + window_bp
            if include_trait:
                handle.write(
                    f"chr{row['chr']}\t{start}\t{end}\t{row['variant_id']}\t{row.get('trait', '')}\n"
                )
            else:
                handle.write(f"chr{row['chr']}\t{start}\t{end}\t{row['variant_id']}\n")


def _parse_bedtools_output(path: Path, project_pos_lookup: dict[str, tuple[str, int]], gwas_pos_lookup: dict[str, tuple[str, int, str]]) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            a_chr, a_start, a_end, a_id = parts[0], parts[1], parts[2], parts[3]
            b_chr, b_start, b_end, b_id, b_trait = parts[4], parts[5], parts[6], parts[7], parts[8]
            project_chr, project_pos = project_pos_lookup[a_id]
            if b_chr == ".":
                results.append(
                    {
                        "project_id": a_id,
                        "chr": project_chr,
                        "pos": project_pos,
                        "nearest_gwas_hit": "",
                        "nearest_gwas_trait": "",
                        "distance_bp": "",
                        "category": "novel",
                    }
                )
                continue
            _, gwas_pos, gwas_trait = gwas_pos_lookup[b_id]
            distance = abs(project_pos - gwas_pos)
            category = "known_exact" if distance == 0 else "near_known"
            results.append(
                {
                    "project_id": a_id,
                    "chr": project_chr,
                    "pos": project_pos,
                    "nearest_gwas_hit": b_id,
                    "nearest_gwas_trait": gwas_trait,
                    "distance_bp": distance,
                    "category": category,
                }
            )
    return results


def _deduplicate_results(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    best: dict[str, dict[str, Any]] = {}
    for row in rows:
        key = row["project_id"]
        current = best.get(key)
        if current is None:
            best[key] = row
            continue
        current_distance = current["distance_bp"]
        row_distance = row["distance_bp"]
        current_numeric = int(current_distance) if str(current_distance).isdigit() else 10**18
        row_numeric = int(row_distance) if str(row_distance).isdigit() else 10**18
        if row_numeric < current_numeric:
            best[key] = row
    return list(best.values())


def run_gwas_overlap(config: dict[str, Any], project_root: Path) -> Path:
    gwas_cfg = config.get("gwas_overlap", {})
    project_result = project_root / get_required(config, "gwas_overlap", "project_result_file")
    window_kb = int(get_required(config, "gwas_overlap", "window_kb"))
    trait_keyword = gwas_cfg.get("trait_keyword")
    bedtools_bin = gwas_cfg.get("bedtools_bin", "bedtools")
    output_root = resolve_output_root(project_root, get_required(config, "output", "root_dir"))
    output_dir = output_root / "gwas_overlap"

    if shutil.which(bedtools_bin) is None:
        raise FileNotFoundError(
            f"bedtools executable not found: {bedtools_bin}. "
            "Install bedtools and ensure it is available on PATH, or set gwas_overlap.bedtools_bin."
        )

    _ensure_exists(project_result, "Project result file")
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = output_dir / "tmp_overlap"
    work_dir.mkdir(parents=True, exist_ok=True)

    gwas_reference_raw, source_mode = _load_gwas_reference(config, project_root, work_dir)

    project_norm = work_dir / "project_hits.normalized.tsv"
    gwas_norm = work_dir / "gwas_reference.normalized.tsv"
    project_rows = _normalize_project_hits(project_result, project_norm)
    gwas_rows = _normalize_gwas_reference(gwas_reference_raw, gwas_norm, trait_keyword)
    project_pos_lookup = {row["variant_id"]: (row["chr"], int(row["pos"])) for row in project_rows}
    gwas_pos_lookup = {
        row["variant_id"]: (row["chr"], int(row["pos"]), row.get("trait", ""))
        for row in gwas_rows
    }

    project_bed = work_dir / "project_hits.bed"
    gwas_bed = work_dir / "gwas_reference.bed"
    intersect_out = work_dir / "bedtools_intersect.tsv"

    _write_bed(project_rows, project_bed, window_bp=window_kb * 1000, include_trait=False)
    _write_bed(gwas_rows, gwas_bed, window_bp=0, include_trait=True)

    command = [
        bedtools_bin,
        "intersect",
        "-a",
        str(project_bed),
        "-b",
        str(gwas_bed),
        "-wao",
    ]
    with intersect_out.open("w", encoding="utf-8", newline="") as handle:
        completed = subprocess.run(command, cwd=project_root, stdout=handle, check=False)
    if completed.returncode != 0:
        raise RuntimeError("bedtools intersect failed. Check whether the BED inputs are valid.")

    raw_results = _parse_bedtools_output(intersect_out, project_pos_lookup, gwas_pos_lookup)
    final_results = _deduplicate_results(raw_results)

    overlap_tsv = output_dir / "gwas_overlap_result.tsv"
    _write_tsv(
        overlap_tsv,
        final_results,
        [
            "project_id",
            "chr",
            "pos",
            "nearest_gwas_hit",
            "nearest_gwas_trait",
            "distance_bp",
            "category",
        ],
    )

    metadata = {
        "method": "gwas_overlap",
        "source_mode": source_mode,
        "project_result_file": str(project_result),
        "gwas_reference_source": str(gwas_reference_raw),
        "trait_keyword": trait_keyword,
        "window_kb": window_kb,
        "bedtools_bin": bedtools_bin,
        "project_bed": str(project_bed),
        "gwas_bed": str(gwas_bed),
        "bedtools_intersect_file": str(intersect_out),
        "result_file": str(overlap_tsv),
    }
    metadata_json = output_dir / "run_metadata.json"
    metadata_json.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return overlap_tsv


def main() -> int:
    raise SystemExit("This module is intended to be imported by main.py, not run directly.")


if __name__ == "__main__":
    sys.exit(main())
