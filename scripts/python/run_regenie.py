from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path
from typing import Any


def ensure_exists(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def _extract_pvalue(row: dict[str, str]) -> float | None:
    for key in ("P", "PVAL", "P_VALUE", "PVALUE"):
        value = row.get(key)
        if value not in (None, "", "NA"):
            return float(value)

    log10p = row.get("LOG10P")
    if log10p not in (None, "", "NA"):
        return 10 ** (-float(log10p))

    return None


def _read_regenie_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        filtered = [
            line for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        ]
    if not filtered:
        return []
    header = filtered[0].strip().split()
    rows: list[dict[str, str]] = []
    for line in filtered[1:]:
        parts = line.strip().split()
        if len(parts) < len(header):
            continue
        rows.append(dict(zip(header, parts)))
    return rows


def _find_result_file(output_prefix: Path, phenotype_name: str) -> Path:
    candidate = output_prefix.parent / f"{output_prefix.name}_{phenotype_name}.regenie"
    if candidate.exists():
        return candidate

    matches = sorted(output_prefix.parent.glob(f"{output_prefix.name}_*.regenie"))
    if not matches:
        raise FileNotFoundError(
            f"No regenie result file found for prefix {output_prefix} and phenotype {phenotype_name}."
        )
    if len(matches) == 1:
        return matches[0]
    raise FileNotFoundError(
        "Multiple regenie result files were found, but none matched the requested phenotype "
        f"'{phenotype_name}': {', '.join(str(path.name) for path in matches)}"
    )


def _match_set_rows(rows: list[dict[str, str]], set_id: str) -> list[dict[str, str]]:
    matched: list[dict[str, str]] = []
    for row in rows:
        variant_id = row.get("ID", "")
        if variant_id == set_id or variant_id.startswith(f"{set_id}."):
            matched.append(row)
    if matched:
        return matched
    available_ids = sorted({row.get("ID", "") for row in rows if row.get("ID")})
    preview = ", ".join(available_ids[:10])
    if len(available_ids) > 10:
        preview += ", ..."
    raise ValueError(
        f"regenie result file does not contain set_id '{set_id}'. "
        f"Available IDs: {preview or 'none'}"
    )


def _pick_best_row(
    rows: list[dict[str, str]],
    *,
    test_names: list[str] | None,
    require_test_match: bool = False,
) -> dict[str, str]:
    candidates = rows
    if test_names:
        filtered = [row for row in rows if row.get("TEST", "").lower() in test_names]
        if filtered:
            candidates = filtered
        elif require_test_match:
            available_tests = sorted({row.get("TEST", "") for row in rows if row.get("TEST")})
            raise ValueError(
                "regenie result file does not contain any expected TEST values. "
                f"Expected one of: {', '.join(test_names)}. "
                f"Available TEST values: {', '.join(available_tests) or 'none'}"
            )

    best_row: dict[str, str] | None = None
    best_p: float | None = None
    for row in candidates:
        p = _extract_pvalue(row)
        if p is None:
            continue
        if best_p is None or p < best_p:
            best_p = p
            best_row = row

    if best_row is not None:
        return best_row
    if candidates:
        return candidates[0]
    raise ValueError("No regenie rows available after filtering.")


def run_regenie_set_tests(
    *,
    regenie_bin: str,
    genotype_prefix: Path,
    genotype_format: str,
    phenotype_file: Path,
    phenotype_name: str,
    phenotype_type: str,
    covariate_file: Path | None,
    covariate_names: list[str],
    anno_file: Path,
    set_list: Path,
    mask_def: Path,
    output_prefix: Path,
    vc_tests: list[str] | None,
    aaf_bins: str,
    build_mask: str,
    ignore_pred: bool,
    pred_file: Path | None,
    bsize: int,
    extra_args: list[str] | None,
    working_directory: Path,
) -> Path:
    genotype_flag = {"pfile": "pgen", "bfile": "bed"}.get(genotype_format, genotype_format)
    command = [
        regenie_bin,
        "--step",
        "2",
        f"--{genotype_flag}",
        str(genotype_prefix),
        "--phenoFile",
        str(phenotype_file),
        "--phenoColList",
        phenotype_name,
        "--anno-file",
        str(anno_file),
        "--set-list",
        str(set_list),
        "--mask-def",
        str(mask_def),
        "--aaf-bins",
        aaf_bins,
        "--build-mask",
        build_mask,
        "--bsize",
        str(bsize),
        "--out",
        str(output_prefix),
    ]

    if phenotype_type == "binary":
        command.append("--bt")

    if covariate_file is not None:
        command.extend(["--covarFile", str(covariate_file)])
        if covariate_names:
            command.extend(["--covarColList", ",".join(covariate_names)])

    if ignore_pred:
        command.append("--ignore-pred")
    elif pred_file is not None:
        command.extend(["--pred", str(pred_file)])
    else:
        raise ValueError(
            "regenie step 2 requires either input.pred_file or regenie_ignore_pred=true."
        )

    if vc_tests:
        command.extend(["--vc-tests", ",".join(vc_tests)])

    if extra_args:
        command.extend(extra_args)

    completed = subprocess.run(
        command,
        cwd=working_directory,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "regenie backend failed.\n"
            f"Command: {' '.join(command)}\n"
            f"Exit code: {completed.returncode}\n"
            f"STDOUT:\n{completed.stdout}\n"
            f"STDERR:\n{completed.stderr}"
        )

    return _find_result_file(output_prefix, phenotype_name)


def parse_regenie_to_burden_result(
    *,
    result_file: Path,
    set_id: str,
    phenotype_name: str,
    phenotype_type: str,
    engine: str,
    output_file: Path,
) -> None:
    rows = _read_regenie_rows(result_file)
    if not rows:
        raise ValueError(f"regenie result file has no data rows: {result_file}")
    matched = _match_set_rows(rows, set_id)
    row = _pick_best_row(
        matched,
        test_names=["add", "joint", "minp", "acat", "acato", "acatv", "sbat"],
        require_test_match=False,
    )
    p = _extract_pvalue(row)

    result = [{
        "set_id": set_id,
        "engine": engine,
        "phenotype_name": phenotype_name,
        "phenotype_type": phenotype_type,
        "n_samples": row.get("N", "NA"),
        "n_variants": row.get("NBURDEN", row.get("NVAR", "NA")),
        "burden_beta_or_logodds": row.get("BETA", "NA"),
        "burden_se": row.get("SE", "NA"),
        "burden_pvalue": "NA" if p is None else f"{p:.12g}",
        "burden_model": f"regenie_{row.get('TEST', 'burden').lower()}",
    }]
    with output_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "set_id",
                "engine",
                "phenotype_name",
                "phenotype_type",
                "n_samples",
                "n_variants",
                "burden_beta_or_logodds",
                "burden_se",
                "burden_pvalue",
                "burden_model",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(result)


def parse_regenie_to_skato_result(
    *,
    result_file: Path,
    set_id: str,
    phenotype_name: str,
    phenotype_type: str,
    engine: str,
    output_file: Path,
) -> None:
    rows = _read_regenie_rows(result_file)
    if not rows:
        raise ValueError(f"regenie result file has no data rows: {result_file}")
    matched = _match_set_rows(rows, set_id)
    row = _pick_best_row(
        matched,
        test_names=["skato", "add-skato", "skato-acat", "add-skat"],
        require_test_match=True,
    )
    p = _extract_pvalue(row)

    result = [{
        "set_id": set_id,
        "engine": engine,
        "phenotype_name": phenotype_name,
        "phenotype_type": phenotype_type,
        "n_samples": row.get("N", "NA"),
        "n_variants": row.get("NBURDEN", row.get("NVAR", "NA")),
        "skato_pvalue": "NA" if p is None else f"{p:.12g}",
        "skato_model": f"regenie_{row.get('TEST', 'skato').lower()}",
    }]
    with output_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "set_id",
                "engine",
                "phenotype_name",
                "phenotype_type",
                "n_samples",
                "n_variants",
                "skato_pvalue",
                "skato_model",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(result)


def build_regenie_metadata(
    *,
    method: str,
    engine: str,
    phenotype_name: str,
    phenotype_type: str,
    genotype_format: str,
    genotype_prefix: Path,
    pred_file: Path | None,
    phenotype_file: Path,
    covariate_file: Path | None,
    anno_file: Path,
    set_list: Path,
    mask_def: Path,
    result_file: Path,
) -> dict[str, Any]:
    return {
        "method": method,
        "engine": engine,
        "phenotype_name": phenotype_name,
        "phenotype_type": phenotype_type,
        "genotype_format": genotype_format,
        "genotype_prefix": str(genotype_prefix),
        "pred_file": None if pred_file is None else str(pred_file),
        "phenotype_file": str(phenotype_file),
        "covariate_file": None if covariate_file is None else str(covariate_file),
        "anno_file": str(anno_file),
        "set_list": str(set_list),
        "mask_def": str(mask_def),
        "regenie_result_file": str(result_file),
    }


def write_metadata(path: Path, metadata: dict[str, Any]) -> None:
    path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
