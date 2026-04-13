from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def dosage_to_alleles(value: str) -> tuple[str, str]:
    normalized = value.strip()
    if normalized in {"", "NA", ".", "nan"}:
        return ("0", "0")
    dosage = int(float(normalized))
    if dosage <= 0:
        return ("A", "A")
    if dosage == 1:
        return ("A", "G")
    return ("G", "G")


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_example(project_root: Path) -> None:
    source_dir = project_root / "data" / "example_strong_signal"
    target_dir = project_root / "data" / "plink2"
    target_dir.mkdir(parents=True, exist_ok=True)

    geno_rows = read_tsv(source_dir / "geno_matrix.tsv")
    pheno_rows = read_tsv(source_dir / "pheno.tsv")
    covar_rows = read_tsv(source_dir / "covar.tsv")

    if not geno_rows:
        raise ValueError("geno_matrix.tsv is empty")

    variant_ids = [key for key in geno_rows[0].keys() if key != "sample_id"]
    ped_path = target_dir / "example_dataset.ped"
    map_path = target_dir / "example_dataset.map"
    pheno_path = target_dir / "example_pheno.tsv"
    covar_path = target_dir / "example_covar.tsv"

    with map_path.open("w", encoding="utf-8", newline="") as handle:
        for index, variant_id in enumerate(variant_ids, start=1):
            handle.write(f"1\t{variant_id}\t0\t{index * 1000}\n")

    with ped_path.open("w", encoding="utf-8", newline="") as handle:
        for row in geno_rows:
            sample_id = row["sample_id"]
            alleles: list[str] = []
            for variant_id in variant_ids:
                a1, a2 = dosage_to_alleles(row[variant_id])
                alleles.extend([a1, a2])
            tokens = ["0", sample_id, "0", "0", "0", "-9", *alleles]
            handle.write(" ".join(tokens) + "\n")

    pheno_out = [{"FID": "0", "IID": row["sample_id"], "phenotype": row["phenotype"]} for row in pheno_rows]
    write_tsv(pheno_path, pheno_out, ["FID", "IID", "phenotype"])

    covar_out = []
    for row in covar_rows:
        covar_out.append(
            {
                "FID": "0",
                "IID": row["sample_id"],
                "age": row["age"],
                "sex": row["sex"],
                "PC1": row["PC1"],
                "PC2": row["PC2"],
            }
        )
    write_tsv(covar_path, covar_out, ["FID", "IID", "age", "sex", "PC1", "PC2"])

    print("Generated PLINK2 example inputs:")
    print(ped_path)
    print(map_path)
    print(pheno_path)
    print(covar_path)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate PLINK2 example inputs from bundled TSV demo data.")
    parser.parse_args()
    project_root = Path(__file__).resolve().parents[2]
    build_example(project_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
