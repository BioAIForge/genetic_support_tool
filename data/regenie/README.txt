REGENIE example files included in this directory:

1. example_dataset.ped / example_dataset.map
   - Toy PLINK PED/MAP genotype source for 8 samples and 4 variants.
   - Convert these files to PGEN before running the example config.

2. example_pheno.tsv
   - Example phenotype file with FID/IID/phenotype columns.

3. example_covar.tsv
   - Example covariate file with FID/IID/age/sex/PC1/PC2 columns.

4. example.annotations
   - Example regenie annotation file.

5. example.setlist
   - Example regenie set-list file.

6. example.masks
   - Example regenie mask definition file.

Suggested server-side preparation:

  cd /home/zj/genos-tools/genetic_support_tool

  # For newer PLINK2 versions:
  plink2 --pedmap data/regenie/example_dataset --make-pgen --out data/regenie/example_dataset

  # If your PLINK2 is older, try:
  # plink2 --file data/regenie/example_dataset --make-pgen --out data/regenie/example_dataset

This repository's regenie_example.yaml is configured to use these files and
defaults to regenie_ignore_pred=true for a lightweight demo run.
The demo also uses aaf-bins=1.0 so that the toy variants are not filtered out
by rare-variant masking.
