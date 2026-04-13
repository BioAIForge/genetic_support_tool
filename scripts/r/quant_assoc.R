suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--geno", type = "character", dest = "geno"),
  make_option("--pheno", type = "character", dest = "pheno"),
  make_option("--covar", type = "character", dest = "covar"),
  make_option("--pheno-name", type = "character", dest = "pheno_name"),
  make_option("--covariates", type = "character", dest = "covariates", default = ""),
  make_option("--result-out", type = "character", dest = "result_out")
)

opt <- parse_args(OptionParser(option_list = option_list))

read_table_checked <- function(path) {
  read.table(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

geno <- read_table_checked(opt$geno)
pheno <- read_table_checked(opt$pheno)
covar <- read_table_checked(opt$covar)

if (!"sample_id" %in% colnames(geno)) stop("geno_matrix.tsv must contain sample_id")
if (!"sample_id" %in% colnames(pheno)) stop("pheno.tsv must contain sample_id")
if (!"sample_id" %in% colnames(covar)) stop("covar.tsv must contain sample_id")
if (!(opt[["pheno_name"]] %in% colnames(pheno))) stop(sprintf("Phenotype column not found: %s", opt[["pheno_name"]]))

merged <- merge(geno, pheno, by = "sample_id", all = FALSE)
merged <- merge(merged, covar, by = "sample_id", all = FALSE)

variant_cols <- setdiff(colnames(geno), "sample_id")
covariate_fields <- character(0)
if (!is.null(opt[["covariates"]]) && nchar(opt[["covariates"]]) > 0) {
  covariate_fields <- trimws(unlist(strsplit(opt[["covariates"]], ",")))
  covariate_fields <- covariate_fields[covariate_fields != ""]
}

result_rows <- lapply(variant_cols, function(variant_col) {
  rhs_terms <- c(variant_col, covariate_fields)
  formula_text <- paste(opt[["pheno_name"]], "~", paste(rhs_terms, collapse = " + "))
  fit <- lm(as.formula(formula_text), data = merged)
  coef_table <- summary(fit)$coefficients
  data.frame(
    variant_id = variant_col,
    n_samples = nrow(merged),
    beta = unname(coef_table[variant_col, 1]),
    se = unname(coef_table[variant_col, 2]),
    p_value = unname(coef_table[variant_col, ncol(coef_table)]),
    stringsAsFactors = FALSE
  )
})

result_df <- do.call(rbind, result_rows)
write.table(result_df, file = opt$result_out, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("Quantitative association analysis finished. Result: %s\n", opt$result_out))
