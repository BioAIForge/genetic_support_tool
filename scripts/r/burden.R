suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--geno", type = "character", dest = "geno"),
  make_option("--pheno", type = "character", dest = "pheno"),
  make_option("--covar", type = "character", dest = "covar"),
  make_option("--pheno-name", type = "character", dest = "pheno_name"),
  make_option("--pheno-type", type = "character", dest = "pheno_type"),
  make_option("--set-id", type = "character", dest = "set_id"),
  make_option("--engine", type = "character", default = "base", dest = "engine"),
  make_option("--covariates", type = "character", default = "", dest = "covariates"),
  make_option("--result-out", type = "character", dest = "result_out"),
  make_option("--score-out", type = "character", dest = "score_out")
)

opt <- parse_args(OptionParser(option_list = option_list))

required_args <- c(
  "geno",
  "pheno",
  "covar",
  "pheno_name",
  "pheno_type",
  "set_id",
  "result_out",
  "score_out"
)

missing_args <- required_args[vapply(required_args, function(x) {
  is.null(opt[[x]]) || length(opt[[x]]) == 0 || identical(opt[[x]], "")
}, logical(1))]

if (length(missing_args) > 0) {
  stop(sprintf(
    "Missing required command line arguments: %s",
    paste(missing_args, collapse = ", ")
  ))
}

read_table_checked <- function(path) {
  read.table(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

load_skat_if_needed <- function() {
  if (!requireNamespace("SKAT", quietly = TRUE)) {
    stop("engine=skat requires the R package 'SKAT'. Please install it first.")
  }
}

geno <- read_table_checked(opt$geno)
pheno <- read_table_checked(opt$pheno)
covar <- read_table_checked(opt$covar)

if (!"sample_id" %in% colnames(geno)) {
  stop("geno_matrix.tsv must contain sample_id")
}
if (!"sample_id" %in% colnames(pheno)) {
  stop("pheno.tsv must contain sample_id")
}
if (!"sample_id" %in% colnames(covar)) {
  stop("covar.tsv must contain sample_id")
}
if (!(opt[["pheno_name"]] %in% colnames(pheno))) {
  stop(sprintf("Phenotype column not found: %s", opt[["pheno_name"]]))
}

merged <- merge(geno, pheno, by = "sample_id", all = FALSE)
merged <- merge(merged, covar, by = "sample_id", all = FALSE)

variant_cols <- setdiff(colnames(geno), "sample_id")
if (length(variant_cols) == 0) {
  stop("No variant columns found in genotype matrix")
}

merged$burden_score <- rowSums(merged[, variant_cols, drop = FALSE], na.rm = TRUE)

covariate_fields <- character(0)
if (!is.null(opt[["covariates"]]) && nchar(opt[["covariates"]]) > 0) {
  covariate_fields <- trimws(unlist(strsplit(opt[["covariates"]], ",")))
  covariate_fields <- covariate_fields[covariate_fields != ""]
}

missing_covars <- setdiff(covariate_fields, colnames(merged))
if (length(missing_covars) > 0) {
  stop(sprintf("Missing covariates in merged table: %s", paste(missing_covars, collapse = ", ")))
}

rhs_terms <- c("burden_score", covariate_fields)
formula_text <- paste(opt[["pheno_name"]], "~", paste(rhs_terms, collapse = " + "))
model_formula <- as.formula(formula_text)

fit <- NULL
model_name <- NULL
beta_value <- NA
se_value <- NA
pvalue_value <- NA

if (opt[["engine"]] == "base") {
  if (opt[["pheno_type"]] == "continuous") {
    fit <- lm(model_formula, data = merged)
    model_name <- "linear_regression"
  } else if (opt[["pheno_type"]] == "binary") {
    fit <- glm(model_formula, data = merged, family = binomial())
    model_name <- "logistic_regression"
  } else {
    stop("pheno-type must be 'continuous' or 'binary'")
  }

  fit_summary <- summary(fit)
  coef_table <- fit_summary$coefficients

  if (!"burden_score" %in% rownames(coef_table)) {
    stop("burden_score coefficient not found in fitted model")
  }

  burden_row <- coef_table["burden_score", ]
  beta_value <- unname(burden_row[1])
  se_value <- unname(burden_row[2])
  pvalue_value <- unname(burden_row[ncol(coef_table)])
} else if (opt[["engine"]] == "skat") {
  load_skat_if_needed()

  outcome_type <- ifelse(opt[["pheno_type"]] == "continuous", "C", "D")
  null_formula_text <- if (length(covariate_fields) > 0) {
    paste(opt[["pheno_name"]], "~", paste(covariate_fields, collapse = " + "))
  } else {
    paste(opt[["pheno_name"]], "~ 1")
  }
  null_formula <- as.formula(null_formula_text)

  z_matrix <- as.matrix(merged[, variant_cols, drop = FALSE])
  storage.mode(z_matrix) <- "numeric"

  null_model <- SKAT::SKAT_Null_Model(
    formula = null_formula,
    data = merged,
    out_type = outcome_type
  )

  skat_res <- SKAT::SKAT(
    Z = z_matrix,
    obj = null_model,
    method = "Burden"
  )

  model_name <- ifelse(
    opt[["pheno_type"]] == "continuous",
    "skat_burden_continuous",
    "skat_burden_binary"
  )
  pvalue_value <- skat_res$p.value
} else {
  stop("Unsupported engine. Supported values: base, skat")
}

result_df <- data.frame(
  set_id = opt[["set_id"]],
  engine = opt[["engine"]],
  phenotype_name = opt[["pheno_name"]],
  phenotype_type = opt[["pheno_type"]],
  n_samples = nrow(merged),
  n_variants = length(variant_cols),
  burden_beta_or_logodds = beta_value,
  burden_se = se_value,
  burden_pvalue = pvalue_value,
  burden_model = model_name,
  stringsAsFactors = FALSE
)

score_df <- merged[, c("sample_id", "burden_score", opt[["pheno_name"]]), drop = FALSE]
colnames(score_df)[3] <- "phenotype"

write.table(
  result_df,
  file = opt$result_out,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  score_df,
  file = opt$score_out,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat(sprintf("Burden analysis finished. Result: %s\n", opt$result_out))
