suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--geno", type = "character", dest = "geno"),
  make_option("--pheno", type = "character", dest = "pheno"),
  make_option("--covar", type = "character", dest = "covar", default = NULL),
  make_option("--pheno-name", type = "character", dest = "pheno_name"),
  make_option("--pheno-type", type = "character", dest = "pheno_type"),
  make_option("--set-id", type = "character", dest = "set_id"),
  make_option("--covariates", type = "character", dest = "covariates", default = ""),
  make_option("--result-out", type = "character", dest = "result_out"),
  make_option("--freq-out", type = "character", dest = "freq_out")
)

opt <- parse_args(OptionParser(option_list = option_list))

read_table_checked <- function(path) {
  read.table(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

if (!requireNamespace("haplo.stats", quietly = TRUE)) {
  stop("Haplotype analysis requires the R package 'haplo.stats'. Please install it first.")
}

required_args <- c("geno", "pheno", "pheno_name", "pheno_type", "set_id", "result_out", "freq_out")
missing_args <- required_args[vapply(required_args, function(x) {
  is.null(opt[[x]]) || length(opt[[x]]) == 0 || identical(opt[[x]], "")
}, logical(1))]
if (length(missing_args) > 0) {
  stop(sprintf("Missing required command line arguments: %s", paste(missing_args, collapse = ", ")))
}

geno <- read_table_checked(opt$geno)
pheno <- read_table_checked(opt$pheno)
covar <- NULL
if (!is.null(opt$covar) && file.exists(opt$covar)) {
  covar <- read_table_checked(opt$covar)
}

if (!"sample_id" %in% colnames(geno)) stop("geno_matrix.tsv must contain sample_id")
if (!"sample_id" %in% colnames(pheno)) stop("pheno.tsv must contain sample_id")
if (!(opt[["pheno_name"]] %in% colnames(pheno))) stop(sprintf("Phenotype column not found: %s", opt[["pheno_name"]]))

merged <- merge(geno, pheno, by = "sample_id", all = FALSE)
if (!is.null(covar)) {
  if (!"sample_id" %in% colnames(covar)) stop("covar.tsv must contain sample_id")
  merged <- merge(merged, covar, by = "sample_id", all = FALSE)
}

variant_cols <- setdiff(colnames(geno), "sample_id")
if (length(variant_cols) == 0) stop("No variant columns found in genotype matrix")

convert_dosage_to_alleles <- function(dosage_vector) {
  a1 <- ifelse(dosage_vector >= 1, 2, 1)
  a2 <- ifelse(dosage_vector == 2, 2, 1)
  cbind(a1, a2)
}

allele_blocks <- lapply(variant_cols, function(variant_col) {
  dosage <- as.numeric(merged[[variant_col]])
  convert_dosage_to_alleles(dosage)
})
allele_matrix <- do.call(cbind, allele_blocks)
colnames(allele_matrix) <- as.vector(rbind(paste0(variant_cols, "_a1"), paste0(variant_cols, "_a2")))
locus_labels <- variant_cols

geno_haplo <- haplo.stats::setupGeno(allele_matrix, miss.val = c(0, NA), locus.label = locus_labels)

covariate_fields <- character(0)
if (!is.null(opt[["covariates"]]) && nchar(opt[["covariates"]]) > 0) {
  covariate_fields <- trimws(unlist(strsplit(opt[["covariates"]], ",")))
  covariate_fields <- covariate_fields[covariate_fields != ""]
}

em_obj <- haplo.stats::haplo.em(geno = allele_matrix, locus.label = locus_labels)
hap_strings <- apply(em_obj$haplotype, 1, function(x) paste(x, collapse = "-"))
freq_table <- data.frame(
  set_id = opt[["set_id"]],
  haplotype = hap_strings,
  frequency = em_obj$hap.prob,
  stringsAsFactors = FALSE
)
freq_table <- freq_table[order(freq_table$frequency, decreasing = TRUE), ]
freq_table$count_estimated <- freq_table$frequency * nrow(merged) * 2

top_haplotype <- freq_table$haplotype[1]
top_haplotype_frequency <- freq_table$frequency[1]

x_adj <- NULL
if (length(covariate_fields) > 0) {
  x_adj <- as.matrix(merged[, covariate_fields, drop = FALSE])
}

trait_type <- ifelse(opt[["pheno_type"]] == "continuous", "gaussian", "binomial")
if (!(opt[["pheno_type"]] %in% c("continuous", "binary"))) {
  stop("pheno-type must be 'continuous' or 'binary'")
}

score_obj <- haplo.stats::haplo.score(
  y = merged[[opt[["pheno_name"]]]],
  geno = allele_matrix,
  trait.type = trait_type,
  x.adj = x_adj,
  locus.label = locus_labels,
  simulate = FALSE
)

global_pvalue <- score_obj$score.global.p
model_name <- ifelse(
  opt[["pheno_type"]] == "continuous",
  "haplo_stats_score_gaussian",
  "haplo_stats_score_binomial"
)

result_df <- data.frame(
  set_id = opt[["set_id"]],
  top_haplotype = top_haplotype,
  top_haplotype_frequency = top_haplotype_frequency,
  n_samples = nrow(merged),
  n_variants = length(variant_cols),
  haplotype_effect = NA,
  haplotype_pvalue = global_pvalue,
  haplotype_model = model_name,
  stringsAsFactors = FALSE
)

write.table(result_df, file = opt$result_out, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(freq_table, file = opt$freq_out, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("Haplotype analysis finished. Result: %s\n", opt$result_out))
