######## Tidy Estimate Algorithm, alternative to Estimate ######
# install.packages("devtools")
devtools::install_github("KaiAragaki/tidy_estimate")
### Load the package ####
library(tidyestimate)
######### Tumor Purity Estimator Using Estimate Algorithm #########

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(tidyestimate)
})

# ---------------------------
# Inputs
# ---------------------------
tpm_file  <- "Input_file/input.csv"
meta_file <- "meatdata/metadata_site_merged.csv"

out_scores_tsv <- "output/Input_tidyestimate_scores.tsv"

# ---------------------------
# Read TPM robustly (auto-detect delimiter)
# ---------------------------
tpm_dt <- fread(tpm_file, sep = "auto", header = TRUE, data.table = FALSE, check.names = FALSE)

# Identify gene column
gene_col <- if ("gene" %in% colnames(tpm_dt)) "gene" else colnames(tpm_dt)[1]

# Drop completely empty columns (safety)
tpm_dt <- tpm_dt[, colSums(!is.na(tpm_dt)) > 0, drop = FALSE]

# Gene column explicit + deduplicate genes
tpm_dt <- tpm_dt %>%
  rename(gene = all_of(gene_col)) %>%
  distinct(gene, .keep_all = TRUE)

# Convert to matrix (genes x samples)
expr_mat <- tpm_dt %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Coerce to numeric safely
expr_mat <- suppressWarnings(apply(expr_mat, 2, function(x) as.numeric(as.character(x))))
rownames(expr_mat) <- rownames(tpm_dt %>% column_to_rownames("gene"))

# Replace NAs with 0 (TPM context)
expr_mat[is.na(expr_mat)] <- 0

# log2(TPM + 1) for tidyestimate
expr_log <- log2(expr_mat + 1)

# ---------------------------
# Run tidyestimate
# ---------------------------
scores <- expr_log |>
  filter_common_genes(id = "hgnc_symbol", tell_missing = FALSE, find_alias = TRUE) |>
  estimate_score(is_affymetrix = FALSE)

# Inspect output columns (optional but helpful)
message("tidyestimate columns: ", paste(colnames(scores), collapse = ", "))

# If purity isn't returned, compute it from ESTIMATE score (ESTIMATE mapping)
# purity = cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)
if (!("purity" %in% colnames(scores))) {
  if (!("estimate" %in% colnames(scores))) {
    stop("Neither 'purity' nor 'estimate' columns found in tidyestimate output. Check tidyestimate version/output.")
  }
  scores <- scores %>%
    mutate(purity = cos(0.6049872018 + 0.0001467884 * estimate))
}

# ---------------------------
# Merge MLIVER metadata
# ---------------------------
meta <- read.csv(meta_file, check.names = FALSE) %>%
  select(trunc_anonymized_sample_ids, MLIVER) %>%
  mutate(MLIVER = factor(MLIVER, levels = c("N","Y")))

final_scores <- scores %>%
  rename(Sample = sample) %>%
  left_join(meta, by = c("Sample" = "trunc_anonymized_sample_ids")) %>%
  relocate(Sample, MLIVER) %>%
  arrange(MLIVER, purity)

print(final_scores)

# Save
write.table(final_scores, out_scores_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nSaved scores:", out_scores_tsv, "\n")
