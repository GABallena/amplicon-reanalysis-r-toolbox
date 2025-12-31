#!/usr/bin/env Rscript

# ---------------------------
# PORTFOLIO / ANONYMIZED TEMPLATE
# ---------------------------
# Portfolio-safe version of a real workflow.
# Client-specific identifiers (paths, sample IDs, locations, coordinates) were removed or generalized.
# Data files are NOT included. Configure paths via environment variables or edit the CONFIG section.
#
# Recommended usage:
#   PROJECT_DIR=/path/to/project Rscript <this_script>.R
# ---------------------------

project_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
setwd(project_dir)

# Load required libraries
library(tidyverse)
library(vegan)

cat("=== CCA recheck: starting ===\n")
perms <- 999  # permutations for anova.cca

# ---------- 1. Read feature table ----------
feat_path <- Sys.getenv("FEATURE_TABLE_TSV", unset = "exported_feature_table/feature-table.tsv")
cat("Reading feature table from:", feat_path, "\n")

feat_raw <- read_tsv(
  feat_path,
  skip = 1,          # skip first comment line
  col_types = cols()
)

# first column = ASV IDs (strip leading "#OTU ID")
colnames(feat_raw)[1] <- "FeatureID"

# make ASVs rows, samples columns
feat_mat <- feat_raw %>%
  column_to_rownames("FeatureID") %>%
  as.matrix()

cat("Feature table dimensions (ASVs x samples):",
    paste(dim(feat_mat), collapse = " x "), "\n")

# ---------- 2. Read metadata (incl. environmental vars) ----------
meta_path <- Sys.getenv("METADATA_CSV", unset = "metadata.csv")
cat("Reading metadata from:", meta_path, "\n")

meta <- read_csv(meta_path, show_col_types = FALSE)

# Optional: if you have known sample-ID inconsistencies between exports and metadata,
# reconcile them here (left as a no-op in the portfolio version).
  )

cat("Metadata columns:\n")
print(colnames(meta))

# Make sure sample_id matches feature-table column names
# (adjust if your column name is different)
if (!"sample_id" %in% colnames(meta)) {
  stop("No 'sample_id' column in map_with_metadata.csv; please adjust script.")
}

# Align samples between feature table and metadata
common_samples <- intersect(colnames(feat_mat), meta$sample_id)
cat("Common samples:", length(common_samples), "\n")

feat_mat <- feat_mat[, common_samples, drop = FALSE]
meta_sub <- meta %>% filter(sample_id %in% common_samples)

# ---------- 3. Define environmental variables ----------
# Adjust this vector to match *your* column names for water quality.
# If none of these exist, we fall back to Source / Province_code / Farming_style.
preferred_env <- c(
  "Atm_temp",       # atmospheric temperature
  "Water_temp",     # water temperature
  "DO",             # dissolved oxygen
  "pH",
  "Salinity",
  "Conductivity",
  "TDS",
  "Nitrite",
  "Total_ammonia",
  "Iron"
)

env_vars <- intersect(preferred_env, colnames(meta_sub))

if (length(env_vars) == 0) {
  cat("WARNING: No preferred environmental variables found.\n")
  fallback_env <- intersect(c("Source", "Province_code", "Farming_style"),
                            colnames(meta_sub))
  if (length(fallback_env) == 0) {
    stop("No usable environmental variables available.")
  }
  env_vars <- fallback_env
  cat("Using fallback env vars:", paste(env_vars, collapse = ", "), "\n")
} else {
  missing_env <- setdiff(preferred_env, colnames(meta_sub))
  if (length(missing_env) > 0) {
    cat("WARNING: These env vars are missing in metadata:\n")
    print(missing_env)
    cat("Using available env vars only.\n")
  }
  cat("Environmental variables used:\n")
  print(env_vars)
}

# Filter metadata to samples with complete env data
env_all <- meta_sub %>%
  select(sample_id, all_of(env_vars)) %>%
  drop_na()

# Re-align feature table to rows in env_all
common_samples_env <- intersect(colnames(feat_mat), env_all$sample_id)
cat("Samples with complete env data:", length(common_samples_env), "\n")

# subset feature table & metadata to those samples
feat_mat_env <- feat_mat[, common_samples_env, drop = FALSE]
env_all <- env_all %>%
  filter(sample_id %in% common_samples_env) %>%
  arrange(match(sample_id, common_samples_env))

# sanity check
stopifnot(identical(env_all$sample_id, colnames(feat_mat_env)))

# ---------- 4. Prepare community matrix for CCA ----------
# transpose: samples as rows, ASVs as columns
comm_all <- t(feat_mat_env)

# relative abundance per sample
comm_rel <- decostand(comm_all, method = "total")  # each row sums to 1

# log10(x + 1) transform to down-weight dominant ASVs
comm_log <- log10(comm_rel + 1)

# environmental matrix: numeric only
env_mat <- env_all %>%
  select(all_of(env_vars))

# Standardize numeric columns; leave factors/characters untouched
env_scaled <- env_mat %>%
  mutate(across(where(is.numeric), scale)) %>%
  mutate(across(where(is.character), as.factor))

cat("Community matrix (samples x ASVs):",
    paste(dim(comm_log), collapse = " x "), "\n")
cat("Env matrix (samples x vars):",
    paste(dim(env_scaled), collapse = " x "), "\n")

# ---------- 5. Run CCA: all samples ----------
cat("\n=== CCA: all samples ===\n")

cca_all <- cca(comm_log ~ ., data = as.data.frame(env_scaled))

cat("Inertia (total):", cca_all$tot.chi, "\n")
cat("Constrained inertia:", cca_all$CCA$tot.chi, "\n")
cat("Unconstrained inertia:", cca_all$CA$tot.chi, "\n")

cat("\nGlobal test (anova.cca, perms =", perms, "):\n")
print(anova.cca(cca_all, permutations = perms))

cat("\nBy-axis test (anova.cca, by = 'axis'):\n")
print(anova.cca(cca_all, by = "axis", permutations = perms))

cat("\nProportion of variance explained by canonical axes:\n")
print(summary(cca_all)$concont$importance)

# ---------- 6. Separate CCA for water-only and gill-only ----------
if (!"Source" %in% colnames(meta_sub)) {
  cat("\nNo 'Source' column in metadata (Gill/Water); skipping split CCA.\n")
} else {
  # attach Source to the env_all table
  src_df <- meta_sub %>%
    filter(sample_id %in% common_samples_env) %>%
    arrange(match(sample_id, env_all$sample_id)) %>%
    select(sample_id, Source)

  stopifnot(identical(src_df$sample_id, env_all$sample_id))

  env_all$Source <- src_df$Source

  # helper function
  run_cca_subset <- function(source_label) {
    cat("\n=== CCA:", source_label, "only ===\n")
    idx <- which(env_all$Source == source_label)
    if (length(idx) < 5) {
      cat("Not enough samples for", source_label, "CCA (n =", length(idx), ")\n")
      return(NULL)
    }
    comm_sub <- comm_log[idx, , drop = FALSE]
    env_sub <- env_scaled[idx, , drop = FALSE]

    cca_sub <- cca(comm_sub ~ ., data = as.data.frame(env_sub))

    cat("Inertia (total):", cca_sub$tot.chi, "\n")
    cat("Constrained inertia:", cca_sub$CCA$tot.chi, "\n")
    cat("Unconstrained inertia:", cca_sub$CA$tot.chi, "\n")

    cat("\nGlobal test (anova.cca, perms =", perms, "):\n")
    print(anova.cca(cca_sub, permutations = perms))

    cat("\nBy-axis test (anova.cca, by = 'axis'):\n")
    print(anova.cca(cca_sub, by = "axis", permutations = perms))

    cat("\nProportion of variance explained by canonical axes:\n")
    print(summary(cca_sub)$concont$importance)

    invisible(cca_sub)
  }

  cca_water <- run_cca_subset("Water")
  cca_gill  <- run_cca_subset("Gill")
}

cat("\n=== CCA recheck: done ===\n")
