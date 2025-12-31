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


# Beta diversity re-check: Bray–Curtis NMDS + PERMANOVA
# Inputs:
#   - exported_bray/distance-matrix.tsv  (from QIIME2)
#   - map_with_metadata.csv              (your cleaned metadata)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(vegan)
})

dist_path <- Sys.getenv("DISTANCE_MATRIX_TSV", unset = "exported_bray/distance-matrix.tsv")
meta_path <- Sys.getenv("METADATA_CSV", unset = "metadata.csv")

cat("Reading Bray–Curtis distance matrix from:", dist_path, "\n\n")
dist_df <- read_tsv(dist_path, show_col_types = FALSE)

# QIIME distance-matrix.tsv typically:
#   col1 = sample IDs, remaining columns = distances with colnames = sample IDs
first_col <- names(dist_df)[1]
sample_ids <- dist_df[[1]]

dist_mat <- dist_df[, -1, drop = FALSE] %>% as.matrix()
rownames(dist_mat) <- sample_ids

# sanity checks
if (nrow(dist_mat) != ncol(dist_mat)) {
  stop("Distance matrix is not square: ", nrow(dist_mat), " x ", ncol(dist_mat))
}

if (!all(rownames(dist_mat) == colnames(dist_mat))) {
  warning("Row and column names differ; attempting to reorder columns to match rows.")
  common_ids <- intersect(rownames(dist_mat), colnames(dist_mat))
  dist_mat <- dist_mat[common_ids, common_ids, drop = FALSE]
}

bray_all <- as.dist(dist_mat)
cat("Distance matrix includes", attr(bray_all, "Size"), "samples\n\n")

# ---- Metadata ----
meta <- read_csv(meta_path, show_col_types = FALSE)

# Optional: if you have known sample-ID inconsistencies between exports and metadata,
# reconcile them here (left as a no-op in the portfolio version).
  )

meta_matched <- meta %>%
  filter(sample_id %in% rownames(dist_mat)) %>%
  arrange(match(sample_id, rownames(dist_mat)))

if (!all(meta_matched$sample_id == rownames(dist_mat))) {
  stop("Metadata sample_id order does not match distance matrix rownames even after match().")
}

cat("Matched metadata rows:", nrow(meta_matched), "\n\n")
cat("Source breakdown:\n")
print(table(meta_matched$Source))
cat("\n")

# ---- NMDS (all samples) ----
set.seed(9999)
nmds_all <- metaMDS(bray_all, k = 2, trymax = 999)

cat("NMDS (all samples):\n")
cat("  Stress =", nmds_all$stress, "\n\n")

# ---- PERMANOVA: all samples ----
cat("PERMANOVA (all samples: Source + Province_code + Farming_style)\n")
perm_all <- adonis2(
  bray_all ~ Source + Province_code + Farming_style,
  data = meta_matched,
  permutations = 9999,
  by = "margin"
)
print(perm_all)
cat("\n")

# ---- PERMANOVA: water only ----
meta_water <- meta_matched %>% filter(Source == "Water")
if (nrow(meta_water) > 2) {
  mat_w <- dist_mat[meta_water$sample_id, meta_water$sample_id, drop = FALSE]
  bray_water <- as.dist(mat_w)

  cat("PERMANOVA (Water only: Province_code + Farming_style)\n")
  perm_water <- adonis2(
    bray_water ~ Province_code + Farming_style,
    data = meta_water,
    permutations = 9999,
    by = "margin"
  )
  print(perm_water)
  cat("\n")
} else {
  cat("Not enough water samples for PERMANOVA.\n\n")
}

# ---- PERMANOVA: gill only ----
meta_gill <- meta_matched %>% filter(Source == "Gill")
if (nrow(meta_gill) > 2) {
  mat_g <- dist_mat[meta_gill$sample_id, meta_gill$sample_id, drop = FALSE]
  bray_gill <- as.dist(mat_g)

  cat("PERMANOVA (Gill only: Province_code + Farming_style)\n")
  perm_gill <- adonis2(
    bray_gill ~ Province_code + Farming_style,
    data = meta_gill,
    permutations = 9999,
    by = "margin"
  )
  print(perm_gill)
  cat("\n")
} else {
  cat("Not enough gill samples for PERMANOVA.\n\n")
}

cat("Done.\n")
