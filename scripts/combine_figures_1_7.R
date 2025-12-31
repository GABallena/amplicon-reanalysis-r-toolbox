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


# Combine Figures 1-7 PDFs into a single file for submission using pdfunite.

suppressPackageStartupMessages(library(fs))
pdfs <- c(
  "figures/Figure1_map.pdf",
  "figures/Figure2_rarefaction.pdf",
  "figures/Figure3_alpha_diversity.pdf",
  "figures/Figure4_NMDS.pdf",
  "figures/Figure5_heatmap.pdf",
  "figures/Figure6_top20_taxa.pdf",
  "figures/Figure7_key_taxa.pdf"
)

if (!all(file_exists(pdfs))) {
  missing <- pdfs[!file_exists(pdfs)]
  stop("Missing PDFs: ", paste(missing, collapse = ", "))
}

output <- "figures/Figures_1_7_combined_portfolio.pdf"
cmd <- c("pdfunite", pdfs, output)
res <- system2(cmd[1], cmd[-1], stdout = TRUE, stderr = TRUE)

if (!file_exists(output)) {
  stop("pdfunite failed: ", paste(res, collapse = "\n"))
}

cat("Combined PDF saved to", output, "\n")
