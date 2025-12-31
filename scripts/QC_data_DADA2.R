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

# --- Libraries ---
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# --- Paths (adjust if needed) ---
stats_file <- Sys.getenv("DADA2_STATS_TSV", unset = "exported_stats/stats.tsv")   # from qiime tools export
metadata_file <- Sys.getenv("MAP_CSV", unset = "map.csv")                # optional; only if you want extra columns

# --- 1. Read DADA2 denoising stats ---
stats_raw <- read_tsv(
  stats_file,
  comment = "#",          # skips the "#q2:types" row
  show_col_types = FALSE
)

# stats_raw columns should look like:
# sample-id, input, filtered, `percentage of input passed filter`,
# denoised, merged, `percentage of input merged`,
# `non-chimeric`, `percentage of input non-chimeric`

# --- 2. Basic parsing from sample IDs (Group, Farming system) ---
stats_parsed <- stats_raw %>%
  rename(
    sample_id = `sample-id`,
    input_reads = input,
    nonchimeric_reads = `non-chimeric`,
    pct_nonchimeric = `percentage of input non-chimeric`
  ) %>%
  mutate(
    # Optional: infer grouping labels from sample IDs (customize patterns below)
    groupA_prefix = Sys.getenv("GROUPA_PREFIX", unset = "A-"),
    groupA_label  = Sys.getenv("GROUPA_LABEL",  unset = "Group A"),
    groupB_label  = Sys.getenv("GROUPB_LABEL",  unset = "Group B"),
    Group = if_else(str_starts(sample_id, groupA_prefix), groupA_label, groupB_label),
    System = case_when(
      str_detect(sample_id, "-OP") ~ "Open cage",
      str_detect(sample_id, "-CP") ~ "Closed pond",
      TRUE ~ NA_character_
    ),
    pct_retained = (nonchimeric_reads / input_reads) * 100
  )

# --- 3. (Optional) Join with metadata if you want site/river/etc. ---
#   Adjust join key based on the column name in map.csv
#   e.g., if map.csv has a column "sample-id" or "sampleid", change below accordingly.

if (file.exists(metadata_file)) {
  meta_raw <- read_csv(metadata_file, show_col_types = FALSE)
  
  # Try to guess a reasonable join column:
  join_col_stats <- "sample_id"
  join_col_meta <- dplyr::case_when(
    "sample-id" %in% names(meta_raw) ~ "sample-id",
    "sampleid"  %in% names(meta_raw) ~ "sampleid",
    "SampleID"  %in% names(meta_raw) ~ "SampleID",
    TRUE ~ NA_character_
  )
  
  if (!is.na(join_col_meta)) {
    stats_parsed <- stats_parsed %>%
      left_join(meta_raw, by = setNames(join_col_meta, join_col_stats))
  } else {
    message("No obvious sample ID column in metadata; skipping join.")
  }
}

# --- 4. Select columns for Supplementary Table S1 ---

supp_S1 <- stats_parsed %>%
  transmute(
    `Sample ID` = sample_id,
    Group,
    `System` = System,
    `Input reads` = input_reads,
    `Non-chimeric reads` = nonchimeric_reads,
    `Reads retained (%)` = round(pct_retained, 2)
  ) %>%
  arrange(`Sample ID`)

# --- 5. Write to CSV (or TSV) ---

write_csv(supp_S1, "Supplementary_Table_S1_sequencing_QC.csv")

# --- 6. (Optional) Print quick summary for plugging into Results ---

summary_stats <- stats_parsed %>%
  summarise(
    n_samples = n(),
    median_input = median(input_reads),
    min_input    = min(input_reads),
    max_input    = max(input_reads),
    median_nonchim = median(nonchimeric_reads),
    min_nonchim    = min(nonchimeric_reads),
    max_nonchim    = max(nonchimeric_reads),
    median_pct_nonchim = median(pct_nonchimeric),
    min_pct_nonchim    = min(pct_nonchimeric),
    max_pct_nonchim    = max(pct_nonchimeric)
  )

print(summary_stats)
