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


# Supplementary Fig. S1. Distribution of primer detection across forward and reverse reads.
# Primer sequences were variably captured within sequenced read windows, consistent with ITS
# length heterogeneity rather than amplification failure.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})
pdf_fonts <- names(grDevices::pdfFonts())
font_family <- if ("Times New Roman" %in% pdf_fonts) {
  "Times New Roman"
} else if ("Times" %in% pdf_fonts) {
  "Times"
} else {
  "serif"
}

summary_path <- Sys.getenv("CUTADAPT_SUMMARY_TSV", unset = "logs/cutadapt_summary.tsv")
meta_path    <- Sys.getenv("METADATA_CSV", unset = "metadata.csv")

# ---- Load cutadapt summary + metadata ----
cutadapt_df <- read_tsv(summary_path, show_col_types = FALSE)

meta <- read_csv(meta_path, show_col_types = FALSE) %>%
  mutate(Source = factor(Source, levels = c("Water", "Gill")))

primer_df <- cutadapt_df %>%
  mutate(
    sample_id = sub("_S\\d+_L\\d+$", "", sample),
    pct_fwd   = with_fwd / total_pairs * 100,
    pct_rev   = with_rev / total_pairs * 100
  ) %>%
  left_join(meta %>% select(sample_id, Source), by = "sample_id")

missing_meta <- primer_df %>%
  filter(is.na(Source)) %>%
  pull(sample_id) %>%
  unique()
if (length(missing_meta) > 0) {
  warning("No metadata match for: ", paste(missing_meta, collapse = ", "))
}

# ---- Long format for distributions ----
long_det <- primer_df %>%
  select(sample_id, Source, pct_fwd, pct_rev) %>%
  pivot_longer(
    cols = c(pct_fwd, pct_rev),
    names_to = "read",
    values_to = "pct_detected"
  ) %>%
  mutate(
    read = recode(
      read,
      pct_fwd = "Forward (R1)",
      pct_rev = "Reverse (R2)"
    ),
    read = factor(read, levels = c("Forward (R1)", "Reverse (R2)"))
  )

cols <- c("Water" = "#1f3c5d", "Gill" = "#d1495b")

# ---- Panel A: distributions by read ----
p_box <- ggplot(long_det, aes(read, pct_detected, fill = Source)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.55, colour = "black") +
  geom_jitter(aes(colour = Source), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = cols, drop = FALSE) +
  scale_colour_manual(values = cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Primer-detected reads (% of pairs)",
    title = "A. Primer detection by read"
  ) +
  theme_minimal(base_size = 12, base_family = font_family) +
  theme(
    plot.background = element_rect(fill = "#fdfdfd", color = NA),
    panel.background = element_rect(fill = "#f7f9fb", color = NA),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# ---- Panel B: concordance between forward and reverse ----
p_scatter <- ggplot(primer_df, aes(pct_fwd, pct_rev, colour = Source)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 2.6, alpha = 0.9) +
  scale_colour_manual(values = cols, name = "Sampling source", drop = FALSE) +
  labs(
    x = "Forward read with primer (%)",
    y = "Reverse read with primer (%)",
    title = "B. Concordance between forward and reverse detection"
  ) +
  coord_equal() +
  theme_minimal(base_size = 12, base_family = font_family) +
  theme(
    plot.background = element_rect(fill = "#fdfdfd", color = NA),
    panel.background = element_rect(fill = "#f7f9fb", color = NA),
    panel.grid.minor = element_blank()
  )

combined <- p_box + p_scatter + plot_layout(widths = c(0.9, 1.1))

dir.create("figures", showWarnings = FALSE)
ggsave(
  filename = "figures/FigureS1_primer_detection.png",
  plot = combined,
  width = 9,
  height = 4.5,
  dpi = 400,
  type = "cairo",
  bg = "white"
)

cat("Saved supplementary figure to figures/FigureS1_primer_detection.png\n")
cat(sprintf(
  "Median primer detection: R1 = %.1f%%, R2 = %.1f%%\n",
  median(primer_df$pct_fwd, na.rm = TRUE),
  median(primer_df$pct_rev, na.rm = TRUE)
))
