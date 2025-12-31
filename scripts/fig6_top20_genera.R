#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# PORTFOLIO-SAFE TEMPLATE
# - This script is anonymized: no client names, no site names, no sample IDs.
# - Provide your own input files (feature tables / metadata / taxonomy exports).
# - Set PROJECT_DIR to the folder containing the expected inputs.
# ------------------------------------------------------------------------------


# Fig. 6. Top 20 most abundant fungal genera in rearing water and gills.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(stringr)
})
# ---- Project directory ----
project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
setwd(project_dir)
pdf_fonts <- names(grDevices::pdfFonts())
font_family <- if ("Times New Roman" %in% pdf_fonts) {
  "Times New Roman"
} else if ("Times" %in% pdf_fonts) {
  "Times"
} else {
  "serif"
}

# ---- Paths ----
feat_path <- "exported_feature_table/feature-table.tsv"
tax_path  <- "taxa/taxonomy.tsv"
meta_path <- "map_with_metadata.csv"

# ---- Feature table ----
feat_raw <- read_tsv(feat_path, skip = 1, col_types = cols())
colnames(feat_raw)[1] <- "FeatureID"
feat_mat <- feat_raw %>%
  column_to_rownames("FeatureID") %>%
  as.matrix()
comm <- t(feat_mat)  # samples x ASVs

# ---- Taxonomy -> genus ----
tax <- read_tsv(tax_path, comment = "", col_types = cols(.default = col_character()))
colnames(tax)[1] <- "FeatureID"
tax <- tax %>%
  mutate(
    genus = str_extract(taxonomy, "g__[^;]+"),
    genus = str_remove(genus, "^g__"),
    genus = if_else(is.na(genus) | genus == "", "Unclassified", genus)
  ) %>%
  select(FeatureID, genus)

tax <- tax %>% filter(FeatureID %in% colnames(comm))
comm <- comm[, tax$FeatureID, drop = FALSE]

# ---- Collapse to genus ----
genus_vec <- tax$genus
genus_abun <- rowsum(t(comm), group = genus_vec)  # genus x samples
genus_abun <- t(genus_abun)  # samples x genus
genus_rel <- sweep(genus_abun, 1, rowSums(genus_abun), "/")

# ---- Metadata ----
meta <- read_csv(meta_path, show_col_types = FALSE) %>%
  mutate(    sample_id = as.character(sample_id),
Source = factor(Source, levels = c("Water", "Gill")),
    Farming_style = factor(Farming_style, levels = c("Open cage", "Closed pond"))
  )

common <- intersect(rownames(genus_rel), meta$sample_id)
genus_rel <- genus_rel[common, , drop = FALSE]
meta <- meta %>% filter(sample_id %in% common)

# ---- Top 20 genera by mean abundance within each Source ----
genus_long <- genus_rel %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "genus", values_to = "rel_abund") %>%
  left_join(meta %>% select(sample_id, Source), by = "sample_id")

top20 <- genus_long %>%
  group_by(Source, genus) %>%
  summarise(mean_abund = mean(rel_abund, na.rm = TRUE), .groups = "drop") %>%
  group_by(Source) %>%
  slice_max(mean_abund, n = 20) %>%
  ungroup()

plot_df <- genus_long %>%
  semi_join(top20, by = c("Source", "genus")) %>%
  group_by(Source, genus) %>%
  summarise(
    mean = mean(rel_abund, na.rm = TRUE),
    se   = sd(rel_abund, na.rm = TRUE) / sqrt(sum(!is.na(rel_abund))),
    .groups = "drop"
  ) %>%
  mutate(
    Source = factor(Source, levels = c("Water", "Gill")),
    genus = reorder(genus, mean)
  )

cols <- c("Water" = "#1f3c5d", "Gill" = "#d1495b")

p <- ggplot(plot_df, aes(x = genus, y = mean, fill = Source)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, color = "gray30") +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    position = position_dodge(width = 0.8),
    width = 0.3
  ) +
  coord_flip() +
  facet_wrap(~ Source, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = cols, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = NULL,
    y = "Mean relative abundance",
    title = "Top 20 most abundant fungal genera",
    subtitle = "Water and gill samples (mean ± SE; relative abundance after rarefaction)"
  ) +
  theme_minimal(base_size = 12, base_family = font_family) +
  theme(
    plot.background = element_rect(fill = "#fdfdfd", color = NA),
    panel.background = element_rect(fill = "#f7f9fb", color = NA),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    strip.text = element_text(family = font_family, face = "plain", size = 12),
    plot.title = element_text(family = font_family, face = "plain", size = 16, color = "#0f1c2d"),
    plot.subtitle = element_text(family = font_family, face = "plain", size = 12, color = "#2c3e50"),
    axis.text = element_text(color = "#1f2d3d")
  )

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure6_top20_genera.png", p, width = 7, height = 8, dpi = 400, type = "cairo", bg = "white")
ggsave("figures/Figure6_top20_genera.tiff", p, width = 7, height = 8, dpi = 400,
       type = "cairo", device = "tiff", compression = "lzw", bg = "white")
ggsave("figures/Figure6_top20_genera.pdf", p, width = 7, height = 8, device = cairo_pdf, bg = "white")

cat("Saved figure to figures/Figure6_top20_genera.png\n")