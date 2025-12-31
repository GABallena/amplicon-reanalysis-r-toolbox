#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# PORTFOLIO-SAFE TEMPLATE
# - This script is anonymized: no client names, no site names, no sample IDs.
# - Provide your own input files (feature tables / metadata / taxonomy exports).
# - Set PROJECT_DIR to the folder containing the expected inputs.
# ------------------------------------------------------------------------------


# Fig. 5. Heatmap of dominant fungal families (>5% in any sample)
# across rearing water and gill samples.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(vegan)
  library(stringr)
  library(ggdendro)
library(patchwork)
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

# samples x ASVs
comm <- t(feat_mat)

# ---- Taxonomy -> family ----
tax <- read_tsv(tax_path, comment = "", col_types = cols(.default = col_character()))
colnames(tax)[1] <- "FeatureID"

tax <- tax %>%
  mutate(
    family = str_extract(taxonomy, "f__[^;]+"),
    family = str_remove(family, "^f__"),
    family = if_else(is.na(family) | family == "", "Unclassified", family)
  ) %>%
  select(FeatureID, family)

# align taxonomy + feature table
tax <- tax %>% filter(FeatureID %in% colnames(comm))
comm <- comm[, tax$FeatureID, drop = FALSE]
families <- tax$family

# ---- Collapse to family-level abundance ----
fam_abun <- rowsum(t(comm), group = families)  # family x samples
fam_abun <- t(fam_abun)                        # samples x family

# relative abundance per sample
fam_rel <- sweep(fam_abun, 1, rowSums(fam_abun), "/")

# ---- Metadata ----
meta <- read_csv(meta_path, show_col_types = FALSE) %>%
  mutate(    sample_id = as.character(sample_id),
Source = factor(Source, levels = c("Water", "Gill")),
    Farming_style = factor(Farming_style, levels = c("Open cage", "Closed pond"))
  )

# keep samples present in both
common <- intersect(rownames(fam_rel), meta$sample_id)
fam_rel <- fam_rel[common, , drop = FALSE]
meta    <- meta %>% filter(sample_id %in% common)

# ---- Filter dominant families (>5% in any sample) ----
keep_fams <- colnames(fam_rel)[apply(fam_rel, 2, max, na.rm = TRUE) > 0.05]
fam_rel_filt <- fam_rel[, keep_fams, drop = FALSE]

# remove unclassified/undefined families
drop_fams <- c("Unclassified", "unidentified", "unassigned", "NA")
fam_rel_filt <- fam_rel_filt[, !colnames(fam_rel_filt) %in% drop_fams, drop = FALSE]

# province code map for concise labels
# Site code → anonymized label (edit as needed)
prov_map <- setNames(
  paste("Site", LETTERS[1:5]),
  c("A1", "A2", "A3", "A4", "A5")
)

# ensure meta and matrix are in same sample order before clustering
fam_rel_filt <- fam_rel_filt[meta$sample_id, , drop = FALSE]

# ---- Cluster samples and families (Bray-Curtis) ----
if (nrow(fam_rel_filt) > 2 && ncol(fam_rel_filt) > 1) {

  # sample dendrogram
  dist_samp  <- vegdist(fam_rel_filt, method = "bray")
  samp_clust <- hclust(dist_samp, method = "average")
  sample_order <- rownames(fam_rel_filt)[samp_clust$order]

  # family dendrogram
  dist_fam  <- vegdist(t(fam_rel_filt), method = "bray")
  fam_clust <- hclust(dist_fam, method = "average")
  family_order <- colnames(fam_rel_filt)[fam_clust$order]

  # reorder matrix by clusters
  fam_rel_ord <- fam_rel_filt[sample_order, family_order, drop = FALSE]

  # order meta to match clustered samples
  meta_ord <- meta %>%
    arrange(match(sample_id, sample_order)) %>%
    mutate(
      sample_label = paste(
        prov_map[Province_code],
        as.character(Source),
        as.character(Farming_style),
        sep = " | "
      )
    )
  rownames(meta_ord) <- meta_ord$sample_id

} else {
  # not enough to cluster, just keep as-is
  samp_clust   <- NULL
  fam_clust    <- NULL
  fam_rel_ord  <- fam_rel_filt
  sample_order <- rownames(fam_rel_ord)
  family_order <- colnames(fam_rel_ord)

  meta_ord <- meta %>%
    mutate(
      sample_label = paste(
        prov_map[Province_code],
        as.character(Source),
        as.character(Farming_style),
        sep = " | "
      )
    )
  rownames(meta_ord) <- meta_ord$sample_id
}

# ---- Long format for plotting ----
plot_df <- fam_rel_ord %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "family", values_to = "rel_abund") %>%
  left_join(
    meta_ord %>%
      select(sample_id, sample_label, Source, Province, Farming_style),
    by = "sample_id"
  ) %>%
  mutate(
    sample_factor = factor(sample_id, levels = sample_order),
    family_factor = factor(family, levels = family_order)
  )

# ---- betadisper diagnostics (Bray on family rel. abundance, all families) ----
# make sure meta_for_beta is aligned with fam_rel rows
meta_for_beta <- meta_ord[match(rownames(fam_rel), meta_ord$sample_id), ]

bray_fam <- vegdist(fam_rel, method = "bray")
bd_source <- betadisper(bray_fam, meta_for_beta$Source)
bd_prov   <- betadisper(bray_fam, meta_for_beta$Province)
bd_style  <- betadisper(bray_fam, meta_for_beta$Farming_style)

cat("Betadisper (Source):\n"); print(anova(bd_source))
cat("\nBetadisper (Province):\n"); print(anova(bd_prov))
cat("\nBetadisper (Farming_style):\n"); print(anova(bd_style))

# ---- Heatmap ----
p_heat <- ggplot(plot_df,
                 aes(x = sample_factor,
                     y = family_factor,
                     fill = rel_abund)) +
  geom_tile(color = "grey60", linewidth = 0.4) +
  scale_fill_gradient2(
    low  = "#f8f8f8",
    mid  = "#9e9e9e",
    high = "#d1495b",
    midpoint = max(plot_df$rel_abund, na.rm = TRUE) / 2,
    name = "Relative abundance",
    limits = c(0, max(plot_df$rel_abund, na.rm = TRUE))
  ) +
  scale_x_discrete(
    labels = meta_ord$sample_label[match(sample_order, meta_ord$sample_id)]
  ) +
  labs(
    x = "Samples (grouped by Source, Province, Farming style)",
    y = "Fungal family",
    title = "Dominant fungal families (>5% in any sample)"
  ) +
  theme_minimal(base_size = 12, base_family = font_family) +
  theme(
    plot.background = element_rect(fill = "#fdfdfd", color = NA),
    panel.background = element_rect(fill = "#f7f9fb", color = NA),
    axis.text.x  = element_text(angle = 55, vjust = 1, hjust = 1, size = 8, family = font_family, face = "plain", color = "#1f2d3d"),
    axis.text.y  = element_text(size = 9, family = font_family, face = "plain", color = "#1f2d3d"),
    panel.grid   = element_blank(),
    legend.position = "right",
    plot.title   = element_text(family = font_family, face = "plain", size = 16, color = "#0f1c2d")
  )

final_plot <- p_heat

# ---- Add dendrograms if available ----
if (!is.null(samp_clust) && !is.null(fam_clust)) {
  final_plot <- p_heat
}

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure5_heatmap.png",
       final_plot,
       width = 10, height = 7.5, dpi = 400, type = "cairo", bg = "white")
ggsave("figures/Figure5_heatmap.tiff",
       final_plot,
       width = 10, height = 7.5, dpi = 400, type = "cairo",
       device = "tiff", compression = "lzw", bg = "white")
ggsave("figures/Figure5_heatmap.pdf",
       final_plot,
       width = 10, height = 7.5, device = cairo_pdf, bg = "white")

cat("Saved figure to figures/Figure5_heatmap.png\n")