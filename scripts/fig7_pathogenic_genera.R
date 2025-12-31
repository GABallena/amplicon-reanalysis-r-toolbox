#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# PORTFOLIO-SAFE TEMPLATE
# - This script is anonymized: no client names, no site names, no sample IDs.
# - Provide your own input files (feature tables / metadata / taxonomy exports).
# - Set PROJECT_DIR to the folder containing the expected inputs.
# ------------------------------------------------------------------------------


# Fig. 7. Relative abundance of pathogenic genera (Aspergillus, Fusarium, Candida, Rhodotorula)
# across sampling sources and farming styles.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(stringr)
  library(RColorBrewer)
  library(ggalluvial)
  library(ggridges)
library(viridis)
library(circlize)
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

target_genera <- c("Aspergillus", "Fusarium", "Candida", "Rhodotorula")

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
    Farming_style = factor(Farming_style, levels = c("Open cage", "Closed pond")),
    Site = factor(Province_code, levels = c("A1", "A2", "A3", "A4", "A5"),
                  labels = c("Site A", "Site B", "Site C", "Site D", "Site E"))
  )

common <- intersect(rownames(genus_rel), meta$sample_id)
genus_rel <- genus_rel[common, , drop = FALSE]
meta <- meta %>% filter(sample_id %in% common)

# ---- Compute combined pathogenic abundance per sample ----
patho_cols <- intersect(colnames(genus_rel), target_genera)
patho_sum <- rowSums(genus_rel[, patho_cols, drop = FALSE], na.rm = TRUE)
patho_df <- tibble(sample_id = rownames(genus_rel), patho_rel_abund = patho_sum) %>%
  left_join(meta, by = "sample_id")

# ---- Plot ----
cols <- c("Water" = "#1f3c5d", "Gill" = "#d1495b")
shapes <- c("Open cage" = 21, "Closed pond" = 24)

p <- ggplot(patho_df, aes(x = Source, y = patho_rel_abund, fill = Source, shape = Farming_style)) +
  geom_boxplot(alpha = 0.4, colour = "gray30", outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6),
              size = 2.2, alpha = 0.8) +
  facet_grid(Farming_style ~ Site, switch = "x") +
  scale_fill_manual(values = cols, name = "Sampling source") +
  scale_shape_manual(values = shapes, name = "Farming style") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Site",
    y = "Combined relative abundance",
    title = "Pathogenic fungal genera across sources and farming styles"
  ) +
  theme_minimal(base_size = 12, base_family = font_family) +
  theme(
    plot.background = element_rect(fill = "#fdfdfd", color = NA),
    panel.background = element_rect(fill = "#f7f9fb", color = NA),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(family = font_family, face = "plain", size = 12),
    strip.placement = "outside",
    axis.text.x = element_text(angle = 36, hjust = 1, family = font_family, face = "plain", color = "#1f2d3d"),
    plot.title = element_text(family = font_family, face = "plain", size = 16, color = "#0f1c2d"),
    axis.title = element_text(color = "#1f2d3d", size = 12)
  )

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure7_pathogenic_genera.png", p, width = 8, height = 5, dpi = 400, type = "cairo", bg = "white")
ggsave("figures/Figure7_pathogenic_genera.tiff", p, width = 8, height = 5, dpi = 400,
       type = "cairo", device = "tiff", compression = "lzw", bg = "white")
ggsave("figures/Figure7_pathogenic_genera.pdf", p, width = 8, height = 5, device = cairo_pdf, bg = "white")

cat("Saved figure to figures/Figure7_pathogenic_genera.png\n")

# ==============================================================================
# Figure 8. Circos plot linking sites to pathogenic genera
# ==============================================================================

# Aggregate relative abundance per site, source, and genus
site_genus_base <- genus_rel %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "genus", values_to = "rel_abund") %>%
  left_join(meta, by = "sample_id") %>%
  mutate(
    genus = case_when(
      str_detect(tolower(genus), "unclassified") ~ "Unclassified",
      str_detect(tolower(genus), "unidentified|unassigned") ~ "Unidentified",
      TRUE ~ genus
    )
  ) %>%
  group_by(Site, Source, genus) %>%
  summarise(weight = sum(rel_abund, na.rm = TRUE), .groups = "drop")

# Top 12 genera across all sites, force-in grey groups
top12_genus <- site_genus_base %>%
  group_by(genus) %>%
  summarise(total = sum(weight), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 12) %>%
  pull(genus)
forced_grey <- intersect(c("Unclassified", "Unidentified"), unique(site_genus_base$genus))
final_genera <- unique(c(top12_genus, forced_grey))

# Separate sectors by Site and Source (Water = W, Gill = G)
site_genus_long <- site_genus_base %>%
  mutate(
    genus_group = if_else(genus %in% final_genera, genus, "Others"),
    source_short = if_else(Source == "Water", "W", "G"),
    site_sector = paste0(Site, " (", source_short, ")")
  ) %>%
  group_by(site_sector, genus_group) %>%
  summarise(weight = sum(weight), .groups = "drop") %>%
  filter(weight > 0)

# Cityscape-inspired palette for sites (matched to map styling)
# Palette by anonymized site label (edit as needed)
site_palette <- c(
  "Site A" = "#0f1c2d",
  "Site B" = "#1f3c5d",
  "Site C" = "#2f5f91",
  "Site D" = "#4f81bd",
  "Site E" = "#7aa6d8"
)

cityscape_pal <- colorRampPalette(site_palette)
# a slightly warmer extension to avoid “all blue” when many colors are needed
cityscape_pal_extended <- function(n) {
  colorRampPalette(c(site_palette, "#d1495b", "#edae49"))(n)
}

# Genus colors + grey for Others
base_genus_cols <- c(
  "Aspergillus" = "#ef6c00",
  "Fusarium"    = "#d32f2f",
  "Candida"     = "#8e24aa",
  "Rhodotorula" = "#00897b"
)
base_genus_cols <- base_genus_cols[names(base_genus_cols) %in% final_genera]
other_needed <- setdiff(final_genera, names(base_genus_cols))
extra_cols <- colorRampPalette(brewer.pal(12, "Set3"))(max(length(other_needed), 0))
genus_palette <- c(
  base_genus_cols,
  setNames(extra_cols[seq_along(other_needed)], other_needed),
  "Unclassified" = "#a0a0a0",
  "Unidentified" = "#c0c0c0",
  "Others" = "#b0b0b0"
)
# ensure unique colors per name (later entries overwrite earlier duplicates)
genus_palette <- genus_palette[!duplicated(names(genus_palette), fromLast = TRUE)]

# Map site_sector back to base site + source for coloring
site_sector_info <- site_genus_base %>%
  distinct(Site, Source) %>%
  mutate(
    source_short = if_else(Source == "Water", "W", "G"),
    site_sector = paste0(Site, " (", source_short, ")")
  )

site_sector_cols <- sapply(seq_len(nrow(site_sector_info)), function(i) {
  s <- site_sector_info$Site[i]
  src <- site_sector_info$Source[i]
  base_col <- site_palette[[as.character(s)]]
  if (src == "Water") {
    adjustcolor(base_col, offset = c(0.12, 0.12, 0.12, 0), alpha.f = 0.9)
  } else {
    adjustcolor(base_col, alpha.f = 0.95)
  }
})
names(site_sector_cols) <- site_sector_info$site_sector

site_order <- site_sector_info$site_sector
sector_order <- c(site_order, final_genera, "Others")
grid_cols <- c(site_sector_cols, genus_palette)
link_cols <- genus_palette[site_genus_long$genus_group]
gap_after <- c(
  rep(4, length(site_order) - 1), 10,
  rep(4, length(final_genera) - 1), 10,
  8
)

dir.create("figures", showWarnings = FALSE)
png("figures/Figure8_circos_sites_taxa.png", width = 7, height = 7, units = "in", res = 400, type = "cairo", bg = "white")
par(family = font_family)
circos.clear()
circos.par(start.degree = 90, gap.after = gap_after)

chordDiagram(
  x = site_genus_long %>% select(site_sector, genus_group, weight),
  grid.col = grid_cols,
  order = sector_order,
  col = link_cols,
  transparency = 0.25,
  annotationTrack = "grid",
  preAllocateTracks = 1
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    circos.text(
      x = get.cell.meta.data("xcenter"),
      y = get.cell.meta.data("ylim")[2] - mm_y(1.5),
      labels = sector_name,
      facing = "clockwise",
      niceFacing = TRUE,
      cex = 0.65
    )
  },
  bg.border = NA
)
dev.off()

cat("Saved figure to figures/Figure8_circos_sites_taxa.png\n")

# ==============================================================================
# Figure 9. Alluvial: Source → Site → Farming style → Genus (top 12)
# ==============================================================================

top_n_alluv <- 12
alluv_df <- genus_rel %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "genus", values_to = "rel_abund") %>%
  left_join(meta, by = "sample_id") %>%
  mutate(
    genus = case_when(
      str_detect(tolower(genus), "unclassified") ~ "Unclassified",
      str_detect(tolower(genus), "unidentified|unassigned") ~ "Unidentified",
      TRUE ~ genus
    )
  ) %>%
  group_by(genus) %>%
  mutate(genus_total = sum(rel_abund, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(genus_group = forcats::fct_lump_n(genus, n = top_n_alluv, w = genus_total, other_level = "Others")) %>%
  group_by(Source, Site, Farming_style, genus_group) %>%
  summarise(weight = sum(rel_abund, na.rm = TRUE), .groups = "drop") %>%
  filter(weight > 0, !genus_group %in% c("Unclassified", "Unidentified", "Others")) %>%
  mutate(genus_group = forcats::fct_drop(genus_group))

# palette (cityscape-inspired gradient)
genus_levels_alluv <- levels(alluv_df$genus_group)
genus_palette_alluv <- setNames(cityscape_pal_extended(length(genus_levels_alluv)), genus_levels_alluv)

p_alluv <- ggplot(
  alluv_df,
  aes(y = weight, axis1 = Source, axis2 = Site, axis3 = Farming_style, axis4 = genus_group, fill = genus_group)
) +
  geom_alluvium(alpha = 0.8, color = "grey15", linewidth = 0.1) +
  geom_stratum(color = "grey10", linewidth = 0.2, fill = "grey95") +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 3,
    color = "black",
    family = font_family
  ) +
  scale_fill_manual(values = genus_palette_alluv, name = "Genus") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Flow of relative abundance across source, site, farming style, and genus",
    y = "Relative abundance", x = NULL
  ) +
  theme_minimal(base_size = 12, base_family = font_family) +
  theme(
    plot.background = element_rect(fill = "#fdfdfd", color = NA),
    panel.background = element_rect(fill = "#f7f9fb", color = NA),
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(family = font_family, face = "plain", size = 16, color = "#0f1c2d"),
    axis.text = element_text(color = "#1f2d3d")
  )

ggsave("figures/Figure9_alluvial_source_site_style_genus.png", p_alluv,
       width = 9, height = 5.8, dpi = 400, type = "cairo", bg = "white")

cat("Saved figure to figures/Figure9_alluvial_source_site_style_genus.png\n")

# ==============================================================================
# Figure 10. Ridgeline of alpha diversity (Shannon, Faith) by Site & Source
# ==============================================================================

alpha_df_full <- read_csv("alpha_div_with_meta.csv", show_col_types = FALSE) %>%
  mutate(
    Site = factor(Province, levels = c("A1", "A2", "A3", "A4", "A5"),
                  labels = c("Site A", "Site B", "Site C", "Site D", "Site E")),
    Source = factor(Source, levels = c("Water", "Gill"))
  ) %>%
  pivot_longer(c(Shannon, Faith_PD), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("Shannon", "Faith_PD"))) %>%
  drop_na(value)

metric_cols <- unname(c(site_palette[2], site_palette[5]))

p_ridge <- ggplot(alpha_df_full, aes(x = value, y = Site, fill = metric)) +
  geom_density_ridges(alpha = 0.75, color = "white", scale = 1.2, linewidth = 0.2) +
  scale_fill_manual(
    values = c("Shannon" = "#1f3c5d", "Faith_PD" = "#972837"),
    name = "Metric",
    breaks = levels(alpha_df_full$metric)
  ) +
  facet_wrap(~Source, scales = "free_x") +
  labs(
    title = "Alpha diversity distributions by source",
    x = "Diversity value", y = "Site"
  ) +
  theme_minimal(base_size = 12, base_family = font_family) +
  theme(
    plot.background = element_rect(fill = "#fdfdfd", color = NA),
    panel.background = element_rect(fill = "#f7f9fb", color = NA),
    legend.position = "top",
    strip.text = element_text(family = font_family, face = "plain", size = 12),
    plot.title = element_text(family = font_family, face = "plain", size = 16, color = "#0f1c2d"),
    axis.text = element_text(color = "#1f2d3d")
  )

ggsave("figures/Figure10_ridgeline_alpha_diversity.png", p_ridge,
       width = 7, height = 7, dpi = 400, type = "cairo", bg = "white")

cat("Saved figure to figures/Figure10_ridgeline_alpha_diversity.png\n")

# ==============================================================================
# Figure 11. Spearman correlation heatmap among top genera
# ==============================================================================

top_corr <- 20
genus_totals <- colSums(genus_rel, na.rm = TRUE)
top_corr_genus <- names(sort(genus_totals, decreasing = TRUE))[seq_len(min(top_corr, length(genus_totals)))]
genus_mat_top <- genus_rel[, top_corr_genus, drop = FALSE]
# Drop grey groups if present
keep_cols_corr <- setdiff(colnames(genus_mat_top), c("Unclassified", "Unidentified"))
genus_mat_top <- genus_mat_top[, keep_cols_corr, drop = FALSE]
corr_mat <- cor(genus_mat_top, method = "spearman")
corr_df <- as.data.frame(as.table(corr_mat))
colnames(corr_df) <- c("Genus1", "Genus2", "rho")

p_corr <- ggplot(corr_df, aes(Genus1, Genus2, fill = rho)) +
  geom_tile(color = "white", size = 0.2) +
  scale_fill_gradient2(low = "#221f1f", mid = "white", high = "#d1495b",
                       limits = c(-1, 1), name = "Spearman") +
  coord_equal() +
  labs(
    title = "Co-variation among dominant genera (Spearman)",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11, base_family = font_family) +
  theme(
    plot.background = element_rect(fill = "#fdfdfd", color = NA),
    panel.background = element_rect(fill = "#f7f9fb", color = NA),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 9, family = font_family, face = "plain", color = "#1f2d3d"),
    axis.text.y = element_text(size = 9, family = font_family, face = "plain", color = "#1f2d3d"),
    plot.title = element_text(family = font_family, face = "plain", size = 16, color = "#0f1c2d"),
    panel.grid = element_blank()
  )

ggsave("figures/Figure11_correlation_top_genera.png", p_corr,
       width = 7, height = 6.2, dpi = 400, type = "cairo", bg = "white")

cat("Saved figure to figures/Figure11_correlation_top_genera.png\n")

# ==============================================================================
# Figure 12. Polar stacked bars: genus composition by Site x Source
# ==============================================================================

top_polar <- 10
polar_df <- genus_rel %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "genus", values_to = "rel_abund") %>%
  left_join(meta, by = "sample_id") %>%
  mutate(
    genus = case_when(
      str_detect(tolower(genus), "unclassified") ~ "Unclassified",
      str_detect(tolower(genus), "unidentified|unassigned") ~ "Unidentified",
      TRUE ~ genus
    )
  ) %>%
  group_by(genus) %>%
  mutate(genus_total = sum(rel_abund, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(genus_group = forcats::fct_lump_n(genus, n = top_polar, w = genus_total, other_level = "Others")) %>%
  group_by(Site, Source, genus_group) %>%
  summarise(weight = sum(rel_abund, na.rm = TRUE), .groups = "drop") %>%
  group_by(Site, Source) %>%
  mutate(weight = weight / sum(weight, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!genus_group %in% c("Unclassified", "Unidentified", "Others")) %>%
  mutate(genus_group = forcats::fct_drop(genus_group))

missing_polar_cols <- setdiff(levels(polar_df$genus_group), names(genus_palette))
genus_palette_polar <- c(
  genus_palette,
  setNames(cityscape_pal(length(missing_polar_cols)), missing_polar_cols)
)
genus_palette_polar <- genus_palette_polar[!duplicated(names(genus_palette_polar), fromLast = TRUE)]

p_polar <- ggplot(polar_df, aes(x = 1, y = weight, fill = genus_group)) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y") +
  facet_grid(Source ~ Site) +
  scale_fill_manual(values = genus_palette_polar, name = "Genus") +
  labs(
    title = "Genus composition by site and source",
    x = NULL, y = NULL
  ) +
  theme_void(base_size = 9, base_family = font_family) +
  theme(
    text = element_text(family = font_family, face = "plain"),
    legend.position = "bottom",
    strip.text = element_text(family = font_family, face = "plain", size = 9),
    plot.title = element_text(family = font_family, face = "plain", size = 12),
    plot.margin = margin(6, 6, 6, 6)
  )

ggsave("figures/Figure12_polar_genus_composition.png", p_polar,
       width = 8.5, height = 6.5, dpi = 400, type = "cairo", bg = "white")

cat("Saved figure to figures/Figure12_polar_genus_composition.png\n")