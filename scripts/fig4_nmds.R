#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# PORTFOLIO-SAFE TEMPLATE
# - This script is anonymized: no client names, no site names, no sample IDs.
# - Provide your own input files (feature tables / metadata / taxonomy exports).
# - Set PROJECT_DIR to the folder containing the expected inputs.
# ------------------------------------------------------------------------------


# Fig. 4. NMDS of fungal community composition (Bray–Curtis on relative abundance).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
library(vegan)
})

set.seed(123)
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
meta_path <- "map_with_metadata.csv"

# ---- Feature table (ASVs) ----
feat_raw <- read_tsv(feat_path, skip = 1, col_types = cols())
colnames(feat_raw)[1] <- "FeatureID"

feat_mat <- feat_raw %>%
  column_to_rownames("FeatureID") %>%
  as.matrix()

# samples as rows
comm <- t(feat_mat)

# ---- Metadata ----
meta <- read_csv(meta_path, show_col_types = FALSE) %>%
  mutate(    sample_id = as.character(sample_id),
Source = factor(Source, levels = c("Water", "Gill")),
    Farming_style = factor(Farming_style, levels = c("Open cage", "Closed pond"))
  )

# Align samples
common <- intersect(rownames(comm), meta$sample_id)
comm   <- comm[common, , drop = FALSE]
meta   <- meta %>% filter(sample_id %in% common)

# Relative abundance
comm_rel <- decostand(comm, method = "total")

# ---- NMDS helper ----
run_nmds <- function(mat, meta_sub, label) {
  if (nrow(mat) < 3) {
    return(list(df = NULL, stress = NA))
  }
  ord <- metaMDS(mat, distance = "bray", k = 2, trymax = 100, trace = FALSE)
  scores_df <- as.data.frame(scores(ord, display = "sites"))
  out <- bind_cols(meta_sub, scores_df) %>%
    mutate(dataset = label)
  list(df = out, stress = ord$stress)
}

all_res   <- run_nmds(comm_rel, meta, "A. All samples")
water_idx <- meta$Source == "Water"
gill_idx  <- meta$Source == "Gill"
water_res <- run_nmds(comm_rel[water_idx, , drop = FALSE], meta[water_idx, ], "B. Water only")
gill_res  <- run_nmds(comm_rel[gill_idx, , drop = FALSE], meta[gill_idx, ], "C. Gill only")

ord_df <- bind_rows(all_res$df, water_res$df, gill_res$df)

cat("Stress (all):", all_res$stress, "\n")
cat("Stress (water):", water_res$stress, "\n")
cat("Stress (gill):", gill_res$stress, "\n")

# ---- PERMANOVA for all samples ----
bray_all <- vegdist(comm_rel, method = "bray")
perma <- adonis2(bray_all ~ Source + Province + Farming_style,
                 data = meta,
                 permutations = 999)
cat("\nPERMANOVA (999 perms):\n")
print(perma)

# ---- Hulls (by Source) for the combined plot only ----
hull_df <- ord_df %>%
  filter(dataset == "A. All samples") %>%
  group_by(dataset, Source) %>%
  filter(n() >= 3) %>%
  slice(chull(NMDS1, NMDS2)) %>%
  ungroup()

# ---- Plot ----
cols <- c("Water" = "#1f3c5d", "Gill" = "#d1495b")
shapes <- c("Open cage" = 21, "Closed pond" = 24)

p <- ggplot() +
  geom_polygon(
    data = hull_df,
    aes(x = NMDS1, y = NMDS2, fill = Source, group = Source),
    alpha = 0.15,
    color = NA
  ) +
  geom_point(
    data = ord_df,
    aes(x = NMDS1, y = NMDS2, color = Source, shape = Farming_style),
    size = 2.8,
    stroke = 0.8,
    alpha = 0.9
  ) +
  facet_wrap(~ dataset, scales = "free") +
  scale_color_manual(values = cols, name = "Sampling source") +
  scale_fill_manual(values = cols, name = "Sampling source") +
  scale_shape_manual(values = shapes, name = "Farming style") +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    title = "NMDS of fungal community composition (Bray–Curtis)"
  ) +
  theme_minimal(base_size = 12, base_family = font_family) +
  theme(
    plot.background = element_rect(fill = "#fdfdfd", color = NA),
    panel.background = element_rect(fill = "#f7f9fb", color = NA),
    legend.position = "bottom",
    strip.text = element_text(family = font_family, face = "plain", size = 12),
    panel.grid.minor = element_blank(),
    plot.title = element_text(color = "#0f1c2d", size = 16),
    axis.text = element_text(color = "#1f2d3d")
  )

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure4_NMDS.png", p, width = 9, height = 5, dpi = 400, type = "cairo", bg="white")
ggsave("figures/Figure4_NMDS.tiff", p, width = 9, height = 5, dpi = 400,
       type = "cairo", device = "tiff", compression = "lzw", bg = "white")
ggsave("figures/Figure4_NMDS.pdf", p, width = 9, height = 5, device = cairo_pdf, bg = "white")

cat("Saved figure to figures/Figure4_NMDS.png\n")