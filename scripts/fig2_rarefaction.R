#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# PORTFOLIO-SAFE TEMPLATE
# - This script is anonymized: no client names, no site names, no sample IDs.
# - Provide your own input files (feature tables / metadata / taxonomy exports).
# - Set PROJECT_DIR to the folder containing the expected inputs.
# ------------------------------------------------------------------------------


# Fig. 2. Rarefaction curves for ITS libraries from rearing water and gill samples.
# Curves for observed ASVs, Shannon, and Faith's PD across sequencing depth.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
library(ggplot2)
library(vegan)
library(ape)
})

set.seed(123)  # reproducible rarefaction
# ---- Project directory ----
project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
setwd(project_dir)
target_depth <- as.integer(Sys.getenv("RAREFACTION_DEPTH", unset = "12000"))
pdf_fonts <- names(grDevices::pdfFonts())
font_family <- if ("Times New Roman" %in% pdf_fonts) {
  "Times New Roman"
} else if ("Times" %in% pdf_fonts) {
  "Times"
} else {
  "serif"
}

# ---- 1. Feature table ----
feat_path <- "exported_feature_table/feature-table.tsv"
feat_raw  <- read_tsv(feat_path, skip = 1, col_types = cols())
colnames(feat_raw)[1] <- "FeatureID"

feat_mat <- feat_raw %>%
  column_to_rownames("FeatureID") %>%
  as.matrix()

# samples as rows
comm      <- t(feat_mat)
lib_sizes <- rowSums(comm)

# ---- 2. Metadata ----
meta <- read_csv("map_with_metadata.csv", show_col_types = FALSE) %>%
  mutate(    sample_id = as.character(sample_id),
Source = factor(Source, levels = c("Water", "Gill"))
  )

# keep only samples present in BOTH feature table and metadata
common_samples <- intersect(rownames(comm), meta$sample_id)
comm <- comm[common_samples, , drop = FALSE]
meta <- meta %>% filter(sample_id %in% common_samples)

# recompute library sizes after subsetting
lib_sizes <- rowSums(comm)
min_depth_all <- min(lib_sizes)
max_depth     <- max(lib_sizes)
cat("Min library size:", min_depth_all, "reads\n")
cat("Max library size:", max_depth, "reads\n")

# ---- 3. Phylogenetic tree (for Faith's PD) ----
tree_qza <- list.files("tree_out", pattern = "tree_rooted\\.qza$", full.names = TRUE)
if (length(tree_qza) == 0) stop("No rooted tree qza found in tree_out/")

tree_contents <- unzip(tree_qza[1], list = TRUE)
nwk_file      <- tree_contents$Name[grepl("tree\\.nwk$", tree_contents$Name)]
if (length(nwk_file) == 0) stop("No tree.nwk found inside tree qza")

nwk_txt <- readLines(unz(tree_qza[1], nwk_file[1]))
tree    <- read.tree(text = nwk_txt)

# keep only tips that are actually in the feature table (helps speed/consistency)
keep_tips <- intersect(tree$tip.label, colnames(comm))
tree      <- keep.tip(tree, keep_tips)

faith_pd_fun <- function(counts) {
  present <- names(counts)[counts > 0]
  present <- intersect(present, tree$tip.label)
  if (length(present) == 0) return(0)
  sub <- keep.tip(tree, present)
  sum(sub$edge.length)
}

# ---- 4. Depth grid for rarefaction ----
depth_grid <- sort(unique(round(
  c(seq(1000, max_depth, length.out = 25), target_depth)
)))
depth_grid <- depth_grid[depth_grid <= max_depth]

# ---- 5. Rarefaction and metrics ----
compute_metrics <- function(depth) {
  elig <- lib_sizes >= depth
  if (!any(elig)) return(NULL)

  rare <- rrarefy(comm[elig, , drop = FALSE], depth)

  data.frame(
    sample_id = rownames(rare),
    depth     = depth,
    richness  = rowSums(rare > 0),
    shannon   = diversity(rare, index = "shannon"),
    faith_pd  = vapply(
      seq_len(nrow(rare)),
      function(i) faith_pd_fun(rare[i, ]),
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
}

metric_list <- lapply(depth_grid, compute_metrics)
metric_df   <- bind_rows(metric_list)

# samples retained at the chosen rarefaction depth (if present)
n_retained <- metric_df %>%
  filter(depth == target_depth) %>%
  distinct(sample_id) %>%
  nrow()
cat("Samples retained at depth", target_depth, ":", n_retained, "\n")

# ---- 6. Long format + summaries (SE ribbons) ----
metric_long <- metric_df %>%
  pivot_longer(
    cols      = c(richness, shannon, faith_pd),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  left_join(meta %>% select(sample_id, Source), by = "sample_id") %>%
  mutate(
    metric = factor(
      metric,
      levels = c("richness", "shannon", "faith_pd"),
      labels = c("A. Observed ASVs",
                 "B. Shannon diversity",
                 "C. Faith's PD")
    )
  )

summary_df <- metric_long %>%
  group_by(metric, Source, depth) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    n    = sum(!is.na(value)),
    se   = sd / sqrt(n),
    .groups = "drop"
  )

# ---- 7. Plot ----
cols <- c("Water" = "#1f3c5d", "Gill" = "#d1495b")
max_plot_depth <- max(depth_grid[depth_grid <= target_depth * 3], na.rm = TRUE)
if (!is.finite(max_plot_depth)) max_plot_depth <- max_depth

p <- ggplot() +
  # SE ribbons by source
  geom_ribbon(
    data = summary_df,
    aes(x = depth,
        ymin = mean - se,
        ymax = mean + se,
        fill = Source),
    alpha = 0.20
  ) +
  # per-sample curves
  geom_line(
    data = metric_long,
    aes(x = depth, y = value,
        group = sample_id, colour = Source),
    alpha = 0.45,
    linewidth = 0.35
  ) +
  # mean curves
  geom_line(
    data = summary_df,
    aes(x = depth, y = mean, colour = Source),
    linewidth = 1.0
  ) +
  # vertical line at chosen rarefaction depth
  geom_vline(
    xintercept = target_depth,
    linetype   = "dashed",
    colour     = "grey30"
  ) +
  scale_colour_manual(values = cols, name = "Sampling source") +
  scale_fill_manual(values = cols,   name = "Sampling source") +
  scale_x_continuous(limits = c(0, max_plot_depth)) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    x = "Sequencing depth (reads per sample)",
    y = "Alpha diversity"
  ) +
  theme_minimal(base_size = 12, base_family = font_family) +
  theme(
    plot.background = element_rect(fill = "#fdfdfd", color = NA),
    panel.background = element_rect(fill = "#f7f9fb", color = NA),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    strip.text       = element_text(family = font_family, face = "plain", size = 12),
    axis.title       = element_text(color = "#1f2d3d", size = 12),
    axis.text        = element_text(color = "#1f2d3d", size = 10),
    plot.title       = element_text(color = "#0f1c2d", size = 15),
    plot.subtitle    = element_text(color = "#2c3e50", size = 12)
  )

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure2_rarefaction.png", p,
       width = 9, height = 5.5, dpi = 400, type = "cairo", bg="white")
ggsave("figures/Figure2_rarefaction.tiff", p,
       width = 9, height = 5.5, dpi = 400, type = "cairo", device = "tiff",
       compression = "lzw", bg = "white")
ggsave("figures/Figure2_rarefaction.pdf", p,
       width = 9, height = 5.5, device = cairo_pdf, bg = "white")

cat("Saved figure to figures/Figure2_rarefaction.png\n")