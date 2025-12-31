#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# PORTFOLIO-SAFE TEMPLATE
# - This script is anonymized: no client names, no site names, no sample IDs.
# - Provide your own input files (feature tables / metadata / taxonomy exports).
# - Set PROJECT_DIR to the folder containing the expected inputs.
# ------------------------------------------------------------------------------


# Fig. 3. Alpha diversity of fungal communities in rearing water and gills
# across provinces and farming styles (Shannon ID and Faith's PD).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
library(ggplot2)
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

# ----- Paths -----
shannon_path <- "exported_shannon_11979/shannon-alpha-diversity.tsv"
faith_path   <- "exported_faith_11979/faith-alpha-diversity.tsv"
meta_path    <- "map_with_metadata.csv"

# ----- Read data -----
shannon_df <- read_tsv(shannon_path, col_types = cols())
faith_df   <- read_tsv(faith_path,   col_types = cols())

colnames(shannon_df)[1] <- "sample_id"
colnames(faith_df)[1]   <- "sample_id"

alpha_df <- shannon_df %>%
  left_join(faith_df, by = "sample_id")

meta <- read_csv(meta_path, show_col_types = FALSE) %>%
  mutate(    sample_id = as.character(sample_id),
Source = factor(Source, levels = c("Water", "Gill")),
    Farming_style = factor(Farming_style, levels = c("Open cage", "Closed pond"))
  )

# ----- Merge and tidy -----
dat <- alpha_df %>%
  left_join(meta, by = "sample_id") %>%
  filter(Source %in% levels(meta$Source)) %>%
  rename(Shannon = shannon_entropy, Faith_PD = faith_pd)

# Anonymized site labels (edit as needed)
province_levels <- c("Site A", "Site B", "Site C", "Site D", "Site E")

dat$Province <- factor(dat$Province, levels = province_levels)

# Long format for two panels: by Province and by Farming style
long_prov <- dat %>%
  pivot_longer(c(Shannon, Faith_PD), names_to = "metric", values_to = "value") %>%
  mutate(panel = "By province", group_var = Province)

long_style <- dat %>%
  pivot_longer(c(Shannon, Faith_PD), names_to = "metric", values_to = "value") %>%
  mutate(panel = "By farming style", group_var = Farming_style)

plot_df <- bind_rows(long_prov, long_style) %>%
  mutate(
    metric = factor(metric, levels = c("Shannon", "Faith_PD"),
                    labels = c("Shannon diversity", "Faith's PD")),
    panel = factor(panel, levels = c("By province", "By farming style"))
  )

cols <- c("Water" = "#1f3c5d", "Gill" = "#d1495b")
shapes <- c("Open cage" = 21, "Closed pond" = 24)

p <- ggplot(plot_df, aes(x = group_var, y = value,
                         fill = Source, colour = Source, shape = Farming_style)) +
  geom_boxplot(alpha = 0.6, colour = "gray40", width = 0.65, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
              size = 2.3, alpha = 0.8) +
  facet_grid(metric ~ panel, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cols, name = "Sampling source") +
  scale_colour_manual(values = cols, name = "Sampling source") +
  scale_shape_manual(values = shapes, name = "Farming style") +
  labs(
    x = NULL,
    y = "Alpha diversity",
    title = "Alpha diversity of fungal communities in rearing water and gills",
    subtitle = "Shannon index and Faith's PD across provinces and farming styles (rarefied to 11,979 reads)"
  ) +
  theme_minimal(base_size = 12, base_family = font_family) +
  theme(
    plot.background = element_rect(fill = "#fdfdfd", color = NA),
    panel.background = element_rect(fill = "#f7f9fb", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(family = font_family, face = "plain", size = 12),
    plot.title = element_text(family = font_family, face = "plain", size = 16, color = "#0f1c2d"),
    plot.subtitle = element_text(family = font_family, face = "plain", size = 12, color = "#2c3e50"),
    axis.text.x = element_text(angle = 28, hjust = 1, family = font_family, face = "plain", color = "#1f2d3d"),
    axis.text.y = element_text(color = "#1f2d3d")
  )

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure3_alpha_diversity.png", p, width = 10, height = 6, dpi = 400, type = "cairo", bg="white")
ggsave("figures/Figure3_alpha_diversity.tiff", p, width = 10, height = 6, dpi = 400,
       type = "cairo", device = "tiff", compression = "lzw", bg = "white")
ggsave("figures/Figure3_alpha_diversity.pdf", p, width = 10, height = 6, device = cairo_pdf, bg = "white")

cat("Saved figure to figures/Figure3_alpha_diversity.png\n")