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


# alpha_div_reanalysis.R
# Re-analysis of alpha diversity (Shannon, Faith PD) for an amplicon dataset (portfolio template)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(car)       # leveneTest
})

## 0. Helper: nice printing ----
cat_rule <- function(x) {
  cat("\n", strrep("=", nchar(x) + 4), "\n  ", x, "\n", strrep("=", nchar(x) + 4), "\n", sep = "")
}

## 1. Parse metadata from map.csv ----
cat_rule("Reading and parsing metadata")

map_raw <- read_csv(Sys.getenv("MAP_CSV", unset = "map.csv"), show_col_types = FALSE)

meta <- map_raw %>%
  rename(sample_id = `sample-id`) %>%
  mutate(
    sample_id    = str_trim(sample_id),
    Source       = if_else(str_starts(sample_id, "W-"), "Water", "Gill"),
    core         = str_remove(sample_id, "^W-"),        # drop W- for water samples
    core_no_ITS  = str_remove(core, "^ITS-"),           # drop ITS-
    # e.g. A1-1-OP1-G1
    Province_code = str_extract(core_no_ITS, "^A[0-9]+"),
    Province      = Province_code,                      # simple mapping; adjust if needed
    Farm_code     = str_extract(core_no_ITS, "(OP[0-9]+|CP[0-9]+)"),
    Farming_style = case_when(
      str_detect(Farm_code, "^OP") ~ "Open cage",
      str_detect(Farm_code, "^CP") ~ "Closed pond",
      TRUE ~ NA_character_
    ),
    Gill_code     = str_extract(core_no_ITS, "G[0-9]+")
  ) %>%
  select(sample_id, Source, Province_code, Province,
         Farming_style, Farm_code, Gill_code)

write_csv(meta, Sys.getenv("METADATA_OUT", unset = "metadata.csv"))
cat("Parsed metadata written to map_with_metadata.csv\n")

## 2. Read alpha diversity vectors exported from QIIME 2 ----
# Before running this script, export from QIIME 2, e.g.:
#   qiime tools export \
#     --input-path diversity/shannon_vector.qza \
#     --output-path exported_shannon
#   qiime tools export \
#     --input-path diversity/faith_pd_vector.qza \
#     --output-path exported_faith

cat_rule("Reading alpha diversity metrics")

shannon <- read_tsv(
  "exported_shannon/alpha-diversity.tsv",
  skip = 1,                                   # first line only has the metric name
  comment = "#",
  col_names = c("sample_id", "Shannon"),
  col_types = cols(
    sample_id = col_character(),
    Shannon   = col_double()
  ),
  show_col_types = FALSE
)

faith <- read_tsv(
  "exported_faith/alpha-diversity.tsv",
  skip = 1,
  comment = "#",
  col_names = c("sample_id", "Faith_PD"),
  col_types = cols(
    sample_id = col_character(),
    Faith_PD  = col_double()
  ),
  show_col_types = FALSE
)

alpha <- shannon %>%
  inner_join(faith, by = "sample_id") %>%
  inner_join(meta, by = "sample_id") %>%
  mutate(
    Source        = factor(Source),
    Province_code = factor(Province_code),
    Farming_style = factor(Farming_style)
  )

# Drop the rare sample that did not reach the rarefaction depth (NA alpha-div)
alpha_valid <- alpha %>%
  filter(!is.na(Shannon), !is.na(Faith_PD))

write_csv(alpha_valid, "alpha_div_with_meta.csv")
cat("Alpha diversity with metadata written to alpha_div_with_meta.csv\n")
cat("Number of samples after filtering:", nrow(alpha_valid), "\n")

## 3. Descriptive summaries ----
cat_rule("Descriptive summaries")

alpha_valid %>%
  group_by(Source) %>%
  summarise(
    n          = n(),
    mean_Shann = mean(Shannon),
    sd_Shann   = sd(Shannon),
    mean_FPD   = mean(Faith_PD),
    sd_FPD     = sd(Faith_PD)
  ) %>%
  print()

## 4. Two-way ANOVA for Shannon (Source × Farming_style) ----
cat_rule("Two-way ANOVA for Shannon")

# Check assumptions quickly
shapiro_shan <- shapiro.test(alpha_valid$Shannon)
cat("Shapiro–Wilk for Shannon: W =", round(shapiro_shan$statistic, 3),
    "p =", shapiro_shan$p.value, "\n")

lev_shan <- leveneTest(Shannon ~ Source * Farming_style, data = alpha_valid)
cat("Levene test for Shannon (Source × Farming_style):\n")
print(lev_shan)

# Fit 2-way ANOVA
fit_shan <- aov(Shannon ~ Source * Farming_style, data = alpha_valid)
summary_shan <- summary(fit_shan)
print(summary_shan)

## 5. Kruskal–Wallis for Shannon & Faith PD by Source / Farming_style ----
cat_rule("Kruskal–Wallis tests")

kw_shan_source <- kruskal.test(Shannon ~ Source, data = alpha_valid)
cat("KW Shannon by Source: H =", kw_shan_source$statistic,
    "p =", kw_shan_source$p.value, "\n")

kw_shan_style <- kruskal.test(Shannon ~ Farming_style, data = alpha_valid)
cat("KW Shannon by Farming_style: H =", kw_shan_style$statistic,
    "p =", kw_shan_style$p.value, "\n")

kw_faith_source <- kruskal.test(Faith_PD ~ Source, data = alpha_valid)
cat("KW Faith PD by Source: H =", kw_faith_source$statistic,
    "p =", kw_faith_source$p.value, "\n")

kw_faith_style <- kruskal.test(Faith_PD ~ Farming_style, data = alpha_valid)
cat("KW Faith PD by Farming_style: H =", kw_faith_style$statistic,
    "p =", kw_faith_style$p.value, "\n")

## 6. Kruskal–Wallis by Province within each Source ----
cat_rule("Kruskal–Wallis by Province (stratified by Source)")

for (src in levels(alpha_valid$Source)) {
  subdat <- alpha_valid %>% filter(Source == src)
  cat("\nSource:", src, "\n")
  
  kw_shan_prov <- kruskal.test(Shannon ~ Province_code, data = subdat)
  cat("  Shannon ~ Province: H =", kw_shan_prov$statistic,
      "p =", kw_shan_prov$p.value, "\n")
  
  kw_faith_prov <- kruskal.test(Faith_PD ~ Province_code, data = subdat)
  cat("  Faith PD ~ Province: H =", kw_faith_prov$statistic,
      "p =", kw_faith_prov$p.value, "\n")
}

cat_rule("Done")
