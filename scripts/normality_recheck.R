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


library(readr)
library(dplyr)
library(car)

## 1. Read alpha diversity tables --------------------------

# Prefer local files; fall back to exported_* folders
shannon_path <- if (file.exists("shannon-alpha-diversity.tsv")) {
  "shannon-alpha-diversity.tsv"
} else {
  "exported_shannon_11979/shannon-alpha-diversity.tsv"
}

faith_path <- if (file.exists("faith-alpha-diversity.tsv")) {
  "faith-alpha-diversity.tsv"
} else {
  "exported_faith_11979/faith-alpha-diversity.tsv"
}

shannon <- read_tsv(
  shannon_path,
  skip = 1,                        # first line only has the metric name
  comment = "#",
  col_names = c("sample_id", "Shannon"),
  col_types = cols(
    sample_id = col_character(),
    Shannon   = col_double()
  ),
  show_col_types = FALSE
)

faith <- read_tsv(
  faith_path,
  skip = 1,
  comment = "#",
  col_names = c("sample_id", "Faith_PD"),
  col_types = cols(
    sample_id = col_character(),
    Faith_PD  = col_double()
  ),
  show_col_types = FALSE
)

cat("Shannon n =", nrow(shannon), "Faith n =", nrow(faith), "\n")

## 2. Read metadata and fix the one mislabel if needed -----

meta <- read_csv(Sys.getenv("METADATA_CSV", unset = "metadata.csv"))

# Optional: reconcile any known sample-ID inconsistencies here (portfolio template leaves as-is).
## 3. Join everything into one table -----------------------

alpha <- shannon %>%
  inner_join(faith, by = "sample_id") %>%
  inner_join(meta,  by = "sample_id") %>%
  mutate(
    Source        = factor(Source),
    Province_code = factor(Province_code),
    Farming_style = factor(Farming_style)
  )

cat("Joined alpha-div + metadata n =", nrow(alpha), "\n")

## 4. Shapiro–Wilk on ANOVA RESIDUALS (correct way) --------

mod_shannon <- lm(Shannon ~ Source + Province_code + Farming_style, data = alpha)
mod_faith   <- lm(Faith_PD ~ Source + Province_code + Farming_style, data = alpha)

shap_shannon <- shapiro.test(residuals(mod_shannon))
shap_faith   <- shapiro.test(residuals(mod_faith))

cat("\nShapiro–Wilk on ANOVA residuals:\n")
cat("Shannon: W =", shap_shannon$statistic,
    "p =", shap_shannon$p.value, "\n")
cat("Faith PD: W =", shap_faith$statistic,
    "p =", shap_faith$p.value, "\n")

## 5. Levene’s test for homogeneity (optional but nice) ----

cat("\nLevene tests (homogeneity of variance):\n")
print(car::leveneTest(Shannon ~ Source * Province_code * Farming_style, data = alpha))
print(car::leveneTest(Faith_PD ~ Source * Province_code * Farming_style, data = alpha))

## 6. ANOVA tables (to cross-check F and p values) ---------

cat("\nTwo-way ANOVA (Type I by default):\n")
cat("\nShannon ID:\n")
print(anova(mod_shannon))

cat("\nFaith PD:\n")
print(anova(mod_faith))
