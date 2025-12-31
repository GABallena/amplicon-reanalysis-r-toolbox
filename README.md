# amplicon-reanalysis-r-toolbox

Portfolio-safe **R + Bash templates** for auditing / re-analyzing **amplicon microbiome** datasets using common exports from **QIIME 2** and/or a **DADA2** workflow.

This repo is designed for “reviewer-style” rechecks:
- QC sanity checks (primer trimming + DADA2 tracking)
- alpha diversity + assumption checks
- beta diversity ordination + PERMANOVA
- constrained ordination (CCA)
- publication-ready figure scripts (rarefaction, NMDS, heatmaps, top taxa, targeted taxa)

**Confidentiality note:** this repository contains **anonymized templates only**. No client data is included, and you should not commit client data to this repo.

---

## What’s included

### Scripts (high level)
- `scripts/cutadapt.sh` — primer trimming + per-sample summary TSV
- `scripts/autotune_truncLen_slidingQ.R` — quick helper for choosing DADA2 truncation parameters from quality profiles
- `scripts/its_dada2_asv.R` — DADA2 inference scaffold for paired-end amplicons
- `scripts/QC_data_DADA2.R` — QC table parsing + basic checks
- `scripts/alpha_div_recheck.R`, `scripts/normality_recheck.R` — alpha diversity + assumption checks
- `scripts/beta_div_recheck.R` — beta diversity ordination + PERMANOVA
- `scripts/cca_recheck.R` — CCA (constrained ordination)
- `scripts/fig2_rarefaction.R` … `scripts/fig7_pathogenic_genera.R` — figure templates
- `scripts/combine_figures_1_7.R` — figure panel bundling (if you build multi-panel figures)

> Script names reflect one “paper-style” figure sequence; you can run them independently.

---

## Requirements

### System
- R (recommended: 4.2+)
- Bash + `cutadapt` (if you use `cutadapt.sh`)
- Enough RAM for your dataset size (DADA2 can be memory-heavy)

### R packages (used across scripts)
A subset is used per script; install what you need:
- `dada2`, `ShortRead`
- `dplyr`, `tidyr`, `readr`, `stringr`, `tibble`, `tidyverse`
- `ggplot2`, `scales`, `patchwork`, `ggdendro`
- `vegan`, `ape`, `car`
- `RColorBrewer`, `viridis`, `ggridges`, `ggalluvial`, `circlize`
- optional mapping helpers: `sf`, `grid`, `fs`

---

## Repo layout (recommended)

Keep data out of git; use a local `data/` folder ignored by `.gitignore`.

```
amplicon-reanalysis-r-toolbox/
  scripts/
  data/                  # local only (ignored)
    raw_fastq/           # optional
    trimmed/             # optional
    exported_feature_table/
    exported_taxonomy/
    exported_bray/
    exported_alpha/
    exported_stats/
    shapes/              # optional, for mapping
  outputs/               # local only (ignored)
  README.md
```

### Suggested `.gitignore` (important)
At minimum, ignore:
- `data/`
- `outputs/`
- `*.fastq.gz`, `*.qza`, `*.qzv`
- large intermediate TSV/CSV that could leak sample IDs

---

## Quick start

### 1) Point scripts at your working directory
Most scripts are written to run from a “project directory” and use relative paths.

Example:

```bash
export PROJECT_DIR="/absolute/path/to/amplicon-reanalysis-r-toolbox"
cd "$PROJECT_DIR"
```

### 2) Configure inputs (either via env vars or by editing the CONFIG section in scripts)
Some scripts read inputs via environment variables. Common ones:

- `PROJECT_DIR` (default: current working directory)
- `METADATA_CSV` (default: `metadata.csv`)
- `MAP_CSV` (default: `map.csv`)
- `FEATURE_TABLE_TSV` (default: `exported_feature_table/feature-table.tsv`)
- `DISTANCE_MATRIX_TSV` (default: `exported_bray/distance-matrix.tsv`)
- `DADA2_STATS_TSV` (default: `exported_stats/stats.tsv`)
- `CUTADAPT_SUMMARY_TSV` (default: `logs/cutadapt_summary.tsv`)
- `TRIMMED_DIR` (default: `data/trimmed`)
- `RAREFACTION_DEPTH` (default in templates: `12000`)

### 3) Run a script
```bash
PROJECT_DIR="." Rscript scripts/beta_div_recheck.R
```

---

## Input formats expected (common cases)

These templates assume you have some combination of:

- **Metadata table** (CSV): sample IDs + grouping variables (e.g., source/site/condition)
- **QIIME 2 exports**:
  - Feature table: `feature-table.tsv` (from a `.qza` export)
  - Taxonomy: `taxonomy.tsv` (from a `.qza` export)
  - Distance matrix: `distance-matrix.tsv` (e.g., Bray–Curtis)
  - Alpha diversity tables (TSV) (e.g., Shannon/Faith PD)
- **DADA2 tracking stats** (TSV/CSV)
- **Cutadapt per-sample summary** (TSV) from `scripts/cutadapt.sh`

If your exports are different, edit the “read_*” lines near the top of each script.

### QIIME 2 export reminder (example)
```bash
qiime tools export --input-path feature-table.qza --output-path exported_feature_table
qiime tools export --input-path taxonomy.qza     --output-path exported_taxonomy
qiime tools export --input-path distance.qza     --output-path exported_bray
```

---

## Notes on anonymity / publishing

If you’re using these as portfolio artifacts:
- Do not include client names, institutions, locations, sample IDs, coordinates, or raw sequences.
- Keep real data in `data/` (ignored) and commit only templates + documentation.
- If you want to demonstrate “proof,” prefer **public datasets** or **synthetic toy data**.

---

## License

Pick a license you’re comfortable with (MIT is common for code templates).  
If you include any client-facing deliverables, ensure you have permission before publishing.

---

## Citation (optional)

If you want to cite this repo in a portfolio or application, use a simple statement like:

> “A portfolio-safe template repo demonstrating my approach to reproducible amplicon QC and re-analysis using QIIME 2 exports and DADA2.”

