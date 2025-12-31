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


# --- make conda R libs available, but don't nuke existing libPaths ---
conda_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib", "R", "library")
.libPaths(c(conda_lib, .libPaths()))

suppressPackageStartupMessages({
  library(dada2)
  library(dplyr)
})

# --- paths ---
base_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
trim_dir <- file.path(base_dir, "data", "trimmed")  # expects paired-end FASTQs after primer trimming
filt_dir <- file.path(trim_dir, "filtered")

dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)

# --- input fastqs ---
fnFs <- sort(list.files(trim_dir, pattern = "_R1_trim\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(trim_dir, pattern = "_R2_trim\\.fastq\\.gz$", full.names = TRUE))

stopifnot(length(fnFs) == length(fnRs), length(fnFs) > 0)

# robust sample IDs: strip lane/index suffix (edit regex for your run naming)
sample.names <- basename(fnFs) |>
  sub("_S[0-9]+_L[0-9]+_R1_trim\\.fastq\\.gz$", "", x = _) |>
  sub("_S[0-9]+_L[0-9]+_R1_trim\\.fastq\\.gz$", "", x = _)  # just in case no leading zeros

cat("Found", length(sample.names), "samples.\n")
print(head(sample.names))

# --- filter & trim ---
filtFs <- file.path(filt_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(0, 0),     # keep full length (no fixed truncation)
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

cat("\nFilterAndTrim summary:\n")
print(out)

# keep only samples with reads after filtering
keep <- out[, "reads.out"] > 0
cat("\nRetaining", sum(keep), "of", length(keep), "samples after filtering.\n")

filtFs <- filtFs[keep]
filtRs <- filtRs[keep]
sample.names <- sample.names[keep]
out <- out[keep, , drop = FALSE]

# --- learn errors ---
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# --- derep + denoise ---
derepFs <- derepFastq(filtFs); names(derepFs) <- sample.names
derepRs <- derepFastq(filtRs); names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# --- merge pairs (ITS often struggles here if overlap is short) ---
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

seqtab <- makeSequenceTable(mergers)
cat("\nASV table dimensions (samples x ASVs):", dim(seqtab), "\n")

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)
cat("After chimera removal:", dim(seqtab.nochim), " (samples x ASVs)\n")

# --- tracking table ---
track <- data.frame(
  SampleID   = sample.names,
  reads.in   = out[, "reads.in"],
  reads.filt = out[, "reads.out"],
  denoisedF  = sapply(dadaFs, function(x) sum(getUniques(x))),
  denoisedR  = sapply(dadaRs, function(x) sum(getUniques(x))),
  merged     = rowSums(seqtab),
  nonchim    = rowSums(seqtab.nochim),
  stringsAsFactors = FALSE
) |>
  mutate(
    filt.prop    = round(reads.filt / reads.in, 3),
    merged.prop  = round(merged / reads.in, 3),
    nonchim.prop = round(nonchim / reads.in, 3)
  )

summary_path <- file.path(base_dir, "its_dada2_denoising_stats.tsv")
write.table(track, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n===== DADA2 READ SUMMARY =====\n")
cat("Per-sample table:", summary_path, "\n")
cat("Median merged.prop:", median(track$merged.prop, na.rm = TRUE), "\n")
cat("Median nonchim.prop:", median(track$nonchim.prop, na.rm = TRUE), "\n")

# --- export ASVs FASTA ---
asv_seqs <- colnames(seqtab.nochim)
asv_ids  <- paste0("ASV", seq_along(asv_seqs))

fasta_path <- file.path(base_dir, "its_asv_sequences.fasta")
con <- file(fasta_path, "w")
for (i in seq_along(asv_seqs)) {
  writeLines(paste0(">", asv_ids[i]), con)
  writeLines(asv_seqs[i], con)
}
close(con)

asv_table_path <- file.path(base_dir, "its_asv_table.tsv")
asv_tab <- cbind(SampleID = rownames(seqtab.nochim), seqtab.nochim)
write.table(asv_tab, asv_table_path, quote = FALSE, sep = "\t", row.names = FALSE)

cat("Wrote FASTA:", fasta_path, "\n")
cat("Wrote ASV table:", asv_table_path, "\n")
