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


# Autotune truncLenF/R deterministically using sliding-window mean Phred
# + robust quantile threshold + DADA2 grid search objective.

conda_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib", "R", "library")
.libPaths(conda_lib)

suppressPackageStartupMessages({
  library(dada2)
  library(jsonlite)
})

# -------- CONFIG --------
trimmed_path <- Sys.getenv("TRIMMED_DIR", unset = "data/trimmed")
q_cut <- 20                 # Q20 threshold
quantile_cut <- 0.10        # 0.10 = 10th percentile (robust); try 0.25 if you want stricter
win <- 10                   # sliding window size (bp); 5 or 10 are common choices
min_overlap <- 20           # stricter than default 12 (more conservative for validation)
max_mismatch <- 0
maxEE <- c(2, 2)
truncQ <- 20
n_reads_per_sample <- 20000 # downsample per sample for speed + determinism
n_samples <- 12             # number of samples used for tuning (sorted deterministic)
grid_step <- 10             # step size for truncLen grid search
grid_span <- 60             # span around initial estimate for grid (bp)
min_len <- 50               # minimum truncLen to consider
max_len <- 300              # MiSeq 2x300
# ------------------------

# Bioc dependency for robust FASTQ quality parsing
if (!requireNamespace("ShortRead", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ShortRead", ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages(library(ShortRead))

# --- locate files ---
fnFs <- sort(list.files(trimmed_path, pattern = "_R1_trim.fastq.gz$", full.names = TRUE))
fnRs <- sort(list.files(trimmed_path, pattern = "_R2_trim.fastq.gz$", full.names = TRUE))

if (length(fnFs) == 0 || length(fnRs) == 0) {
  stop("No trimmed FASTQs found. Expected trimmed/*_R1_trim.fastq.gz and trimmed/*_R2_trim.fastq.gz")
}
if (length(fnFs) != length(fnRs)) stop("R1/R2 file count mismatch.")

use_idx <- seq_len(min(n_samples, length(fnFs)))
fnFs_sub <- fnFs[use_idx]
fnRs_sub <- fnRs[use_idx]

cat("Autotuning using", length(fnFs_sub), "samples (deterministic subset).\n")

# --- helper: sliding window mean along positions, per read ---
# Input: integer matrix (reads x positions)
# Output: numeric matrix (reads x positions), where value at pos j = mean of last win bases ending at j
sliding_mean_mat <- function(mat, win = 10) {
  nr <- nrow(mat); nc <- ncol(mat)
  out <- matrix(NA_real_, nrow = nr, ncol = nc)

  # cumulative sum for fast rolling mean
  cs <- t(apply(mat, 1, cumsum))
  # For j < win, use mean over 1..j
  for (j in seq_len(min(win - 1, nc))) {
    out[, j] <- cs[, j] / j
  }
  if (nc >= win) {
    for (j in win:nc) {
      out[, j] <- (cs[, j] - cs[, j - win]) / win
    }
  }
  out
}

# --- helper: get pooled per-position quantile of sliding-window mean quality ---
per_pos_sliding_quantile <- function(fns, q = 0.10, n_reads = 20000, win = 10) {
  mats <- list()
  maxlen <- 0

  for (fp in fns) {
    fq <- readFastq(fp)
    if (length(fq) > n_reads) fq <- fq[seq_len(n_reads)]
    qmat <- as(quality(fq), "matrix")  # reads x positions
    maxlen <- max(maxlen, ncol(qmat))
    mats[[fp]] <- qmat
  }

  # pad to maxlen so rbind works
  padded <- lapply(mats, function(m) {
    if (ncol(m) < maxlen) {
      pad <- matrix(NA_integer_, nrow = nrow(m), ncol = maxlen - ncol(m))
      cbind(m, pad)
    } else m
  })

  pooled <- do.call(rbind, padded)

  # Convert to sliding-window mean; treat NA as missing
  # Replace NA with 0 for cumsum, but track valid lengths; easiest is per-row apply with NA handling
  # We'll compute sliding mean ignoring NA by truncating at last non-NA column per read.
  # Approach: compute sliding mean on each row's valid prefix, then pad with NA.
  nr <- nrow(pooled); nc <- ncol(pooled)
  sm <- matrix(NA_real_, nrow = nr, ncol = nc)

  # chunk to avoid huge memory spikes
  chunk_size <- 5000
  starts <- seq(1, nr, by = chunk_size)
  for (s in starts) {
    e <- min(nr, s + chunk_size - 1)
    block <- pooled[s:e, , drop = FALSE]

    for (i in seq_len(nrow(block))) {
      row <- block[i, ]
      last <- max(which(!is.na(row)), 0)
      if (last > 0) {
        rowv <- row[seq_len(last)]
        # rolling mean on prefix
        cs <- cumsum(rowv)
        outv <- numeric(last)
        if (last >= 1) outv[1:min(win - 1, last)] <- cs[1:min(win - 1, last)] / seq_len(min(win - 1, last))
        if (last >= win) {
          for (j in win:last) outv[j] <- (cs[j] - cs[j - win]) / win
        }
        sm[s + i - 1, seq_len(last)] <- outv
      }
    }
  }

  # quantile per position
  qvec <- apply(sm, 2, function(col) {
    col <- col[!is.na(col)]
    if (length(col) == 0) return(NA_real_)
    as.numeric(quantile(col, probs = q, names = FALSE, type = 7))
  })
  qvec
}

# compute profiles
qF <- per_pos_sliding_quantile(fnFs_sub, q = quantile_cut, n_reads = n_reads_per_sample, win = win)
qR <- per_pos_sliding_quantile(fnRs_sub, q = quantile_cut, n_reads = n_reads_per_sample, win = win)

# save profiles for audit trail
profile_df <- data.frame(
  position = seq_along(qF),
  slidingQ_quantile_F = qF,
  slidingQ_quantile_R = qR
)
write.table(profile_df, file = "autotune_quality_profiles.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote autotune_quality_profiles.tsv\n")

# propose trunc: largest j such that slidingQ_quantile >= q_cut for all positions up to j
propose_trunc <- function(qvec, q_cut) {
  if (all(is.na(qvec))) return(0)
  fail <- which(qvec < q_cut)
  if (length(fail) == 0) return(length(qvec))
  max(0, min(fail) - 1)
}

truncF0 <- propose_trunc(qF, q_cut)
truncR0 <- propose_trunc(qR, q_cut)

cat("Initial trunc proposals from slidingQ rule:\n",
    "  truncF0=", truncF0, " truncR0=", truncR0, "\n", sep = "")

grid_around <- function(x, minx = 50, maxx = 300, step = 10, span = 60) {
  lo <- max(minx, x - span)
  hi <- min(maxx, x + span)
  unique(seq(floor(lo / step) * step, ceiling(hi / step) * step, by = step))
}

candF <- grid_around(truncF0, minx = min_len, maxx = max_len, step = grid_step, span = grid_span)
candR <- grid_around(truncR0, minx = min_len, maxx = max_len, step = grid_step, span = grid_span)

# ---------- evaluate candidate pairs via lightweight DADA2 ----------
eval_pair <- function(trF, trR) {
  filt_path <- file.path(trimmed_path, "autotune_filtered")
  dir.create(filt_path, showWarnings = FALSE)

  sample.names <- sub("_R1_trim.fastq.gz$", "", basename(fnFs_sub))
  filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

  out <- filterAndTrim(
    fnFs_sub, filtFs, fnRs_sub, filtRs,
    truncLen = c(trF, trR),
    maxN = 0,
    maxEE = maxEE,
    truncQ = truncQ,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = TRUE
  )
  keep <- out[, "reads.out"] > 0
  if (sum(keep) < 3) return(NULL)

  filtFs <- filtFs[keep]; filtRs <- filtRs[keep]
  sn <- sample.names[keep]

  errF <- learnErrors(filtFs, multithread = TRUE)
  errR <- learnErrors(filtRs, multithread = TRUE)

  derepFs <- derepFastq(filtFs); derepRs <- derepFastq(filtRs)
  names(derepFs) <- sn; names(derepRs) <- sn

  dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
  dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

  mergers <- mergePairs(
    dadaFs, derepFs, dadaRs, derepRs,
    minOverlap = min_overlap,
    maxMismatch = max_mismatch,
    verbose = FALSE
  )

  seqtab <- makeSequenceTable(mergers)
  if (ncol(seqtab) == 0) return(NULL)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                      multithread = TRUE, verbose = FALSE)

  data.frame(
    truncLenF = trF,
    truncLenR = trR,
    n_samples = nrow(seqtab.nochim),
    merged_reads = sum(rowSums(seqtab)),
    nonchim_reads = sum(rowSums(seqtab.nochim)),
    n_asv = ncol(seqtab.nochim)
  )
}

results <- list()
k <- 1

for (trF in candF) {
  for (trR in candR) {
    cat("Testing truncLenF=", trF, " truncLenR=", trR, "\n", sep = "")
    res <- tryCatch(eval_pair(trF, trR), error = function(e) NULL)
    if (!is.null(res)) { results[[k]] <- res; k <- k + 1 }
  }
}

if (length(results) == 0) stop("No candidate truncLen pairs produced usable output. Loosen constraints or inspect trimming.")

df <- do.call(rbind, results)

# Objective: maximize nonchim_reads
# Tie-breakers: maximize n_samples, then maximize merged_reads, then prefer longer truncs (less info loss)
df <- df[order(-df$nonchim_reads, -df$n_samples, -df$merged_reads, -df$truncLenF, -df$truncLenR), ]
best <- df[1, ]

write.table(df, file = "autotune_trunclen_grid.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote autotune_trunclen_grid.tsv\n")

params <- list(
  tuning = list(
    q_cut = q_cut,
    quantile_cut = quantile_cut,
    window_bp = win,
    n_reads_per_sample = n_reads_per_sample,
    n_samples_used = length(fnFs_sub)
  ),
  merge = list(
    min_overlap = min_overlap,
    max_mismatch = max_mismatch
  ),
  chosen = list(
    truncLenF = as.integer(best$truncLenF),
    truncLenR = as.integer(best$truncLenR),
    maxEE = maxEE,
    truncQ = truncQ
  )
)

write(toJSON(params, pretty = TRUE, auto_unbox = TRUE), file = "chosen_params.json")
cat("\nChosen params:\n")
print(params$chosen)
cat("\nWrote chosen_params.json\n")
