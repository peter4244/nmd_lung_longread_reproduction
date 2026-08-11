#!/usr/bin/env Rscript
# =============================================================================
# 02_score_nmdetective_b.R
#
# Apply the published NMDetective-B 4-rule decision tree (Lindeboom et al.,
# Nat Genet 2019, Figure 1c) to each isoform's pre-computed feature row,
# producing a continuous NMD-efficacy estimate and a binary trigger/escape
# call.
#
# Decision tree (thresholds 50 / 150 / 407 nt):
#   1. In last exon?                       → escape (leaf "Last exon", score 0.00)
#   2. Distance to start codon < 150 nt?   → escape (leaf "Start-proximal", 0.12)
#   3. Exon containing PTC > 407 nt long?  → escape (leaf "Long exon", 0.41)
#   4. Within 50 nt of the last EJ?        → escape (leaf "50 nt rule", 0.20)
#   else                                   → trigger NMD (leaf "Trigger NMD", 0.65)
#
# Leaf NMD-efficacy scores are read off the published Figure 1c percentages
# (lower = stronger escape; higher = stronger trigger):
#   last_exon = 0.00,  start_proximal = 0.12,  long_exon = 0.41,
#   50nt_rule = 0.20,  trigger        = 0.65
#
# Binary call: score >= 0.5 → "trigger NMD", else "escape NMD" (per the
# Lindeboom 2019 convention).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
# Intermediates go to .P$OUT, never the source tree -- writing them beside the scripts is
# how six generated TSVs ended up tracked in a code-only repository.
PC_OUT <- file.path(.P$OUT, "predictor_comparison")
dir.create(PC_OUT, recursive = TRUE, showWarnings = FALSE)

HERE <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grep("--file=", args)]
  if (length(fa) > 0) return(dirname(normalizePath(sub("^--file=", "", fa[1]))))
  normalizePath(getwd())
})()
setwd(HERE)   # scripts still resolve siblings relative to themselves

# 2026.8.6, not 2026.7.11: these outputs are the DEPOSITED 1000nt model rescored.
# Leaving the July stamp would put new data in a file named for the old vintage.
DATESTAMP <- "2026.8.6"
# 2026-07-26: was HERE (the source dir). 01 writes to PC_OUT, so the chain was broken --
# 01 -> PC_OUT, 02/03/04 read HERE, 04 wrote HERE, and SF43 read a legacy code/ path that
# does not exist in this repo. Three locations for one chain. Everything is PC_OUT now.
IN_FILE   <- file.path(PC_OUT, sprintf("isoforms_%s.tsv",     DATESTAMP))
OUT_FILE  <- file.path(PC_OUT, sprintf("nmdetective_b_scores_%s.tsv", DATESTAMP))

# Thresholds — Lindeboom 2019 NMDetective-B (Fig 1c)
TH_START      <- 150L   # start-proximal: distance to start < TH_START
TH_LONG_EXON  <- 407L   # long exon:      exon length > TH_LONG_EXON
TH_50NT       <- 50L    # 50 nt rule:     within TH_50NT of last EJ

# Leaf NMD-efficacy scores
SCORE_LAST_EXON    <- 0.00
SCORE_START_PROX   <- 0.12
SCORE_LONG_EXON    <- 0.41
SCORE_50NT_RULE    <- 0.20
SCORE_TRIGGER      <- 0.65

# Binary threshold
TH_BINARY <- 0.5

# ── Decision-tree scorer ──
score_one <- function(in_last_exon, distance_to_start_nt,
                      exon_length_at_stop_nt, within_50nt_of_last_ej) {
  # NB: `within_50nt_of_last_ej == TRUE` means stop is close to last EJ (≤50 nt),
  # i.e., NOT a PTC by the 50-nt rule → "escape NMD"
  if (is.na(in_last_exon) || is.na(distance_to_start_nt) ||
      is.na(exon_length_at_stop_nt) || is.na(within_50nt_of_last_ej)) {
    return(list(score = NA_real_, leaf = NA_character_,
                call  = NA_character_))
  }
  if (isTRUE(in_last_exon)) {
    return(list(score = SCORE_LAST_EXON, leaf = "last_exon",
                call  = "escape"))
  }
  if (distance_to_start_nt < TH_START) {
    return(list(score = SCORE_START_PROX, leaf = "start_proximal",
                call  = "escape"))
  }
  if (exon_length_at_stop_nt > TH_LONG_EXON) {
    return(list(score = SCORE_LONG_EXON, leaf = "long_exon",
                call  = "escape"))
  }
  if (isTRUE(within_50nt_of_last_ej)) {
    return(list(score = SCORE_50NT_RULE, leaf = "50nt_rule",
                call  = "escape"))
  }
  list(score = SCORE_TRIGGER, leaf = "trigger_nmd", call = "trigger")
}

# ── Apply ──
d <- fread(IN_FILE)
cat(sprintf("Loaded %d isoforms from %s\n", nrow(d), basename(IN_FILE)))

results <- vector("list", nrow(d))
for (i in seq_len(nrow(d))) {
  r <- score_one(d$in_last_exon[i],
                 d$distance_to_start_nt[i],
                 d$exon_length_at_stop_nt[i],
                 d$within_50nt_of_last_ej[i])
  results[[i]] <- data.table(
    comparator_isoform_id = d$comparator_isoform_id[i],
    nmdetective_b_score   = r$score,
    nmdetective_b_leaf    = r$leaf,
    nmdetective_b_call    = r$call
  )
}
out <- rbindlist(results)

fwrite(out, OUT_FILE, sep = "\t")
cat(sprintf("Wrote %s (%d rows)\n", basename(OUT_FILE), nrow(out)))

# Summary
cat("\n=== NMDetective-B leaf assignments by subclass ===\n")
joined <- merge(d[, .(comparator_isoform_id, subclass)], out,
                by = "comparator_isoform_id")
print(joined[, .N, by = .(subclass, nmdetective_b_leaf)][
        order(subclass, -N)])

cat("\n=== Binary trigger/escape calls vs gold-standard subclass ===\n")
print(joined[, .N, by = .(subclass, nmdetective_b_call)])
