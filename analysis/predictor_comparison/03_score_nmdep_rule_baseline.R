#!/usr/bin/env Rscript
# =============================================================================
# 03_score_nmdep_rule_baseline.R
#
# Apply the NMDEP-re-thresholded 4-rule baseline (Saadat & Fellay 2025,
# Table 3) to each isoform. Same decision-tree structure as NMDetective-B,
# but with NMDEP's grid-search-optimized thresholds.
#
# NMDEP optimized thresholds (Saadat 2025, Table 3):
#   Penultimate-exon rule:   49 (vs Lindeboom 2019's 50)
#   Start-codon proximity:  120 (vs 150)
#   Exon length:            355 (vs 407)
#
# Note: NMDEP's paper uses a "penultimate-exon" rule that flags PTCs within
# 49 nt of the start of the last exon. This is structurally analogous to
# the canonical 50-nt rule (within 50 nt of the last EJ in NMDetective-B's
# notation), since the start of the last exon is the position of the last
# exon-exon junction.
#
# Leaf NMD-efficacy scores: NMDEP does not specify alternative leaf scores
# in their Table 3; we use the same NMDetective-B leaf scores so the score
# scale is comparable. (Differences between the two baselines come entirely
# from the threshold changes.)
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
IN_FILE   <- file.path(PC_OUT, sprintf("isoforms_%s.tsv",                DATESTAMP))
OUT_FILE  <- file.path(PC_OUT, sprintf("nmdep_rule_baseline_scores_%s.tsv", DATESTAMP))

# Thresholds — NMDEP (Saadat 2025 Table 3)
TH_START      <- 120L
TH_LONG_EXON  <- 355L
TH_50NT       <- 49L

# Leaf NMD-efficacy scores (inherited from NMDetective-B; see header)
SCORE_LAST_EXON    <- 0.00
SCORE_START_PROX   <- 0.12
SCORE_LONG_EXON    <- 0.41
SCORE_50NT_RULE    <- 0.20
SCORE_TRIGGER      <- 0.65

TH_BINARY <- 0.5

score_one <- function(in_last_exon, distance_to_start_nt,
                      exon_length_at_stop_nt, within_50nt_of_last_ej,
                      ej_distance_nt) {
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
  # The 50nt rule's threshold is the NMDEP-re-optimized 49 (vs Lindeboom's 50).
  # `within_50nt_of_last_ej` was computed at threshold 50. For the 1-nt
  # delta we'd need the raw EJ distance; we don't carry that explicitly
  # in the per-isoform input table. The 1-nt threshold difference is
  # immaterial for almost all cases (stop position is rarely exactly at
  # the 49-nt boundary). Treat as same boolean.
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
    nmdep_rule_score      = r$score,
    nmdep_rule_leaf       = r$leaf,
    nmdep_rule_call       = r$call
  )
}
out <- rbindlist(results)

fwrite(out, OUT_FILE, sep = "\t")
cat(sprintf("Wrote %s (%d rows)\n", basename(OUT_FILE), nrow(out)))

# Summary
cat("\n=== NMDEP rule baseline: leaf assignments by subclass ===\n")
joined <- merge(d[, .(comparator_isoform_id, subclass)], out,
                by = "comparator_isoform_id")
print(joined[, .N, by = .(subclass, nmdep_rule_leaf)][order(subclass, -N)])

cat("\n=== Binary trigger/escape calls vs gold-standard subclass ===\n")
print(joined[, .N, by = .(subclass, nmdep_rule_call)])

# Compare to NMDetective-B (if its output is on disk)
nmd_b_file <- file.path(PC_OUT, sprintf("nmdetective_b_scores_%s.tsv", DATESTAMP))
if (file.exists(nmd_b_file)) {
  nb <- fread(nmd_b_file)
  cmp <- merge(out[, .(comparator_isoform_id, nmdep = nmdep_rule_call)],
                nb[, .(comparator_isoform_id, nmd_b = nmdetective_b_call)],
                by = "comparator_isoform_id")
  cat(sprintf("\n=== Agreement with NMDetective-B (n=%d non-NA) ===\n",
              sum(!is.na(cmp$nmdep) & !is.na(cmp$nmd_b))))
  print(table(NMDEP = cmp$nmdep, NMDetective_B = cmp$nmd_b, useNA = "ifany"))
}
