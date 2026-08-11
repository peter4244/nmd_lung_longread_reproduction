#!/usr/bin/env Rscript
# =============================================================================
# 04_compute_metrics.R
#
# Join the three scoring runs (NMDetective-B, NMDEP rule baseline, our model)
# against the per-isoform input table, compute Spearman / Pearson / R² / MAE /
# RMSE against the gold standard (mashr posterior-mean log-FC under SMG1i),
# stratified by subclass.
#
# Outputs:
#   per_isoform_scores_2026.7.11.tsv  — combined per-isoform score table
#   metrics_summary_2026.7.11.tsv     — pooled and stratified metrics
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
IN_FILE     <- file.path(PC_OUT, sprintf("isoforms_%s.tsv",                  DATESTAMP))
NMD_B_FILE  <- file.path(PC_OUT, sprintf("nmdetective_b_scores_%s.tsv",      DATESTAMP))
NMDEP_FILE  <- file.path(PC_OUT, sprintf("nmdep_rule_baseline_scores_%s.tsv",DATESTAMP))
OUT_SCORES  <- file.path(PC_OUT, sprintf("per_isoform_scores_%s.tsv",        DATESTAMP))
OUT_METRICS <- file.path(PC_OUT, sprintf("metrics_summary_%s.tsv",           DATESTAMP))

# ── Load + join ──
d   <- fread(IN_FILE)
nb  <- fread(NMD_B_FILE)
nep <- fread(NMDEP_FILE)

scores <- merge(d[, .(comparator_isoform_id, gene_id, reference_isoform_id,
                       subclass, chr, h5_split, mashr_posterior_mean_logfc,
                       our_model_prob)],
                 nb[, .(comparator_isoform_id, nmdetective_b_score,
                        nmdetective_b_leaf, nmdetective_b_call)],
                 by = "comparator_isoform_id")
scores <- merge(scores,
                 nep[, .(comparator_isoform_id, nmdep_rule_score,
                         nmdep_rule_leaf, nmdep_rule_call)],
                 by = "comparator_isoform_id")
fwrite(scores, OUT_SCORES, sep = "\t")
cat(sprintf("Wrote %s (%d rows)\n", basename(OUT_SCORES), nrow(scores)))

# ── Metric helpers ──
metrics_one <- function(model, gold) {
  ok <- !is.na(model) & !is.na(gold)
  m <- model[ok]; g <- gold[ok]
  n <- length(m)
  if (n < 5) {
    return(data.table(n = n, spearman = NA_real_, pearson = NA_real_,
                       r2 = NA_real_, mae = NA_real_, rmse = NA_real_))
  }
  sp <- suppressWarnings(cor(m, g, method = "spearman"))
  pe <- suppressWarnings(cor(m, g, method = "pearson"))
  data.table(
    n        = n,
    spearman = round(sp, 4),
    pearson  = round(pe, 4),
    r2       = round(pe^2, 4),
    mae      = round(mean(abs(m - g)), 4),
    rmse     = round(sqrt(mean((m - g)^2)), 4)
  )
}

# ── Compute per-model × per-stratum metrics ──
models <- list(
  nmdetective_b = scores$nmdetective_b_score,
  nmdep_rule    = scores$nmdep_rule_score,
  our_model     = scores$our_model_prob
)
strata <- c("pooled", "NMD+/PTC+", "NMD+/PTC-", "Control")

mk_metrics <- function(model_name, model_vec, stratum) {
  if (stratum == "pooled") {
    mask <- rep(TRUE, nrow(scores))
  } else {
    mask <- scores$subclass == stratum
  }
  m <- metrics_one(model_vec[mask], scores$mashr_posterior_mean_logfc[mask])
  cbind(data.table(model = model_name, stratum = stratum), m)
}

results <- rbindlist(unlist(lapply(names(models), function(nm) {
  lapply(strata, function(s) mk_metrics(nm, models[[nm]], s))
}), recursive = FALSE))

# ── Also report a "head-to-head" set: isoforms with all three scores ──
hh_mask <- !is.na(scores$nmdetective_b_score) &
           !is.na(scores$nmdep_rule_score) &
           !is.na(scores$our_model_prob)
cat(sprintf("\nHead-to-head intersection (all three models score the same isoform): %d\n",
            sum(hh_mask)))
hh_results <- rbindlist(lapply(names(models), function(nm) {
  m <- metrics_one(models[[nm]][hh_mask],
                    scores$mashr_posterior_mean_logfc[hh_mask])
  cbind(data.table(model = nm, stratum = "head-to-head"), m)
}))
results <- rbind(results, hh_results)

# Per-subclass within head-to-head set
for (s in c("NMD+/PTC+", "NMD+/PTC-", "Control")) {
  s_mask <- hh_mask & scores$subclass == s
  if (sum(s_mask) < 5) next
  sub <- rbindlist(lapply(names(models), function(nm) {
    m <- metrics_one(models[[nm]][s_mask],
                      scores$mashr_posterior_mean_logfc[s_mask])
    cbind(data.table(model = nm, stratum = paste0("head-to-head:", s)), m)
  }))
  results <- rbind(results, sub)
}

# ── Test-only sensitivity within head-to-head set ──
# Our model has training-data overlap on train/val splits; rerun head-to-head
# pooled + per-subclass metrics restricted to H5 split == "test" (chr1/3/5/7
# holdout, never seen at training time).
test_mask <- hh_mask & !is.na(scores$h5_split) & scores$h5_split == "test"
cat(sprintf("\nTest-only head-to-head (H5 split == 'test'): %d\n", sum(test_mask)))
test_pooled <- rbindlist(lapply(names(models), function(nm) {
  m <- metrics_one(models[[nm]][test_mask],
                    scores$mashr_posterior_mean_logfc[test_mask])
  cbind(data.table(model = nm, stratum = "head-to-head:test"), m)
}))
results <- rbind(results, test_pooled)
for (s in c("NMD+/PTC+", "NMD+/PTC-", "Control")) {
  ts_mask <- test_mask & scores$subclass == s
  if (sum(ts_mask) < 5) next
  sub <- rbindlist(lapply(names(models), function(nm) {
    m <- metrics_one(models[[nm]][ts_mask],
                      scores$mashr_posterior_mean_logfc[ts_mask])
    cbind(data.table(model = nm, stratum = paste0("head-to-head:test:", s)), m)
  }))
  results <- rbind(results, sub)
}

fwrite(results, OUT_METRICS, sep = "\t")
cat(sprintf("\nWrote %s\n", basename(OUT_METRICS)))

cat("\n=== Metrics (model × stratum) ===\n")
print(results)
