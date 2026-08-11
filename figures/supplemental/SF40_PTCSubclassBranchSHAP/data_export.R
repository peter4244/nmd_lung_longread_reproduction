#!/usr/bin/env Rscript
# ==============================================================================
# Data export — per-isoform branch-level KernelSHAP joined with the n=1,166
# Subset 2 group labels, plus per-isoform NMD-prediction probability, then
# stratified summaries by PTC subclass.
#
# Outputs (shared with the sibling PTCSubclassPerformance/ supplement):
#
#   data/branch_shap_by_subclass_refaug.tsv
#       Long form, one row per (isoform, branch). Columns:
#         arm, group, gene_id, reference_isoform_id, comparator_isoform_id,
#         category, branch ∈ {ATG, Stop, Structural}, shap_signed, abs_shap,
#         prediction, label, prob.
#
#   data/branch_shap_by_subclass_refaug_descriptives.tsv
#       Per-group × per-branch mean |SHAP|, mean signed SHAP, n_isoforms.
#
#   data/branch_shap_by_subclass_refaug_pairwise.tsv
#       Wilcoxon per-branch contrast (PTC+ vs PTC−, PTC− vs Control, PTC+ vs Control).
#
#   data/predprob_by_subclass_refaug.tsv  (consumed by PTCSubclassPerformance/)
#       One row per isoform: arm, group, comparator_isoform_id, label, prob.
#
# Sources:
#   ~/claude_projects/NMD_orf_model_v5_4ct/results_4ct/kernel_shap_branch_atg1000_stop1000_seed42_all.tsv
#       isoform_id, label, prediction, expected_value, shap_atg, shap_stop,
#       shap_structural, shap_sum, residual
#   ~/claude_projects/NMD_orf_model_v5_4ct/results_4ct/predictions_atg1000_stop1000_seed42_test_clean.tsv
#       isoform_id, chr, label, logit, prob
#   ~/claude_projects/nmd/figures/multipanel/figure5_dl_model/data/subset2_refaug_isoforms.tsv
#       arm, group, gene_id, reference_isoform_id, comparator_isoform_id, category
# ==============================================================================

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({ library(data.table) })

HERE <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grep("--file=", args)]
  if (length(fa) > 0)
    return(dirname(normalizePath(sub("^--file=", "", fa[1]))))
  normalizePath(getwd())
})()
OUT_DIR <- file.path(HERE, "data")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Shared output also written next to SF41's figure (canonical dir has the
# SF41_ prefix; the unprefixed PTCSubclassPerformance is a rename leftover).
PERF_OUT <- normalizePath(file.path(HERE, "..", "SF41_PTCSubclassPerformance", "data"))
dir.create(PERF_OUT, showWarnings = FALSE, recursive = TRUE)

# REPOINTED 2026-08-06 onto the DEPOSITED 1000nt model (seed 42). This read
# .P$MODEL/results_4ct -- the superseded 500nt vintage -- so SF40 and SF41 would have been rebuilt
# from the model the paper is replacing. No fallback between the two roots: a search that tries one
# then the other silently takes whichever exists.
MODEL_RES  <- file.path(.P$ROOT, "data_deposit", "source_data", "model")
TAG        <- "atg1000_stop1000_seed42"
# Repointed 2026-07-26 (W31 step 2). NOTE: this artifact is MISSING -- figure5_dl_model/
# data/ does not exist and nothing writes it (W32). Pointed in-repo so the failure is honest
# rather than silently satisfied from the legacy tree.
N1166_LIST <- file.path(.P$ROOT, "figures", "multipanel", "figure5_dl_model", "data",
                        "subset2_refaug_isoforms.tsv")

# ── Load ──
# FULL-COHORT KernelSHAP (39,939 isoforms; matches Figure 5 Panel C scope).
# Use `_all.tsv`, NOT the test-set-only `kernel_shap_branch_atg1000_stop1000_seed42_all.tsv`,
# so the n=1,166 join is not test-set-restricted (which dropped PTC- to n=30).
for (f in c(N1166_LIST,
            file.path(MODEL_RES, sprintf("kernel_shap_branch_%s_all.tsv", TAG)),
            file.path(MODEL_RES, sprintf("predictions_%s_test_clean.tsv", TAG))))
  if (!file.exists(f)) stop("missing input: ", f,
                            "\n  SF40/SF41 read the DEPOSITED 1000nt model. Regenerate the refaug",
                            "\n  list with figures/multipanel/figure5_dl_model/data_export_refaug.R",
                            "\n  (NMD_FEATURES must point at a stage C results dir).")

refaug  <- fread(N1166_LIST)
kshap  <- fread(file.path(MODEL_RES, sprintf("kernel_shap_branch_%s_all.tsv", TAG)))
# TEST-SET predictions only, and deliberately: SF41 is a PERFORMANCE figure, and performance is the
# one thing that stays test-only (D-row on analysis scope). Interpretability -- SF40's branch SHAP
# above -- pools train and held-out per D77. The two SFs share a producer and NOT a population,
# which is worth saying out loud in a file whose job is to write both.
preds  <- fread(file.path(MODEL_RES, sprintf("predictions_%s_test_clean.tsv", TAG)),
                select = c("isoform_id", "prob"))

cat(sprintf("refaug: %d rows\nkshap (full cohort): %d rows\npreds (test set): %d rows\n",
            nrow(refaug), nrow(kshap), nrow(preds)))

# ── Branch SHAP join — full cohort scope ──
# Inner-join n=1,166 list with the full-cohort branch SHAP. Do NOT join with
# preds here (preds is test-set only — joining would silently restrict the
# SHAP-side scope to the test split, which is what happened in the first
# build of this SF).
have <- merge(refaug, kshap,
              by.x = "comparator_isoform_id", by.y = "isoform_id",
              all.x = FALSE)
cat(sprintf("\nBranch-SHAP join: %d / %d (%.1f%% coverage; full cohort)\n",
            nrow(have), nrow(refaug), 100 * nrow(have) / nrow(refaug)))

cat("\nBranch-SHAP subgroup counts:\n")
print(have[, .N, by = group])

# ── Per-isoform predicted-prob table — test-set only ──
# Separate inner-join branch — used by the sibling PTCSubclassPerformance/
# supplement (test-set scope per the test-only-for-performance-metrics policy).
have_preds <- merge(refaug, preds,
                    by.x = "comparator_isoform_id", by.y = "isoform_id",
                    all.x = FALSE)
cat(sprintf("\nPredprob join (test set): %d / %d (%.1f%% coverage)\n",
            nrow(have_preds), nrow(refaug),
            100 * nrow(have_preds) / nrow(refaug)))

predprob <- have_preds[, .(arm, group, gene_id, reference_isoform_id,
                            comparator_isoform_id, category, prob)]
fwrite(predprob, file.path(OUT_DIR,  "predprob_by_subclass_refaug.tsv"), sep = "\t")
fwrite(predprob, file.path(PERF_OUT, "predprob_by_subclass_refaug.tsv"), sep = "\t")

# ── Long-form branch SHAP for SF1 (full cohort) ──
long <- melt(
  have,
  id.vars      = c("arm", "group", "gene_id", "reference_isoform_id",
                   "comparator_isoform_id", "category",
                   "prediction", "label"),
  measure.vars = c("shap_atg", "shap_stop", "shap_structural"),
  variable.name = "branch_raw", value.name = "shap_signed"
)
long[, branch := fcase(
  branch_raw == "shap_atg",        "ATG",
  branch_raw == "shap_stop",       "Stop",
  branch_raw == "shap_structural", "Structural"
)]
long[, abs_shap := abs(shap_signed)]
long[, branch_raw := NULL]

fwrite(long, file.path(OUT_DIR, "branch_shap_by_subclass_refaug.tsv"), sep = "\t")

# ── Per-group × per-branch descriptives ──
gord <- c("NMD+/PTC+", "NMD+/PTC- retained", "Control")
desc <- long[, .(
  n_isoforms      = .N,
  mean_abs_shap   = round(mean(abs_shap), 4),
  mean_signed     = round(mean(shap_signed), 4),
  median_abs      = round(median(abs_shap), 4)
), by = .(group, branch)]
desc <- desc[order(factor(group, levels = gord), branch)]
fwrite(desc, file.path(OUT_DIR, "branch_shap_by_subclass_refaug_descriptives.tsv"), sep = "\t")

cat("\n=== Per-group × per-branch mean |SHAP| ===\n")
print(dcast(desc, group ~ branch, value.var = "mean_abs_shap")[
        order(factor(group, levels = gord))])

cat("\n=== START ratio: PTC− vs PTC+ ===\n")
atg_pct_pos  <- desc[group == "NMD+/PTC+",          mean_abs_shap[branch == "ATG"]]
atg_pct_neg  <- desc[group == "NMD+/PTC- retained", mean_abs_shap[branch == "ATG"]]
atg_pct_ctrl <- desc[group == "Control",            mean_abs_shap[branch == "ATG"]]
cat(sprintf("  PTC+ mean |SHAP_ATG|       = %.4f\n", atg_pct_pos))
cat(sprintf("  PTC− mean |SHAP_ATG|       = %.4f\n", atg_pct_neg))
cat(sprintf("  Control mean |SHAP_ATG|    = %.4f\n", atg_pct_ctrl))
cat(sprintf("  Ratio  PTC−/PTC+           = %.2f×\n", atg_pct_neg / atg_pct_pos))

# ── Pairwise Wilcoxon per branch ──
branches <- c("ATG", "Stop", "Structural")
pairs <- list(c("NMD+/PTC+", "NMD+/PTC- retained"),
              c("NMD+/PTC- retained", "Control"),
              c("NMD+/PTC+", "Control"))
pw <- rbindlist(lapply(branches, function(b) {
  rbindlist(lapply(pairs, function(p) {
    x <- long[group == p[1] & branch == b, abs_shap]
    y <- long[group == p[2] & branch == b, abs_shap]
    if (length(x) < 3 || length(y) < 3)
      return(data.table(branch = b, group_x = p[1], group_y = p[2],
                        nx = length(x), ny = length(y), wilcox_p = NA_real_))
    wt <- wilcox.test(x, y, exact = FALSE)
    data.table(branch = b, group_x = p[1], group_y = p[2],
               nx = length(x), ny = length(y),
               wilcox_p = signif(wt$p.value, 3))
  }))
}))
fwrite(pw, file.path(OUT_DIR, "branch_shap_by_subclass_refaug_pairwise.tsv"), sep = "\t")

cat("\n=== Pairwise Wilcoxon (|SHAP|) per branch × contrast ===\n")
print(pw)

cat(sprintf("\nWrote:\n  %s\n  %s\n  %s\n  %s (shared with PTCSubclassPerformance/)\n",
            file.path(OUT_DIR, "branch_shap_by_subclass_refaug.tsv"),
            file.path(OUT_DIR, "branch_shap_by_subclass_refaug_descriptives.tsv"),
            file.path(OUT_DIR, "branch_shap_by_subclass_refaug_pairwise.tsv"),
            file.path(OUT_DIR, "predprob_by_subclass_refaug.tsv")))
