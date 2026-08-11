#!/usr/bin/env Rscript
#
# NMD Isoform Transitions — Rebuild Analysis Cache
#
# Purpose: Regenerate all analysis_cache/ objects from gene-matched profiles.
#          Standalone script that replicates the caching logic from
#          03_nmd_analysis_mashr.Rmd without rendering the full exploratory report.
#
# Prerequisite: 01_prepare_data_mashr.R and 02_build_profiles_mashr.R must have
#               been run first to produce data_mashr/*.rds files.
#
# Output: data_mashr/analysis_cache/{div,fc,fw,ptc,ri,er,cooc,cps}_*.rds
#
# Usage:
#   Rscript 03b_rebuild_cache.R [--force]

library(Isopair)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args

# --- data_dir: config-resolved (2026-07-26) -------------------------------------------
# This was a bare relative "data_mashr", with NO config bootstrap anywhere in the file, so
# the directory it read AND WROTE depended entirely on the caller's working directory. Run
# from the legacy wrapper dir -- the natural place to run it -- it writes analysis_cache/
# INTO AN INPUT ROOT, which is the one thing tools/check_output_roots.R exists to prevent.
#
# Default is .P$OUT/data_mashr. 01b and 02 write the same filenames there and 02 overwrites
# 01b; only 02's pairs are floored, and only the floored ones are the published quantity, so
# pass --data-dir explicitly rather than trusting the default.
.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)

data_dir <- file.path(.P$OUT, "data_mashr")
if ("--data-dir" %in% args) {
  .i <- which(args == "--data-dir")
  if (.i < length(args)) data_dir <- args[.i + 1]
}
data_dir <- normalizePath(data_dir, mustWork = FALSE)
if (!dir.exists(data_dir)) {
  stop(sprintf("data_dir does not exist: %s\n  Run 01_prepare_data_mashr.R -> 01b -> 02 first, or pass --data-dir.", data_dir))
}
cat("data dir: ", data_dir, "\n", sep = "")
# Rule in code: this script WRITES here, so it must not be an input root. A --data-dir
# pointing at a reference tree would overwrite the copies verification compares against.
nmd_assert_writable(data_dir, .P, "--data-dir")

cache_dir <- file.path(data_dir, "analysis_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

ct_to_label <- function(ct) tolower(gsub("_", "", ct))

cached_compute <- function(key, expr, force_recompute = force) {
  cache_file <- file.path(cache_dir, paste0(key, ".rds"))
  if (!force_recompute && file.exists(cache_file)) {
    message(sprintf("  [cache hit] %s", key))
    return(readRDS(cache_file))
  }
  message(sprintf("  [computing] %s ...", key))
  result <- eval(expr)
  # A FAILED COMPUTATION MUST NOT BE CACHED. Every call site below wraps its expression in
  # tryCatch(..., error = function(e) { message(e$message); NULL }), so a throw arrives here
  # as NULL with its reason on stderr and nothing else. Saving that NULL made the failure
  # PERMANENT: the guard above is existence-only, so every later run found the file, reported
  # [cache hit] and returned the NULL, and the computation was never retried even after the
  # cause was fixed. Exit status was 0 throughout and the consumers read a NULL where they
  # expect a data frame. P5 -- a guard whose failure branch writes a plausible-looking
  # artifact is not a guard. Measured before this changed: 83 cache files, 0 NULLs, so this
  # has not yet fired; it is the re-run loop that makes it reachable.
  if (is.null(result))
    stop(sprintf(paste0("cached_compute: '%s' computed to NULL and was NOT cached.\n",
                        "  The underlying error is on stderr above -- the tryCatch at the call ",
                        "site turned it into a message.\n",
                        "  Fix the cause and re-run; nothing was written, so there is no stale ",
                        "cache entry to clear."), key))
  saveRDS(result, cache_file)
  result
}

cat("=== Rebuilding Analysis Cache ===\n\n")

# ==============================================================================
# 1. Load data
# ==============================================================================

cat("Loading data...\n")
structures   <- readRDS(file.path(data_dir, "structures.rds"))
cds          <- readRDS(file.path(data_dir, "cds.rds"))
ptc          <- readRDS(file.path(data_dir, "ptc.rds"))
union_exons  <- readRDS(file.path(data_dir, "union_exons.rds"))
ue_mapping   <- readRDS(file.path(data_dir, "isoform_union_mapping.rds"))
annotated_ue <- readRDS(file.path(data_dir, "annotated_ue.rds"))
cds_exons    <- readRDS(file.path(data_dir, "cds_exons.rds"))

cell_types <- c("all_samples", "AT", "DD", "FB", "MV")  # paper scope (2026-04-29)
comparisons <- c("c2", "c4")

# ==============================================================================
# 2. Load profiles and apply gene-matching
# ==============================================================================

cat("Loading profiles...\n")
# ABSENCE MUST NOT BE SILENT -- AND MUST NOT BE OVER-ASSERTED EITHER.
#
# This loop had no else branch, and all three consumers below skip quietly: the gene-matching
# loop and the comparePairSets loop `next` on a missing key, and the compute loop iterates
# names(profiles_list) rather than cell_types, so a short list cannot be noticed by the thing
# that reads it. One profile 02 failed to write therefore cost 7 cache objects plus a cps_
# entry, printed nothing about it, and exited 0.
#
# A bare count assertion would be wrong. 02 writes each profile inside `if (!is.null(pairs_cN))`
# (02_build_profiles_mashr.R:326,:335), so a cell type with no pairs for a comparison
# legitimately has no profile, and a guard that fires on a valid tree is a guard that gets
# disabled. 02 writes pairs_<key>.rds BEFORE the profile in the same branch, which separates
# the two cases exactly:
#     pairs present, profile absent -> 02 started that key and did not finish -> STOP
#     both absent                   -> no pairs for that key at all           -> report, continue
profiles_list <- list()
absent_with_pairs <- character()
absent_no_pairs   <- character()
for (ct in cell_types) {
  ct_label <- ct_to_label(ct)
  for (comp in comparisons) {
    key <- paste(comp, ct_label, sep = "_")
    prof_file  <- file.path(data_dir, sprintf("profiles_%s.rds", key))
    pairs_file <- file.path(data_dir, sprintf("pairs_%s.rds", key))
    if (file.exists(prof_file)) {
      profiles_list[[key]] <- readRDS(prof_file)
    } else if (file.exists(pairs_file)) {
      absent_with_pairs <- c(absent_with_pairs, key)
    } else {
      absent_no_pairs <- c(absent_no_pairs, key)
    }
  }
}
if (length(absent_with_pairs))
  stop(sprintf(paste0("%d profile(s) missing from %s while their pairs_*.rds ARE present, ",
                      "so 02 started these keys and did not finish:\n  %s\n",
                      "  Re-run 02_build_profiles_mashr.R, or point --data-dir at the ",
                      "directory it actually wrote."),
               length(absent_with_pairs), data_dir,
               paste(absent_with_pairs, collapse = "\n  ")))
cat(sprintf("  %d of %d profiles loaded%s\n",
            length(profiles_list), length(cell_types) * length(comparisons),
            if (length(absent_no_pairs))
              sprintf("; no pairs at all for %s", paste(absent_no_pairs, collapse = ", "))
            else ""))

cat("Applying gene-matching...\n")
for (ct in cell_types) {
  ct_label <- ct_to_label(ct)
  c2_key <- paste0("c2_", ct_label)
  c4_key <- paste0("c4_", ct_label)
  if (!c2_key %in% names(profiles_list) || !c4_key %in% names(profiles_list)) next

  n_c2_before <- nrow(profiles_list[[c2_key]])
  n_c4_before <- nrow(profiles_list[[c4_key]])

  c2_gene_refs <- unique(profiles_list[[c2_key]][, c("gene_id", "reference_isoform_id")])
  c4_gene_refs <- unique(profiles_list[[c4_key]][, c("gene_id", "reference_isoform_id")])
  shared_refs  <- inner_join(c2_gene_refs, c4_gene_refs,
                              by = c("gene_id", "reference_isoform_id"))

  profiles_list[[c2_key]] <- semi_join(profiles_list[[c2_key]], shared_refs,
                                        by = c("gene_id", "reference_isoform_id"))
  profiles_list[[c4_key]] <- semi_join(profiles_list[[c4_key]], shared_refs,
                                        by = c("gene_id", "reference_isoform_id"))

  cat(sprintf("  %s: C2 %d→%d, C4 %d→%d (matched genes: %d)\n",
    ct, n_c2_before, nrow(profiles_list[[c2_key]]),
    n_c4_before, nrow(profiles_list[[c4_key]]), nrow(shared_refs)))
}

# ==============================================================================
# 3. Compute and cache all analysis objects
# ==============================================================================

cat("\nComputing analyses...\n")
for (key in names(profiles_list)) {
  profiles <- profiles_list[[key]]
  cat(sprintf("\n--- %s (%d pairs) ---\n", key, nrow(profiles)))

  cached_compute(paste0("div_", key),
    tryCatch(quantifyPairDivergence(profiles, structures, cds), error = function(e) { message(e$message); NULL }))
  cached_compute(paste0("ri_", key),
    tryCatch(quantifyRegionalImpact(profiles, annotated_ue, cds), error = function(e) { message(e$message); NULL }))
  cached_compute(paste0("cooc_", key),
    tryCatch(testCooccurrence(profiles), error = function(e) { message(e$message); NULL }))
  cached_compute(paste0("er_", key),
    tryCatch(mapEventsToRegions(profiles, annotated_ue, cds), error = function(e) { message(e$message); NULL }))
  cached_compute(paste0("ptc_", key),
    tryCatch(testPtcAssociation(ptc, profiles), error = function(e) { message(e$message); NULL }))
  cached_compute(paste0("fw_", key),
    tryCatch(analyzeFrameWalk(profiles, cds, ptc), error = function(e) { message(e$message); NULL }))
  cached_compute(paste0("fc_", key),
    tryCatch(compareIsoformFrames(profiles, cds_exons, cds), error = function(e) { message(e$message); NULL }))
}

# ==============================================================================
# 4. comparePairSets per cell type
# ==============================================================================

cat("\nComputing NMD vs Control cross-comparisons...\n")
for (ct in cell_types) {
  ct_label <- ct_to_label(ct)
  c2_key <- paste0("c2_", ct_label)
  c4_key <- paste0("c4_", ct_label)
  if (!c2_key %in% names(profiles_list) || !c4_key %in% names(profiles_list)) next

  ptc_c2 <- ptc[ptc$isoform_id %in% profiles_list[[c2_key]]$comparator_isoform_id, ]
  ptc_c4 <- ptc[ptc$isoform_id %in% profiles_list[[c4_key]]$comparator_isoform_id, ]

  cached_compute(paste0("cps_", ct_label),
    tryCatch(
      comparePairSets(profiles_list[[c2_key]], profiles_list[[c4_key]],
                       ptc_a = ptc_c2, ptc_b = ptc_c4,
                       paired = TRUE, verbose = FALSE),
      error = function(e) { message(e$message); NULL }))
}

cat("\n=== Cache rebuild complete ===\n")
