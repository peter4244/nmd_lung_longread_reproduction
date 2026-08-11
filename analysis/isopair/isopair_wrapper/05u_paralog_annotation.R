#!/usr/bin/env Rscript
# ==============================================================================
# 05u_paralog_annotation.R
#
# Query Ensembl for human gene paralogs via biomaRt.
#
# For predicting NMD from isoform structure, genes with high-similarity
# expressed paralogs pose a data leakage risk when one gene is in the training
# set (non-holdout chromosomes) and its paralog is in the holdout test set.
# Structurally similar isoforms could inflate test performance.
#
# We flag genes that meet ALL three criteria:
#   1. Paralog protein sequence identity >= 80%
#   2. Both the gene and its paralog are expressed in our dataset
#   3. The gene and its paralog are on DIFFERENT sides of the train/test split
#      (one on a holdout chromosome, one on a training chromosome)
#
# Only these cross-split high-similarity pairs represent actual leakage risk.
#
# Inputs:
#   - data_mashr/structures.rds (gene IDs + chromosomes)
#   - data_mashr/nmd_classification.rds (to identify expressed genes)
#
# Output:
#   - data_mashr/analysis_cache/paralog_genes.rds
#     $leakage_genes: versioned gene_ids to remove from test set
#     $leakage_pairs: the cross-split high-similarity pairs (for reporting)
#     $all_expressed_pairs: all expressed paralog pairs (before filtering)
#     $paralog_pairs: full biomaRt result with percent identity
#     $metadata: timestamp, Ensembl version, thresholds, counts
#
# Usage:
#   Rscript 05u_paralog_annotation.R
# ==============================================================================

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
})

# REPOINTED 2026-07-27, same fix and same reason as 05s_orfik_scan.R and 05t.
#
# This did setwd(.P$ISOPAIR) and read "data_mashr/..." RELATIVE to it -- the LEGACY tree.
# paralog_genes.rds feeds export_rds.R -> data_prep.py -> the model, so a legacy read here
# defines the model's paralog feature on the wrong isoform universe.
#
# No fallback to the legacy tree: it is always populated, so a fallback could only ever
# silently serve the vintage the repoint exists to stop using.
DATA_MASHR <- .P$ISOPAIR_OUT
if (!dir.exists(DATA_MASHR))
  stop("intermediates root not found: ", DATA_MASHR,
       "\n  This script reads the DEPOSIT-REBUILT tree and does NOT fall back to the legacy one.")
cat("  data_mashr root:", DATA_MASHR, "\n")

OUTPUT_PATH <- file.path(.P$OUT, "data_mashr/analysis_cache/paralog_genes.rds")
# The parent must exist before the write. 05s_orfik_scan.R lost 75 MINUTES of completed,
# fully-verified work to its absence on 2026-07-28: it scanned 1,540,674 ORFs, passed every
# self-check, then died on saveRDS with "cannot open compressed file ... No such file or
# directory" because tmp/out/data_mashr/analysis_cache/ had never been created. Seven
# necessary files write into that directory and the only one that CREATES it, 03b_rebuild_
# cache.R, runs at index 42 -- after five of them. Each writer now makes its own parent, so
# the guarantee is order-independent rather than an accident of run_order. W81.
dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)

# Thresholds
MIN_IDENTITY <- 80     # minimum % protein sequence identity
HOLDOUT_CHRS <- c("chr1", "chr3", "chr5", "chr7")

# THE VALIDATION SPLIT IS SCREENED TOO (2026-07-29, D35). The screen below was built for the
# train/test split only, so `filter(gene_in_holdout != paralog_in_holdout)` said nothing about
# chr2/chr4 -- no validation gene could enter leakage_genes, and val_paralog came out 0 by
# construction rather than because validation is clean. That was harmless while validation only
# drove early stopping; D31 makes it the SELECTION set for the window sweep, at which point a
# paralog shared with train or test biases the choice of configuration and flatters the final
# test estimate. Must match data_prep.py:48.
VAL_CHRS <- c("chr2", "chr4")

# --from-cache reuses paralog_pairs from an existing paralog_genes.rds instead of querying
# Ensembl. C59: the query result is ALREADY CACHED in that object, so extending the screen is a
# re-filter rather than a network operation. The filters are recomputed from scratch either way,
# which is what lets the run below prove it reproduces the existing leakage_genes exactly.
FROM_CACHE <- "--from-cache" %in% commandArgs(trailingOnly = TRUE)

cat("=== Paralog Annotation ===\n")
cat("Started:", format(Sys.time()), "\n\n")
cat("Thresholds:\n")
cat("  Min protein identity:", MIN_IDENTITY, "%\n")
cat("  Holdout chromosomes:", paste(HOLDOUT_CHRS, collapse = ", "), "\n\n")

# ==============================================================================
# 1. Extract gene IDs and chromosome assignments
# ==============================================================================
cat("Loading gene IDs...\n")

structures <- readRDS(file.path(DATA_MASHR, "structures.rds"))
nmd_class <- readRDS(file.path(DATA_MASHR, "nmd_classification.rds"))

# All classified genes (NMD + non-NMD from all_samples)
classified_isoforms <- union(nmd_class$all_samples$nmd, nmd_class$all_samples$non_nmd)
classified_structures <- structures[structures$isoform_id %in% classified_isoforms, ]

# Gene → chromosome mapping (one gene can only be on one chromosome)
gene_chr <- classified_structures %>%
  distinct(gene_id, chr) %>%
  mutate(in_holdout = chr %in% HOLDOUT_CHRS,
         in_val     = chr %in% VAL_CHRS)

classified_genes_versioned <- unique(gene_chr$gene_id)

# Strip version for biomaRt query
gene_id_map <- data.frame(
  gene_id_versioned = classified_genes_versioned,
  gene_id_unversioned = sub("\\.[0-9]+$", "", classified_genes_versioned),
  stringsAsFactors = FALSE
)

cat("  Classified genes:", length(classified_genes_versioned), "\n")
cat("  On holdout chromosomes:", sum(gene_chr$in_holdout), "\n")
cat("  On training chromosomes:", sum(!gene_chr$in_holdout), "\n")

# ==============================================================================
# 2. Query Ensembl for paralogs WITH percent identity
# ==============================================================================
# THE DEPOSITED TABLE IS THE DEFAULT SOURCE (2026-08-03). This step used to query biomaRt live,
# which made the entire chain require internet -- and this table defines the paralog screen behind
# the train/test split, so a reader who cannot reach Ensembl loses the LEAKAGE SCREEN, not a
# figure. Homology is not derivable from the GTFs in the record (a GTF says where exons are, never
# that gene A is a paralog of gene B at 82% identity), so the answer is deposited rather than
# recomputed. --live-query re-fetches, which is how the deposited copy is regenerated.
#
# NOTE this REPLACES --from-cache as the offline route. That flag read paralog_pairs back out of a
# previous OUTPUT of this script, which a fresh reader does not have -- offline it could only work
# for someone who had already run it online.
LIVE_QUERY  <- "--live-query" %in% commandArgs(trailingOnly = TRUE)
DEPOSITED   <- file.path(.P$ANNOT, "ensembl115_human_paralogs.tsv.gz")

if (!LIVE_QUERY && file.exists(DEPOSITED)) {
  cat("\n=== reading deposited paralogs, NOT querying Ensembl ===\n")
  paralog_pairs <- read.delim(gzfile(DEPOSITED), stringsAsFactors = FALSE)
  ens_version   <- paste0("deposited:", basename(DEPOSITED))
  stopifnot("deposited paralog table has unexpected columns" =
              all(c("ensembl_gene_id", "hsapiens_paralog_ensembl_gene",
                    "hsapiens_paralog_perc_id") %in% names(paralog_pairs)))
  cat("  from:", DEPOSITED, "\n")
  cat("  paralog_pairs:", nrow(paralog_pairs), "rows (Ensembl v115, deposited)\n")
} else if (FROM_CACHE) {
  stopifnot("--from-cache needs an existing paralog_genes.rds" = file.exists(OUTPUT_PATH))
  cat("\n=== --from-cache: reusing paralog_pairs, NOT querying Ensembl ===\n")
  .cached <- readRDS(OUTPUT_PATH)
  paralog_pairs <- .cached$paralog_pairs
  ens_version <- .cached$metadata$ensembl_host
  cat("  from:", OUTPUT_PATH, "\n")
  cat("  paralog_pairs:", nrow(paralog_pairs), "rows (cached, Ensembl v115)\n")
} else {
cat("\nConnecting to Ensembl...\n")

# Pinned to Ensembl v115 to match the deposited annotation maps (gmap_*_ENSGv115_*.rds).
# Unpinned, this silently tracks whatever Ensembl is current. Still a LIVE network call.
ensembl <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl", version = 115)
ens_version <- ensembl@host

cat("  Ensembl host:", ens_version, "\n")
cat("  Querying paralogs for", length(unique(gene_id_map$gene_id_unversioned)), "genes...\n")

# Query in batches to avoid timeout
batch_size <- 2000L
uniq_ids <- unique(gene_id_map$gene_id_unversioned)
n_batches <- ceiling(length(uniq_ids) / batch_size)

paralog_list <- vector("list", n_batches)
t0 <- proc.time()

for (b in seq_len(n_batches)) {
  start_idx <- (b - 1L) * batch_size + 1L
  end_idx <- min(b * batch_size, length(uniq_ids))
  batch_ids <- uniq_ids[start_idx:end_idx]

  # Retry up to 3 times on transient biomaRt failures
  for (attempt in 1:3) {
    paralog_list[[b]] <- tryCatch(
      getBM(
        attributes = c("ensembl_gene_id",
                        "hsapiens_paralog_ensembl_gene",
                        "hsapiens_paralog_perc_id"),
        filters = "ensembl_gene_id",
        values = batch_ids,
        mart = ensembl
      ),
      error = function(e) {
        if (attempt < 3) {
          cat(sprintf("  Batch %d attempt %d failed: %s. Retrying in 10s...\n",
                      b, attempt, conditionMessage(e)))
          Sys.sleep(10)
          NULL
        } else {
          stop(sprintf("biomaRt query failed after 3 attempts (batch %d): %s",
                       b, conditionMessage(e)))
        }
      }
    )
    if (!is.null(paralog_list[[b]])) break
  }

  cat(sprintf("  Batch %d/%d: %d genes, %d paralog rows (%.0f sec)\n",
              b, n_batches, length(batch_ids), nrow(paralog_list[[b]]),
              (proc.time() - t0)[3]))
}

paralog_pairs <- do.call(rbind, paralog_list)
stopifnot("biomaRt returned no paralog pairs" = nrow(paralog_pairs) > 0)

# Remove empty paralog entries
paralog_pairs <- paralog_pairs[paralog_pairs$hsapiens_paralog_ensembl_gene != "", ]
}   # end of the live-query branch (see FROM_CACHE above)

cat("\n  Total paralog pairs:", nrow(paralog_pairs), "\n")
cat("  Genes with >= 1 paralog:", length(unique(paralog_pairs$ensembl_gene_id)), "\n")

# Percent identity distribution
cat("\n  Percent identity distribution:\n")
print(summary(paralog_pairs$hsapiens_paralog_perc_id))

# ==============================================================================
# 3. Filter to high-similarity expressed cross-split pairs
# ==============================================================================
cat("\n=== Filtering paralog pairs ===\n")

expressed_unversioned <- unique(gene_id_map$gene_id_unversioned)

# Step 1: Both genes expressed in our dataset
expressed_pairs <- paralog_pairs %>%
  filter(
    ensembl_gene_id %in% expressed_unversioned,
    hsapiens_paralog_ensembl_gene %in% expressed_unversioned
  )
cat("  Step 1 — Both expressed:", nrow(expressed_pairs), "pairs,",
    length(unique(expressed_pairs$ensembl_gene_id)), "genes\n")

# Step 2: High protein sequence identity
high_sim_pairs <- expressed_pairs %>%
  filter(hsapiens_paralog_perc_id >= MIN_IDENTITY)
cat("  Step 2 — Identity >=", MIN_IDENTITY, "%:", nrow(high_sim_pairs), "pairs,",
    length(unique(high_sim_pairs$ensembl_gene_id)), "genes\n")

# Step 3: Cross-split (one in holdout, one in training)
# Map unversioned IDs to chromosomes via versioned → gene_chr
gene_chr_unver <- gene_chr %>%
  mutate(gene_id_unversioned = sub("\\.[0-9]+$", "", gene_id)) %>%
  distinct(gene_id_unversioned, in_holdout)

cross_split_pairs <- high_sim_pairs %>%
  inner_join(gene_chr_unver, by = c("ensembl_gene_id" = "gene_id_unversioned")) %>%
  rename(gene_in_holdout = in_holdout) %>%
  inner_join(gene_chr_unver, by = c("hsapiens_paralog_ensembl_gene" = "gene_id_unversioned")) %>%
  rename(paralog_in_holdout = in_holdout) %>%
  filter(gene_in_holdout != paralog_in_holdout)

cat("  Step 3 — Cross-split:", nrow(cross_split_pairs), "pairs,",
    length(unique(cross_split_pairs$ensembl_gene_id)), "genes\n")

# Genes to remove from test set: holdout-side genes in cross-split pairs
# (Remove the holdout gene, not the training gene)
leakage_genes_unversioned <- unique(c(
  cross_split_pairs$ensembl_gene_id[cross_split_pairs$gene_in_holdout],
  cross_split_pairs$hsapiens_paralog_ensembl_gene[cross_split_pairs$paralog_in_holdout]
))

# Map back to versioned IDs
leakage_genes_versioned <- gene_id_map$gene_id_versioned[
  gene_id_map$gene_id_unversioned %in% leakage_genes_unversioned
]
# Keep only those actually on holdout chromosomes
leakage_genes_versioned <- leakage_genes_versioned[
  leakage_genes_versioned %in% gene_chr$gene_id[gene_chr$in_holdout]
]

cat("\n  Leakage genes to remove from test set:", length(leakage_genes_versioned), "\n")

# ==============================================================================
# 3b. The SAME screen, for the VALIDATION split (2026-07-29, D35)
# ==============================================================================
# Deliberately a SECOND, SEPARATE screen rather than a widening of the one above. Everything
# feeding leakage_genes is untouched, so test_paralog = 122 and claim 5.6.4 stand exactly as
# baselined (provenance/PARALOG_RDS_BASELINE.md). Widening the existing screen would have moved
# a published number for no scientific gain.
#
# A validation gene leaks if it has a >=80%-identity expressed paralog OUTSIDE the validation
# chromosomes -- i.e. in train OR test, which is Pete's ruling under D35. Both matter and they
# protect different things: a paralog in train biases the CONFIGURATION CHOICE; a paralog in
# test means a config chosen on validation is also flattered on test, so the FINAL estimate is
# optimistic. `!in_val` covers both at once, since anything not on chr2/chr4 is one or the other.
#
# Note the existing screen already covers the reverse direction: "not in holdout" includes
# chr2/chr4, so a test gene whose paralog sits in validation is already caught above.
cat("\n=== Filtering paralog pairs (VALIDATION split) ===\n")

gene_val_unver <- gene_chr %>%
  mutate(gene_id_unversioned = sub("\\.[0-9]+$", "", gene_id)) %>%
  distinct(gene_id_unversioned, in_val)

cross_val_pairs <- high_sim_pairs %>%
  inner_join(gene_val_unver, by = c("ensembl_gene_id" = "gene_id_unversioned")) %>%
  rename(gene_in_val = in_val) %>%
  inner_join(gene_val_unver, by = c("hsapiens_paralog_ensembl_gene" = "gene_id_unversioned")) %>%
  rename(paralog_in_val = in_val) %>%
  filter(gene_in_val != paralog_in_val)

cat("  Cross-validation-boundary pairs:", nrow(cross_val_pairs), "\n")

val_leakage_unversioned <- unique(c(
  cross_val_pairs$ensembl_gene_id[cross_val_pairs$gene_in_val],
  cross_val_pairs$hsapiens_paralog_ensembl_gene[cross_val_pairs$paralog_in_val]
))

val_leakage_genes_versioned <- gene_id_map$gene_id_versioned[
  gene_id_map$gene_id_unversioned %in% val_leakage_unversioned
]
# Keep only genes actually ON validation chromosomes -- the mirror of the holdout step.
val_leakage_genes_versioned <- val_leakage_genes_versioned[
  val_leakage_genes_versioned %in% gene_chr$gene_id[gene_chr$in_val]
]

cat("  Leakage genes to remove from the VALIDATION set:",
    length(val_leakage_genes_versioned), "\n")
cat("  (validation genes total:", sum(gene_chr$in_val), ")\n")

# ==============================================================================
# 4. Summary
# ==============================================================================
cat("\n=== SUMMARY ===\n")
cat("  Total classified genes:", length(classified_genes_versioned), "\n")
cat("  Genes with any paralog:", length(unique(paralog_pairs$ensembl_gene_id)), "\n")
cat("  High-similarity (>=", MIN_IDENTITY, "%) expressed cross-split pairs:",
    nrow(cross_split_pairs), "\n")
cat("  Holdout genes to remove (leakage risk):", length(leakage_genes_versioned),
    sprintf("(%.1f%% of holdout genes)\n",
            100 * length(leakage_genes_versioned) / sum(gene_chr$in_holdout)))

# Identity distribution of the leakage pairs
if (nrow(cross_split_pairs) > 0) {
  cat("\n  Identity distribution of leakage pairs:\n")
  print(summary(cross_split_pairs$hsapiens_paralog_perc_id))
}

# ==============================================================================
# 5. Save
# ==============================================================================
results <- list(
  # Primary output: genes to remove from test set
  leakage_genes = leakage_genes_versioned,

  # For reporting
  leakage_pairs = cross_split_pairs,
  all_expressed_pairs = expressed_pairs,
  paralog_pairs = paralog_pairs,

  # ADDED 2026-07-29 (D35). New elements only; every pre-existing element above is unchanged.
  val_leakage_genes = val_leakage_genes_versioned,
  val_leakage_pairs = cross_val_pairs,

  metadata = list(
    run_timestamp = Sys.time(),
    ensembl_host = ens_version,
    min_identity = MIN_IDENTITY,
    holdout_chrs = HOLDOUT_CHRS,
    val_chrs = VAL_CHRS,
    n_classified_genes = length(classified_genes_versioned),
    n_holdout_genes = sum(gene_chr$in_holdout),
    n_leakage_genes = length(leakage_genes_versioned),
    n_all_paralog_pairs = nrow(paralog_pairs),
    n_expressed_pairs = nrow(expressed_pairs),
    n_high_sim_pairs = nrow(high_sim_pairs),
    n_cross_split_pairs = nrow(cross_split_pairs),
    n_val_genes = sum(gene_chr$in_val),
    n_val_leakage_genes = length(val_leakage_genes_versioned),
    n_cross_val_pairs = nrow(cross_val_pairs)
  )
)

saveRDS(results, OUTPUT_PATH)
cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Completed:", format(Sys.time()), "\n")
