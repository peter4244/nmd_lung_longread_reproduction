#!/usr/bin/env Rscript
# ==============================================================================
# 05k_utr5_all_isoforms.R
#
# Compute 5'UTR features for ALL classified coding isoforms (~62K), not just
# the paired subset (~10K). This provides the feature set needed for the
# unified predictive model.
#
# Input files:
#   - data_mashr/structures.rds
#   - data_mashr/cds.rds
#   - data_mashr/nmd_classification.rds
#   - FASTA: sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta
#
# Output:
#   - data_mashr/analysis_cache/utr5_features_all.rds
#
# Usage:
#   Rscript 05k_utr5_all_isoforms.R
# ==============================================================================

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(Isopair)
  library(dplyr)
})

# setwd(.P$ISOPAIR) REMOVED 2026-07-28. It made this script UNRUNNABLE for a reader: .P$ISOPAIR
# is the legacy authoring tree, which does not exist outside this machine, and setwd() to a
# missing directory is a hard error. The clean-room rerun died here in 0.9s with
# "cannot change working directory".
#
# THE COMMENT IT REPLACES CLAIMED A REVIEW THAT HAD ASKED THE WRONG QUESTION. It said "the only
# remaining relative read is the helper source()" -- and this file has NO relative read at all,
# measured: 0 relative source/read calls. The review checked what the setwd was USED FOR and
# concluded correctly that nothing needed it; it never asked whether the directory would EXIST.
# So the call was pure dead weight and the sole reason the step could not run.
#
# INVISIBLE TO THE GRAPH, WHICH IS THE TRANSFERABLE PART. setwd() is neither a read nor a write,
# so extract.py emits no edge for it (0 occurrences in extract_edges_G2.tsv) and views.py scored
# this file needs_update = FALSE -- the view whose whole job is "which files cannot run from the
# deposit" declared it ready. A dependency on a directory's EXISTENCE, rather than on a file
# inside it, is outside what the graph models. W80.

# --- DATA_DIR: the half-repoint this script used to have ------------------------------
# Every data read below was a bare relative path under the setwd() above, i.e. the LEGACY
# tree, while every write already went to .P$OUT. That combination silently builds the
# analysis cache from the PUBLISHED pairs and files it beside a rebuild, so a downstream
# report mixes vintages and looks fine. Anchor the reads too. (2026-07-26)
#
# Default is .P$OUT/data_mashr. Note that 01b and 02 write the SAME filenames there and 02
# overwrites 01b -- only 02's are floored, and only the floored ones are the published
# quantity. Point --data-dir at the floored tree explicitly rather than trusting the default.
DATA_DIR <- file.path(.P$OUT, "data_mashr")
local({
  a <- commandArgs(trailingOnly = TRUE)
  if ("--data-dir" %in% a) { i <- which(a == "--data-dir"); if (i < length(a)) DATA_DIR <<- a[i + 1] }
})
DATA_DIR <- normalizePath(DATA_DIR, mustWork = FALSE)
if (!dir.exists(DATA_DIR)) {
  stop(sprintf("DATA_DIR does not exist: %s\n  Run 01_prepare_data_mashr.R -> 01b -> 02 first, or pass --data-dir.", DATA_DIR))
}
cat("data dir: ", DATA_DIR, "\n", sep = "")
# Rule in code: this script WRITES here, so it must not be an input root. A --data-dir
# pointing at a reference tree would overwrite the copies verification compares against.
nmd_assert_writable(DATA_DIR, .P, "--data-dir")


FASTA_PATH <- file.path(.P$SQANTI, "nmd_lungcells_corrected.fasta")
OUTPUT_PATH <- file.path(DATA_DIR, "analysis_cache", "utr5_features_all.rds")   # the cache belongs beside the data it was built FROM   # anchored: setwd() is into an INPUT root
# The parent must exist before the write. 05s_orfik_scan.R lost 75 MINUTES of completed,
# fully-verified work to its absence on 2026-07-28: it scanned 1,540,674 ORFs, passed every
# self-check, then died on saveRDS with "cannot open compressed file ... No such file or
# directory" because tmp/out/data_mashr/analysis_cache/ had never been created. Seven
# necessary files write into that directory and the only one that CREATES it, 03b_rebuild_
# cache.R, runs at index 42 -- after five of them. Each writer now makes its own parent, so
# the guarantee is order-independent rather than an accident of run_order. W81.
dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)

cat("=== 5'UTR Feature Scan: All Classified Coding Isoforms ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Load data and identify target isoforms
# ==============================================================================
cat("Loading data files...\n")

structures <- readRDS(file.path(DATA_DIR, "structures.rds"))
cds <- readRDS(file.path(DATA_DIR, "cds.rds"))
nmd <- readRDS(file.path(DATA_DIR, "nmd_classification.rds"))

# All classified isoform IDs from all_samples
nmd_ids <- nmd$all_samples$nmd
non_nmd_ids <- nmd$all_samples$non_nmd
classified_ids <- union(nmd_ids, non_nmd_ids)

cat("  Classified isoforms (all_samples): ", length(classified_ids), "\n")
cat("    NMD: ", length(nmd_ids), "\n")
cat("    non-NMD: ", length(non_nmd_ids), "\n")

# Filter to coding only
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]
target_ids <- intersect(classified_ids, coding_ids)

cat("  Coding isoforms in CDS metadata: ", length(coding_ids), "\n")
cat("  Target isoforms (classified & coding): ", length(target_ids), "\n\n")

# Filter structures and cds to target isoforms
structures_filtered <- structures[structures$isoform_id %in% target_ids, ]
cds_filtered <- cds[cds$isoform_id %in% target_ids, ]

n_in_structures <- length(unique(structures_filtered$isoform_id))
n_in_cds <- length(unique(cds_filtered$isoform_id))
cat("  Isoforms in filtered structures: ", n_in_structures, "\n")
cat("  Isoforms in filtered cds: ", n_in_cds, "\n")

# Final target = intersection of all three
target_ids <- intersect(target_ids, unique(structures_filtered$isoform_id))
target_ids <- intersect(target_ids, unique(cds_filtered$isoform_id))
cat("  Final target isoforms (in structures & cds): ", length(target_ids), "\n\n")

# ==============================================================================
# 2. Extract sequences from FASTA
# ==============================================================================
cat("Extracting sequences from FASTA...\n")

id_file <- tempfile(fileext = ".txt")
writeLines(target_ids, id_file)

awk_cmd <- sprintf(
  "awk 'BEGIN{while((getline id < \"%s\") > 0) want[id]=1} /^>/{if(seq && found) print hdr \"\\t\" seq; hdr=substr($1,2); found=(hdr in want); seq=\"\"; next} found{seq=seq$0} END{if(seq && found) print hdr \"\\t\" seq}' '%s'",
  id_file, FASTA_PATH
)

t_awk <- proc.time()
raw_output <- system(awk_cmd, intern = TRUE)
cat("  awk completed in", round((proc.time() - t_awk)[3], 1), "seconds\n")
unlink(id_file)

seq_vec <- character(0)
if (length(raw_output) > 0) {
  split_lines <- strsplit(raw_output, "\t", fixed = TRUE)
  seq_names <- vapply(split_lines, `[`, character(1), 1)
  seq_seqs <- vapply(split_lines, `[`, character(1), 2)
  seq_vec <- setNames(seq_seqs, seq_names)
}

n_with_seq <- sum(target_ids %in% names(seq_vec))
n_missing_seq <- sum(!target_ids %in% names(seq_vec))
cat("  Found", length(seq_vec), "sequences\n")
cat("  Target IDs with sequence:", n_with_seq, "\n")
cat("  Target IDs missing sequence:", n_missing_seq, "\n\n")

# ASSERT WHAT THE EXTRACTION DID, NOT THAT THE CALL RETURNED (W122).
#
# system(intern = TRUE) raises a WARNING and never an error on a non-zero exit, so a failed
# awk -- a missing binary, an unreadable FASTA, a quoting slip -- returns character(0), the
# length() guard above leaves seq_vec empty, and the scan below then runs over NO sequences
# and writes a plausible-looking empty feature table at exit 0. Exit status is not evidence
# about what the script did; the counts are. Same rule TRACK_A_PLAN states for 05s_b's
# extended_n.
#
# A PARTIAL shortfall is legitimate -- an isoform can be absent from the FASTA -- and is
# reported above rather than asserted, because a guard that fires on a valid tree is a guard
# that gets disabled. Zero is never legitimate when there are targets.
.awk_status <- attr(raw_output, "status")
if (!is.null(.awk_status) && .awk_status != 0)
  stop(sprintf("awk sequence extraction exited %s; no sequences were extracted from %s",
               .awk_status, FASTA_PATH))
if (length(target_ids) > 0 && n_with_seq == 0)
  stop(sprintf("awk returned no sequence for ANY of the %d target isoforms in %s",
               length(target_ids), FASTA_PATH))

# ==============================================================================
# 3. Scan 5'UTR features
# ==============================================================================
cat("Scanning 5'UTR features for", length(target_ids), "isoforms...\n")

t_scan <- proc.time()
isoform_features <- scan5UtrFeatures(
  structures = structures_filtered,
  cds_metadata = cds_filtered,
  sequences = seq_vec,
  min_orf_nt = 9L,
  verbose = TRUE
)
elapsed <- round((proc.time() - t_scan)[3], 1)
cat("  Feature scan completed in", elapsed, "seconds\n")
cat("  Features computed for", nrow(isoform_features), "isoforms\n")

# ==============================================================================
# 4. ATG validation and exclusion flag
# ==============================================================================
n_validated <- sum(!is.na(isoform_features$atg_validated))
n_pass <- sum(isoform_features$atg_validated, na.rm = TRUE)
n_fail <- sum(!is.na(isoform_features$atg_validated) & !isoform_features$atg_validated)
cat("  ATG validation:", n_pass, "/", n_validated,
    "(", round(100 * n_pass / n_validated, 1), "%)\n")

isoform_features$excluded <- isoform_features$utr5_length == 0L |
  (!is.na(isoform_features$atg_validated) & !isoform_features$atg_validated)

n_zero_utr <- sum(isoform_features$utr5_length == 0L)
n_atg_fail <- sum(!is.na(isoform_features$atg_validated) & !isoform_features$atg_validated)
n_excluded <- sum(isoform_features$excluded)
n_included <- sum(!isoform_features$excluded)

cat("  Excluded:", n_excluded,
    "(", n_zero_utr, "zero-length,",
    n_atg_fail, "ATG validation failed)\n")
cat("  Included (usable):", n_included, "\n\n")

# ==============================================================================
# 5. Add NMD group labels
# ==============================================================================
cat("Adding group labels...\n")

isoform_features$group <- ifelse(
  isoform_features$isoform_id %in% nmd_ids, "NMD",
  ifelse(isoform_features$isoform_id %in% non_nmd_ids, "non_NMD", NA_character_)
)

group_counts <- table(isoform_features$group, useNA = "ifany")
cat("  Group counts:\n")
for (g in names(group_counts)) {
  cat("    ", g, ":", group_counts[g], "\n")
}

# Usable by group
usable <- isoform_features[!isoform_features$excluded, ]
usable_counts <- table(usable$group, useNA = "ifany")
cat("  Usable (non-excluded) by group:\n")
for (g in names(usable_counts)) {
  cat("    ", g, ":", usable_counts[g], "\n")
}
cat("\n")

# ==============================================================================
# 6. Save results
# ==============================================================================
output <- list(
  isoform_features = isoform_features,
  metadata = list(
    run_timestamp = Sys.time(),
    script_path = "05k_utr5_all_isoforms.R",
    fasta_path = FASTA_PATH,
    n_classified = length(classified_ids),
    n_coding = length(coding_ids),
    n_target = length(target_ids),
    n_with_sequence = n_with_seq,
    n_missing_sequence = n_missing_seq,
    n_features_computed = nrow(isoform_features),
    n_excluded = n_excluded,
    n_zero_utr = n_zero_utr,
    n_atg_fail = n_atg_fail,
    n_included = n_included,
    group_counts = group_counts,
    usable_group_counts = usable_counts,
    elapsed_scan_seconds = elapsed
  )
)

saveRDS(output, OUTPUT_PATH)
cat("Saved results to:", OUTPUT_PATH, "\n")
cat("Finished:", format(Sys.time()), "\n")
