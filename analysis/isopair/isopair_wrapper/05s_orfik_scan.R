#!/usr/bin/env Rscript
# ==============================================================================
# 05s_orfik_scan.R
#
# Comprehensive ORF detection using ORFik on all expressed isoforms.
# Computes per-ORF features: Kozak context, ORF length, position in transcript,
# upstream ATG count, downstream EJC count.
#
# Inputs:
#   - SQANTI corrected FASTA (spliced transcript sequences)
#   - data_mashr/structures.rds (exon coordinates for EJC positions)
#   - data_mashr/nmd_classification.rds (NMD/non-NMD labels)
#   - data_mashr/cds.rds (for population filtering only)
#
# Output:
#   - data_mashr/analysis_cache/orfik_scan.rds
#     $orfs_per_tx: GRangesList of all ORFs per transcript
#     $orf_features: data.frame with per-ORF features
#     $tx_summary: data.frame with per-transcript summary
#
# Usage:
#   Rscript 05s_orfik_scan.R
# ==============================================================================

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(ORFik)
  library(Biostrings)
  library(GenomicRanges)
  library(dplyr)
})

# REPOINTED 2026-07-27 (Pete's decision) off the LEGACY tree onto the deposit-rebuilt one.
#
# This file used to do setwd(.P$ISOPAIR) and read "data_mashr/..." RELATIVE to it, which made it
# the last un-repointed legacy-root reader (the W11 class). For a rebuild-from-deposit it must
# read the tree the deposit pipeline produces.
#
# The repoint was measured before it was made, because it changes which isoforms the model's
# universe is built from. It is STRICTLY ADDITIVE: legacy eligible (41,634) is a strict SUBSET of
# deposit eligible (42,063), setdiff(legacy, deposit) = 0, and the 429 gained are absent from the
# legacy cds entirely. The 1,450 legacy-only SCANNED isoforms that this appeared to drop are 6-CT
# era scan residue -- 0 of them are eligible under EITHER tree's current 4-CT labels.
# See docs/FINDINGS.md 2026-07-27 and verification/phaseC/model_universe_repoint_impact.R.
#
# .P$ISOPAIR_OUT is the declared intermediates root, NOT .P$OUT: 05s_b_orfik_scan_extend.R
# re-invokes this script with NMD_OUT overridden to a temp dir, so keying the INPUT off the output
# root would make the extension read its own empty output tree. Input and output roots are
# separately overridable (NMD_ISOPAIR_OUT / NMD_OUT) for exactly that reason.
#
# NO FALLBACK TO THE LEGACY TREE, deliberately. A fallback that always succeeds converts a loud
# failure into a silent one: the legacy tree is always populated, so a fallback would silently
# serve the vintage this repoint exists to stop using.
DATA_MASHR <- .P$ISOPAIR_OUT
if (!dir.exists(DATA_MASHR))
  stop("intermediates root not found: ", DATA_MASHR,
       "\n  This script reads the DEPOSIT-REBUILT tree. Build it first, or override",
       "\n  NMD_ISOPAIR_OUT. It deliberately does NOT fall back to the legacy tree.")

FASTA_PATH <- file.path(.P$SQANTI, "nmd_lungcells_corrected.fasta")
OUTPUT_PATH <- file.path(.P$OUT, "data_mashr/analysis_cache/orfik_scan.rds")
# The parent must exist before the write. 05s_orfik_scan.R lost 75 MINUTES of completed,
# fully-verified work to its absence on 2026-07-28: it scanned 1,540,674 ORFs, passed every
# self-check, then died on saveRDS with "cannot open compressed file ... No such file or
# directory" because tmp/out/data_mashr/analysis_cache/ had never been created. Seven
# necessary files write into that directory and the only one that CREATES it, 03b_rebuild_
# cache.R, runs at index 42 -- after five of them. Each writer now makes its own parent, so
# the guarantee is order-independent rather than an accident of run_order. W81.
dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)

cat("=== ORFik Comprehensive ORF Scan ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Define population
# ==============================================================================
cat("Loading metadata...\n")

cat("  data_mashr root:", DATA_MASHR, "\n")   # print WHICH tree; the vintage is the whole question
nmd_class <- readRDS(file.path(DATA_MASHR, "nmd_classification.rds"))
cds <- readRDS(file.path(DATA_MASHR, "cds.rds"))
structures <- readRDS(file.path(DATA_MASHR, "structures.rds"))

nmd_ids <- nmd_class$all_samples$nmd
non_nmd_ids <- nmd_class$all_samples$non_nmd
coding_ids <- cds$isoform_id[cds$coding_status == "coding"]
target_ids <- intersect(c(nmd_ids, non_nmd_ids), coding_ids)
cat("  Target isoforms:", length(target_ids), "\n")

# ==============================================================================
# 2. Load transcript sequences
# ==============================================================================
cat("\nLoading FASTA...\n")
t0 <- proc.time()
all_seqs <- readDNAStringSet(FASTA_PATH)
cat("  Loaded", length(all_seqs), "sequences in",
    round((proc.time() - t0)[3], 1), "seconds\n")

# Filter to target isoforms
target_seqs <- all_seqs[names(all_seqs) %in% target_ids]
cat("  Target sequences:", length(target_seqs), "\n")
cat("  Missing:", sum(!target_ids %in% names(target_seqs)), "\n")
rm(all_seqs)

# ==============================================================================
# 3. Find all ORFs with ORFik
# ==============================================================================
cat("\nFinding ORFs with ORFik...\n")
t0 <- proc.time()
# MINIMUMLENGTH IS IN CODONS, EXCLUDING THE START AND STOP CODONS -- not in
# nucleotides. ORFik's own documentation: "Minimum length of ORF, without
# counting 3bps for START and STOP codons. For example minimumLength = 8 will
# result in size of ORFs to be at least START + 8*3 (bp) + STOP = 30 bases."
# So the previous value of 9 imposed a floor of 3 + 27 + 3 = 33 NUCLEOTIDES,
# which is exactly what the scan produced: of 1,540,674 ORFs the minimum
# orf_length was 33, none below, 59,425 sitting on it. That floor was never
# intended and was read as "9 nt" everywhere downstream, including in the
# metadata field below, which is named min_orf_length.
#
# 0 IS THE FLOOR PETE ASKED FOR (2026-08-01, and in the retrain plan: "we
# shouldn't have a minimum ORF length, just go with the Kozak scoring threshold
# idea"). ORFik documents minimumLength = 0 as START + STOP = 6 bp, i.e. a
# start-stop element, which is what admission-by-start-codon-score requires.
# ATF4's textbook uORF1 is 3 codons and is inadmissible above minimumLength = 1.
#
# POOL SIZE MEASURED BEFORE CHANGING IT, as the plan requires. On a 3,000
# transcript sample (median length 3,002 nt), ORFs per transcript:
#   minimumLength  9 -> 36.9   (33 nt floor, the old behaviour)
#                  3 -> 47.8   (15 nt)
#                  1 -> 52.4   ( 9 nt)
#                  0 -> 55.1   ( 6 nt)   = 1.50x the old pool
# So removing the floor grows the candidate pool by half, not by orders of
# magnitude. Selection still takes MAX_ORFS per transcript downstream; what
# changes is what is available to select from, which is the point.
ORF_MIN_LENGTH_CODONS <- 0L   # codons EXCLUDING start and stop; 0 => 6 nt floor
orfs <- findORFs(target_seqs, startCodon = "ATG",
                  longestORF = FALSE, minimumLength = ORF_MIN_LENGTH_CODONS)
cat("  Completed in", round((proc.time() - t0)[3], 1), "seconds\n")
cat("  Transcripts with ORFs:", length(orfs), "\n")
cat("  Total ORFs:", sum(lengths(orfs)), "\n")
cat("  Median ORFs per transcript:", median(lengths(orfs)), "\n")

# ==============================================================================
# 4. Compute per-ORF features
# ==============================================================================
cat("\nComputing per-ORF features...\n")

# Precompute junction positions in transcript space
cat("  Computing junction positions...\n")
target_structs <- structures[structures$isoform_id %in% names(target_seqs), ]
junction_map <- setNames(vector("list", nrow(target_structs)), target_structs$isoform_id)
for (i in seq_len(nrow(target_structs))) {
  starts <- target_structs$exon_starts[[i]]
  ends <- target_structs$exon_ends[[i]]
  strand <- target_structs$strand[i]
  n_exons <- length(starts)
  if (n_exons <= 1) { junction_map[[i]] <- integer(0); next }
  exon_lengths <- ends - starts + 1L
  if (strand == "-") exon_lengths <- rev(exon_lengths)
  junction_map[[i]] <- cumsum(exon_lengths)[-n_exons]
}

# Build feature table
cat("  Computing features for each ORF...\n")
t0 <- proc.time()

# ORFik findORFs returns a CompressedIRangesList that DROPS transcripts with
# zero ORFs. The names are the original numeric indices (as strings) into the
# input DNAStringSet. Use these indices to map back to transcript names.
# BUG FIX: seq_along(orfs) gives 1:length(orfs), which misaligns when
# transcripts are dropped. Use as.integer(names(orfs)) instead.
names(orfs) <- names(target_seqs)[as.integer(names(orfs))]

feature_rows <- vector("list", length(orfs))
tx_names <- names(orfs)

for (tx_idx in seq_along(orfs)) {
  tx_id <- tx_names[tx_idx]
  tx_orfs <- orfs[[tx_idx]]
  n_orfs <- length(tx_orfs)
  if (n_orfs == 0) next

  if (!tx_id %in% names(target_seqs)) next
  tx_seq <- target_seqs[[tx_id]]
  if (is.null(tx_seq) || length(tx_seq) == 0) next
  tx_len <- length(tx_seq)
  junctions <- junction_map[[tx_id]]
  if (is.null(junctions)) junctions <- integer(0)

  orf_starts <- start(tx_orfs)
  orf_ends <- end(tx_orfs)
  orf_widths <- width(tx_orfs)

  # Per-ORF features
  kozak_m3 <- character(n_orfs)
  kozak_p4 <- character(n_orfs)
  kozak_score <- integer(n_orfs)
  n_upstream_atgs <- integer(n_orfs)
  n_downstream_ejc <- integer(n_orfs)
  start_codon <- character(n_orfs)
  stop_codon <- character(n_orfs)

  # Find all ATG positions in transcript for upstream ATG counting
  atg_matches <- tryCatch(matchPattern("ATG", tx_seq), error = function(e) NULL)
  if (is.null(atg_matches)) next
  all_atg_pos <- start(atg_matches)

  for (j in seq_len(n_orfs)) {
    os <- orf_starts[j]
    oe <- orf_ends[j]

    # Bounds check
    if (os < 1L || oe > tx_len || os > oe) next

    # Start codon (should be ATG)
    if (os + 2L <= tx_len) {
      start_codon[j] <- as.character(subseq(tx_seq, os, os + 2L))
    } else {
      start_codon[j] <- NA_character_
    }

    # Stop codon: last 3 nt of the ORF (ORFik includes the stop codon in the range)
    if (oe >= 3L && oe <= tx_len) {
      stop_codon[j] <- as.character(subseq(tx_seq, oe - 2L, oe))
    } else {
      stop_codon[j] <- NA_character_
    }

    # Kozak context: -3 position (A/G = strong) and +4 position (G = strong)
    if (os >= 4L) {
      kozak_m3[j] <- as.character(subseq(tx_seq, os - 3L, os - 3L))
    } else {
      kozak_m3[j] <- NA_character_
    }
    if (os + 3L >= 1L && os + 3L <= tx_len) {
      kozak_p4[j] <- as.character(subseq(tx_seq, os + 3L, os + 3L))
    } else {
      kozak_p4[j] <- NA_character_
    }

    # Kozak score: 0=weak, 1=moderate, 2=strong
    has_m3 <- !is.na(kozak_m3[j]) && kozak_m3[j] %in% c("A", "G")
    has_p4 <- !is.na(kozak_p4[j]) && kozak_p4[j] == "G"
    kozak_score[j] <- as.integer(has_m3) + as.integer(has_p4)

    # Upstream ATGs (ATG positions before this ORF's start)
    n_upstream_atgs[j] <- sum(all_atg_pos < os)

    # Downstream EJCs: junctions after the stop codon
    # ORFik ORF end = last nt of ORF (including stop codon). Stop ends at oe.
    n_downstream_ejc[j] <- sum(junctions > oe)
  }

  feature_rows[[tx_idx]] <- data.frame(
    isoform_id = tx_id,
    orf_idx = seq_len(n_orfs),
    orf_start = orf_starts,
    orf_end = orf_ends,
    orf_length = orf_widths,
    tx_length = tx_len,
    frac_position = orf_starts / tx_len,
    frac_tx_covered = orf_widths / tx_len,
    start_codon = start_codon,
    stop_codon = stop_codon,
    kozak_m3 = kozak_m3,
    kozak_p4 = kozak_p4,
    kozak_score = kozak_score,
    n_upstream_atgs = n_upstream_atgs,
    n_downstream_ejc = n_downstream_ejc,
    has_downstream_ejc = n_downstream_ejc > 0L,
    n_junctions = length(junctions),
    stringsAsFactors = FALSE
  )

  if (tx_idx %% 10000 == 0) {
    cat(sprintf("    %d / %d transcripts (%.0f sec)\n",
                tx_idx, length(orfs), (proc.time() - t0)[3]))
  }
}

orf_features <- do.call(rbind, feature_rows)

# AN ALL-EMPTY BATCH IS A LEGITIMATE OUTCOME, AND IT USED TO BE A NULL DEREFERENCE.
# findORFs drops transcripts with no ORF, so a batch in which EVERY transcript lacks an
# ATG-initiated ORF >= 9 nt makes do.call(rbind, <all NULL>) return NULL, and the group_by
# below then failed with "no applicable method for 'group_by' applied to an object of class
# NULL" -- a message that says nothing about the cause. Measured 2026-07-29: the 20 isoforms
# eligible-but-unscanned are exactly this case (0 ORFs across all 20, sequences 335-940 nt).
# Raise a NAMED sentinel the caller can recognise instead. 05s_b keys on this string to
# record the transcripts as legitimately ORF-less rather than retrying them forever.
if (is.null(orf_features) || nrow(orf_features) == 0L) {
  stop(sprintf("NO_ORFS_IN_BATCH: all %d scanned transcript(s) yielded zero ORFs under ",
               length(target_seqs)),
       "startCodon=ATG, minimumLength=9. This is a valid result for a small batch and a ",
       "red flag for a full scan -- the caller must decide which.")
}
cat("  Completed in", round((proc.time() - t0)[3], 1), "seconds\n")
cat("  Total ORF feature rows:", nrow(orf_features), "\n")

# ==============================================================================
# 5. Per-transcript summary
# ==============================================================================
cat("\nComputing per-transcript summary...\n")

tx_summary <- orf_features %>%
  group_by(isoform_id) %>%
  summarise(
    n_orfs = n(),
    n_nmd_competent = sum(has_downstream_ejc),
    max_orf_length = max(orf_length),
    max_downstream_ejc = max(n_downstream_ejc),
    n_strong_kozak = sum(kozak_score == 2L),
    n_strong_kozak_nmd_comp = sum(kozak_score == 2L & has_downstream_ejc),
    best_kozak_nmd_comp = max(kozak_score[has_downstream_ejc], 0L),
    longest_nmd_comp_orf = ifelse(any(has_downstream_ejc),
                                   max(orf_length[has_downstream_ejc]), 0L),
    tx_length = first(tx_length),
    n_junctions = first(n_junctions),
    .groups = "drop"
  )

# Add NMD status
tx_summary$is_nmd <- as.integer(tx_summary$isoform_id %in% nmd_ids)

# Add chromosome
chr_map <- structures[, c("isoform_id", "chr")]
tx_summary <- merge(tx_summary, chr_map, by = "isoform_id", all.x = TRUE)

# ==============================================================================
# 6. Summary statistics
# ==============================================================================
cat("\n=== SUMMARY ===\n\n")
cat("Transcripts scanned:", nrow(tx_summary), "\n")
cat("  NMD:", sum(tx_summary$is_nmd), "  Non-NMD:", sum(!tx_summary$is_nmd), "\n\n")

cat("ORF features:\n")
cat("  Total ORFs:", nrow(orf_features), "\n")
cat("  Median ORFs per transcript:", median(tx_summary$n_orfs), "\n")
cat("  Median NMD-competent ORFs:", median(tx_summary$n_nmd_competent), "\n")
cat("  Median max ORF length:", median(tx_summary$max_orf_length), "nt\n")
cat("  ORFs with strong Kozak (score=2):",
    sum(orf_features$kozak_score == 2),
    sprintf("(%.1f%%)\n", 100 * mean(orf_features$kozak_score == 2)))
cat("  ORFs with downstream EJC:",
    sum(orf_features$has_downstream_ejc),
    sprintf("(%.1f%%)\n", 100 * mean(orf_features$has_downstream_ejc)))

cat("\nBy NMD status:\n")
for (status in c(1, 0)) {
  sub <- tx_summary[tx_summary$is_nmd == status, ]
  label <- if (status == 1) "NMD" else "Non-NMD"
  # %.0f, not %d: median() returns a DOUBLE and lands on a half-value whenever the
  # subgroup has an even count, which aborts the script HERE -- before saveRDS at the
  # bottom, so the whole scan is lost to a printf in a summary. The 43,039-transcript
  # run happened to give whole-number medians and never hit it; a 474-isoform extension
  # did on its first try (Non-NMD n=22).
  cat(sprintf("  %s (n=%d): median ORFs=%.0f, median NMD-comp=%.0f, median max_orf=%.0f nt\n",
              label, nrow(sub), median(sub$n_orfs), median(sub$n_nmd_competent),
              median(sub$max_orf_length)))
}

# ==============================================================================
# 7. Save
# ==============================================================================
results <- list(
  orf_features = orf_features,
  tx_summary = tx_summary,
  metadata = list(
    n_transcripts = nrow(tx_summary),
    n_orfs = nrow(orf_features),
    # BOTH THE PARAMETER AND ITS CONSEQUENCE, in their own units, because
    # recording only the parameter is what hid the 33 nt floor: the field was
    # called min_orf_length and carried 9, and every reader took that for
    # nucleotides. Whoever checks this next should check the OBSERVED floor.
    min_orf_length_codons = ORF_MIN_LENGTH_CODONS,
    min_orf_length_nt_implied = 3L + ORF_MIN_LENGTH_CODONS * 3L + 3L,
    min_orf_length_nt_observed = min(orf_features$orf_length),
    start_codon = "ATG",
    # WHICH STRUCTURES PRODUCED THIS SCAN. Added 2026-07-27 to close a guard gap: a scan is a
    # function of the exon structures as much as of the parameters, and 05s_b would merge two
    # scans agreeing on min_orf_length and start_codon while disagreeing on the underlying
    # structures -- silently mixing vintages. Recording it here is what lets the merge refuse.
    structures_path = file.path(DATA_MASHR, "structures.rds"),
    structures_md5  = unname(tools::md5sum(file.path(DATA_MASHR, "structures.rds"))),
    structures_n    = nrow(structures),
    run_timestamp = Sys.time()
  )
)

# ==============================================================================
# 8. Verification: confirm all ORFs start with ATG
# ==============================================================================
cat("\n=== Verification ===\n")
n_atg <- sum(orf_features$start_codon == "ATG", na.rm = TRUE)
n_non_atg <- sum(orf_features$start_codon != "ATG" & !is.na(orf_features$start_codon))
n_na <- sum(is.na(orf_features$start_codon))
cat("  ATG starts:", n_atg, sprintf("(%.1f%%)\n", 100 * n_atg / nrow(orf_features)))
cat("  Non-ATG starts:", n_non_atg, "\n")
cat("  NA starts:", n_na, "\n")
if (n_non_atg > 0 || n_na > 0) {
  stop("VERIFICATION FAILED: Found ", n_non_atg, " non-ATG and ", n_na,
       " NA start codons. Name mapping or coordinate extraction is likely wrong.")
}
cat("  PASSED: All ORFs start with ATG\n")

# Also verify transcript IDs are valid
n_valid_tx <- sum(orf_features$isoform_id %in% names(target_seqs))
cat("  Valid transcript IDs:", n_valid_tx, "/", nrow(orf_features), "\n")
if (n_valid_tx != nrow(orf_features)) {
  stop("VERIFICATION FAILED: ", nrow(orf_features) - n_valid_tx,
       " ORFs have isoform_ids not in the input sequences.")
}
cat("  PASSED: All isoform_ids match input sequences\n")

saveRDS(results, OUTPUT_PATH)
cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Completed:", format(Sys.time()), "\n")
