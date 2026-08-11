#!/usr/bin/env Rscript
# ==============================================================================
# 05r_ref_atg_analysis.R
#
# Reference-ATG-anchored NMD analysis. Backend: Isopair::traceReferenceAtg().
#
# For each gene-matched NMD pair (C2) and Control pair (C4):
#   1. Trace the reference ATG through the comparator (Isopair::traceReferenceAtg)
#   2. For effectively-PTC cases, attribute the responsible splice event
#      using the shared attribute_ptc_events() / attribute_3utr_splice()
#      helpers in analysis_functions.R
#
# Categories produced (from Isopair::traceReferenceAtg):
#   effectively_ptc, truncated_no_ejc, no_downstream_ejc, ref_atg_lost,
#   no_ref_cds, mapping_failed
#
# Inputs:
#   - data_mashr/cds.rds, structures.rds, profiles_c2/c4_allsamples.rds
#   - data_mashr/analysis_cache/ptc_c2_allsamples.rds
#   - data_mashr/analysis_cache/fw_c2_allsamples.rds
#   - SQANTI corrected FASTA
#
# Output:
#   - data_mashr/analysis_cache/ref_atg_analysis.rds
#
# Usage:
#   Rscript 05r_ref_atg_analysis.R
# ==============================================================================

.nmd_root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config", "paths.yml"))) {
  p <- dirname(d); if (identical(p, d)) stop("repo root not found"); d <- p }; d })
source(file.path(.nmd_root, "R", "load_config.R")); .P <- nmd_paths(.nmd_root)
suppressPackageStartupMessages({
  library(dplyr)
  library(Isopair)
})

# setwd(.P$ISOPAIR) REMOVED 2026-07-28. Unlike its two siblings this one served a REAL
# dependency -- the analysis_functions.R helper sourced below -- but it still made the
# script unrunnable for a reader, because .P$ISOPAIR is the legacy authoring tree. The
# clean-room rerun dies on it with "cannot change working directory".
#
# The old comment recorded the fix without taking it: the legacy copy was "verified
# 2026-07-26 byte-identical to the in-repo copy". So the helper is now sourced from the
# in-repo copy by absolute path and the setwd is unnecessary. Equivalence is not assumed --
# it was measured on 2026-07-26 and that comment is the record of it. W80.

# --- DATA_DIR: the half-repoint this script used to have ------------------------------
# Every data read below was a bare relative path under the setwd() above, i.e. the LEGACY
# tree, while every write already went to .P$OUT. That combination silently builds the
# analysis cache from the PUBLISHED pairs and files it beside a rebuild, so a downstream
# report mixes vintages and looks fine. Anchor the reads too.
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
OUTPUT_PATH <- file.path(DATA_DIR, "analysis_cache", "ref_atg_analysis.rds")   # the cache belongs beside the data it was built FROM   # anchored: setwd() is into an INPUT root
# The parent must exist before the write. 05s_orfik_scan.R lost 75 MINUTES of completed,
# fully-verified work to its absence on 2026-07-28: it scanned 1,540,674 ORFs, passed every
# self-check, then died on saveRDS with "cannot open compressed file ... No such file or
# directory" because tmp/out/data_mashr/analysis_cache/ had never been created. Seven
# necessary files write into that directory and the only one that CREATES it, 03b_rebuild_
# cache.R, runs at index 42 -- after five of them. Each writer now makes its own parent, so
# the guarantee is order-independent rather than an accident of run_order. W81.
dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)

cat("=== Reference-ATG-Anchored NMD Analysis (Isopair::traceReferenceAtg) ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Load data and define populations
# ==============================================================================
cat("Loading data...\n")

cds <- readRDS(file.path(DATA_DIR, "cds.rds"))
structures <- readRDS(file.path(DATA_DIR, "structures.rds"))
profiles_c2 <- readRDS(file.path(DATA_DIR, "profiles_c2_allsamples.rds"))
profiles_c4 <- readRDS(file.path(DATA_DIR, "profiles_c4_allsamples.rds"))
ptc_assoc <- readRDS(file.path(DATA_DIR, "analysis_cache/ptc_c2_allsamples.rds"))

# Gene-match
shared <- inner_join(
  unique(profiles_c2[, c("gene_id", "reference_isoform_id")]),
  unique(profiles_c4[, c("gene_id", "reference_isoform_id")]),
  by = c("gene_id", "reference_isoform_id"))

gm_c2 <- semi_join(profiles_c2, shared, by = c("gene_id", "reference_isoform_id"))
gm_c4 <- semi_join(profiles_c4, shared, by = c("gene_id", "reference_isoform_id"))

coding_ids <- cds$isoform_id[cds$coding_status == "coding"]
c2_coding <- gm_c2 %>%
  filter(reference_isoform_id %in% coding_ids, comparator_isoform_id %in% coding_ids)
c4_coding <- gm_c4 %>%
  filter(reference_isoform_id %in% coding_ids, comparator_isoform_id %in% coding_ids)

# Add original PTC status
ptc_lookup <- ptc_assoc$pair_summary[, c("comparator_isoform_id", "comp_has_ptc")]
c2_coding <- merge(c2_coding, ptc_lookup, by = "comparator_isoform_id", all.x = TRUE)
c2_coding$comp_has_ptc[is.na(c2_coding$comp_has_ptc)] <- FALSE

cat("  Gene-matched coding pairs: C2 =", nrow(c2_coding), " C4 =", nrow(c4_coding), "\n")
cat("  Original PTC+:", sum(c2_coding$comp_has_ptc),
    " PTC-:", sum(!c2_coding$comp_has_ptc), "\n")

# ==============================================================================
# 2. Extract transcript sequences from FASTA
# ==============================================================================
cat("\nExtracting sequences...\n")

all_isoforms <- unique(c(c2_coding$reference_isoform_id,
                          c2_coding$comparator_isoform_id,
                          c4_coding$reference_isoform_id,
                          c4_coding$comparator_isoform_id))

id_file <- tempfile(fileext = ".txt")
writeLines(all_isoforms, id_file)
awk_cmd <- sprintf(
  "awk 'BEGIN{while((getline id < \"%s\") > 0) want[id]=1} /^>/{if(seq && found) print hdr \"\\t\" seq; hdr=substr($1,2); found=(hdr in want); seq=\"\"; next} found{seq=seq$0} END{if(seq && found) print hdr \"\\t\" seq}' '%s'",
  id_file, FASTA_PATH)
t0 <- proc.time()
raw <- system(awk_cmd, intern = TRUE)
cat("  awk completed in", round((proc.time() - t0)[3], 1), "seconds\n")
unlink(id_file)

split_lines <- strsplit(raw, "\t", fixed = TRUE)
seq_vec <- setNames(
  toupper(vapply(split_lines, `[`, character(1), 2)),
  vapply(split_lines, `[`, character(1), 1))
n_with_seq <- sum(all_isoforms %in% names(seq_vec))
cat("  Sequences loaded:", length(seq_vec),
    sprintf(" (targets with sequence: %d / %d)\n", n_with_seq, length(all_isoforms)))

# ASSERT WHAT THE EXTRACTION DID, NOT THAT THE CALL RETURNED (W122). system(intern = TRUE)
# raises a WARNING and never an error on non-zero exit, so a failed awk returns character(0),
# strsplit/vapply reduce that to an EMPTY seq_vec without complaint, and the ref-ATG trace
# below then runs over no sequences at exit 0. Nothing here printed the target overlap at all
# before this, so the shortfall was not even visible. A PARTIAL shortfall is legitimate -- an
# isoform can be absent from the FASTA -- and is reported rather than asserted; zero is not.
.awk_status <- attr(raw, "status")
if (!is.null(.awk_status) && .awk_status != 0)
  stop(sprintf("awk sequence extraction exited %s; no sequences were extracted from %s",
               .awk_status, FASTA_PATH))
if (length(all_isoforms) > 0 && n_with_seq == 0)
  stop(sprintf("awk returned no sequence for ANY of the %d target isoforms in %s",
               length(all_isoforms), FASTA_PATH))

# Structures lookup used by Section 4 (transcript-to-genomic mapping for PTC stops)
structs_lookup <- setNames(
  lapply(seq_len(nrow(structures)), function(i) {
    list(starts = structures$exon_starts[[i]],
         ends = structures$exon_ends[[i]],
         strand = structures$strand[i])
  }), structures$isoform_id)

# ==============================================================================
# 3. Trace ref ATG through comparator (Isopair::traceReferenceAtg)
# ==============================================================================
analyze_pairs <- function(pairs, label) {
  cat(sprintf("\n=== Analyzing %s pairs (n=%d) ===\n", label, nrow(pairs)))

  t0 <- proc.time()
  trace <- Isopair::traceReferenceAtg(
    pairs = pairs[, c("reference_isoform_id", "comparator_isoform_id")],
    structures = structures,
    cds_metadata = cds,
    sequences = seq_vec,
    ejc_threshold = 50L,
    resolve_alt_start = FALSE  # legacy behavior; alt-start handled separately if needed
  )
  cat(sprintf("  traceReferenceAtg completed in %.0f seconds\n",
              (proc.time() - t0)[3]))

  # Augment with wrapper-specific columns expected by Section 4 and downstream consumers
  pair_key <- paste(pairs$reference_isoform_id, pairs$comparator_isoform_id, sep = "|")
  trace_key <- paste(trace$reference_isoform_id, trace$comparator_isoform_id, sep = "|")
  trace$gene_id <- pairs$gene_id[match(trace_key, pair_key)]
  trace$strand <- cds$strand[match(trace$comparator_isoform_id, cds$isoform_id)]

  # Backward-compatible column name used by 05_final_report_mashr.Rmd
  trace$comp_orf_length_from_ref_atg <- trace$comp_orf_length

  # Summary
  cat("\n  Results:\n")
  tab <- sort(table(trace$category), decreasing = TRUE)
  for (nm in names(tab)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", nm, tab[nm], 100 * tab[nm] / nrow(trace)))
  }

  trace
}

# Run on C2 (NMD) pairs
c2_results <- analyze_pairs(c2_coding, "NMD (C2)")
c2_results$original_ptc <- c2_coding$comp_has_ptc[
  match(c2_results$comparator_isoform_id, c2_coding$comparator_isoform_id)]

# Run on C4 (Control) pairs
c4_results <- analyze_pairs(c4_coding, "Control (C4)")

# ==============================================================================
# 4. Splice event attribution for ALL effectively-PTC cases
#    (a) Truncated ORF → frameshift / in-frame stop (attribute_ptc_events)
#    (b) Same/longer ORF with downstream EJCs → 3'UTR splice (attribute_3utr_splice)
#
# 2026-06-15: Extended the eff_ptc_mask to include ALL effectively_ptc cases
# (both original_ptc=FALSE and original_ptc=TRUE). This makes ref-AUG-derived
# attribution the canonical attribution for the entire PTC+ class under the
# §4 framework — no TD2-PTC-stop-based path is needed downstream. The 1,012
# original_ptc=TRUE effectively_ptc cases now also get attr_mechanism /
# attr_event populated. The original 900 original_ptc=FALSE attributions are
# unchanged (deterministic algorithm, same inputs).
# See figures/multipanel/figure4_ptcneg_and_model/RATIONALE.md §3 for rationale.
# ==============================================================================
cat("\n=== Splice Event Attribution (Shared Function) ===\n")

source(file.path(.nmd_root, "analysis/isopair/isopair_wrapper/analysis_functions.R"))

fw_c2 <- readRDS(file.path(DATA_DIR, "analysis_cache/fw_c2_allsamples.rds"))

# All effectively-PTC cases in C2 (ref-AUG-derived PTC determination)
eff_ptc_mask <- c2_results$category == "effectively_ptc"
eff_ptc_idx <- which(eff_ptc_mask)
cat("  All effectively-PTC (ref-AUG-determined):", length(eff_ptc_idx), "\n")
cat("    of which original_ptc=TRUE: ",
    sum(c2_results$original_ptc[eff_ptc_idx], na.rm=TRUE), "\n")
cat("    of which original_ptc=FALSE:",
    sum(!c2_results$original_ptc[eff_ptc_idx], na.rm=TRUE), "\n")

# Map premature stop from transcript-space to genomic coordinates
# (Isopair::transcriptToGenomic — strand-aware)
comp_stop_genomic <- rep(NA_integer_, nrow(c2_results))
for (k in eff_ptc_idx) {
  comp_id <- c2_results$comparator_isoform_id[k]
  strand <- c2_results$strand[k]
  comp_stop_tx <- c2_results$comp_stop_tx_pos[k]
  comp_s <- structs_lookup[[comp_id]]
  if (!is.null(comp_s) && !is.na(comp_stop_tx)) {
    comp_stop_genomic[k] <- Isopair::transcriptToGenomic(
      comp_stop_tx, comp_s$starts, comp_s$ends, strand)
  }
}
c2_results$comp_stop_genomic <- comp_stop_genomic

# Split reclassified pairs into:
#   (a) Truncated ORF: comp_orf < ref_orf → frameshift/in-frame stop attribution
#   (b) Same/longer ORF with downstream EJCs → 3'UTR splice mechanism
eff_ptc_all <- c2_results[eff_ptc_mask, ]
is_truncated <- !is.na(eff_ptc_all$ref_orf_length) &
  eff_ptc_all$comp_orf_length_from_ref_atg < eff_ptc_all$ref_orf_length
is_truncated[is.na(is_truncated)] <- FALSE  # treat NA as non-truncated

cat("  Truncated ORF (frameshift/in-frame stop):", sum(is_truncated), "\n")
cat("  Same/longer ORF (3'UTR splice):", sum(!is_truncated), "\n")

# --- (a) Truncated ORF attribution via shared function ---
trunc_pairs <- eff_ptc_all[is_truncated, ]
attr_mechanism <- rep(NA_character_, nrow(c2_results))
attr_event <- rep(NA_character_, nrow(c2_results))

if (nrow(trunc_pairs) > 0) {
  ptc_stop_vec <- setNames(trunc_pairs$comp_stop_genomic, trunc_pairs$comparator_isoform_id)
  atg_vec <- setNames(trunc_pairs$ref_atg_genomic, trunc_pairs$comparator_isoform_id)
  strand_vec_t <- setNames(trunc_pairs$strand, trunc_pairs$comparator_isoform_id)

  attr_result <- attribute_ptc_events(
    pairs = trunc_pairs[, c("comparator_isoform_id", "reference_isoform_id")],
    fw_events = fw_c2$events,
    profiles = c2_coding,
    ptc_genomic_pos = ptc_stop_vec,
    atg_genomic_pos = atg_vec,
    strand_vec = strand_vec_t,
    is_frameshift_vec = NULL
  )

  attr_lookup <- setNames(seq_len(nrow(attr_result)), attr_result$comparator_isoform_id)
  trunc_idx <- eff_ptc_idx[is_truncated]
  for (k in seq_along(trunc_idx)) {
    comp_id <- c2_results$comparator_isoform_id[trunc_idx[k]]
    j <- attr_lookup[comp_id]
    if (!is.na(j)) {
      attr_mechanism[trunc_idx[k]] <- attr_result$mechanism[j]
      attr_event[trunc_idx[k]] <- attr_result$ptc_causing_event[j]
    }
  }
}

# --- (b) Same/longer ORF attribution: 3'UTR splice ---
nontrunc_pairs <- eff_ptc_all[!is_truncated, ]
nontrunc_idx <- eff_ptc_idx[!is_truncated]

if (nrow(nontrunc_pairs) > 0) {
  utr3_stop_vec <- setNames(nontrunc_pairs$comp_stop_genomic, nontrunc_pairs$comparator_isoform_id)
  utr3_strand_vec <- setNames(nontrunc_pairs$strand, nontrunc_pairs$comparator_isoform_id)

  utr3_result <- attribute_3utr_splice(
    pairs = nontrunc_pairs[, c("comparator_isoform_id", "reference_isoform_id")],
    profiles = c2_coding,
    stop_genomic_pos = utr3_stop_vec,
    strand_vec = utr3_strand_vec
  )

  utr3_lookup <- setNames(seq_len(nrow(utr3_result)), utr3_result$comparator_isoform_id)
  for (k in seq_along(nontrunc_idx)) {
    comp_id <- c2_results$comparator_isoform_id[nontrunc_idx[k]]
    j <- utr3_lookup[comp_id]
    if (!is.na(j)) {
      attr_mechanism[nontrunc_idx[k]] <- utr3_result$mechanism[j]
      attr_event[nontrunc_idx[k]] <- utr3_result$ptc_causing_event[j]
    }
  }
}

c2_results$attr_mechanism <- attr_mechanism
c2_results$attr_event <- attr_event

# Summary of attribution
reclassified <- c2_results[eff_ptc_mask, ]
cat("\n  PTC-causing mechanism for reclassified pairs:\n")
mech_tab <- sort(table(reclassified$attr_mechanism), decreasing = TRUE)
for (nm in names(mech_tab)) {
  cat(sprintf("    %s: %d (%.1f%%)\n", nm, mech_tab[nm],
              100 * mech_tab[nm] / nrow(reclassified)))
}

cat("\n  Event types within Frameshift mechanism:\n")
fs_reclass <- reclassified[reclassified$attr_mechanism == "Frameshift", ]
if (nrow(fs_reclass) > 0) {
  evt_tab <- sort(table(fs_reclass$attr_event), decreasing = TRUE)
  for (nm in names(evt_tab)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", nm, evt_tab[nm],
                100 * evt_tab[nm] / nrow(fs_reclass)))
  }
}

cat("\n  Event types within In-frame stop mechanism:\n")
ifs_reclass <- reclassified[reclassified$attr_mechanism == "In-frame stop", ]
if (nrow(ifs_reclass) > 0) {
  evt_tab2 <- sort(table(ifs_reclass$attr_event), decreasing = TRUE)
  for (nm in names(evt_tab2)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", nm, evt_tab2[nm],
                100 * evt_tab2[nm] / nrow(ifs_reclass)))
  }
}

cat("\n  Event types within 3'UTR splice mechanism:\n")
utr3_reclass <- reclassified[reclassified$attr_mechanism == "3'UTR splice", ]
if (nrow(utr3_reclass) > 0) {
  evt_tab3 <- sort(table(utr3_reclass$attr_event), decreasing = TRUE)
  for (nm in names(evt_tab3)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", nm, evt_tab3[nm],
                100 * evt_tab3[nm] / nrow(utr3_reclass)))
  }
}

# ==============================================================================
# 5. Summary and save
# ==============================================================================
# ==============================================================================
# 4.5 Additive fields for §4 / Figure 4 (RATIONALE.md §8 step 1).
# Strictly additive — existing columns must remain byte-identical.
#   ref_source              — "GENCODE" (ENST-prefixed) / "TD2" (novel) per pair
#   comp_tx_len_nt          — transcript length of comparator isoform (nt)
#   ref_tx_len_nt           — transcript length of reference isoform (nt)
#   tx_len_delta_nt         — comp_tx_len_nt − ref_tx_len_nt (TES-variation flag)
#   utr3_to_tx_end_nt       — tx_len(comp) − comp_stop_tx_pos. Translation-based
#                             3'UTR length (the region the cell processes after
#                             the comparator's ORF ends — IS the PTC for
#                             effectively_ptc pairs, hence definitionally
#                             inflated for that group; see RATIONALE §4.2).
#   non_ptc_stop_tx_pos     — position of the first stop codon downstream of
#                             ref-AUG that is NOT a PTC by the 50-nt rule
#                             (i.e., has no downstream EJC > 50 nt past its
#                             end). For non-effectively_ptc traceable categories
#                             this equals comp_stop_tx_pos by construction.
#                             For effectively_ptc, walks past the PTC to find
#                             the next in-frame stop in the last exon. NA if
#                             no such stop exists within the transcript or if
#                             ref AUG is not traceable.
#   utr3_via_non_ptc_stop_nt — tx_len(comp) − non_ptc_stop_tx_pos. Bias-free
#                             3'UTR measure: defined uniformly across all
#                             groups as "post-natural-stop region".
# Note: `mechanism_class` is intentionally NOT added here. It is a derived
# field computed at use-site by figures/lib/mechanism_class.R, so the
# classification can evolve without invalidating this cache.
# ==============================================================================

# Helper: junction positions in transcript coordinates (strand-aware).
.junction_positions <- function(exon_starts, exon_ends, strand) {
  exon_lens <- exon_ends - exon_starts + 1L
  if (strand == "-") exon_lens <- rev(exon_lens)
  if (length(exon_lens) <= 1L) return(integer(0))
  cumsum(exon_lens[-length(exon_lens)])
}

# Helper: walking downstream of a position (e.g., a PTC), find the next
# in-frame stop codon that is NOT a PTC by the 50-nt rule. Returns the
# 1-based transcript position of that stop, or NA if none exists.
.next_non_ptc_stop <- function(seq, last_stop_pos, exon_starts, exon_ends,
                               strand, ejc_threshold = 50L) {
  if (is.na(last_stop_pos) || nchar(seq) < last_stop_pos + 5L) {
    return(NA_integer_)
  }
  stops <- c("TAA", "TAG", "TGA")
  junctions <- .junction_positions(exon_starts, exon_ends, strand)
  tx_len <- nchar(seq)
  current_pos <- last_stop_pos + 3L
  while (current_pos + 2L <= tx_len) {
    codon <- substr(seq, current_pos, current_pos + 2L)
    if (codon %in% stops) {
      stop_end <- current_pos + 2L
      n_dejc <- sum(junctions > (stop_end + ejc_threshold - 1L))
      if (n_dejc == 0L) return(as.integer(current_pos))
    }
    current_pos <- current_pos + 3L
  }
  NA_integer_
}

augment_with_additive_fields <- function(trace_df) {
  trace_df$ref_source <- ifelse(grepl("^ENST", trace_df$reference_isoform_id),
                                "GENCODE", "TD2")
  trace_df$comp_tx_len_nt  <- tx_len_lookup[trace_df$comparator_isoform_id]
  trace_df$ref_tx_len_nt   <- tx_len_lookup[trace_df$reference_isoform_id]
  trace_df$tx_len_delta_nt <- trace_df$comp_tx_len_nt - trace_df$ref_tx_len_nt
  utr3 <- trace_df$comp_tx_len_nt - trace_df$comp_stop_tx_pos
  utr3[!is.na(utr3) & utr3 < 0] <- NA_integer_
  trace_df$utr3_to_tx_end_nt <- utr3

  # First non-PTC stop downstream of ref AUG.
  # For non-effectively_ptc traceable categories: comp_stop_tx_pos IS already
  # a non-PTC stop (by category definition, n_downstream_ejc == 0 → not a PTC).
  # For effectively_ptc: walk past the PTC to find the next in-frame non-PTC stop.
  # For ref_atg_lost / mapping_failed / NA: NA.
  TRACEABLE_NON_PTC <- c("no_downstream_ejc", "truncated_no_ejc")
  non_ptc_stop <- rep(NA_integer_, nrow(trace_df))

  # Cheap path: copy comp_stop_tx_pos for categories already non-PTC
  cheap_mask <- !is.na(trace_df$category) &
                trace_df$category %in% TRACEABLE_NON_PTC &
                !is.na(trace_df$comp_stop_tx_pos)
  non_ptc_stop[cheap_mask] <- trace_df$comp_stop_tx_pos[cheap_mask]

  # Walk path: effectively_ptc cases — find next non-PTC stop downstream of PTC
  walk_mask <- !is.na(trace_df$category) & trace_df$category == "effectively_ptc" &
               !is.na(trace_df$comp_stop_tx_pos)
  walk_idx <- which(walk_mask)
  for (i in walk_idx) {
    cid <- trace_df$comparator_isoform_id[i]
    cs <- structs_lookup[[cid]]
    if (is.null(cs)) next
    strand_i <- trace_df$strand[i]
    if (is.na(strand_i) || strand_i == "") next
    if (!cid %in% names(seq_vec)) next
    non_ptc_stop[i] <- .next_non_ptc_stop(
      seq_vec[cid], trace_df$comp_stop_tx_pos[i],
      cs$starts, cs$ends, strand_i, ejc_threshold = 50L)
  }
  trace_df$non_ptc_stop_tx_pos <- non_ptc_stop

  utr3_via_npp <- trace_df$comp_tx_len_nt - non_ptc_stop
  utr3_via_npp[!is.na(utr3_via_npp) & utr3_via_npp < 0] <- NA_integer_
  trace_df$utr3_via_non_ptc_stop_nt <- utr3_via_npp

  trace_df
}

tx_len_lookup <- setNames(
  vapply(seq_len(nrow(structures)),
         function(i) sum(structures$exon_ends[[i]] -
                         structures$exon_starts[[i]] + 1L),
         integer(1)),
  structures$isoform_id)

c2_results <- augment_with_additive_fields(c2_results)
c4_results <- augment_with_additive_fields(c4_results)

cat("  Added additive fields: ref_source, comp_tx_len_nt, ref_tx_len_nt,",
    "tx_len_delta_nt, utr3_to_tx_end_nt, non_ptc_stop_tx_pos,",
    "utr3_via_non_ptc_stop_nt\n")

cat("\n=== FINAL SUMMARY ===\n\n")

# Three groups
group2 <- sum(c2_results$ref_atg_exonic_in_comp == TRUE, na.rm = TRUE)
group3 <- sum(c2_results$ref_atg_exonic_in_comp == FALSE, na.rm = TRUE)
cat("Three groups:\n")
cat("  Group 2 (ref ATG available):", group2, sprintf("(%.1f%%)\n", 100 * group2 / nrow(c2_results)))
cat("  Group 3 (ref ATG lost):", group3, sprintf("(%.1f%%)\n", 100 * group3 / nrow(c2_results)))

# Reclassification
cat("\nAmong originally PTC- pairs (n =", sum(!c2_results$original_ptc), "):\n")
ptc_neg <- c2_results[!c2_results$original_ptc, ]
tab <- sort(table(ptc_neg$category), decreasing = TRUE)
for (nm in names(tab)) {
  cat(sprintf("  %s: %d (%.1f%%)\n", nm, tab[nm], 100 * tab[nm] / nrow(ptc_neg)))
}

# Updated accounting
orig_ptc_pos <- sum(c2_results$original_ptc)
reclass_ptc <- sum(c2_results$category == "effectively_ptc" & !c2_results$original_ptc, na.rm = TRUE)
total_ptc <- orig_ptc_pos + reclass_ptc
cat(sprintf("\nUpdated PTC accounting:\n"))
cat(sprintf("  Original PTC+: %d\n", orig_ptc_pos))
cat(sprintf("  Reclassified effectively-PTC: %d\n", reclass_ptc))
cat(sprintf("  Total PTC-mediated: %d of %d (%.1f%%)\n",
            total_ptc, nrow(c2_results), 100 * total_ptc / nrow(c2_results)))

results <- list(
  c2 = c2_results,
  c4 = c4_results,
  summary = list(
    n_c2 = nrow(c2_results),
    n_c4 = nrow(c4_results),
    group2_n = group2,
    group3_n = group3,
    original_ptc_pos = orig_ptc_pos,
    reclassified_ptc = reclass_ptc,
    total_ptc = total_ptc
  ),
  metadata = list(
    run_timestamp = Sys.time(),
    fasta_path = FASTA_PATH,
    isopair_version = as.character(packageVersion("Isopair")),
    backend = "Isopair::traceReferenceAtg"
  )
)

saveRDS(results, OUTPUT_PATH)
cat("\nSaved to:", OUTPUT_PATH, "\n")
cat("Completed:", format(Sys.time()), "\n")
